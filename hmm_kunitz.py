#!/usr/bin/env python3
"""
hmm_kunitz.py
A pipeline to build a profile HMM for Kunitz-type protease inhibitor domain,
search SwissProt sequences, and evaluate predictions.

Steps:
1. Select human Kunitz sequences for training.
2. Build multiple sequence alignment (MSA) with Clustal Omega.
3. Convert alignment to Stockholm format for HMMER.
4. Build profile HMM from MSA.
5. Prepare validation sets: non-human Kunitz positives and negatives.
6. Run hmmsearch on validation set.
7. Evaluate predictions: confusion matrix, precision, recall, F1, accuracy.
8. Plot performance metrics.


Usage:
    python hmm_kunitz.py --swissprot uniprot_sprot.fasta --outdir results

Required Tools:
- HMMER (hmmbuild, hmmsearch)
- Clustal Omega (clustalo)
- Biopython, scikit-learn, matplotlib

Vanessa EL DEBS
Lab of Bioinformatics - Prof. Capriotti - University of Bologna
"""

import argparse
import subprocess
import os
import sys
import random
import logging

from Bio import SeqIO, AlignIO
from sklearn.metrics import confusion_matrix, precision_score, recall_score, f1_score, accuracy_score
import matplotlib.pyplot as plt

# Setup logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="HMM-based Kunitz domain predictor")
    parser.add_argument("--swissprot", required=True, help="SwissProt FASTA file (uncompressed)")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value cutoff for hmmsearch hits")
    parser.add_argument("--n_neg", type=int, default=50, help="Number of negative sequences to sample")
    parser.add_argument("--outdir", default="results", help="Output directory")
    parser.add_argument("--clustalo", default="clustalo", help="Path to Clustal Omega executable")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    return parser.parse_args()


def run_command(command, description):
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error during {description}: {e}")
        sys.exit(1)


def select_human_kunitz(fasta_file, outdir):
    logger.info("Selecting human Kunitz sequences for training...")
    human_kunitz = [r for r in SeqIO.parse(fasta_file, "fasta")
                    if "kunitz" in r.description.lower() and "homo sapiens" in r.description.lower()]
    if not human_kunitz:
        raise RuntimeError("No human Kunitz sequences found in SwissProt.")
    out_file = os.path.join(outdir, "train", "training_human_kunitz.fasta")
    SeqIO.write(human_kunitz, out_file, "fasta")
    return out_file


def build_msa(input_fasta, clustal_file, clustalo_path):
    logger.info("Building MSA with Clustal Omega...")
    run_command([clustalo_path, "-i", input_fasta, "-o", clustal_file, "--force", "--outfmt", "clustal"],
                "Clustal Omega alignment")


def convert_to_stockholm(clustal_file, stockholm_file):
    alignment = AlignIO.read(clustal_file, "clustal")
    AlignIO.write(alignment, stockholm_file, "stockholm")


def build_hmm(stockholm_file, hmm_file):
    logger.info("Building profile HMM...")
    run_command(["hmmbuild", hmm_file, stockholm_file], "HMM building")


def create_validation_sets(fasta_file, outdir, n_negatives, seed):
    logger.info("Creating validation sets...")
    positives = []
    others = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        desc = record.description.lower()
        if "kunitz" in desc and "homo sapiens" not in desc:
            positives.append(record)
        elif "kunitz" not in desc:
            others.append(record)

    random.seed(seed)
    negatives = random.sample(others, min(n_negatives, len(others)))

    val_dir = os.path.join(outdir, "validation")
    os.makedirs(val_dir, exist_ok=True)

    pos_file = os.path.join(val_dir, "positives.fasta")
    neg_file = os.path.join(val_dir, "negatives.fasta")
    test_file = os.path.join(val_dir, "test_set.fasta")

    SeqIO.write(positives, pos_file, "fasta")
    SeqIO.write(negatives, neg_file, "fasta")
    SeqIO.write(positives + negatives, test_file, "fasta")

    return pos_file, neg_file, test_file


def run_hmmsearch(hmm_file, fasta_file, tblout_file):
    logger.info("Running hmmsearch...")
    run_command(["hmmsearch", "--tblout", tblout_file, hmm_file, fasta_file], "hmmsearch")


def parse_tblout(tblout_file, evalue_cutoff):
    hits = set()
    with open(tblout_file) as f:
        for line in f:
            if not line.startswith("#"):
                parts = line.strip().split()
                if len(parts) > 4:
                    try:
                        evalue = float(parts[4])
                        if evalue <= evalue_cutoff:
                            hits.add(parts[0])
                    except ValueError:
                        continue
    return hits


def evaluate(hits, pos_fa, neg_fa, outdir):
    logger.info("Evaluating predictions...")
    pos_ids = {r.id for r in SeqIO.parse(pos_fa, "fasta")}
    neg_ids = {r.id for r in SeqIO.parse(neg_fa, "fasta")}

    y_true = [1] * len(pos_ids) + [0] * len(neg_ids)
    y_pred = [1 if i in hits else 0 for i in pos_ids] + [1 if i in hits else 0 for i in neg_ids]

    logger.info("\nEvaluation Metrics:")
    logger.info("Confusion Matrix:")
    logger.info(confusion_matrix(y_true, y_pred))
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)
    acc = accuracy_score(y_true, y_pred)
    logger.info(f"Precision: {prec:.3f}\nRecall: {rec:.3f}\nF1 Score: {f1:.3f}\nAccuracy: {acc:.3f}")

    plot_metrics(prec, rec, f1, acc, outdir)


def plot_metrics(precision, recall, f1, accuracy, outdir):
    metrics = [precision, recall, f1, accuracy]
    names = ['Precision', 'Recall', 'F1 Score', 'Accuracy']
    plt.figure(figsize=(6, 4))
    bars = plt.bar(names, metrics, color='skyblue')
    plt.ylim(0, 1)
    plt.ylabel('Score')
    plt.title('Model Performance Metrics')
    for bar, value in zip(bars, metrics):
        plt.text(bar.get_x() + bar.get_width() / 2, value + 0.02, f"{value:.2f}", ha='center')
    plt.tight_layout()
    plot_path = os.path.join(outdir, "metrics", "performance_metrics.png")
    os.makedirs(os.path.dirname(plot_path), exist_ok=True)
    plt.savefig(plot_path)
    plt.close()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(os.path.join(args.outdir, "train"), exist_ok=True)

    training_fasta = select_human_kunitz(args.swissprot, args.outdir)
    clustal_file = os.path.join(args.outdir, "train", "training.aln")
    stockholm_file = os.path.join(args.outdir, "train", "training.sto")
    hmm_file = os.path.join(args.outdir, "train", "kunitz.hmm")
    tblout_file = os.path.join(args.outdir, "validation", "hits.tbl")

    build_msa(training_fasta, clustal_file, args.clustalo)
    convert_to_stockholm(clustal_file, stockholm_file)
    build_hmm(stockholm_file, hmm_file)

    pos_fa, neg_fa, test_fa = create_validation_sets(args.swissprot, args.outdir, args.n_neg, args.seed)
    run_hmmsearch(hmm_file, test_fa, tblout_file)
    hits = parse_tblout(tblout_file, args.evalue)
    evaluate(hits, pos_fa, neg_fa, args.outdir)


if __name__ == "__main__":
    main()
