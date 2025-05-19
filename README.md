# hmm-pipeline

## Kunitz Domain Prediction with Profile Hidden Markov Models

A bioinformatics pipeline for building and validating a Profile Hidden Markov Model (HMM) to predict Kunitz-type protease inhibitor domains in protein sequences.

---

## Table of Contents
- [Project Overview](#project-overview)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Methodology](#methodology)
- [License](#license)

---

## Project Overview

This project implements a computational workflow to:
1. Build a profile HMM from human Kunitz domain sequences.
2. Validate the model against non-human Kunitz proteins.
3. Evaluate prediction performance using standard metrics.

Additionally, this pipeline supports generating entropy plots to analyze the sequence conservation in multiple sequence alignments (MSA) of Kunitz-type protease inhibitors.

Key features:
- Processes SwissProt-curated sequences.
- Uses Clustal Omega for multiple sequence alignment.
- Implements HMMER3 for model building and searching.
- Generates comprehensive performance metrics.
- **Entropy Calculation**: Calculates entropy values for each position in the MSA to visualize sequence conservation.
- **Plot Generation**: Generates a plot visualizing the entropy values using `matplotlib`.

---

## Installation

### Prerequisites
- Python 3.8+
- HMMER (v3.3+)
- Clustal Omega (v1.2.4+)

### Setup
```bash
# Clone repository
git clone https://github.com/YOUR-USERNAME/kunitz-hmm.git
cd kunitz-hmm

# Install Python dependencies
pip install -r requirements.txt

# Install system tools (Linux/macOS)
sudo apt-get install hmmer clustalo  # Debian/Ubuntu
brew install hmmer clustal-omega    # macOS
Usage
Quick Start
bash
Copy code
# Download SwissProt data
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Run pipeline
python hmm_kunitz.py --swissprot uniprot_sprot.fasta --outdir results
Entropy Plot Usage
Run the entropy calculation script:

bash
Copy code
python3 entropy_plot.py
This will generate a plot (entropy_plot.png) showing the entropy values for each position in the MSA.

Command Line Options
Parameter	Description	Default
--swissprot	Path to SwissProt FASTA	Required
--outdir	Output directory	results
--evalue	E-value threshold	1e-5
--n_neg	Number of negative samples	50
--seed	Random seed	42

Output
The pipeline generates:

Copy code
results/
├── train/
│   ├── training_human_kunitz.fasta
│   ├── training.aln
│   ├── training.sto
│   └── kunitz.hmm
├── validation/
│   ├── positives.fasta
│   ├── negatives.fasta
│   ├── test_set.fasta
│   └── hits.tbl
└── metrics/
    └── performance_metrics.png
└── entropy_plot.png
Methodology
Data Collection: Human Kunitz sequences are extracted from SwissProt.

Alignment: Multiple sequence alignment with Clustal Omega.

Model Building: HMM construction using hmmbuild.

Validation: Performance evaluation against curated test sets.

Metrics: Precision, recall, F1-score, and accuracy calculations.

Entropy Plot: Visualize sequence conservation using entropy values calculated from the MSA.

License
This project is licensed under the MIT License.

For questions, please contact: vanessaxdebs@gmail.com# hmm-pipeline

```markdown
# Kunitz Domain Prediction with Profile Hidden Markov Models

A bioinformatics pipeline for building and validating a Profile Hidden Markov Model (HMM) to predict Kunitz-type protease inhibitor domains in protein sequences.

## Table of Contents
- [Project Overview](#project-overview)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [Methodology](#methodology)
- [License](#license)

## Project Overview

This project implements a computational workflow to:
1. Build a profile HMM from human Kunitz domain sequences
2. Validate the model against non-human Kunitz proteins
3. Evaluate prediction performance using standard metrics

Key features:
- Processes SwissProt-curated sequences
- Uses Clustal Omega for multiple sequence alignment
- Implements HMMER3 for model building and searching
- Generates comprehensive performance metrics

## Installation

### Prerequisites
- Python 3.8+
- HMMER (v3.3+)
- Clustal Omega (v1.2.4+)

### Setup
```bash
# Clone repository
git clone https://github.com/YOUR-USERNAME/kunitz-hmm.git
cd kunitz-hmm

# Install Python dependencies
pip install -r requirements.txt

# Install system tools (Linux/macOS)
sudo apt-get install hmmer clustalo  # Debian/Ubuntu
brew install hmmer clustal-omega    # macOS
```

## Usage

### Quick Start
```bash
# Download SwissProt data
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Run pipeline
python hmm_kunitz.py --swissprot uniprot_sprot.fasta --outdir results
```

### Command Line Options
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--swissprot` | Path to SwissProt FASTA | Required |
| `--outdir` | Output directory | `results` |
| `--evalue` | E-value threshold | 1e-5 |
| `--n_neg` | Number of negative samples | 50 |
| `--seed` | Random seed | 42 |

## Output

The pipeline generates:
```
results/
├── train/
│   ├── training_human_kunitz.fasta
│   ├── training.aln
│   ├── training.sto
│   └── kunitz.hmm
├── validation/
│   ├── positives.fasta
│   ├── negatives.fasta
│   ├── test_set.fasta
│   └── hits.tbl
└── metrics/
    └── performance_metrics.png
```

## Methodology

1. **Data Collection**: Human Kunitz sequences extracted from SwissProt
2. **Alignment**: Multiple sequence alignment with Clustal Omega
3. **Model Building**: HMM construction using hmmbuild
4. **Validation**: Performance evaluation against curated test sets
5. **Metrics**: Precision, recall, F1-score, and accuracy calculations

## License

This project is licensed under the MIT License

---

*For questions, please contact: vanessaxdebs@gmail.com*
```

