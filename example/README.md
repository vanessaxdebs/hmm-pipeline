# Example: Kunitz Domain Pipeline Test Case

This folder provides a minimal working example to test the pipeline.

## Files

- `kunitz_example.fasta`: 5 protein sequences (3 with Kunitz domains, 2 without)
- `msa.aln`: Multiple sequence alignment of the 3 known Kunitz-domain proteins
- `kunitz_example.hmm`: Profile HMM built from the alignment
- `hmmsearch_results.txt`: Results from searching the HMM against the 5 sequences

## Instructions

To regenerate the HMM and search results:

```bash
hmmbuild example/kunitz_example.hmm example/msa.aln
hmmsearch --tblout example/hmmsearch_results.txt example/kunitz_example.hmm example/kunitz_example.fasta

