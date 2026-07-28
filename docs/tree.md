# Tree construction

This workflow is a modular Snakemake pipeline for building phylogenetic trees from mitochondrial genomes or mitochondrial protein-coding genes. It supports two major modes:
- mitogenome mode
  Builds from mitogenome nucl seuqences. 
- Mito proteins mode
  Builds concatenated protein alignments and phylogenetic trees from annotated mitochondrial genomes.
- Gene mode
  Builds DNA-based phylogenetic trees directly from user-provided gene FASTA files.

- The pipeline uses:
    - MITOS2 for mitochondrial genome annotation
    - MAFFT for multiple sequence alignment
    - IQ-TREE for phylogenetic inference
    - SeqKit for FASTA processing

