#!/usr/bin/env python3
import os
import sys
from Bio import SeqIO
from collections import defaultdict

# -----------------------
# Arguments
# -----------------------
# 1: folder containing aligned genes
# 2: output concatenated alignment file
# -----------------------
aligned_dir = sys.argv[1]
output_file = sys.argv[2]

#3+: optional fixed gene order
#below I have the mitogenomes order
genes = ["atp6", "atp8", "cob", "cox1", "cox2", "cox3",
         "nad1", "nad2", "nad3", "nad4", "nad4L", "nad5", "nad6"]

# -----------------------
# Read alignments
# -----------------------
seqs_by_sample = defaultdict(list)
gene_lengths = {}  # length of each gene alignment

for gene in genes:
    aln_file = os.path.join(aligned_dir, f"{gene}_aligned.faa")
    if os.path.exists(aln_file):
        print(f"Processing {gene}")
        records = list(SeqIO.parse(aln_file, "fasta"))
        gene_lengths[gene] = len(records[0].seq)  # assume all sequences same length
        for record in records:
            seqs_by_sample[record.id].append(str(record.seq))
    else:
        print(f"Warning: {gene}_aligned.faa not found, adding gaps")
        # store the gene length for gaps; use max length so all samples match
        gene_lengths[gene] = max(gene_lengths.values(), default=50)
        for sample in seqs_by_sample:
            seqs_by_sample[sample].append("-" * gene_lengths[gene])

# -----------------------
# Fill in missing samples for genes that appear later
# -----------------------
for gene in genes:
    if gene not in gene_lengths:
        continue
    for sample in seqs_by_sample:
        if len(seqs_by_sample[sample]) < genes.index(gene) + 1:
            seqs_by_sample[sample].append("-" * gene_lengths[gene])

# -----------------------
# Write concatenated alignment
# -----------------------
with open(output_file, "w") as outfh:
    for sample, seqs in seqs_by_sample.items():
        concatenated = "".join(seqs)
        outfh.write(f">{sample}\n{concatenated}\n")

print(f"Concatenated alignment written to {output_file}")
