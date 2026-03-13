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

for gene in genes:
    aln_file = os.path.join(aligned_dir, f"{gene}_aligned.faa")
    if os.path.exists(aln_file):
        print(f"Processing {gene}")
        for record in SeqIO.parse(aln_file, "fasta"):
            seqs_by_sample[record.id].append(str(record.seq))
    else:
        # gene missing, skip it
        print(f"Warning: {gene}_aligned.faa not found, skipping this gene")

# -----------------------
# Write concatenated alignment
# -----------------------
with open(output_file, "w") as outfh:
    for sample, seqs in seqs_by_sample.items():
        concatenated = "".join(seqs)
        outfh.write(f">{sample}\n{concatenated}\n")

print(f"Concatenated alignment written to {output_file}")
