[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![](https://img.shields.io/static/v1?label=CLI&message=Snaketool&color=blueviolet)](https://github.com/beardymcjohnface/Snaketool)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/npbhavya/MitoBee?branch=main)

<img src="mitobee/mitobeelogo.png" width="220" alt="MitoBee logo">

# MitoBee

## Snakemake workflow to get mitogenomes from metagenomic data
Still under development! Stable release out as a version, but only if there is a closely related mitogenome available.

Documentation: [Read the Docs](https://mitobee.readthedocs.io/en/latest/)

### Install 

**Source install**
Run the below commands:

    git clone https://github.com/npbhavya/MitoBee.git
    cd MitBee
    mamba create -y -n mitobee python=3.13
    conda activate mitobee
    pip install -e . 

### Input files
Input files:
- Input directory with metagenomes
- Reference directory
  - If running `run`  or `tree` module, provide a **(one)** reference genome. 
  - If running `gene` module, provide a reference gene set

### Output files
Output files: Provide the output folder, contains subdirectories
- PROCESSING: Folder containing intermediate files
- REPORTS: Final results including the mitogenome fasta files from (hopefully) each metagenome sample \
      Also inlcudes the QC reports, to include stats on how many reads were processed, and not
