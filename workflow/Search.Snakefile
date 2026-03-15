"""
Separate Snakefile for gene search workflow
- User must provide a directory with reference genes
- Run mmseqs2 to generate a table of hits for each metagenome in the reference genes

This would allow the user to determine which is the most closely related genome. Help make an informed decison about which reference genome \
    to use for the assembly and binning steps. It also allows the user to determine if there are any closely related reference genomes \
    that could be used for downstream analyses (e.g. pangenome analysis, phylogenetics, etc.).
"""

import yaml
import os
import glob

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""Rules"""
include: os.path.join("rules", "0.preflight.smk")
include: os.path.join("rules", "1.qc.smk")
include: os.path.join("rules", "search.2.mapping-search.smk")
include: os.path.join("rules", "search.2.mappingGene.smk")

"""Mark target rules"""
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn

"""
Defining the targets dictionary
"""
targets ={'qc':[], 'hostsearch':[]}


if config['args']['mode'] == "mitogenome":
    for sample in sample_names:
        targets['qc'].append(expand(os.path.join(dir_fastp,"{sample}.stats.html"), sample=sample)),
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_temp.bam"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch, "{sample}_all_idxstats.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_primary_idxstats.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_primary_coverage.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_strict_idxstats.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_reports, "Reference_search", "{sample}_host_ranking.tsv"), sample=sample))
elif config['args']['mode'] == "gene":
    for sample in sample_names:
        targets['qc'].append(expand(os.path.join(dir_fastp,"{sample}.stats.html"), sample=sample)),
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_temp_gene.bam"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch, "{sample}_all_idxstats_gene.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_primary_idxstats_gene.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_primary_coverage_gene.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_hostsearch,"{sample}_strict_idxstats_gene.txt"), sample=sample))
        targets['hostsearch'].append(expand(os.path.join(dir_reports, "Reference_search", "{sample}_host_ranking_gene.tsv"), sample=sample))

@targetRule
rule all:
    input:
        targets['qc'],
        targets['hostsearch']