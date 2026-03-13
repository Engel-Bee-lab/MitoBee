"""
Separate Snakefile for building phylogenetic trees
- download and add any reference mitogenomes to "REPORTS/mitogenome/"
- build a phylogenetic tree using FastTree
"""
import yaml
import os
import glob

"""Parse config"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
Input files and directories
"""
input_dir = config['args']['input']
extn=config['args']['extn']

# Pattern: all files with this extension
pattern = os.path.join(input_dir, f"*.{extn}")
file_paths = sorted(glob.glob(pattern))
# Extract sample names (strip extension only)
sample_names = [re.sub(rf"\.{re.escape(extn)}$", '', os.path.basename(fp)) for fp in file_paths]

print(f"Samples found: {sample_names}")

"""
Defining directories
"""
dir = {}
#declaring output file
try:
    if config['args']['output'] is None:
        dir_out = os.path.join('output')
    else:
	    dir_out = config['args']['output']
except KeyError:
    dir_out = os.path.join('output')

# temp dir
if config['args']['temp_dir'] is None:
    dir_temp = os.path.join(dir_out, "temp")
else:
    dir_temp = config['args']['temp_dir']


#declaring some the base directories
dir_env = os.path.join(workflow.basedir,"envs")
dir_script = os.path.join(workflow.basedir,"scripts")

dir_hostcleaned = os.path.join(dir_out, 'PROCESSING' ,'Host_cleaned')
dir_mitos = os.path.join(dir_out, 'PROCESSING', "mitogenome")
dir_reports = os.path.join(dir_out, 'REPORTS')


"""
Rules
"""
include: os.path.join("rules", "mitos.smk")
include: os.path.join("rules", "alignment.smk")

rule build_alignment_fasta:
    input:
        fasta=expand(os.path.join(input_dir, "{sample}.{extn}"), sample=sample_names, extn=extn)
    output:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome.aln")
    params:
        folder=os.path.join(input_dir),
        concat=os.path.join(dir_hostcleaned, "mitogenome", "all_samples.fasta"),
        concat_clean=os.path.join(dir_hostcleaned, "mitogenome", "all_samples_clean.fasta"),
        host= config['args']['host_seq']
    conda:
        os.path.join(dir_env, "mafft.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.final_fasta} ]; then
            echo "Final alignment fasta already exists. Skipping..."
            exit 0
        else
            cat {params.folder}/*.fasta > {params.concat}
            cat {params.host} >> {params.concat}
            awk '/^>/{{print $1; next}} {{print}}' {params.concat} > {params.concat_clean}
            mafft --auto {params.concat_clean} > {output.final_fasta}
        fi
        """

rule phylo_tree:
    input:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome.aln")
    output:
        tree = os.path.join(dir_reports, "mitogenome_phylo_tree.nwk")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    shell:
        """
            iqtree -s {input.final_fasta} -m GTR+G -bb 3000 -nt AUTO -pre tmp_prefix
            mv tmp_prefix.treefile {output.tree}
            rm tmp_prefix.*
        """

"""Mark target rules"""
rule all:
    input:
        os.path.join(dir_out, "database", "mitos_db" "mitos_downlaoded.txt"),
        expand(os.path.join(dir_mitos, "{sample}_mitogenome", "{sample}_result.faa"), sample=sample_names, extn=extn),
        expand(os.path.join(dir_mitos, "{sample}_mitogenome", "done.txt"), sample=sample_names),
        expand(os.path.join(dir_out, "temp", "{sample}_merged.txt"), sample=sample_names),
        #os.path.join(dir_out, "temp", "aligned_done.txt"),
        #os.path.join(dir_mitos, "mafft", "concatenated_alignment.faa"),
        #os.path.join(dir_mitos, "mafft", "Partitions.txt"),
        #os.path.join(dir_mitos, "mafft", "Partitions.nex")
        #os.path.join(dir_reports, "mitogenome_phylo_tree.nwk")