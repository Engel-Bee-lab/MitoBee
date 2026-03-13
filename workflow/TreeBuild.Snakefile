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
include: os.path.join("rules", "tree.1.mitos.smk")
include: os.path.join("rules", "tree.2.alignment.smk")
include: os.path.join("rules", "tree.3.genetrees.smk")

"""Mark target rules"""
rule all:
    input:
        os.path.join(dir_out, "database", "mitos_db" "mitos_downlaoded.txt"),
        expand(os.path.join(dir_mitos, "{sample}_mitogenome", "{sample}_result.faa"), sample=sample_names, extn=extn),
        expand(os.path.join(dir_mitos, "{sample}_mitogenome", "done.txt"), sample=sample_names),
        os.path.join(dir_out, "temp", "genes_merged.txt"),
        os.path.join(dir_out, "temp", "aligned_done.txt"),
        os.path.join(dir_mitos, "mafft", "concatenated_alignment.faa"),
        os.path.join(dir_mitos, "mafft", "Partitions.txt"),
        os.path.join(dir_mitos, "mafft", "Partitions.nex"),
        os.path.join(dir_reports, "mitogenome_phylo_tree.nwk"),
        os.path.join(dir_reports, "mitogenome_phylo_tree.log"),
        os.path.join(dir_out, "temp", "gene_trees_done")