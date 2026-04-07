import os
import glob
import yaml

# -------------------------
# CONFIG
# -------------------------
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

# -------------------------
# DIRECTORIES
# -------------------------
dir_out = config.get('args', {}).get('output', 'output')
dir_fastp = os.path.join(dir_out, 'PROCESSING', '1_fastp')
dir_env = os.path.join(workflow.basedir, "envs")

os.makedirs(dir_fastp, exist_ok=True)

# -------------------------
# INPUT FILES
# -------------------------
input_dir = config['args']['input']
extn = config['args']['extn']
pattern_r1 = config['args']['pattern_r1']
pattern_r2 = config['args']['pattern_r2']

# Find all R1 and R2 files
r1_files = glob.glob(os.path.join(input_dir, f"*{pattern_r1}*.{extn}"))
r2_files = glob.glob(os.path.join(input_dir, f"*{pattern_r2}*.{extn}"))

# Extract sample names
def extract_sample_names(file_list, pattern, ext):
    samples = []
    for f in file_list:
        name = os.path.basename(f)
        sample = name.replace(pattern, "").replace(f".{ext}", "")
        samples.append(sample)
    return set(samples)

samples_r1 = extract_sample_names(r1_files, pattern_r1, extn)
samples_r2 = extract_sample_names(r2_files, pattern_r2, extn)
samples = sorted(samples_r1 & samples_r2)  # only keep pairs

# Create sample dictionary
config["sample_names"] = {
    sample: {
        "r1": os.path.join(input_dir, f"{sample}{pattern_r1}.{extn}"),
        "r2": os.path.join(input_dir, f"{sample}{pattern_r2}.{extn}")
    } for sample in samples
}

# -------------------------
# WILDCARD CONSTRAINTS
# -------------------------
wildcard_constraints:
    sample="|".join(samples)

# -------------------------
# FASTP RULE
# -------------------------
rule fastp:
    input:
        r1 = lambda wc: config["sample_names"][wc.sample]["r1"],
        r2 = lambda wc: config["sample_names"][wc.sample]["r2"]
    output:
        r1 = os.path.join(dir_fastp,"{sample}_R1.fastq.gz"),
        r2 = os.path.join(dir_fastp,"{sample}_R2.fastq.gz"),
        stats = os.path.join(dir_fastp,"{sample}.stats.json"),
        html = os.path.join(dir_fastp,"{sample}.stats.html")
    conda:
        os.path.join(dir_env, "fastp.yaml")
    threads: config['resources']['smalljob']['threads']
    resources:
        mem_mb = config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --thread {threads} \
              --json {output.stats} \
              --html {output.html}
        """