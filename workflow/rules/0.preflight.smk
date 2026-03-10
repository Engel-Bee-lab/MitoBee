#Snakemake rule to amke sure the input files are correct and all the params are provided to start the workflow

#To start with the input files have to be in a folders
#can be compressed or decompressed
#For now there can be no sub folders, the files have to be saved in the same folder

import os
import glob
import yaml
import re
import sys
import shutil
from metasnek import fastq_finder

"""
CONFIG FILE
"""
configfile: os.path.join(workflow.basedir, "..", "config", "config.yaml")

"""
DIRECTORIES
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

#making directories for each step
#Saving most of the files to PROCESSING, sine they are intermediate files
dir_fastp = os.path.join(dir_out, 'PROCESSING' ,'1_fastp')
dir_hostcleaned = os.path.join(dir_out, 'PROCESSING' ,'Host_cleaned')
dir_hostsearch = os.path.join(dir_out, 'PROCESSING' ,'Reference_search')
dir_reports = os.path.join(dir_out, 'REPORTS')


"""
CHECK INPUT FILES
"""
#get user inputs 
input_dir = config['args']['input']
extn=config['args']['extn']
pattern_r1 = config['args']['pattern_r1']
pattern_r2 = config['args']['pattern_r2']

# Convert {sample} into wildcard *
glob_pattern = f"*{pattern_r1}.{extn}"
file_paths = sorted(glob.glob(os.path.join(input_dir, glob_pattern)))

sample_names = []
for fp in file_paths:
    filename = os.path.basename(fp)

    # Extract sample name by reversing the pattern
    sample = filename.replace(f"{pattern_r1}.{extn}", "")
    sample_names.append(sample)

sample_names = list(dict.fromkeys(sample_names))

print(f"Samples are {sample_names}")
FQEXTN = extn[0]
PATTERN_R1 = f'{{sample}}{pattern_r1}.{extn}'
PATTERN_R2 = f'{{sample}}{pattern_r2}.{extn}'

"""ONSTART/END/ERROR
Tasks to perform at various stages the start and end of a run.
"""
def copy_log_file():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    target_log = os.path.join(dir['log'], os.path.basename(current_log))
    shutil.copy(current_log, target_log)

dir = {'log': os.path.join(dir_out, 'logs')}
onstart:
    """Cleanup old log files before starting"""
    if os.path.isdir(dir["log"]):
        oldLogs = filter(re.compile(r'.*.log').match, os.listdir(dir["log"]))
        for logfile in oldLogs:
            os.unlink(os.path.join(dir["log"], logfile))
onsuccess:
    """Print a success message"""
    sys.stderr.write('\n\n Workflow ran successfully!\n\n')
onerror:
    """Print an error message"""
    sys.stderr.write('\n\n Workflow run failed\n\n')
