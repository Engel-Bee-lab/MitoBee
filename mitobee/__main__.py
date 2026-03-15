"""
Entrypoint for this workflow

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click
from shutil import copyfile
from snaketool_utils.cli_utils import OrderedCommands, run_snakemake, copy_config, echo_click

"""Get the filepath to a Snaketool system file (relative to __main__.py)"""
PACKAGE_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.dirname(PACKAGE_DIR)  # one level up

def snake_base(rel_path):
    """Get the filepath to a project file relative to the repo root."""
    return os.path.join(PROJECT_ROOT, rel_path)

def get_version():
    try:
        from mitobee._version import version
    except Exception:
        version = "0.1.0"
    return version


def print_citation():
    with open(snake_base("../MitoBee.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value

# This is to ensure that the output directory gets generated if it hasnt already been generated
def ensure_directory_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option('--input', '_input', help='Directory of reads', type=click.Path(), required=False, default='testReads/paired', show_default=True),
        click.option('--extn', 'extn',  help='Reads extension; fastq, fq, fastq.gz', type=click.Path(), required=False, default='fastq', show_default=True),
        click.option('--pattern_r1', 'r1', help='Pattern to identify R1 reads (for paired-end data)', default='_R1', show_default=True),
        click.option('--pattern_r2', 'r2', help='Pattern to identify R2 reads (for paired-end data)', default='_R2', show_default=True),
        click.option('--host_seq', 'host_seq', help='Path to host genome index for host read removal', type=click.Path(), show_default=True),
        click.option('--mode', type=click.Choice(['mitogenome', 'gene'], case_sensitive=False), help='Reference search mode either gene or mitochondrial genome', required=False),
        click.option('--ref_set', 'refdb', help='Path to reference database of mitochondrial genomes or mitochondrial genes', type=click.Path(), required=False, show_default=True),
        click.option('--output', 'output', help='Output directory', type=click.Path(),
                     default='output', show_default=True),
        click.option("--configfile", default="config.yaml", show_default=False, callback=default_to_output,
                     help="Custom config file [default: config.yaml]"),
        click.option('--threads', help='Number of threads to use', default=1, show_default=True),
        click.option('--profile', help='Snakemake profile', default=None, show_default=False),
        click.option('--db_dir', 'db_dir', help='Custom database directory', type=click.Path(), required=False),
        click.option('--temp-dir', help='Temp directory', required=False),
        click.option('--snake-default', multiple=True,
                     default=['--rerun-incomplete', '--printshellcmds', '--nolock', '--show-failed-logs'],
                     help="Customise Snakemake runtime args", show_default=True),
        click.option("--log", default="mitobee.log", callback=default_to_output, hidden=True,),
        click.option('--use-conda', default=True, help='Use conda for Snakemake rules',
                     show_default=True),
        click.option('--conda-frontend', default='mamba', help='Use mamba for Snakemake rules',
                     show_default=True),
        click.option('--conda-prefix', default=snake_base(os.path.join('workflow', 'conda')),
                     help='Custom conda env directory', type=click.Path(), show_default=True),
        click.option("--system-config", default=snake_base(os.path.join("..","config", "config.yaml")),hidden=True,),
        click.argument("snake_args", nargs=-1),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """ Extract and assemble host mitochondrial genomes from metagenomic sequencing data.
    \b
    For more options, run:
    mitobee --help"""
    pass


help_msg_run = """
\b
RUN EXAMPLES 
mitobee run --input <input directory with metagenome reads> --pattern_r1 R1 --patern_r2 R2 --extn fq --host_seq <path to host mitochondrial genome> --sequencing paired --output <output directory> -k
"""
@click.command(epilog=help_msg_run, 
    context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True)
    )

@common_options
def run(_input, extn, r1, r2, host_seq, output, temp_dir, configfile, conda_frontend, **kwargs):
    """Run mitobee workflow"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    if not host_seq:
        raise click.UsageError("Option --host_seq is required for 'mitobee run'.")
    
    merge_config = {
        "args": {
            "input": _input, 
            "output": output, 
            "extn": extn,
            "pattern_r1": r1,
            "pattern_r2": r2,
            "host_seq": host_seq,
            "configfile": configfile,
            "temp_dir": temp_dir,
        }
    }

    snake_default = list(kwargs.get('snake_default', []))
    if conda_frontend and not any('--conda-frontend' in str(arg) for arg in snake_default):
        snake_default.extend(['--conda-frontend', conda_frontend])
    kwargs['snake_default'] = tuple(snake_default)

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'Snakefile')),
        configfile=configfile,
        merge_config=merge_config,
        **kwargs
    )

help_msg_run = """
\b
TREE EXAMPLES 
mitobee tree --input output/REPORTS/mitogenome --extn fasta --host_seq test-files/am-dh4.fasta --output output -k
"""
@click.command(epilog=help_msg_run, 
    context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True)
    )

@common_options
def tree(_input, extn, output, temp_dir, configfile, conda_frontend, **kwargs):
    """Run mitobee workflow"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))

    merge_config = {
        "args": {
            "input": _input, 
            "output": output, 
            "extn": extn,
            "configfile": configfile,
            "temp_dir": temp_dir,
        }
    }

    snake_default = list(kwargs.get('snake_default', []))
    if conda_frontend and not any('--conda-frontend' in str(arg) for arg in snake_default):
        snake_default.extend(['--conda-frontend', conda_frontend])
    kwargs['snake_default'] = tuple(snake_default)

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'TreeBuild.Snakefile')),
        configfile=configfile,
        merge_config=merge_config,
        targets=["all"], 
        **kwargs
    )

help_msg_run = """
\b
SEARCH EXAMPLES 
mitobee search --input <input directory with metagenome reads> --extn fastq.gz --pattern_r1 R1 --patern_r2 R2 --ref_set test-files/ref-db --output output -k
"""
@click.command(epilog=help_msg_run, 
    context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True)
    )

@common_options
def search(_input, extn, r1, r2, refdb, mode, output, temp_dir, configfile, conda_frontend, **kwargs):
    """Run mitobee workflow"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))
    
    merge_config = {
        "args": {
            "input": _input, 
            "output": output, 
            "extn": extn,
            "pattern_r1": r1,
            "pattern_r2": r2,
            "ref_set": refdb,
            "configfile": configfile,
            "temp_dir": temp_dir,
            "mode": mode,
        }
    }

    snake_default = list(kwargs.get('snake_default', []))
    if conda_frontend and not any('--conda-frontend' in str(arg) for arg in snake_default):
        snake_default.extend(['--conda-frontend', conda_frontend])
    kwargs['snake_default'] = tuple(snake_default)

    # run!
    run_snakemake(
        snakefile_path=snake_base(os.path.join('workflow', 'Search.Snakefile')),
        configfile=configfile,
        merge_config=merge_config,
        targets=["all"], 
        **kwargs
    )

@click.command()
@click.option('--configfile', default='config.yaml', help='Copy template config to file', show_default=True)
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile, system_config=snake_base(os.path.join('config', 'config.yaml')))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(tree)
cli.add_command(search)
cli.add_command(config)
cli.add_command(citation)

def main():
    cli()

if __name__ == '__main__':
    main()
