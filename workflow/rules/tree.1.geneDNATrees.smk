"""
Rules to build gene DNA trees.
"""
rule build_alignment_fasta:
    input:
        fasta=expand(os.path.join(input_dir, "{sample}.{extn}"), sample=sample_names, extn=extn)
    output:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenes", "final_mitogenome.aln")
    params:
        folder=os.path.join(input_dir),
        concat=os.path.join(dir_hostcleaned, "mitogenes", "all_samples.fasta"),
        concat_clean=os.path.join(dir_hostcleaned, "mitogenes", "all_samples_clean.fasta"),
    conda:
        os.path.join(dir_env, "mafft.yaml")
    shell:
        """
        set -euo pipefail
        if [ -f {output.final_fasta} ]; then
            echo "Final alignment fasta already exists. Skipping..."
            exit 0
        else
            cat {input.fasta} > {params.concat}
            awk '/^>/{{print $1; next}} {{print}}' {params.concat} > {params.concat_clean}
            mafft --auto {params.concat_clean} > {output.final_fasta}
        fi
        """
    
rule phylo_tree:
    input:
        final_fasta = os.path.join(dir_hostcleaned, "mitogenome", "final_mitogenome.aln")
    output:
        tree = os.path.join(dir_reports, "gene_trees", "mitogenome_phylo_tree.treefile")
    params:
        odir=os.path.join(dir_reports, "gene_trees"),
        prefix="gene_tree"
    conda:
        os.path.join(dir_env, "mafft.yaml")
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime']
    threads:
        config['resources']['smalljob']['threads']
    shell:
        """
            iqtree -s {input.final_fasta} -m MFP -bb 1000 -nt {threads} -pre {params.prefix}
            mv {params.prefix}* {params.odir}/.
            touch {output.tree}
        """