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
