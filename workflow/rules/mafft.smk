rule merge_proteins:
    input:
        os.path.join(dir_mitos, "{sample}_mitogenome", "done.txt")
    output:
        os.path.join(dir_out, "temp", "{sample}_merged.txt")
    conda:
        os.path.join(envs_dir, "mafft.yaml")
    params:
        indir=os.path.join(dir_mitos, "{sample}_mitogenome", "{sample}_consensus_updated_result.faa.split"),
        folder=os.path.join(dir_mitos, "mafft"),
        sample="{sample}"
    shell:
        """
        mkdir -p {params.folder}

        for gene in atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6
        do
            cat {params.indir}/{params.sample}_consensus_updated_result.part_${{gene}}.faa >> {params.folder}${{gene}}.faa
        done
        touch {output}
        """