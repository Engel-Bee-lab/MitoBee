"""
Merges all the proteins commonly found in all of the samples.
So if a mitocgenome is partial and missing a gene, then that gene will not be included in the alignment and tree
"""
rule merge_proteins:
    input:
        os.path.join(dir_mitos, "{sample}_mitogenome", "done.txt")
    output:
        os.path.join(dir_out, "temp", "{sample}_merged.txt")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    params:
        indir=os.path.join(dir_mitos, "{sample}_mitogenome"),
        folder=os.path.join(dir_mitos, "mafft"),
        sample="{sample}"
    shell:
        """
        mkdir -p {params.folder}
        rm -rf {params.folder}/*.faa
        for gene in atp6 atp8 cob cox1 cox2 cox3 nad1 nad2 nad3 nad4 nad4l nad5 nad6
        do
            f={params.indir}/{params.sample}_updated_result.faa.split/{params.sample}_updated_result.part_${{gene}}.faa

            if [ -f "$f" ]; then
                cat "$f" >> {params.folder}/${{gene}}.faa
            else
                echo "File $f does not exist."
            fi
        done

        touch {output}
        """

import os 
import glob

"""Aligns the merged proteins using mafft"""
rule mafft:
    input:
        expand(os.path.join(dir_out, "temp", "{sample}_merged.txt"), sample=sample_names)
    output:
        folder=os.path.join(dir_out, "temp", "aligned_done.txt")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    params:
        indir=os.path.join(dir_mitos, "mafft"),
    resources:
        mem_mb =config['resources']['smalljob']['mem_mb'],
        runtime = config['resources']['smalljob']['runtime']
    threads: 
        config['resources']['smalljob']['threads']
    shell:
        """
        for f in {params.indir}/*.faa
        do
            gene=$(basename "$f" .faa)
            mafft --auto --thread {threads} "$f" > {params.indir}/${{gene}}_aligned.faa
        done
        touch {output.folder}
        """

"""Concatenates the aligned proteins into a single file for tree building"""
rule concat_alignments:
    input:
        folder=os.path.join(dir_out, "temp", "aligned_done.txt") 
    output:
        fasta=os.path.join(dir_mitos, "mafft", "concatenated_alignment.faa"),
        partitions_txt=os.path.join(dir_mitos, "mafft", "Partitions.txt"),
        partitions_nex=os.path.join(dir_mitos, "mafft", "Partitions.nex")
    params:
        script="workflow/scripts/concat_alignments.py",
        indir=os.path.join(dir_mitos, "mafft"), # folder with *_aligned.faa,
        tmp_dir=os.path.join(dir_out, "temp")
    conda:
        os.path.join(dir_env, "seqkit.yaml")  # make sure Biopython is installed
    shell:
        """
        ./{params.script} \
            {params.tmp_dir} \
            {output.fasta} \
            {output.partitions_txt} \
            {output.partitions_nex}
        """