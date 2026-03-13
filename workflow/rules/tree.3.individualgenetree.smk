"""
Buidling from tree.2.aligment.smk, building gene trees 
Only of the genes present in all the genomes 
"""

"""Build a tree with the provided concatenated alignment and partition files using IQ-TREE
User provided which genes, otherwise defaults to cox1, and cob 
To change this, edit the config.yaml file and add the genes to be used for tree building under tree:genes
"""
rule build_tree:
    """
    Build trees for individual genes using IQ-TREE.
    Skips genes that were not aligned (e.g., missing in some genomes).
    """
    input:
        folder=os.path.join(dir_out, "temp", "aligned_done.txt"),
        gene_list=config["tree"]["genes"],
    output:
        tree_dir=os.path.join(dir_out, "temp", "gene_trees_done")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    params:
        iqtree_dir=os.path.join(dir_mitos, "gene_trees_tmp"),
        output_dir=os.path.join(dir_reports, "gene_trees")
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime']
    threads:
        config['resources']['smalljob']['threads']
    run:
        import os

        # make output directories
        os.makedirs(output.tree_dir, exist_ok=True)
        os.makedirs(params.iqtree_dir, exist_ok=True)

        # iterate over each gene
        for gene in input.gene_list:
            aln_file = os.path.join(dir_mitos, "mafft", f"{gene}_aligned.faa")
            
            # check if alignment exists
            if not os.path.exists(aln_file):
                print(f"[WARN] Alignment file for gene '{gene}' not found. Skipping tree for this gene.")
                continue
            
            print(f"[INFO] Building tree for gene: {gene}")

            tmp_prefix = os.path.join(params.iqtree_dir, f"{gene}_tmp")

            # run IQ-TREE
            shell(
                f"iqtree -s {aln_file} -m MFP -bb 1000 -nt {threads} -pre {tmp_prefix}"
            )

            # move outputs to final directory
            tree_out = os.path.join(params.tree_dir, f"{gene}.treefile")
            log_out = os.path.join(params.tree_dir, f"{gene}.log")
            shell(
                f"mv {tmp_prefix}.treefile {tree_out}; "
                f"mv {tmp_prefix}.log {log_out}; "
                f"mv {tmp_prefix}.iqtree {params.iqtree_dir}/."
            )
        touch {output}
        """