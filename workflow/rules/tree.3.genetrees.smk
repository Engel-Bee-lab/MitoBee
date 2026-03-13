"""
Build a tree with the provided concatenated alignment and partition files using IQ-TREE
User provided which genes, otherwise defaults to cox1, and cob 
To change this, edit the config.yaml file and add the genes to be used for tree building under tree:genes
"""
gene_list = config["tree"]["genes"] 

rule build_gene_tree:
    input:
        folder=os.path.join(dir_out, "temp", "aligned_done.txt"),
        faa=os.path.join(dir_mitos, "mafft", "{gene}.faa")
    output:
        tree_file=os.path.join(dir_reports, "gene_trees", "{gene}.treefile"),
        log_file=os.path.join(dir_reports, "gene_trees", "{gene}.log")
    conda:
        os.path.join(dir_env, "mafft.yaml")
    params:
        aln_file=os.path.join(dir_mitos, "mafft", "{gene}_aligned.faa"),
        iqtree_dir=os.path.join(dir_mitos, "gene_trees_tmp"),
        output_dir=os.path.join(dir_reports, "gene_trees"),
        gene="{gene}"
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime']
    threads:
        config['resources']['smalljob']['threads']
    shell:
        """
        mkdir {params.iqtree_dir} 
        mkdir {params.output_dir}

        if [ -f {params.aln_file} ]; then
            echo "Alignment file for gene {params.gene} found. Proceeding with tree building."
            iqtree -s {params.aln_file} -m MFP -bb 1000 -nt {threads} -pre "{params.gene}"
            mv {params.gene}.treefile {params.output_dir}/.
            mv {params.gene}.log {params.output_dir}/.
            mv {params.gene}.* {params.iqtree_dir}/.
        else
            echo "Alignment file for gene {params.gene} not found. Skipping tree building for this gene."
            touch {output.tree_file}
            touch {output.log_file}
        fi
        """