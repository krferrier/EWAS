rule run_combined_ewas:
    input:
        script = "scripts/ewas.R",
        pheno_file = CW.pheno,
        methyl_file = CW.mvals
    params:
        assoc_var = CW.assoc_var,
        stratified = CW.stratified,
        cs = CW.chunk_size,
        pt = CW.processing_type,
        n_workers = CW.n_workers,
        o_dir = CW.out_dir,
        o_type = CW.out_type
    output: 
        CW.raw_results
    log:
        CW.log_path("run_combined_ewas")
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        # Everything below -- stdout and stderr of every command -- goes to the
        # rule's log. Snakemake names the log in its error report if a job fails.
        exec >{log} 2>&1
        export R_PROGRESSR_ENABLE=TRUE 
        Rscript {input.script} \
        --pheno {input.pheno_file} \
        --methyl {input.methyl_file} \
        --assoc {params.assoc_var} \
        --stratified {params.stratified} \
        --chunk-size {params.cs} \
        --processing-type {params.pt} \
        --workers {params.n_workers} \
        --out-dir {params.o_dir} \
        --out-type {params.o_type}
        """


rule run_bacon:
    input:
        in_file = CW.raw_results,
        script = "scripts/run_bacon.R"
    params:
        o_dir = CW.out_dir,
        o_type = CW.out_type,
        o_prefix = ""
    output: 
        CW.bacon_results,
        CW.bacon_plot_files()
    log:
        CW.log_path("run_bacon")
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        # Everything below -- stdout and stderr of every command -- goes to the
        # rule's log. Snakemake names the log in its error report if a job fails.
        exec >{log} 2>&1
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {params.o_dir} \
        --out-prefix {params.o_prefix} \
        --out-type {params.o_type} \
        """
    
    