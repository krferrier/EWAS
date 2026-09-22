rule plot_results:
    input: 
        in_file = CW.annotated_results,
        script = "scripts/plots.R"
    params:
        o_dir = CW.out_dir,
        strat = CW.stratified,
        assoc = CW.assoc_var
    output: 
        CW.manhattan_qq_plot
    log:
        CW.log_path("plot_results")
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
        --stratified {params.strat} \
        --assoc {params.assoc}
        """