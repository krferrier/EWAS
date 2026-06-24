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
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc}
        """