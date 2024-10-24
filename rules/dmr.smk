rule make_bed:
    input: 
        in_file = annotated_results,
        script = "scripts/make_bed.R"
    params:
        o_dir = OUT_DIR,
        strat = STRATIFIED,
        assoc = ASSOC
    output: 
        results_bed
    conda:
        "../envs/dmr.yaml"
    shell:
        f"""
        Rscript {{input.script}} \
        --input-file {{input.in_file}} \
        --out-dir {{params.o_dir}} \
        --stratified {{params.strat}} \
        --assoc {{params.assoc}} 
        """

rule run_dmr:
    input:
        in_file = results_bed
    params:
        o_prefix = OUT_DIR + ASSOC,
        min_p = MIN_P,
        win_sz = WIN_SZ,
        region_filter = REGION_FILTER
        anno = ANNO
    output: 
        dmr_files
    conda:
        "../envs/dmr.yaml"
    shell:
        f"""
        comb-p pipeline \
		-c 4 \ 
		--seed {{params.min_p}} \
		--dist {{params.win_sz}} \
		-p {{params.assoc}} \
		--region-filter-p {{params.region_filter}} \
		--anno {{params.anno}} \
		{{input.in_file}} 
        """