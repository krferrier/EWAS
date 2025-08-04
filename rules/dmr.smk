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
        """
        Rscript {input.script} \
        --results {input.in_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} 
        """

rule run_dmr:
    input:
        in_file = results_bed
    params:
        o_prefix = OUT_DIR +  "dmr/" + ASSOC + "_ewas",
        min_p = MIN_P,
        win_sz = WIN_SZ,
        region_filter = REGION_FILTER,
        anno = ANNO
    output: 
        dmr_outfiles
    conda:
        "../envs/dmr.yaml"
    shell:
        """
        comb-p pipeline \
		--seed {params.min_p} \
		--dist {params.win_sz} \
		-p {params.o_prefix} \
		--region-filter-p {params.region_filter} \
		--anno {params.anno} \
		{input.in_file} 
        """

rule annotate_dmr:
    input:
        ewas_bed = rules.run_dmr.output,
        dmr_bed = dmr_anno
    output:
        dmr_cpg_anno
    conda:
        "../envs/dmr.yaml"
    shell:
        """
        Rscript scripts/dmr_annotation.R {input.dmr_bed} {input.ewas_bed}
        """