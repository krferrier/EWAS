localrules: fetch_ewas_atlas

rule fetch_ewas_atlas:
    output:
        atlas = protected(CW.ewas_atlas_txt),
        manifest = CW.anno_cache_dir.joinpath("ewas_atlas", "ewas_atlas_manifest.tsv")
    params:
        url = CW.ewas_atlas_url
    shell:
        """
        wget -O {output.atlas}.tmp "{params.url}"
        mv {output.atlas}.tmp {output.atlas}

        {{
          echo -e "resource\tsource\trows\tcreated"
          echo -e "ewas_atlas_associations\t{params.url}\t$(wc -l < {output.atlas})\t$(date -Iseconds)"
        }} > {output.manifest}
        """

rule enrich_features:
    input:
        results = CW.annotated_results,
        kycg_dir = rules.fetch_kycg_features.output.kycg,
        ordering = rules.fetch_kycg_features.output.ordering,
        script = "scripts/enrich_features.R",
        common = "scripts/enrichment_common.R"
    output:
        CW.enrichment_feature_results
    params:
        strat = CW.stratified,
        significance = CW.enrich_significance,
        threshold = CW.enrich_threshold,
        min_set = CW.enrich_min_set_size,
        sets = ",".join(CW.kycg_tested_sets),
        qc_sets = ",".join(CW.kycg_qc_sets) or "NONE",
    log:
        CW.log_path("enrich_features")
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        exec >{log} 2>&1
        Rscript {input.script} \
        --input-file {input.results} \
        --kycg-dir {input.kycg_dir} \
        --sets {params.sets} \
        --qc-sets {params.qc_sets} \
        --probe-order {input.ordering} \
        --stratified {params.strat} \
        --significance {params.significance} \
        --threshold {params.threshold} \
        --min-set-size {params.min_set} \
        --output {output}
        """

rule enrich_pathways:
    input:
        results = CW.annotated_results,
        script = "scripts/enrich_pathways.R",
        common = "scripts/enrichment_common.R"
    output:
        CW.enrichment_pathway_results
    params:
        platform = CW.array_platform,
        strat = CW.stratified,
        significance = CW.enrich_significance,
        threshold = CW.enrich_threshold,
        min_set = CW.enrich_min_set_size
    log:
        CW.log_path("enrich_pathways")
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        exec >{log} 2>&1
        Rscript {input.script} \
        --input-file {input.results} \
        --platform {params.platform} \
        --stratified {params.strat} \
        --significance {params.significance} \
        --threshold {params.threshold} \
        --min-set-size {params.min_set} \
        --output {output}
        """

def _enrichment_table(wildcards):
    return CW.enrichment_table(wildcards.kind)

rule plot_enrichment:
    input:
        table = _enrichment_table,
        script = "scripts/plot_enrichment.R"
    output:
        CW.enrichment_plot("{kind}")
    params:
        assoc = CW.assoc_var,
        top_n = CW.enrich_plot_top_n,
        threshold = CW.enrich_threshold
    wildcard_constraints:
        kind = "|".join(CW.ENRICHMENT_KINDS)
    log:
        CW.log_path("plot_enrichment", "{kind}")
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        exec >{log} 2>&1
        Rscript {input.script} \
        --input-file {input.table} \
        --kind {wildcards.kind} \
        --assoc {params.assoc} \
        --top-n {params.top_n} \
        --fdr-threshold {params.threshold} \
        --output {output}
        """

rule enrich_traits:
    input:
        results = CW.annotated_results,
        atlas = rules.fetch_ewas_atlas.output.atlas,
        script = "scripts/enrich_traits.R",
        common = "scripts/enrichment_common.R"
    output:
        CW.enrichment_trait_results
    params:
        strat = CW.stratified,
        significance = CW.enrich_significance,
        threshold = CW.enrich_threshold,
        min_set = CW.enrich_min_set_size
    log:
        CW.log_path("enrich_traits")
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        exec >{log} 2>&1
        Rscript {input.script} \
        --input-file {input.results} \
        --ewas-atlas {input.atlas} \
        --stratified {params.strat} \
        --significance {params.significance} \
        --threshold {params.threshold} \
        --min-set-size {params.min_set} \
        --output {output}
        """
