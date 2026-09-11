localrules: fetch_ewas_atlas

# Functional assessment of the significant CpG set. Three independent tests,
# each writing its own table; a test that cannot run writes an empty table
# naming the reason rather than failing the workflow.
#
# All three share one background: the CpGs actually tested, never the full
# array manifest. See scripts/enrichment_common.R for why.

rule fetch_ewas_atlas:
    # EWAS Atlas association table (CpG to published trait association).
    # Roughly 100 MB; cached alongside the other annotation resources.
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
    # KYCG knowledgebases: chromatin states, histone marks, TF binding,
    # repeats, PMDs, A/B compartments, metagene position.
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
        min_set = CW.enrich_min_set_size
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.results} \
        --kycg-dir {input.kycg_dir} \
        --probe-order {input.ordering} \
        --stratified {params.strat} \
        --significance {params.significance} \
        --threshold {params.threshold} \
        --min-set-size {params.min_set} \
        --output {output}
        """

rule enrich_pathways:
    # GO and KEGG via missMethyl, which corrects for probes-per-gene bias.
    # Skips with an explanatory empty table on platforms it cannot map.
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
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.results} \
        --platform {params.platform} \
        --stratified {params.strat} \
        --significance {params.significance} \
        --threshold {params.threshold} \
        --min-set-size {params.min_set} \
        --output {output}
        """

rule enrich_traits:
    # EWAS Atlas trait over-representation.
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
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.results} \
        --ewas-atlas {input.atlas} \
        --stratified {params.strat} \
        --significance {params.significance} \
        --threshold {params.threshold} \
        --min-set-size {params.min_set} \
        --output {output}
        """
