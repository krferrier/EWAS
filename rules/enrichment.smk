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
    # KYCG knowledgebases. Not everything the platform publishes: a short
    # biological selection, set by enrichment: knowledgebases: in the config
    # and defaulting to ConfigWizard.KYCG_ENRICHMENT_SETS, which explains the
    # reasoning for each inclusion and exclusion.
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
        # Everything below -- stdout and stderr of every command -- goes to the
        # rule's log. Snakemake names the log in its error report if a job fails.
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
    log:
        CW.log_path("enrich_pathways")
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        # Everything below -- stdout and stderr of every command -- goes to the
        # rule's log. Snakemake names the log in its error report if a job fails.
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
    # One dot plot per enrichment analysis: top N by p-value, point size and
    # transparency showing how many CpGs or genes drive each result, solid
    # versus hollow separating what passes FDR from what does not.
    #
    # The plot is always written, so it stays a tracked output. When the table
    # is empty or nothing is significant, the plot says that instead -- the
    # same convention the enrichment tables use. That avoids the alternative
    # of an untracked output, which would need a sentinel file for Snakemake
    # to know the rule had run.
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
        # Everything below -- stdout and stderr of every command -- goes to the
        # rule's log. Snakemake names the log in its error report if a job fails.
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
    log:
        CW.log_path("enrich_traits")
    conda:
        "../envs/enrichment.yaml"
    shell:
        """
        # Everything below -- stdout and stderr of every command -- goes to the
        # rule's log. Snakemake names the log in its error report if a job fails.
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
