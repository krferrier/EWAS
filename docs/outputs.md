# Outputs

## Directory layout

Everything is written under the `out_directory` set in the config file. In the
paths below, `<assoc>` is `association_variable` and `<stratum>` is one level of
`stratify_variables` (so a stratified run repeats those files once per stratum).

```
<out_directory>/
|
|-- provenance/                                  how this run was invoked
|   |-- runs.tsv                                 append-only index: run_id, times, status, command
|   `-- <run_id>/
|       |-- command.txt                          the exact command line
|       |-- config_snapshot.yml                  verbatim copy of the config file used
|       |-- config_resolved.yml                  config after --config overrides
|       `-- run_info.yml                         versions, git state, host, timings, status
|
|-- logs/                                        what the analysis scripts printed (stdout + stderr)
|   |-- run_combined_ewas.log, run_bacon.log     STANDARD ONLY: EWAS and BACON
|   |-- run_ewas_group/<stratum>.log             STRATIFIED ONLY: EWAS, one per stratum
|   |-- run_bacon_group/<stratum>.log            STRATIFIED ONLY: BACON, one per stratum
|   |-- run_metal.log                            STRATIFIED ONLY: METAL meta-analysis
|   |-- plot_results.log                         Manhattan and QQ plots
|   |-- make_bed.log, run_dmr.log,               DMR: BED export, comb-p, annotation, plots
|   |   annotate_dmrs.log, plot_dmrs.log
|   |-- enrich_features.log, enrich_pathways.log, enrich_traits.log
|   `-- plot_enrichment/<kind>.log               one per enrichment plot
|
|-- <assoc>_ewas_results.csv.gz                  STANDARD ONLY: raw per-CpG regression
|-- <assoc>_ewas_bacon_results.csv.gz            STANDARD ONLY: bias/inflation adjusted
|-- bacon_plots/                                 STANDARD ONLY: BACON diagnostics
|   |-- <assoc>_traces.jpg
|   |-- <assoc>_posteriors.jpg
|   |-- <assoc>_fit.jpg
|   `-- <assoc>_qqs.jpg
|
|-- <stratum>/                                   STRATIFIED ONLY: one directory per stratum
|   |-- <stratum>_pheno.fst                      temporary, removed when the stratum finishes
|   |-- <stratum>_mvals.fst                      temporary, removed when the stratum finishes
|   |-- <stratum>_<assoc>_ewas_results.csv.gz    raw per-CpG regression for this stratum
|   |-- <stratum>_<assoc>_ewas_bacon_results.csv.gz
|   `-- bacon_plots/
|       |-- <stratum>_<assoc>_traces.jpg
|       |-- <stratum>_<assoc>_posteriors.jpg
|       |-- <stratum>_<assoc>_fit.jpg
|       `-- <stratum>_<assoc>_qqs.jpg
|-- meta_analysis/                               STRATIFIED ONLY: METAL meta-analysis
|   |-- <assoc>_metal_commands.sh                generated METAL script
|   |-- <assoc>_ewas_meta_analysis_results_1.txt       METAL output
|   `-- <assoc>_ewas_meta_analysis_results_1.txt.info  METAL's run summary
|
|-- <assoc>_ewas_annotated_results.csv.gz        final results, annotated (both modes)
|-- <assoc>_ewas_manhattan_qq_plots.jpg          manhattan + QQ of the final results
|
|-- enrichment/                                  only when enrichment.run: "yes"
|   |-- <assoc>_enrichment_features.tsv          KYCG feature over-representation
|   |-- <assoc>_enrichment_pathways.tsv          GO and KEGG terms
|   |-- <assoc>_enrichment_traits.tsv            EWAS Atlas trait associations
|   |-- <assoc>_enrichment_features.jpg          top hits per analysis;
|   |-- <assoc>_enrichment_features_qc.jpg        the design/QC sets, plotted
|   |                                              apart from the biology
|   |-- <assoc>_enrichment_pathways.jpg            all four only when
|   `-- <assoc>_enrichment_traits.jpg              enrichment.make_plots: "yes"
|
|-- <assoc>_ewas_annotated_results.bed           only when dmr_analysis: "yes": comb-p input
`-- dmr/                                         only when dmr_analysis: "yes"
    |-- <assoc>_ewas.args.txt                    comb-p arguments used
    |-- <assoc>_ewas.acf.txt                     autocorrelation by distance lag
    |-- <assoc>_ewas.slk.bed.gz                  per-CpG p-values after Stouffer-Liptak-Kechris
    |-- <assoc>_ewas.fdr.bed.gz                  the above with Benjamini-Hochberg FDR
    |-- <assoc>_ewas.regions.bed.gz              called regions
    |-- <assoc>_ewas.regions-p.bed.gz            called regions with region-level p-values
    |-- <assoc>_dmr_annotated_results.tsv        regions annotated with genes and CpG islands
    |-- <assoc>_dmr_manhattan.jpg                manhattan plot of the annotated regions
    |-- <assoc>_dmr_zoom_cluster_<i>.jpg         only when dmr_plots.make_zoom: "yes";
    |-- <assoc>_dmr_zoom_cluster_<i>.refGene_genes.tsv   one pair per region cluster,
    |                                              count unknown in advance, untracked
    `-- <assoc>_dmr_combined.pdf                 only when dmr_plots.make_combined: "yes"
```

A note on the DMR files: `run_dmr` declares six of comb-p's outputs as rule
outputs. comb-p writes some additional intermediates (for example
`<assoc>_ewas.regions-t.bed.gz`) that Snakemake does not track because if no region
reaches significance they are not created and the workflow breaks thinking there are
missing files.

The zoom plots are untracked for a similar reason. If and how
many appear depends on how many clusters of significant regions comb-p finds,
which is not knowable before it runs. They are side outputs of `plot_dmrs`,
whose manhattan plot *is* tracked, so the rule still reruns when the region
calls change. Because Snakemake will not clean them up, the `plot_dmrs` rule deletes any
zoom files from a previous run before writing new ones.

You can also choose to have a pdf with all the different DMR plots (manhattan and regional)
combined into one output pdf.

Annotation resources are cached outside `out_directory`, under
`annotation.cache_dir` (default `resources/annotation/`) so they are downloaded
once and shared across runs. METAL is built once into `software/metal/`. See
[Annotation Resources](annotation.md#annotation-resources).

### Descriptions of outputs

**EWAS** (`run_combined_ewas` / `stratify_data` + `run_ewas_group`)

* Raw per-CpG linear regression results: effect size, standard error, test
  statistic and p-value for the association variable, one row per CpG.
* A stratified run writes one of these per stratum, into `<stratum>/`. The
  per-stratum `.fst` phenotype and M-value subsets are temporary files, so
  Snakemake deletes them once the stratum's results exist.

**BACON** (`run_bacon` / `run_bacon_group`)

* Bias- and inflation-adjusted results, adding `bacon.es`, `bacon.se`,
  `bacon.statistic`, `bacon.pval` and the genomic inflation factors `lambda`
  (before) and `b.lambda` (after).
* Four diagnostic plots per run or per stratum -- traces, posteriors, fit and
  QQ. Check these before trusting the final results; see
  [Interpretation of BACON Performance Plots](bacon-plots.md).

**Meta-analysis** (`make_metal_script`, `run_metal`, stratified only)

* A generated METAL command script, kept so the meta-analysis is reproducible
  and inspectable.
* METAL's inverse-variance weighted results across strata: `MarkerName`,
  `Effect`, `StdErr`, `P-value` and `Direction`, the last giving one character
  per stratum so you can see whether strata agree in direction of effect.

**Annotation** (`add_annotation`)

* The final results table: the EWAS or meta-analysis statistics with annotation
  columns, sorted by p-value.
* Gene and transcript assignments, CpG-island context, whole-blood eQTM genes
  and KYCG feature columns. Column descriptions are listed in
  [Annotated Results Columns](annotation.md#annotated-results-columns).
* The rule reports the fraction of result CpGs it could annotate and stops if
  that falls below a threshold. This catches mismatches between
  `annotation.array_platform` and the array the M-values came from.

**Enrichment** (`enrich_features`, `enrich_pathways`, `enrich_traits`)

* Three tables testing the significant CpG set for enrichment: KYCG
  features, GO/KEGG terms, and published EWAS Atlas traits. All enrichment tests are against the CpGs
  actually tested rather than the whole array. Details and caveats in
  [Functional Enrichment](enrichment.md).

**DMR** (`make_bed`, `run_dmr`, `annotate_dmrs`, `plot_dmrs`)

* A BED format of the annotated results, which is comb-p's input.
* comb-p's intermediates: the arguments used, the autocorrelation function by
  distance lag, per-CpG p-values corrected for local correlation, and the same
  with FDR. The `.acf.txt` shows the distance over which
  neighbouring CpGs are correlated, which sets how regions get built.
* The called regions, with and without region-level p-values. comb-p also draws
  its own manhattan plot, but `run_dmr` deletes it; the workflow's DMR manhattan
  plot comes from `plot_dmrs` (below).
* The final annotated region table, with genes and CpG-island context attached,
  and a manhattan plot of those regions.

**Provenance** (`onstart` / `onsuccess` / `onerror` handlers)

* A per-run record of the command and configuration at runtime. See [Run Provenance](outputs.md#run-provenance).

## Run provenance

Every run that executes at least one job writes a record of how it was run to
`<out_directory>/provenance/`. It is written by the Snakefile's `onstart`,
`onsuccess` and `onerror` handlers rather than by a rule, so it is never
skipped as "already up to date".

Each run gets its own directory, `provenance/<run_id>/`, where `<run_id>` is
the start time (`YYYYMMDDTHHMMSS`):

| File | Contents |
|------|----------|
| `command.txt` | The exact command line, with each argument also on its own line |
| `config_snapshot.yml` | A verbatim copy of every configuration file used |
| `config_resolved.yml` | The merged configuration Snakemake ran with, including any `--config` overrides |
| `run_info.yml` | Snakemake and Python versions, the workflow's git state (commit, branch, tag, uncommitted changes), host, user, working directory, start and end time, and final status |

`provenance/runs.tsv` is an append-only index with one row per run: run ID,
times, status and command.

A dry run, or a run that reports "Nothing to be done", writes no record. No
results were produced, and the record of the run that did produce them is
already on disk. A failure to write provenance only prints a warning; it never
stops the analysis.
