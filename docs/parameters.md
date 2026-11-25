# Parameters

All settings live in `config.yml` and are validated against
`config.schema.yml` when the workflow starts. The settings for annotation and
enrichment are described on the [Annotation](annotation.md#config-parameters-to-check)
and [Functional enrichment](enrichment.md) pages.

## Parallelization Parameters

The EWAS rules (`run_combined_ewas`, or `run_ewas_group` for each stratum) first split the methylation dataset into chunks of CpGs, sized by the required `chunk_size` setting in the config file (5000 in the provided `config.yml`; larger chunks use more memory per worker). Then, linear regressions are run for each chunk and the results combined back into one dataframe once all chunks have been processed. Run sequentially, this step could take several hours. However, each chunk can be processed in parallel to reduce the total computation time.

The main step of the EWAS utilizes the [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) R package to perform the linear regressions in parallel. You can specify the type of parallelization to be performed at this step in the config file. The options include: sequential, multisession (threads), multicore, or cluster. You must also specify the number of workers that you want to be used for the parallelization in the `workers` parameter of the config file. For more details on the different types of parallelization you can reference the BiocParallel documentation.

If you are performing a stratified analysis, you can have each strata processed in parallel using the `-j` parameter in the snakemake command. For resource allocation, keep in mind that each EWAS performed will use the number of workers specified in the config file, so you will need n=(workers x jobs) resources available. For example, if you use multicore parallelization with 2 workers and `-j 2` in the snakemake command, you will need to have 4 cores available to run the analysis.

## DMR Parameters

If you do not want DMR results, you can simply change the `dmr_analysis` parameter to "no".

Note that `genome_build` is not a DMR-only setting: it selects both the UCSC tracks and the Zhou manifest, so it applies whether or not DMR analysis is run.

| Parameter | Default value | Description                                                          |
|-----------|---------------|----------------------------------------------------------------------|
|min_pvalue | 1e-04         | P-value threshold for beginning a region                             |
| window_size | 200           | Maximum distance to search for another CpG < min_pvalue              |
| region_filter | 0.05          | max adjusted region-level p-value to be reported                     |

### DMR plots

`plot_dmrs` always exports a region manhattan plot. The `dmr_plots:` block
controls two optional additions.

| Parameter | Default | Description |
|-----------|---------|-------------|
| make_zoom | "no" | Write a locus-zoom plot plus a refGene gene table per cluster of significant regions |
| make_combined | "no" | Combine the manhattan and zoom plots into one multi-panel figure; needs `make_zoom: "yes"` |
| combined_formats | "pdf" | Format for the combined figure |
| min_probes | 2 | Minimum CpGs in a region. Passed to comb-p as `--region-filter-n`, so regions with fewer CpGs are dropped from the region-level results (`regions-p`), and also used as the threshold for highlighting and zoom plots |
| max_y | -1 | Manhattan y-axis maximum; -1 lets the data set it |
| zoom_padding | 2000 | bp of context on each side of a zoomed cluster |
| cluster_gap | 3000 | Largest gap in bp between regions still grouped into one zoom window |

`make_zoom` is off by default because it is unknown before you run how many significant regions might be detected, which could mean a large number of files could be created that you might not need or want. For each region, two files are generate: a plot of the CpGs in the regional window with a gene track underneath and a tsv listing the genes in that window.  Turn it on when
you want to look at individual regions.
