# EWAS workflow

!!! info "Version 2.0"
    The previous stable release is available as tag [`v1.2.0`](https://github.com/krferrier/EWAS/tree/v1.2.0).

This repository is a snakemake workflow for performing Epigenome-Wide Association Studies (EWAS) of methylation measured by an Illumina Infinium methylation array. It was developed and tested against the EPIC (850k) array. The annotation step is platform-agnostic — every array Wanding Zhou has published annotation data for (EPIC, EPICv2, MSA, HM450, HM27, Mammal40, MM285) is available and can be specified via `annotation.array_platform`  in the `config.yaml`. Only EPIC has been executed and evaluated end to end, so treat the others as untested. This workflow can perform a standard EWAS or a stratified EWAS based on the variables provided for stratification. In either a standard or stratified EWAS, a linear regression model is used where the outcome is the methylation value (Beta or M-values) and the trait/phenotype you want to perform association testing with is the main predictor. All other variables included in the phenotype dataframe will be added as covariates to the model.

If your main predictor is categorical (with 2 levels), the values must be dummy coded to numerics. All other categorical covariates can remain as characters. The R glm() function is unable to use unordered categorical data as the main predictor because it cannot determine which of the values to use as the reference. The pipeline has not yet been tested with a main predictor that has more than two categorical levels. Automatic handling of predictors that are categorical/character values will be added in the next major update.

Results from the linear regression analyses will be adjusted for bias and inflation using a Bayesian approach implemented by the [BACON](https://www.bioconductor.org/packages/release/bioc/html/bacon.html) R package. BACON >= 1.32.0 (Bioconductor 3.19) is required and pinned in `envs/ewas.yaml`: earlier versions lack the `globalSeed` and `parallelSeed` arguments this workflow relies on for reproducible Gibbs sampling. The ggplot2 diagnostic plots in `scripts/updated_bacon/modified_bacon_plots.R` are local modifications to the package, so the plots are created with ggplot instead of base R. It is strongly recommended to assess the performance of the bias and inflation adjustment by viewing the traces, posteriors, fit, and qq plots output at this step. See [Interpreting BACON plots](bacon-plots.md) for how to assess them. If the EWAS was stratified, after adjustment for bias and inflation, the results from all the strata will be combined using an inverse-variance weighted meta-analysis approach with the command line tool [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation).

The final EWAS results are annotated from three sources: gene and promoter assignments from the [Wanding Zhou](https://zwdzwd.github.io/InfiniumAnnotation) Infinium manifest, CpG-island context derived from the UCSC `cpgIslandExt` track, and whole-blood *cis*-eQTM genes from BIOS. A manhattan and qq plot of the final EWAS results is also output. See [Annotation Resources](annotation.md#annotation-resources) for what each contributes and how it is versioned, and [Annotated Results Columns](annotation.md#annotated-results-columns) for the output file layout.

The significant CpG set is then assessed automatically for functional, pathway and trait enrichment, chromatin states and other [Know Your CpG (KYCG)](https://zhou-lab.github.io/kycg/) feature sets, GO and KEGG terms, and published EWAS Atlas trait associations. All tests are compared against the CpGs actually analyzed rather than the whole array. See [Functional Enrichment](enrichment.md).

This workflow can also perform differentially methylated region (DMR) analysis using the EWAS summary statistics. The DMR is conducted with the [comb-p](https://github.com/brentp/combined-pvalues) command-line tool. Parameters for performing the DMR can be modified in the `config.yaml`.

## Dependencies

* Conda/Mamba
* Snakemake, version >= 9.3.3

## Quick start

```shell
git clone https://github.com/krferrier/EWAS.git
cd EWAS
snakemake --dry-run      # check the configuration
snakemake -j 1           # run the bundled simulated test data
```

`config.yml` is set up to run the simulated test data in `data/` with no
changes. See [Running the workflow](usage.md) to point it at your own data.

## Citation

If you use this workflow, please cite [Ferrier *et al.* 2025, *HGG Advances*](https://doi.org/10.1016/j.xhgg.2025.100552). See [Citing](citing.md) for the full reference and BibTeX.

## License

This workflow is free for non-commercial use under the
[PolyForm Noncommercial License 1.0.0](license.md). Commercial use is not
permitted.
