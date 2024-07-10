# Overview
This repository is a snakemake workflow for performing Epigenome-Wide Association Studies (EWAS) of methylation measured by the Illumina EPIC (850k) methylation array. This workflow can perform a standard EWAS or a stratified EWAS based on the variables provided for stratification. In either the standard or stratified EWAS, a linear regression model is used where the outcome is the methylation value (Beta or M-values) and the trait/phenotype you want to perform association testing with is the main predictor. All other variables included in the phenotype dataframe will be added as covariates to the model.

Results from the linear regression analyses will be adjusted for bias and inflation using a Bayesian approach implemented by the [BACON](https://www.bioconductor.org/packages/release/bioc/html/bacon.html) R package. It is strongly recommended to assess the performance of the bias and inflation adjustment by viewing the traces, posteriors, fit, and qq plots output at this step. Further details on how to assess the performance plots are provided below. If the EWAS was stratified, after adjustment for bias and inflation, the results from all the strata will be combined using an inverse-variance weighted meta-analysis approach with the command line tool [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation).

The final results will be annotated with hg38/GRCh38 human genome build information collated by [Wanding Zhou](https://zwdzwd.github.io/InfiniumAnnotation). A manhattan and qq plot of the final results will also be output.

# Dependencies
* Conda/Mamba
* Snakemake
* Singularity (only for stratified EWAS)

# Using this Workflow

## Modifying the Configuration File
This workflow uses a configuration file, `config.yml`, to specify the paths for input and output files, what kind of EWAS to perform (standard or stratified), and parallelization parameters.
### Input Files
This workflow is intended to be used with phenotype and methylation data that has already been cleaned and has no missing/NA data. The first column of the phenotype data and the methylation data must be the sample IDs. Examples of what the phenotype and methylation data should look like can be found in the `data/` directory. Accepted input filetypes:
- regular delimited (.csv, .tsv, .txt, etc.)
    - file path or URL
    - can be `.gz` or `.bz2` compressed
- fast-storate (.fst):
    - file path
### Output Files
Output files will be generated in the 'out_directory' specified in the config file.
Standard EWAS:
- raw linear regression results (.csv or .csv.gz)
- BACON bias- and inflation-adjusted results (.csv or .csv.gz) and plots (.jpg)
- final annotated results (.csv or .csv.gz)
- manhattan and qq plots (.jpg)
Stratified EWAS:
- raw linear regression results (.csv or .csv.gz) for each stratum
- BACON bias- and inflation-adjusted results (.csv or .csv.gz) and plots (.jpg) for each stratum
- final meta-analyzed and annotated results (.csv or .csv.gz)
- manhattan and qq plots (.jpg)
### Parallelization Parameters
The rule 'ewas' first chunks the methylation dataset into sets of CpGs where the length of each set is specified by the parameter `chunks` in the config file (default of 1000 CpGs if a number is not provided). Then, linear regressions are run for each chunk and the results combined back into one dataframe once all chunks have been processed. Run sequentially, this step could take several hours. However, each chunk can be processed in parallel to reduce the total computation time. The main step of the EWAS utilizes the [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) R package to perform the linear regressions in parallel. You can specify the type of parallelization to be performed at this step in the config file. The options include: sequential, multisession (threads), multicore, or cluster. You must also specify the number of workers that you want to be used for the parallelization in the 'workers' parameter of the config file. For more details on the different types of parallelization you can reference the BiocParallel documentation.

If you are performing a stratified analysis, you can have each strata processed in parallel using the `-j` parameter in the snakemake command. For resource allocation, keep in mind that each EWAS performed will use the number of workers specified in the config file, so you will need n=(workers x jobs) resources available. For example, if you use multicore parallelization with 2 workers and `-j 2` in the snakemake command, you will need to have 4 cores available to run the analysis.

## Run EWAS

Once you have modified the `config.yml` file to match your projects specifications and computational resources you are ready to run the EWAS.

For a standard (non-stratified) EWAS:
```shell
snakemake -j 1 --use-conda
```

For a stratified EWAS:
```shell
snakemake -j 1 --use-conda --use-singularity
```
## Interpretation of BACON Performance Plots
There are four plot types output from the bias- and inflation- adjustment analysis step to assess the performance of the Gibbs Sampler algorithm including:
    - Traces: traces-plots of all estimates (sigma, p, mu)
    - Posteriors: scatter plot of the Gibbs Sampler posterior probabilities for inflation (sigma.0) and proportion of null features (p.0). Elliptical curves correspond to 75%, 90%, and 95% probability regions.
    - Fit: fit of the three component mixture over-layed on a histogram of z-scores. The black curve shows the overall fit, red shows the fit of the null distribution, and the blue and green curves show the alternatives.
    - QQ: qq plots before and after adjustment with BACON.

Below are examples of what these plots should look like when the Gibbs Sampler algorithm has performed well. Poor performance can be an indication that there is something problematic with the underlying data; perhaps there are extreme outliers that were not removed prior to analysis, or maybe there is a confounding factor that isn't properly accounted for in the EWAS model. If the plots look like those shown below, then the EWAS results are likely to be reliable. If the plots look different, then further investigation is needed to determine the cause of the poor performance.

### Traces Plot

