# Overview

---
This repository is a snakemake workflow for performing Epigenome-Wide Association Studies (EWAS) of methylation measured by the Illumina EPIC (850k) methylation array. This workflow can perform a standard EWAS or a stratified EWAS based on the variables provided for stratification. In either the standard or stratified EWAS, a linear regression model is used where the outcome is the methylation value (Beta or M-values) and the trait/phenotype you want to perform association testing with is the main predictor. All other variables included in the phenotype dataframe will be added as covariates to the model.

Results from the linear regression analyses will be adjusted for bias and inflation using a Bayesian approach implemented by the [BACON](https://www.bioconductor.org/packages/release/bioc/html/bacon.html) R package. It is strongly recommended to assess the performance of the bias and inflation adjustment by viewing the traces, posteriors, fit, and qq plots output at this step. Further details on how to assess the performance plots are provided below. If the EWAS was stratified, after adjustment for bias and inflation, the results from all the strata will be combined using an inverse-variance weighted meta-analysis approach with the command line tool [METAL](https://genome.sph.umich.edu/wiki/METAL_Documentation).

The final results will be annotated with hg38/GRCh38 human genome build information collated by [Wanding Zhou](https://zwdzwd.github.io/InfiniumAnnotation). A manhattan and qq plot of the final results will also be output.

# Dependencies

---

* Conda/Mamba
* Snakemake
* Singularity (only for stratified EWAS)

# Using this Workflow

---

## Modifying the Configuration File

This workflow uses a configuration file, `config.yml`, to specify the paths for input and output files, what kind of EWAS to perform (standard or stratified), and parallelization parameters.

### Input Files

This workflow is intended to be used with phenotype and methylation data that has already been cleaned and has no missing/`NA` data. The first column of the phenotype data and the methylation data must be the sample IDs. Examples of what the phenotype and methylation data should look like can be found in the `data/` directory. Accepted input file types:

* regular delimited (.csv, .tsv, .txt, etc.)
  * file path or URL
  * can be `.gz` or `.bz2` compressed
* fast-storate (.fst):
  * file path

### Output Files

Output files will be generated in the 'out_directory' specified in the config file.

**Standard EWAS**:

* raw linear regression results (.csv or .csv.gz)
* BACON bias- and inflation-adjusted results (.csv or .csv.gz) and plots (.jpg)
* final annotated results (.csv or .csv.gz)
* manhattan and qq plots (.jpg)

**Stratified EWAS**:

* raw linear regression results (.csv or .csv.gz) for each stratum
* BACON bias- and inflation-adjusted results (.csv or .csv.gz) and plots (.jpg) for each stratum
* final meta-analyzed and annotated results (.csv or .csv.gz)
* manhattan and qq plots (.jpg)

### Parallelization Parameters

The rule 'ewas' first chunks the methylation dataset into sets of CpGs where the length of each set is specified by the parameter `chunks` in the config file (default of 1000 CpGs if a number is not provided). Then, linear regressions are run for each chunk and the results combined back into one dataframe once all chunks have been processed. Run sequentially, this step could take several hours. However, each chunk can be processed in parallel to reduce the total computation time.

The main step of the EWAS utilizes the [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) R package to perform the linear regressions in parallel. You can specify the type of parallelization to be performed at this step in the config file. The options include: sequential, multisession (threads), multicore, or cluster. You must also specify the number of workers that you want to be used for the parallelization in the `workers` parameter of the config file. For more details on the different types of parallelization you can reference the BiocParallel documentation.

If you are performing a stratified analysis, you can have each strata processed in parallel using the `-j` parameter in the snakemake command. For resource allocation, keep in mind that each EWAS performed will use the number of workers specified in the config file, so you will need n=(workers x jobs) resources available. For example, if you use multicore parallelization with 2 workers and `-j 2` in the snakemake command, you will need to have 4 cores available to run the analysis.

## Dry Run & DAG

Before running the full pipeline, you can perform a dry run to make sure that your configuration file is set up correctly. 

```shell
snakemake --dry-run
```
It can also be helpful to export a Directed Acyclic Graph (DAG) of the workflow, especially in the case of a stratified EWAS so that you can check that the data is subset how you expect it to be. If not already installed, you will need Graphviz to create a png of the DAG (`conda install graphviz`).

```shell
snakemake --dag --rule plot_results | dot -Tpng > dag.png
```

Examples of DAGs for a standard and stratified EWAS can be found in `data/example_dags`.

## Run EWAS

Once you have modified the `config.yml` file to match your projects specifications and computational resources you are ready to run the EWAS. The first time you run the pipeline it may take a while to begin as snakemake will need to download the annotation files and build the conda environment for performing the EWAS. These steps will be cached and will not need to be performed again in subsequent runs. 

**Standard (non-stratified) EWAS**:

```shell
snakemake -j 1 --use-conda
```

**Stratified EWAS**:

```shell
snakemake -j 1 --use-conda --use-singularity
```

## Interpretation of BACON Performance Plots

There are four plot types output from the bias- and inflation- adjustment analysis step to assess the performance of the Gibbs Sampler algorithm including:

* **Traces**: traces-plots of all estimates (sigma, p, mu)
* **Posteriors**: scatter plot of the Gibbs Sampler posterior probabilities for inflation (sigma.0) and proportion of null features (p.0). Elliptical curves correspond to 75%, 90%, and 95% probability regions.
* **Fit**: fit of the three component mixture over-layed on a histogram of z-scores. The black curve shows the overall fit, red shows the fit of the null distribution, and the blue and green curves show the alternatives.
* **QQ**: qq plots before and after adjustment with BACON.

Below are examples of what these plots should look like when the Gibbs Sampler algorithm has performed well. Poor performance can be an indication that there is something problematic with the underlying data; perhaps there are extreme outliers that were not removed prior to analysis, or maybe there is a confounding factor that isn't properly accounted for in the EWAS model. If the plots look like those shown below, then the EWAS results are likely to be reliable. If the plots look different, then further investigation is needed to determine the cause of the poor performance.

### Traces Plot

The traces plot shows the evolution of each parameter over time during the Gibbs Sampler algorithm. The algorithm begins with the provided priors and iterates until it converges on a set of estimates for each parameter. If the initial prior for a parameter is not close to the final estimate, there will be a period of divergence before convergence occurs. This can be seen in the traces plot as a sharp change in the value of the parameter at the beginning of the algorithm and then a gradual approach towards the final estimate. If the initial prior is close to the final estimate, then the algorithm will converge quickly and the plot will appear like a 'hairy caterpillar' (not a smooth line).

![Example of the output traces plot from performing bias- and inflation-adjustment with BACON.](/data/example_plots/traces_plot.jpg "Example Traces Plot")

In the above plot, you can see that the traces for converged on their respective estimates after a period of divergence for all parameters except mu.1. The initial prior estimate for mu.1 was close to the final estimate (the range of estimates is relatively small), so the algorithm converged quickly and resulted in a 'hairy caterpillar' plot.

### Posteriors Plot

The posteriors plot is a scatterplot of the estimated posterior values of the parameters mu.0 and p.0. The distribution of posterior values should be normally distributed for each parameter, thus the scatterplot should appear like a dense cloud of points. A sparse cloud indicates the algorithm was not able to converge on an estimate quickly (if at all).

![Example of the output posteriors plot from performing bias- and inflation-adjustment with BACON.](/data/example_plots/posteriors_plot.jpg "Example Posteriors Plot")

### Fit Plot

The fit plot shows the distribution of the z-scores for each observation in the dataset as a histogram. Overlayed is a black density line of the overall fit of the model estimated by the Gibb's Sampling algorithm based on three components: a null component (red density line) and two alternate components (blue and green density lines). In the case of an EWAS:

* **Null Component**: The background noise or null hypothesis that a CpG is not significantly associated with the phenotype.
* **Alternate Components**: The alternative hypotheses; the CpG is significantly associated with the disease in a negative or positive way.

If the estimated distribution (black line) is vastly different from the observed distribution, then the model may be mis-specified and/or the Gibb's Sampling algorithm did not converge on an estimate. If the estimated distribution is similar to the observed distribution, then the model is likely correctly specified.

![Example of the output fit plot from performing bias- and inflation-adjustment with BACON.](/data/example_plots/fit_plot.jpg "Example Fit Plot")

### QQ Plots

QQ plots are a scatterplot of the observed p-values against the expected p-values under the null hypothesis. The diagonal line represents the perfect fit between the observed and expected values. In the case of an EWAS, the points should follow the diagonal line relatively closely, with a spike of points at the end of the plot which represent the significant associations. Deviation from the diagonal line prior to the spike (if there is one) indicates there is global bias and/or inflation of the observed data. After adjustment  with BACON ('corrected' panel), the points prior to the spike should follow the diagonal line more closely. If there is still evidence of global bias and/or inflation after adjustment, then the model may be mis-specified.

![Example QQ plot from performing bias- and inflation-adjustment with BACON.](/data/example_plots/qqs_plot.jpg "Example QQ Plot")

Often, the BACON adjustment does not drastically change the observed p-values compared with the unadjusted p-values and it can be difficult to assess visually whether the BACON adjustment 'worked'. The effect of the adjustment can more easily be discerned with a metric often referred to as the 'genomic inflation factor', or lambda (λ), which is a calculation of the deviation of the observed p-values from the expected p-values. A value of 1 indicates no inflation, while values greater than 1 indicate inflation and less than 1 indicate deflation. In general, the BACON adjustment should bring the λ closer to 1. The unadjusted λ ('lambda') and BACON-adjusted ('b-lambda') values can be found in the output results file ending in 'ewas_bacon_results'.
