# Running the workflow

## Modifying the Configuration File

This workflow uses a configuration file, `config.yml`, to specify the paths for input and output files, what kind of EWAS to perform (standard or stratified), parallelization parameters, whether to run DMR, and whether to perform the functional analyses.

### Input Files

This workflow is intended to be used with phenotype and methylation data that has already been cleaned and has no missing/`NA` data. The first column of the phenotype data and the methylation data must be the sample IDs. Examples of what the phenotype and methylation data should look like can be found in the `data/` directory. Accepted input file types:

* regular delimited (.csv, .tsv, .txt, etc.)
  * file path or URL
  * can be `.gz` or `.bz2` compressed
* fast-storage (.fst):
  * file path

#### Test data

`data/pheno.csv` and `data/mvals.csv.gz` are a small simulated test set, sized
to run the whole workflow quickly: 120 samples (60 F, 60 M) and 4,988 real EPIC CpGs across all chromosomes,
including X and Y. The phenotype columns are `sampleID`, `sex`, `re` and `BMI`.
The `config.yaml` is set to run the dummy dataset with no changes needed.

The dummy dataset is seeded with the following signals:

* 60 single CpGs with a BMI effect, two thirds of them at active promoters
  (ChromHMM `TssA`), which gives the KYCG enrichment a known signal to detect;
* 10 differentially methylated regions of 6-15 neighbouring probes each, plus 10 correlated regions with no effect as a negative
  control that should not be called;
* X-inactivation in females, so chrX behaves as it does in real data.
* a few CpGs on unplaced/alt contigs and a few with no hg38 position, as real
  EPIC data has.

`data/test_truth.tsv` lists all planted CpGs with its true effect so you can
check what you found. Effects are moderate (strongest hits around p = 1e-24), so the plots
look like a real study rather than a stress test.

The set is regenerated, byte-identically, by `data/make_test_data.py` (seeded;
needs numpy and pandas, plus `yame` and the KYCG cache for the promoter bias).
Its options change the sample size, number of CpGs and planted signal.

## Dry Run & DAG

Before running the full pipeline, you can perform a dry run to make sure that your configuration file is set up correctly.

```shell
snakemake --dry-run
```
It can also be helpful to export a Directed Acyclic Graph (DAG) of the workflow, especially in the case of a stratified EWAS so that you can check that the data is subset how you expect it to be. If not already installed, you will need Graphviz to create a png of the DAG (`conda install graphviz`).

```shell
snakemake --dag --rule plot_results | dot -Tpng > dag.png
```

Examples of DAGs for a standard and stratified EWAS can be found in
`data/example_dags`. That directory also holds rule-level graphs, which are
much easier to read than the job DAG on a stratified run because they show one
node per rule rather than one per job:

```shell
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

* `data/example_dags/standard_ewas_rulegraph.png`
* `data/example_dags/stratified_ewas_rulegraph.png`

## Run EWAS

Once you have modified the `config.yml` file to match your projects specifications and computational resources you are ready to run the EWAS. The first time you run the pipeline it may take a while to begin as Snakemake will need to download the annotation files and build the Conda environments. These steps will be cached and will not need to be performed again in subsequent runs.

**Standard (non-stratified) or Stratified EWAS**:

```shell
snakemake -j <n_jobs>
```

## Track Progress

The EWAS, DMR, functional-enrichment and plotting rules write what their
scripts print -- messages, warnings, errors and the EWAS progress bar -- to
a log file under `<out_directory>/logs/`, not to the terminal. Snakemake's
own log in `.snakemake/log/` records only which jobs ran; the scripts' output is
in these files. If one of these jobs fails, Snakemake's error report tells you which
log to read. The download, annotation-preparation and METAL-setup rules print
little and still write to the terminal.

To see the EWAS progress bar while it runs, open a second terminal in the
workflow directory and run:

```shell
tail -F <out_directory>/logs/run_combined_ewas.log          # standard EWAS
tail -F <out_directory>/logs/run_ewas_group/<stratum>.log   # stratified EWAS
```

Each log is rewritten when its job reruns, so it always shows output/progress from the current run.
