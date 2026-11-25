# Workflow overview

Optional stages are labelled in the config;
rule names are in parentheses.

```mermaid
flowchart TD
    IN["M-value matrix + phenotype table"]
    RES["Annotation resources: Zhou manifest, KYCG,<br/>UCSC tracks, BIOS eQTM, EWAS Atlas"]

    IN -- "stratified_ewas: no" --> REG["Linear regression per CpG<br/>(run_combined_ewas)"]
    IN -- "stratified_ewas: yes" --> ST["Split into strata<br/>(stratify_data)"]
    ST --> REGG["Linear regression per CpG, per stratum<br/>(run_ewas_group)"]

    REG --> BAC["Bias and inflation adjustment<br/>+ diagnostic plots<br/>(run_bacon)"]
    REGG --> BACG["Bias and inflation adjustment<br/>+ diagnostic plots, per stratum<br/>(run_bacon_group)"]
    BACG --> META["Inverse-variance weighted meta-analysis<br/>(make_metal_script, run_metal)"]

    BAC --> ANN["Annotate CpGs: genes, CpG islands,<br/>eQTM, KYCG features<br/>(add_annotation)"]
    META --> ANN
    RES --> ANN

    ANN --> PLT["Manhattan and QQ plots<br/>(plot_results)"]
    ANN -- "enrichment.run: yes" --> ENR["Feature, GO/KEGG and trait enrichment<br/>(enrich_features, enrich_pathways, enrich_traits)"]
    RES --> ENR
    ANN -- "dmr_analysis: yes" --> BED["Annotated results to BED<br/>(make_bed)"]
    BED --> DMR["Region calling with comb-p<br/>(run_dmr)"]
    DMR --> DANN["Annotate regions, manhattan plot<br/>(annotate_dmrs, plot_dmrs)"]
    RES --> DANN
```

The run command and configuration is recorded for each run that executes at least one job-- see [Run Provenance](outputs.md#run-provenance). Rule-level graphs of the complete workflow for both standard and stratified EWAS are shown below.

## Rule graphs

=== "Standard EWAS"

    ![Rule graph of a standard EWAS](assets/dags/standard_ewas_rulegraph.png)

=== "Stratified EWAS"

    ![Rule graph of a stratified EWAS](assets/dags/stratified_ewas_rulegraph.png)
