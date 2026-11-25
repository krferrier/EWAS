# Annotation

## Annotation Resources

Everything the workflow downloads for annotation is configured in the single
`annotation:` block of `config.yml` and cached under one root
(`annotation.cache_dir`, default `resources/annotation/`). The cache is
subdivided by source because each source is versioned differently:

```
resources/annotation/
  zhou/<array_platform>/<zhou_release>/   Zhou Infinium manifest + provenance
  ucsc/<genome_build>/<cache_tag>/        UCSC tracks + provenance
  eqtm/                                   BIOS eQTM and HGNC tables
```

Each fetch rule writes a small `*_manifest.tsv`
recording the source URL, version and timestamp, so a results directory can be
traced back to the exact annotation used.

| Source | Supplies | Fetched by | Needed when |
|--------|----------|------------|-------------|
| Zhou Infinium manifest (`<platform>.<genome>.manifest.gencode.<release>.tsv.gz`) | `genesUniq`, `geneNames`, `transcriptTypes`, `transcriptIDs`, `distToTSS`, probe coordinates | `get_annotation_data` | always |
| UCSC `cpgIslandExt` | `CGI`, `CGIposition` (Island / N_Shore / S_Shore / N_Shelf / S_Shelf) | `fetch_cpg_island_cache` | always |
| BIOS *cis*-eQTM + HGNC | `BIOS_eQTM_genes` | `get_annotation_data`, `prep_bios_eqtm_annotation` | always |
| UCSC `refGene` + HGNC BigBed | DMR gene annotation | `fetch_dmr_annotation_cache` | only when `dmr_analysis: "yes"` |

### Config parameters to check

| Setting | Default | Notes |
|---------|---------|-------|
| `array_platform` | `EPIC` | **Must match the probe IDs in your M-value matrix.** EPIC, HM450, HM27 and Mammal40 use bare IDs (`cg00000029`); EPICv2 and MSA use design-suffixed IDs (`cg00000029_TC21`). All platforms share the same manifest columns, so a mismatch does not obviously fail. `annotation.R` therefore reports the match rate and aborts below 50%. |
| `zhou_release` | `v8.1` | Git tag in `zhou-lab/InfiniumAnnotationData`. Use `main` to track the newest release; pin a tag for a reproducible run. |
| `gencode_release` | `v41` | Tied to `genome_build`: hg38 &rarr; `v41`, hg19 &rarr; `v26lift37`, mm10 &rarr; `vM25`, mm39 &rarr; `vM31`. |
| `cache_tag` | `latest` | UCSC track cache version. Set an ISO date (e.g. `2026-06-08`) for a manuscript run so the cache path records when the tracks were pulled. |

The remaining keys are base URLs and rarely need changing. The whole block may
be omitted; `helper_fxns.ConfigWizard` supplies the same values as defaults.

### Notes on the annotation sources

Zhou's resources moved off `zhouserver.research.chop.edu` and were reorganised
around a versioned "coherent" release (currently v8.1). This has caused two notable
changes for this pipeline. First, GENCODE v41 no longer carries the `CGI` and `CGIposition` columns
that v36 supplied, so island context is recomputed from the UCSC track — which
is where Zhou's own CGI annotation was derived from. Re-deriving reproduces the
KYCG v8.1 CGI bitset for 99.90% of probes and the legacy island coordinate
strings for 99.99%, and the residual differences are probes the frozen legacy
table missed. Shores extend 2 kb from an island and shelves a further 2 kb.

Second, the old `EPIC.hg38.commonsnp.tsv.gz` SNP annotation has no successor in
v8.1 and is no longer used. Probes affected by common variants are
expected to be masked before the M-value matrix reaches this workflow.

### Network access

The fetch rules require outbound access to `github.com`, `api.github.com`,
`ngdc.cncb.ac.cn`,
`hgdownload.soe.ucsc.edu`, `molgenis26.gcc.rug.nl` and
`storage.googleapis.com`. All downloads are cached, so only the first run needs
them.

## Annotated Results Columns

`<out_directory>/<association_variable>_ewas_annotated_results<out_type>`, one row
per tested CpG, sorted by adjusted p-value.

The identifier and statistic columns differ by EWAS mode: an unstratified run
keys on `cpgid` with BACON columns, a stratified run keys on `MarkerName` with
METAL columns.

| Column | Source | Meaning |
|--------|--------|---------|
| `cpgid` *(unstratified)* | EWAS | Probe identifier. |
| `MarkerName` *(stratified)* | METAL | Probe identifier. |
| `bacon.es`, `bacon.se`, `bacon.statistic`, `bacon.pval` *(unstratified)* | BACON | Bias- and inflation-adjusted effect size, standard error, test statistic and p-value. |
| `lambda`, `b.lambda` *(unstratified)* | `run_bacon.R` | Genomic inflation before and after BACON adjustment: median observed chi-square over its null expectation, computed as `QCEWAS::P_lambda` does. |
| `Effect`, `StdErr`, `P-value`, `Direction` *(stratified)* | METAL | Inverse-variance weighted meta-analysis across strata. `Direction` gives one character per stratum. |
| `CpG_chrm`, `CpG_beg`, `CpG_end` | Zhou manifest | Probe coordinates on `genome_build`, 0-based half-open. |
| `probe_strand` | Zhou manifest | Strand the probe interrogates. |
| `genesUniq` | Zhou manifest | Unique gene symbols the probe maps to, semicolon separated. |
| `geneNames` | Zhou manifest | Gene symbol per overlapping transcript, aligned with `transcriptIDs`. |
| `transcriptTypes` | Zhou manifest | Biotype per transcript (`protein_coding`, `lncRNA`, ...). |
| `transcriptIDs` | Zhou manifest | Ensembl transcript identifiers. |
| `distToTSS` | Zhou manifest | Signed distance in bp to each transcript's TSS; negative is upstream. |
| `CGI` | UCSC `cpgIslandExt` | Coordinates of the nearest CpG island, e.g. `CGI:chr1:28735-29737`. Empty in open sea. |
| `CGIposition` | derived | `Island`, `N_Shore`, `S_Shore`, `N_Shelf`, `S_Shelf`, or empty for open sea. Shores span 2 kb from the island, shelves a further 2 kb. |
| `BIOS_eQTM_genes` | BIOS | Genes whose expression correlates with methylation at this CpG in whole blood (FDR 0.05), HGNC-resolved. |
| `chromHMM_state` | KYCG | Roadmap/ENCODE 18-state chromatin state: `TssA`, `TssFlnk`, `TssFlnkU`, `TssFlnkD`, `Tx`, `TxWk`, `EnhG1`, `EnhG2`, `EnhA1`, `EnhA2`, `EnhWk`, `ZNF/Rpts`, `Het`, `TssBiv`, `EnhBiv`, `ReprPC`, `ReprPCWk`, `Quies`. |
| `PMD` | KYCG | `commonPMD` (partially methylated domain) or `commonHMD` (highly methylated domain). |
| `AB_compartment` | KYCG | Hi-C compartment: `A1`, `A2` (open/active) or `B1`-`B4` (closed/inactive). |
| `repeat_class` | KYCG | RepeatMasker class: `LINE`, `SINE`, `LTR`, `Satellite`, `Simple_repeat` and others. Semicolon separated where a probe falls in more than one. |
| `imprinting_DMR` | KYCG | `ImprintingDMR` when the probe sits in a known imprinting control region. |
| `CTCF_binding` | KYCG | `CTCFbind` when the probe overlaps a CTCF binding site. |
| `ENCODE_blacklist` | KYCG | `Blacklist` when the probe falls in an ENCODE blacklist region. Treat such hits with suspicion. |

KYCG columns appear only when the configured `array_platform` publishes those
sets; a missing set is skipped rather than erroring. An empty value means the
probe is not annotated for that feature, not that the feature is absent.

The seven KYCG columns are described in
[What the KYCG knowledgebases are](enrichment.md#what-the-kycg-knowledgebases-are).
