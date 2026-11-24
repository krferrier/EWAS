# Functional enrichment

Three independent tests of the significant CpG set, controlled by the
`enrichment:` block in `config.yml` and written to
`<out_directory>/enrichment/`. Set `run: "no"` to skip all three.

**The background is the set of CpGs actually tested, not the full array.**
Probes are dropped before an EWAS for mapping quality, masking and QC, and
those exclusions are not uniform across the genome. Testing against the whole
manifest would show associations based on what was excluded. The counts here will
therefore not match tools that default to the full array as the background set.

A test that cannot run writes an empty table with the reason it wasn't run recorded in its log
(otherwise the workflow would fail for a 'missing' file).

| Setting | Default | Notes |
|---------|---------|-------|
| `run` | `yes` | `no` skips all three rules. |
| `significance` | `fdr` | `fdr` (Benjamini-Hochberg), `bonferroni`, or `nominal` for the raw p-value. |
| `threshold` | `0.05` | Cutoff on the adjusted (or raw) p-value that defines the significant set. |
| `min_set_size` | `20` | Features, terms and traits with fewer probes in the background are not tested. |
| `ewas_atlas_url` | NGDC | EWAS Atlas association table, about 100 MB, cached on first use. |

## `<assoc>_enrichment_features.tsv`

The significant set against the KYCG knowledgebases from the same pinned v8.1
release the annotation uses: chromatin states, 82 histone marks, 1188
transcription-factor binding sets, repeats, PMDs, A/B compartments, metagene
position and CpG islands.

For what each knowledgebase is and what its feature labels mean, see
[What the KYCG knowledgebases are](enrichment.md#what-the-kycg-knowledgebases-are).

Implemented with `yame`: the query is packed as a format-6 record with two
bits per probe, one for the background and one for the test set, so the restricted
background is applied inside the overlap counting. `yame summary` then reports
the 2x2 per feature. The hypergeometric p-value and FDR are computed in R.

Columns: `knowledgebase`, `role`, `feature`, `n_universe`, `n_significant`,
`n_in_feature`, `n_overlap`, `expected`, `fold_enrichment`, `neg_log10_p`,
`neg_log10_fdr`, `n_tested_in_kb`, `log2_odds_ratio`, `p_value`, `fdr`.

`role` is `biological` or `qc`, which is how the two feature plots split this
table. `n_tested_in_kb` is the size of the testing family the row's FDR was
computed in. `neg_log10_p` and `neg_log10_fdr` are exact where `p_value` and
`fdr` have underflowed to zero.

The background here is **CpGs**, *not genes* -- see
[KYCG tests CpGs; GO and KEGG test genes](enrichment.md#kycg-tests-cpgs-go-and-kegg-test-genes).

## `<assoc>_enrichment_pathways.tsv`

GO and KEGG via `missMethyl::gometh`, which corrects for the number of probes
per gene. This matters because a gene with 80 probes is far more likely to have a significant CpG than a gene covered by 3, and an uncorrected gene-set test
reflects that difference in coverage and is misconstrued as biological significance.

Of note, missMethyl maps probes with Illumina's own annotation
packages, so **only 450K, EPIC and EPICv2 are supported**. Any other
`array_platform` will output an empty table for GO and KEGG enrichment. Feature and trait
enrichment still run.

Columns: `collection`, `term_id`, `term`, `ontology`, `n_genes_in_term`,
`n_significant_genes`, `p_value`, `fdr`.

## `<assoc>_enrichment_traits.tsv`

Over-representation of published trait associations from the
[EWAS Atlas](https://ngdc.cncb.ac.cn/ewas/atlas).

Columns: `trait`, `n_universe`, `n_significant`, `n_trait_probes`,
`n_overlap`, `expected`, `fold_enrichment`, `p_value`, `fdr`, `n_studies`,
`pmids`, `overlapping_probes`.

The EWAS Atlas is a catalogue of what has been published, so its
coverage reflects study volume. So, commonly studied traits/phenotypes, like smoking, ageing and sex are overrepresented compared to rare traits/phenotypes, and an
overlap partly reflects how often they have been measured. As a reminder, each
trait's probe list is restricted here to probes in your background, so counts
will not match the Atlas website if you were to use their online webtool, which reports across all arrays at once.

The Atlas is keyed on `cg` identifiers, so an EPICv2 or MSA run that keeps
the design suffix will not match; the rule reports this rather than silently returning
empty results.

## What the KYCG knowledgebases are

The names of the knowledgebases referenced in the Know Your CpG (KYCG) database
are not always intuitive or clear as to what that knowledgebase represents for
biological interpretation. And, the feature labels within a knowledgebase annotation
can be even more unclear unless you happen to already be familiar with the knowledgebase.
For example the `ABCompartment` is an enrichment knowledgebase of the Hi-C subcompartments
as defined by Rao et al. 2014. The `ABCompartment` set is comprised of six subcompartment labels,
which generally tell you how open/active or closed/inactive the chromatin is.

### Which sets are tested

Zhou publishes 17 sets for EPIC and 32 for MSA. This workflow tests two
named groups, both written to the same table (differentiated by the `role`
column), but plotted separately:

`enrichment.knowledgebases` -- the biological sets (13 by default):

| Set | Description |
|---|---|
| `ChromHMM`, `REMCChromHMM` | chromatin state (ENCODE and Roadmap models) |
| `HM` | histone mark or variant |
| `TFBSrm` | transcription factor binding sites |
| `CTCFbind` | CTCF sites, whose binding is methylation-sensitive |
| `CGI` | islands, shores, shelves or open sea |
| `MetagenePC` | particular position relative to genes |
| `ABCompartment` | Hi-C A/B subcompartment |
| `PMD` | partially methylated domains |
| `rmsk1`, `rmsk2` | repeat class or family |
| `Tetranuc2` | share a WCGW/SCGS sequence context |
| `ImprintingDMR` | allele-specific biology |

There are four knowledgebases which are more related to the array design:
`ProbeType`, `InfiniumChemistry`, `Blacklist` and `nFlankCG`. Zhou's registry frames the first three as post-hoc *controls*.
You *don't* want to see your results enriched for these features. `nFlankCG` tells
you if your results are enriched for areas with high or low density of SNPs. Typically,
you will see hits enriched in CpG Islands and/or promotor regions where there are a higher
density of CpGs, so this is more of a positive-control check.

On EPIC those 17 happen to be everything the platform publishes. Other
platforms publish more. MSA has 32, including tissue signatures, CoRSIVs,
evolutionary conservation and G-quadruplex peaks. Any published set can
be added by name. If you want to see all options, you can find them in:

* `data/kycg_set_coverage.tsv` in this repository: all 32 sets across the six
  array platforms, with which platforms publish each, its role, upstream
  source and citation.
* the cached registry, `<cache>/zhou/<platform>/<release>/KYCG/knowledgebases.tsv`,
  which is Zhou's own provenance for every set -- see above. Upstream removed
  this file from `zhou-lab/kycg` in September 2026 (the definitions are now
  compiled into the `kycg` tool, where `kycg info` shows them), so the workflow
  fetches it from the last commit that still has it, `d6df6f36`.

The FDR is computed within each knowledgebase. Some knowledge bases have only one or
two features while others have thousands, so the FDR for each knowledgebase is reflective of
the multiple testing burden for each set.

Further reading: the KnowYourCG framework paper is
[Goldberg et al. 2025, *Science Advances*](https://doi.org/10.1126/sciadv.adw3027),
and the sets are browsable at
[zwdzwd.github.io/InfiniumAnnotation](https://zwdzwd.github.io/InfiniumAnnotation).

### KYCG tests CpGs; GO and KEGG test genes

The KYCG knowledgebases are based on **CpGs**: the background is the probes
actually tested and each feature is a set of probes. A significant result says *your
top CpGs fall in this annotation more often than the full set of CpGs you tested*.
GO and KEGG are based on **genes**: `enrich_pathways.R` maps significant CpGs to
genes first, and a hit says *the genes your CpGs map to are over-represented in
this term*.

Two points to keep in mind:

* A CpG-level result does not rely on a gene mapping. CpGs occur in both coding and
  non-coding regions of the genome. The KYCG annotation information is particularly useful for providing
  biological information on CpGs that are non-coding or intergenic.
* A gene-level result uses the CpG-to-gene mapping. Not all CpGs *have* gene mappings. And the CpG-to-gene mapping assumes that the CpG affects that gene's expression (cis-effect) when it possible a CpG affects a more distant (by a linear measurement of basepairs) gene (trans-effect).

## Enrichment plots

With `enrichment.make_plots: "yes"` four plots are produced.

| Plot | Description |
|---|---|
| `<assoc>_enrichment_features.jpg` | KYCG enrichment blocks |
| `<assoc>_enrichment_features_qc.jpg` | KYCG enrichment blocks, QC sets only |
| `<assoc>_enrichment_pathways.jpg` | GO/KEGG dot plot, one panel per collection |
| `<assoc>_enrichment_traits.jpg` | Trait Enrichment dot plot |

**KYCG plots** follow the shape of knowYourCG's
`KYCG_plotEnrichAll`: all tested knowledgebase are along the x
axis, labelled underneath with the number of features tested, and its
enriched features directly above it. The y axis is -log10(FDR), point size
is the log2 odds ratio, and the strongest hit in each knowledgebase is
labelled. Related knowledgebases are adjacent, grouped by what they describe,
with alternating background bands marking the groups

**The pathway and trait plots** dot plots of the top
`enrichment.plot_top_n` by FDR (that setting does not affect the feature
plots, which draw every enriched feature). Point size and transparency show how many
CpGs or genes drive each result; solid points pass the FDR threshold and
hollow ones do not.
