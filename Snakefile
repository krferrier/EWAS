import sys

from helper_fxns import (ConfigWizard, record_run_provenance,
                         finalize_run_provenance)
from snakemake.utils import validate

configfile: "config.yml"
CW = ConfigWizard(config)

validate(config, "config.schema.yml")

dmr_targets = [CW.dmr_anno_final, CW.dmr_manhattan_plot]
enrichment_targets = [CW.enrichment_feature_results,
                      CW.enrichment_pathway_results,
                      CW.enrichment_trait_results]
enrichment_targets.extend(CW.enrichment_plot_files())

#---- DETERMINE INPUT FILES FOR RULE ALL ----#
if CW.stratified:
    in_files = [CW.pheno, CW.mvals, CW.strat_raw_results(), CW.strat_bacon_results(),
                CW.strat_bacon_plot_files(), CW.meta_analysis_results,
                CW.annotated_results, CW.manhattan_qq_plot]
else:
    in_files = [CW.pheno, CW.mvals, CW.raw_results, CW.bacon_results,
                CW.bacon_plot_files(), CW.annotated_results,
                CW.manhattan_qq_plot]

if CW.dmr:
    in_files.extend(dmr_targets)

if CW.enrichment:
    in_files.extend(enrichment_targets)

#---- BEGIN WORKFLOW ----#
rule all:
    input:
        in_files


include: "rules/combined_ewas.smk"
include: "rules/stratified_ewas.smk"
include: "rules/annotate.smk"
include: "rules/plots.smk"
include: "rules/dmr.smk"
# after annotate.smk: references rules.fetch_kycg_features
include: "rules/enrichment.smk"


#---- RUN PROVENANCE ----#
# Saves the command line and the configuration file into
# <out_directory>/provenance/<run_id>/ so a results directory always carries a
# record of how it was produced. See helper_fxns.record_run_provenance for why
# this is a handler rather than a rule.
_run_provenance = None

onstart:
    global _run_provenance
    _run_provenance = record_run_provenance(CW, workflow, config, sys.argv)

onsuccess:
    finalize_run_provenance(_run_provenance, "success")

onerror:
    finalize_run_provenance(_run_provenance, "error")