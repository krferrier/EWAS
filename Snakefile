from helper_fxns import ConfigWizard
from snakemake.utils import validate

configfile: "config.yml"
CW = ConfigWizard(config)

validate(config, "config.schema.yml")

#---- DETERMINE INPUT FILES FOR RULE ALL ----#
if CW.stratified:
    in_files = [CW.pheno, CW.mvals, CW.strat_raw_results(), CW.strat_bacon_results(),
                CW.strat_bacon_plot_files(), CW.meta_analysis_results,
                CW.annotated_results, CW.manhattan_qq_plot]
else:
    in_files = [CW.pheno, CW.mvals, CW.raw_results, CW.bacon_results,
                CW.bacon_plot_files(), CW.annotated_results,
                CW.manhattan_qq_plot]

dmr_infile = [CW.dmr_results_bed]
dmr_outfiles = [CW.dmr_acf, CW.dmr_args, CW.dmr_fdr, CW.dmr_regions, CW.dmr_regions_p, CW.dmr_slk]
dmr_anno_files = [CW.dmr_anno_final, CW.dmr_manhattan_plot]

if CW.dmr:
    in_files = in_files + dmr_infile + dmr_outfiles + dmr_anno_files
else:
    in_files = in_files

#---- BEGIN WORKFLOW ----#
rule all:
    input:
        in_files


include: "rules/combined_ewas.smk"
include: "rules/stratified_ewas.smk"
include: "rules/annotate.smk"
include: "rules/plots.smk"
include: "rules/dmr.smk"