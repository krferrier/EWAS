from helper_fxns import ConfigWizard

configfile: "config.yml"
CW = ConfigWizard(config)


#----SET VARIABLES----#
PHENO = str(CW.pheno)
MVALS = str(CW.mvals)
ASSOC = CW.assoc_var
STRATIFIED = CW.stratified
STRAT_VARS = CW.strat_vars
DMR = CW.dmr
CHUNK_SIZE = CW.chunk_size
PROCESSING_TYPE = CW.processing_type
N_WORKERS = CW.n_workers
OUT_DIR = str(CW.out_dir)
OUT_TYPE = CW.out_type

# Plots / groups
PLOTS = CW.bacon_plot_kinds
GROUPS = CW.groups  # e.g., ["female"] or ["female","male","..."] or ["all"]

# Unstratified outputs
raw_results      = str(CW.raw_results)
bacon_results    = str(CW.bacon_results)
bacon_plots      = CW.bacon_plot_files()

# Stratified outputs
strat_raw_results    = CW.strat_raw_results()
strat_bacon_results  = CW.strat_bacon_results()
strat_bacon_plots    = CW.strat_bacon_plot_files()

# Final Outputs
annotated_results = str(CW.annotated_results)
manhattan_qq_plot = str(CW.manhattan_qq_plot)
meta_analysis_results = str(CW.meta_analysis_results)

# DMR outputs
results_bed = str(CW.dmr_results_bed)
dmr_acf     = str(CW.dmr_acf)
dmr_args    = str(CW.dmr_args)
dmr_fdr     = str(CW.dmr_fdr)
dmr_regions = str(CW.dmr_regions)
dmr_slk = str(CW.dmr_slk)
dmr_anno_final = str(CW.dmr_anno_final)
#dmr_anno = str(CW.dmr_anno)
#dmr_cpg_anno = [str(CW.dmr_cpg_anno)]

# Local DMR annotation-cache resources
dmr_refgene_txt = str(CW.dmr_refgene_txt)
dmr_cpg_island_txt = str(CW.dmr_cpg_island_txt)
dmr_hgnc_bb = str(CW.dmr_hgnc_bb)
#dmr_hgnc_complete_set = str(CW.dmr_hgnc_complete_set)
dmr_refgene_bed = str(CW.dmr_refgene_bed)
dmr_cpg_island_bed = str(CW.dmr_cpg_island_bed)
#dmr_hgnc_raw_bed = str(CW.dmr_hgnc_raw_bed)
dmr_hgnc_bed = str(CW.dmr_hgnc_bed)
dmr_annotation_manifest = str(CW.dmr_annotation_manifest)

dmr_infile = [results_bed]
dmr_outfiles = [dmr_acf, dmr_args, dmr_fdr, dmr_regions, dmr_slk]
dmr_anno_files = [dmr_anno_final]
#---- DETERMINE INPUT FILES FOR RULE ALL ----#
if STRATIFIED:
    in_files = [PHENO, MVALS, strat_raw_results, strat_bacon_results,
                strat_bacon_plots, meta_analysis_results,
                annotated_results, manhattan_qq_plot]
else:
    in_files = [PHENO, MVALS, raw_results, bacon_results, 
                bacon_plots, annotated_results, 
                manhattan_qq_plot]

if DMR:
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