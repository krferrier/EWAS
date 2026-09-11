# GO and KEGG enrichment for the significant CpG set, via missMethyl.
#
# gometh/gsameth are used rather than a plain gene-set test on the mapped genes
# because array probes are not distributed evenly across genes: a gene covered
# by 80 probes is far more likely to pick up a significant CpG than one covered
# by 3, and an uncorrected test reports that coverage as biology. missMethyl
# models the probe-per-gene bias explicitly.
#
# The cost is that missMethyl maps probes with Illumina's own annotation
# packages, so it supports only 450K, EPIC and EPICv2. On any other platform
# this script writes an empty table naming the reason rather than failing.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
       "enrichment_common.R"))

parser <- argparse::ArgumentParser(
    description = "GO/KEGG enrichment for the significant CpG set")
parser$add_argument('--input-file', required = TRUE)
parser$add_argument('--platform', required = TRUE,
                    help = "annotation.array_platform from config.yml")
parser$add_argument('--stratified', required = TRUE)
parser$add_argument('--significance', default = "fdr")
parser$add_argument('--threshold', type = "double", default = 0.05)
parser$add_argument('--min-set-size', type = "integer", default = 20L)
parser$add_argument('--output', required = TRUE)
args <- parser$parse_args()

OUT_COLS <- c("collection", "term_id", "term", "ontology", "n_genes_in_term",
              "n_significant_genes", "p_value", "fdr")

# missMethyl's array codes and the annotation package each needs. gometh
# resolves the annotation through the search list rather than the namespace,
# so the package must be attached with library() -- having it installed is not
# enough, and the failure surfaces as an opaque "no item called package:..."
# error from updateObject().
ARRAY_TYPES <- list(
    HM450  = list(code = "450K",
                  anno = "IlluminaHumanMethylation450kanno.ilmn12.hg19"),
    EPIC   = list(code = "EPIC",
                  anno = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"),
    EPICv2 = list(code = "EPICv2",
                  anno = "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
)

sets <- read_cpg_sets(args$input_file, args$stratified,
                      args$significance, args$threshold)

if (!args$platform %in% names(ARRAY_TYPES)) {
    write_empty_result(args$output, OUT_COLS, sprintf(
        paste("missMethyl supports only 450K, EPIC and EPICv2;",
              "annotation.array_platform is '%s'. Feature and trait enrichment",
              "still run -- only GO/KEGG is skipped."), args$platform))
    quit(save = "no", status = 0)
}
if (length(sets$foreground) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       "no CpGs passed the significance threshold")
    quit(save = "no", status = 0)
}
if (!requireNamespace("missMethyl", quietly = TRUE)) {
    write_empty_result(args$output, OUT_COLS, "missMethyl is not installed")
    quit(save = "no", status = 0)
}

spec <- ARRAY_TYPES[[args$platform]]
array_type <- spec$code
attached <- suppressWarnings(suppressPackageStartupMessages(
    require(spec$anno, character.only = TRUE, quietly = TRUE)))
if (!attached) {
    write_empty_result(args$output, OUT_COLS, sprintf(
        "annotation package %s could not be attached; add it to envs/enrichment.yaml",
        spec$anno))
    quit(save = "no", status = 0)
}
message(sprintf("missMethyl array type: %s (annotation: %s)", array_type, spec$anno))

run_collection <- function(collection) {
    tryCatch({
        res <- missMethyl::gometh(
            sig.cpg  = sets$foreground,
            all.cpg  = sets$background,
            collection = collection,
            array.type = array_type,
            plot.bias = FALSE)
        res <- as.data.table(res, keep.rownames = "term_id")
        res[, collection := collection]
        res
    }, error = function(e) {
        warning(sprintf("%s enrichment failed: %s", collection, conditionMessage(e)))
        NULL
    })
}

collected <- Filter(Negate(is.null), lapply(c("GO", "KEGG"), run_collection))
if (length(collected) == 0L) {
    write_empty_result(args$output, OUT_COLS, "gometh returned no results")
    quit(save = "no", status = 0)
}

all <- rbindlist(collected, fill = TRUE)
# gometh names differ a little between collections
if (!"ONTOLOGY" %in% names(all)) all[, ONTOLOGY := NA_character_]
if (!"TERM" %in% names(all) && "Description" %in% names(all)) {
    setnames(all, "Description", "TERM")
}
setnames(all,
         c("TERM", "ONTOLOGY", "N", "DE", "P.DE", "FDR"),
         c("term", "ontology", "n_genes_in_term", "n_significant_genes",
           "p_value", "fdr"),
         skip_absent = TRUE)

all <- all[n_genes_in_term >= args$min_set_size]
if (nrow(all) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       sprintf("no term reached min_set_size=%d", args$min_set_size))
    quit(save = "no", status = 0)
}

setorder(all, p_value)
missing <- setdiff(OUT_COLS, names(all))
for (nm in missing) all[[nm]] <- NA
fwrite(all[, ..OUT_COLS], args$output, sep = "\t")
message(sprintf("Tested %d terms; %d at FDR < 0.05",
                nrow(all), sum(all$fdr < 0.05, na.rm = TRUE)))
