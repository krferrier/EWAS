# Feature enrichment of the significant CpG set against the v8.1 KYCG
# knowledgebases (chromatin states, histone marks, transcription-factor binding,
# repeats, PMDs, A/B compartments, metagene position and the rest).
#
# The sets are bit-packed and row-aligned to <platform>.ordering.tsv.gz, so the
# query is built as a yame format-6 record: two bits per probe, one marking the
# universe (tested) and one marking the set (significant). yame summary then
# reports the 2x2 per knowledgebase record against that universe, and the
# hypergeometric p-value and FDR are computed here.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
       "enrichment_common.R"))

parser <- argparse::ArgumentParser(
    description = "KYCG feature enrichment for the significant CpG set")
parser$add_argument('--input-file', required = TRUE,
                    help = "Annotated EWAS results")
parser$add_argument('--kycg-dir', required = TRUE,
                    help = "Directory of normalised KYCG .cm feature sets")
parser$add_argument('--probe-order', required = TRUE,
                    help = "Path to <platform>.ordering.tsv.gz")
parser$add_argument('--stratified', required = TRUE)
parser$add_argument('--significance', default = "fdr")
parser$add_argument('--threshold', type = "double", default = 0.05)
parser$add_argument('--min-set-size', type = "integer", default = 20L)
parser$add_argument('--output', required = TRUE)
args <- parser$parse_args()

OUT_COLS <- c("knowledgebase", "feature", "n_universe", "n_significant",
              "n_in_feature", "n_overlap", "expected", "fold_enrichment",
              "log2_odds_ratio", "p_value", "fdr")

sets <- read_cpg_sets(args$input_file, args$stratified,
                      args$significance, args$threshold)

if (length(sets$foreground) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       "no CpGs passed the significance threshold")
    quit(save = "no", status = 0)
}
if (Sys.which("yame") == "") {
    write_empty_result(args$output, OUT_COLS, "yame not found on PATH")
    quit(save = "no", status = 0)
}

cm_files <- list.files(args$kycg_dir, pattern = "\\.cm$", full.names = TRUE)
if (length(cm_files) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       sprintf("no .cm knowledgebases in %s", args$kycg_dir))
    quit(save = "no", status = 0)
}

# --- build the format-6 query over the platform's row order ---
probe_order <- fread(args$probe_order, select = 1L)
setnames(probe_order, 1L, "cpgid")
in_bg <- probe_order$cpgid %in% sets$background
in_fg <- probe_order$cpgid %in% sets$foreground
if (!any(in_bg)) {
    stop(paste("no tested CpG matched the probe ordering file --",
               "check that annotation.array_platform matches the results"))
}
message(sprintf("Query: %d significant within %d tested, over %d array probes",
                sum(in_fg), sum(in_bg), nrow(probe_order)))

tmp_txt <- tempfile(fileext = ".txt")
tmp_cx  <- tempfile(fileext = ".cx")
on.exit(unlink(c(tmp_txt, tmp_cx)), add = TRUE)
fwrite(data.table(S = as.integer(in_fg), U = as.integer(in_bg)),
       tmp_txt, sep = "\t", col.names = FALSE)
system2("yame", c("pack", "-f", "d", shQuote(tmp_txt), shQuote(tmp_cx)))

# --- summarize against every knowledgebase ---
collected <- list()
for (cm in cm_files) {
    kb <- sub("\\.cm$", "", basename(cm))
    res <- tryCatch({
        txt <- system2("yame", c("summary", "-m", shQuote(cm), shQuote(tmp_cx)),
                       stdout = TRUE, stderr = FALSE)
        if (length(txt) < 2L) NULL else fread(text = paste(txt, collapse = "\n"))
    }, error = function(e) {
        warning(sprintf("knowledgebase '%s' skipped: %s", kb, conditionMessage(e)))
        NULL
    })
    if (is.null(res) || nrow(res) == 0L) next
    res[, knowledgebase := kb]
    collected[[kb]] <- res
    message(sprintf("  %-22s %d features", kb, nrow(res)))
}

if (length(collected) == 0L) {
    write_empty_result(args$output, OUT_COLS, "yame summary returned no rows")
    quit(save = "no", status = 0)
}

all <- rbindlist(collected, fill = TRUE)
setnames(all, c("Mask", "N_univ", "N_query", "N_mask", "N_overlap"),
         c("feature", "n_universe", "n_significant", "n_in_feature", "n_overlap"))

# Drop features too small to test, and the NA/unassigned bucket the state
# knowledgebases carry, which is an absence of annotation rather than a feature.
# fread parses a literal "NA" mask name as missing, so both forms are dropped.
all <- all[n_in_feature >= args$min_set_size &
           !is.na(feature) & !feature %in% c("NA", "", ".")]
if (nrow(all) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       sprintf("no feature reached min_set_size=%d", args$min_set_size))
    quit(save = "no", status = 0)
}

all[, expected := n_significant * (n_in_feature / n_universe)]
all[, fold_enrichment := fifelse(expected > 0, n_overlap / expected, NA_real_)]
all[, log2_odds_ratio := as.numeric(Log2OddsRatio)]
all[, p_value := hyper_test(n_overlap, n_significant, n_in_feature, n_universe)]
all[, fdr := p.adjust(p_value, method = "BH")]

setorder(all, p_value, -fold_enrichment)
fwrite(all[, ..OUT_COLS], args$output, sep = "\t")
message(sprintf("Tested %d features across %d knowledgebases; %d at FDR < 0.05",
                nrow(all), length(collected), sum(all$fdr < 0.05, na.rm = TRUE)))
