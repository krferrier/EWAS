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
parser$add_argument('--fdr-by-knowledgebase', default = "yes",
                    help = paste("Adjust p-values within each knowledgebase",
                                 "(yes, the default and what",
                                 "knowYourCG::testEnrichment does) rather than",
                                 "pooled across all of them (no)."))
parser$add_argument('--qc-sets', default = "",
                    help = paste("Comma-separated subset of --sets that are",
                                 "design/QC knowledgebases rather than",
                                 "biological ones. They are tested and",
                                 "reported in the same table, tagged",
                                 "role = 'qc', and plotted separately:",
                                 "enrichment in one of these says the hit list",
                                 "tracks array design rather than biology."))
parser$add_argument('--sets', default = "all",
                    help = paste("Comma-separated KYCG set names to test, or",
                                 "'all' for every .cm in --kycg-dir. Named",
                                 "explicitly rather than globbed so the tested",
                                 "set is recorded in the run, and does not",
                                 "change when the cache gains a file."))
parser$add_argument('--output', required = TRUE)
args <- parser$parse_args()

OUT_COLS <- c("knowledgebase", "role", "feature", "n_universe", "n_significant",
              "n_in_feature", "n_overlap", "expected", "fold_enrichment",
              "neg_log10_p", "neg_log10_fdr", "n_tested_in_kb",
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

# Restrict to the requested sets. A requested set that is absent is reported
# and skipped rather than fatal: a platform need not publish every set.
if (!identical(tolower(trimws(args$sets)), "all")) {
    want <- trimws(strsplit(args$sets, ",", fixed = TRUE)[[1]])
    want <- want[nzchar(want)]
    have <- sub("\\.cm$", "", basename(cm_files))
    missing <- setdiff(want, have)
    if (length(missing)) {
        message(sprintf("not published for this platform, skipping: %s",
                        paste(missing, collapse = ", ")))
    }
    cm_files <- cm_files[have %in% want]
    if (length(cm_files) == 0L) {
        write_empty_result(args$output, OUT_COLS,
                           sprintf("none of the requested knowledgebases (%s) are in %s",
                                   args$sets, args$kycg_dir))
        quit(save = "no", status = 0)
    }
    message(sprintf("testing %d knowledgebase(s): %s", length(cm_files),
                    paste(sub("\\.cm$", "", basename(cm_files)), collapse = ", ")))
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

# Tag each row biological or qc. The role does not change the testing family --
# that is the knowledgebase -- it only decides which figure a row appears in.
# "NONE" is the sentinel the rule passes when no QC sets are configured; an
# empty string cannot be handed through the shell command cleanly.
.qc <- if (identical(toupper(trimws(args$qc_sets)), "NONE")) character(0) else
    trimws(strsplit(args$qc_sets, ",", fixed = TRUE)[[1]])
.qc <- .qc[nzchar(.qc)]
all[, role := fifelse(knowledgebase %in% .qc, "qc", "biological")]
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
.ht <- hyper_test(n_overlap = all$n_overlap, n_query = all$n_significant,
                  n_mask = all$n_in_feature, n_univ = all$n_universe)
all[, log_p := .ht$log_p]
all[, p_value := .ht$p_value]
# Exact magnitudes even where p_value/fdr have underflowed to 0; the plots put
# -log10(FDR) on the x axis. See hyper_test() and bh_log().
all[, neg_log10_p := .ht$neg_log10_p]
# Each knowledgebase is a separate enrichment analysis, so the FDR is computed
# within one by default rather than pooled across all of them: pooling makes
# the FDR for a chromatin state depend on how many TF motifs were tested
# alongside it, and lets one large set (TFBSrm is 1188 motifs) spend the power
# the others need. knowYourCG::testEnrichment does the same -- its
# mtc_by_group defaults to TRUE and splits on the knowledgebase group.
#
# The trade is that q-values from families of very different sizes are no
# longer on one ranking: a motif needs stronger evidence than a chromatin
# state to reach the same FDR. n_tested_in_kb records each family's size so
# that is visible in the table.
if (identical(tolower(trimws(args$fdr_by_knowledgebase)), "yes")) {
    all[, n_tested_in_kb := .N, by = knowledgebase]
    all[, c("fdr", "neg_log10_fdr") := {
            b <- bh_log(log_p)
            list(b$fdr, b$neg_log10_fdr)
        }, by = knowledgebase]
} else {
    all[, n_tested_in_kb := .N]
    .bh <- bh_log(all$log_p)
    all[, fdr := .bh$fdr]
    all[, neg_log10_fdr := .bh$neg_log10_fdr]
}

setorder(all, p_value, -fold_enrichment)
fwrite(all[, ..OUT_COLS], args$output, sep = "\t")
message(sprintf("Tested %d features across %d knowledgebases; %d at FDR < 0.05",
                nrow(all), length(collected), sum(all$fdr < 0.05, na.rm = TRUE)))
