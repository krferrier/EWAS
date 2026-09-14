# Shared helpers for the enrichment rules.
#
# One rule about the background set is enforced here rather than left to each
# script: the universe is always the CpGs that were actually tested, never the
# full array manifest. Probes are dropped before an EWAS for mapping quality,
# masking and QC, and those exclusions are not uniform across the genome --
# testing against the whole array would score that removal pattern as
# enrichment.

suppressPackageStartupMessages({
    library(data.table)
})

#' Column names in the annotated results, which differ by EWAS mode.
result_columns <- function(stratified) {
    if (stratified == "no" || stratified == "False") {
        list(id = "cpgid", p = "bacon.pval")
    } else {
        list(id = "MarkerName", p = "P-value")
    }
}

#' Split the annotated results into a significant set and the tested background.
#'
#' @param significance "fdr" (Benjamini-Hochberg), "bonferroni", or "nominal".
#' @return list(foreground, background, adjusted) of probe IDs plus the
#'   adjusted p-values, or NULL when the results carry no usable p-value.
read_cpg_sets <- function(results_file, stratified, significance, threshold) {
    cols <- result_columns(stratified)
    res <- fread(results_file)
    for (nm in c(cols$id, cols$p)) {
        if (!nm %in% names(res)) {
            stop(sprintf("annotated results have no '%s' column (stratified=%s)",
                         nm, stratified))
        }
    }

    ids <- as.character(res[[cols$id]])
    pvals <- suppressWarnings(as.numeric(res[[cols$p]]))
    keep <- !is.na(ids) & !is.na(pvals)
    ids <- ids[keep]; pvals <- pvals[keep]
    if (length(ids) == 0L) stop("no CpGs with a usable p-value in the results")

    adjusted <- switch(
        significance,
        fdr        = p.adjust(pvals, method = "BH"),
        bonferroni = p.adjust(pvals, method = "bonferroni"),
        nominal    = pvals,
        stop(sprintf("unknown significance method '%s'", significance))
    )

    fg <- ids[adjusted < threshold]
    message(sprintf(
        "Significant CpGs: %d of %d tested (%s < %g)",
        length(fg), length(ids), significance, threshold))

    list(foreground = fg, background = ids, adjusted = adjusted, pvals = pvals,
         ids = ids, method = significance, threshold = threshold)
}

#' Write an empty results table with the expected header.
#'
#' A skipped analysis still produces its output file so the workflow does not
#' fail and the reason is visible in the file itself.
write_empty_result <- function(path, columns, reason) {
    empty <- setNames(data.table(matrix(character(0), ncol = length(columns))),
                      columns)
    fwrite(empty, path, sep = "\t")
    message(sprintf("Wrote empty %s: %s", basename(path), reason))
}

#' Hypergeometric over-representation test from 2x2 counts.
#'
#' Tests over-representation only (upper tail): a set depleted in the
#' significant CpGs is not evidence of biological depletion here, because the
#' foreground is defined by a p-value cutoff rather than by sampling.
#' Returns a list of `p_value` and `neg_log10_p`.
#'
#' The test is evaluated with log.p = TRUE and exponentiated afterwards, rather
#' than read off the linear scale directly. This is not a micro-optimisation: a
#' strong enrichment over a large feature genuinely produces p-values below the
#' smallest positive double (~1e-308), where the linear result collapses to
#' exactly 0 and the ranking among the top hits is destroyed. For example, 1400
#' of 2000 significant CpGs falling in a feature covering 6507 of 50000 tested
#' CpGs gives p = 1e-802: representable as a log, not as a double.
#'
#' `p_value` is therefore allowed to underflow to 0 in the output table -- that
#' is what a p-value column can express -- while `neg_log10_p` carries the exact
#' magnitude and is what the plots should use.
hyper_test <- function(n_overlap, n_query, n_mask, n_univ) {
    log_p <- stats::phyper(n_overlap - 1L, n_mask, n_univ - n_mask, n_query,
                           lower.tail = FALSE, log.p = TRUE)
    list(log_p = log_p, p_value = exp(log_p), neg_log10_p = -log_p / log(10))
}

#' Benjamini-Hochberg FDR on the log scale.
#'
#' `p.adjust(method = "BH")` cannot be used here: its input is linear p-values,
#' so any p that has underflowed to 0 gives an FDR of exactly 0, and -log10(0)
#' is Inf. Since the plots put -log10(FDR) on the x axis (the knowYourCG
#' convention), the FDR needs the same treatment the p-value already gets.
#'
#' BH is order statistics and one multiplication, both of which are exact in
#' logs: adjusted_(i) = min over j >= i of (n / j) * p_(j), which becomes
#' log n - log j + log p_(j), followed by a running minimum from the largest
#' p-value down and a clamp at log(1) = 0.
#'
#' Returns a list of `log_fdr`, `fdr` (exponentiated, may underflow to 0) and
#' `neg_log10_fdr` (exact). Verified against p.adjust() on values that do not
#' underflow.
bh_log <- function(log_p) {
    n <- length(log_p)
    if (n == 0L) {
        return(list(log_fdr = numeric(0), fdr = numeric(0),
                    neg_log10_fdr = numeric(0)))
    }
    ord <- order(log_p)
    scaled <- log(n) - log(seq_len(n)) + log_p[ord]
    # Running minimum from the largest p-value downwards, then enforce fdr <= 1.
    stepped <- pmin(rev(cummin(rev(scaled))), 0)
    log_fdr <- numeric(n)
    log_fdr[ord] <- stepped
    list(log_fdr = log_fdr, fdr = exp(log_fdr),
         neg_log10_fdr = -log_fdr / log(10))
}
