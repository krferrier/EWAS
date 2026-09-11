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
hyper_test <- function(n_overlap, n_query, n_mask, n_univ) {
    stats::phyper(n_overlap - 1L, n_mask, n_univ - n_mask, n_query,
                  lower.tail = FALSE)
}
