#!/usr/bin/env Rscript
#
# One summary plot per enrichment analysis: the top N features, GO/KEGG terms
# or EWAS Atlas traits as a dot plot, with -log10(p) on the x axis, point size
# and transparency carrying how many CpGs or genes drive the result, and solid
# versus hollow points separating what passes FDR from what does not.
#
# The plot is ALWAYS written, including when the enrichment table is empty or
# nothing reached significance -- it then carries a panel saying so. This keeps
# the output trackable by Snakemake, and matches how the enrichment rules write
# an empty table naming the reason rather than failing the workflow.
#
# Non-significant terms are shown rather than dropped, because "the top hit was
# only p = 0.2" is a useful thing to see. Hollow points mark them, so a plot
# with no solid points is a plot with no findings.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(dplyr)
    library(ggplot2)
})

parser <- argparse::ArgumentParser(
    description = "Plot the top enrichment results for one analysis")
parser$add_argument('--input-file', required = TRUE,
                    help = "Enrichment results table from one of the enrich_* rules")
parser$add_argument('--kind', required = TRUE,
                    choices = c("features", "features_by_knowledgebase",
                                "pathways", "traits"),
                    help = "Which enrichment analysis this table came from")
parser$add_argument('--assoc', required = TRUE,
                    help = "Association variable, used in the plot title")
parser$add_argument('--top-n', type = "integer", default = 10L,
                    help = "Number of results to show [default %(default)s]")
parser$add_argument('--fdr-threshold', type = "double", default = 0.05,
                    help = "FDR below which a result counts as significant [default %(default)s]")
parser$add_argument('--label-width', type = "integer", default = 55L,
                    help = "Wrap y-axis labels at this many characters [default %(default)s]")
parser$add_argument('--output', required = TRUE)
args <- parser$parse_args()

# How each table names its label, grouping and count columns, and what one
# count actually is. Keep in step with the OUT_COLS in the enrich_* scripts.
#
# `facet` matters statistically, not just visually. enrich_features.R adjusts
# one pooled BH family across all knowledgebases, so its results form a single
# coherent ranking and belong in one panel. enrich_pathways.R instead calls
# gometh once per collection and keeps each collection's own FDR, so GO and
# KEGG are SEPARATE testing families -- and GO carries some 22,000 terms
# against KEGG's ~350, so a pooled top N is GO-dominated and hides KEGG
# entirely. Those get one panel each, with top N taken within each.
SPEC <- list(
    features = list(
        label = "feature", group = "knowledgebase", facet = NULL,
        count = "n_overlap",
        unit = "CpGs", axis = "Feature",
        title = "KYCG feature enrichment", colour = "#1B9E77"),
    pathways = list(
        label = "term", group = NULL, facet = "collection",
        count = "n_significant_genes",
        unit = "Genes", axis = "Term",
        title = "GO and KEGG enrichment", colour = "#D95F02"),
    traits = list(
        label = "trait", group = NULL, facet = NULL, count = "n_overlap",
        unit = "CpGs", axis = "Trait",
        title = "EWAS Atlas trait enrichment", colour = "#7570B3"),
    # One row per knowledgebase instead of a pooled top N, and x is fold
    # enrichment rather than -log10(p). Effect size is the comparable quantity
    # across knowledgebases of wildly different sizes: a large feature earns a
    # tiny p-value at the same fold enrichment as a small one, so ranking these
    # by p would just re-sort them by feature size.
    features_by_knowledgebase = list(
        label = "feature", group = NULL, facet = NULL, count = "n_overlap",
        unit = "CpGs", axis = "Knowledgebase (strongest feature)",
        title = "KYCG best hit per knowledgebase", colour = "#1B9E77",
        per_group = "knowledgebase", force_effect = TRUE)
)
spec <- SPEC[[args$kind]]

#' Wrap long labels without pulling in stringr; base strwrap is equivalent and
#' the other plotting scripts in this workflow do not depend on it either.
wrap_labels <- function(x, width) {
    vapply(as.character(x), function(s) {
        if (is.na(s) || !nzchar(s)) return("(unnamed)")
        paste(strwrap(s, width = width), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}

#' Placeholder used when there is nothing to plot, so the output always exists.
write_placeholder <- function(path, title, reason) {
    p <- ggplot() +
        annotate("text", x = 0, y = 0.15, label = title,
                 size = 5.5, fontface = "bold") +
        annotate("text", x = 0, y = -0.15, label = paste(strwrap(reason, 60),
                                                         collapse = "\n"),
                 size = 4, colour = "grey30") +
        scale_x_continuous(limits = c(-1, 1)) +
        scale_y_continuous(limits = c(-1, 1)) +
        theme_void() +
        theme(plot.background = element_rect(fill = "white", colour = NA))
    ggsave(path, plot = p, width = 7, height = 3.2, dpi = 300, bg = "white")
    message("Wrote placeholder plot: ", reason)
}

plot_title <- sprintf("%s: %s", args$assoc, spec$title)

res <- tryCatch(fread(args$input_file), error = function(e) NULL)
if (is.null(res) || nrow(res) == 0L) {
    write_placeholder(args$output, plot_title,
                      paste("No results to plot. The enrichment table is empty --",
                            "see the run log for why the test was skipped."))
    quit(save = "no", status = 0)
}
for (nm in c(spec$label, spec$count, "p_value", "fdr")) {
    if (!nm %in% names(res)) {
        stop(sprintf("column '%s' missing from %s", nm, args$input_file))
    }
}

dt <- as.data.table(res)
dt <- dt[!is.na(p_value) & is.finite(p_value)]
if (nrow(dt) == 0L) {
    write_placeholder(args$output, plot_title,
                      "No results carried a usable p-value.")
    quit(save = "no", status = 0)
}

# Break ties on effect size where there is one, so a saturated p-value column
# still produces a stable, meaningful top N.
if ("fold_enrichment" %in% names(dt)) {
    setorder(dt, p_value, -fold_enrichment)
} else {
    setorder(dt, p_value)
}
per_group <- !is.null(spec$per_group)
if (per_group) {
    if (!spec$per_group %in% names(dt)) {
        stop(sprintf("column '%s' missing from %s", spec$per_group,
                     args$input_file))
    }
    dt[, .facet := NA_character_]
    # Strongest hit per knowledgebase, "strongest" meaning smallest p-value;
    # the table is already ordered by p with effect size as the tiebreak.
    top <- dt[, head(.SD, 1L), by = c(spec$per_group)]
    faceted <- FALSE
} else {

faceted <- !is.null(spec$facet) && spec$facet %in% names(dt) &&
           uniqueN(dt[[spec$facet]]) > 0L
if (faceted) {
    dt[, .facet := as.character(get(spec$facet))]
    # Top N per testing family, not top N overall.
    top <- dt[, head(.SD, args$top_n), by = .facet]
} else {
    dt[, .facet := NA_character_]
    top <- head(dt, args$top_n)
}
}

top[, .label := as.character(get(spec$label))]
top[, .count := suppressWarnings(as.numeric(get(spec$count)))]
top[, .count := fifelse(is.na(.count), 0, .count)]
top[, .signif := !is.na(fdr) & fdr < args$fdr_threshold]

# A strong enrichment over a large feature can produce a p-value below the
# smallest positive double, where p_value underflows to exactly 0 and
# -log10(0) is Inf -- which ggplot drops silently, giving an empty panel with
# a full legend.
#
# enrich_features.R and enrich_traits.R therefore also emit neg_log10_p,
# computed on the log scale, which is exact; use it when present. gometh
# returns only a linear p-value, so the pathways table has no such column and
# still needs the cap below.
# The axis is -log10(FDR), following knowYourCG, whose KYCG_plotDot and
# KYCG_plotBar both default to y = "-log10(FDR)". It also makes the figure
# self-consistent: the threshold that decides solid vs hollow points sits on
# the axis, marked with a reference line.
LOGP_CAP <- -log10(.Machine$double.xmin)
if ("neg_log10_fdr" %in% names(top)) {
    top[, .logp := as.numeric(neg_log10_fdr)]
    n_capped <- 0L
} else {
    # gometh returns a linear FDR only, so the pathways table still needs the
    # cap: an FDR below the smallest positive double reads as exactly 0.
    top[, .logp := -log10(fdr)]
    n_capped <- sum(!is.finite(top$.logp) | top$.logp > LOGP_CAP)
    top[, .logp := pmin(fifelse(is.finite(.logp), .logp, LOGP_CAP), LOGP_CAP)]
}

if (per_group) {
    top[, .label := paste0(get(spec$per_group), " (", get(spec$label), ")")]
    top[, .group := NA_character_]
} else if (!is.null(spec$group) && spec$group %in% names(top)) {
# Two knowledgebases can use the same feature name, and a duplicated label
# would collapse two rows onto one y position.
    top[, .group := as.character(get(spec$group))]
    top[, .label := fifelse(duplicated(.label) | duplicated(.label, fromLast = TRUE),
                            paste0(.label, " (", .group, ")"), .label)]
} else {
    top[, .label := as.character(get(spec$label))]
    top[, .group := NA_character_]
}
top[, .wrapped := wrap_labels(.label, args$label_width)]
top[, .wrapped := factor(.wrapped, levels = .wrapped[order(.logp)])]

n_sig <- sum(top$.signif)

# Fallback for a collapsed p-value axis: if the p-values available to this plot
# have underflowed, capping them puts every point on one vertical line and
# hides the ranking, so plot effect size instead. Both the axis title and the
# subtitle say which is shown, so the two cases are never confused.
#
# With neg_log10_p present this is unreachable; it remains for the pathways
# table, which carries only gometh's linear p-value.
use_effect <- (isTRUE(spec$force_effect) ||
               n_capped > nrow(top) / 2) && "fold_enrichment" %in% names(top)
if (use_effect) {
    top[, .x := suppressWarnings(as.numeric(fold_enrichment))]
    top <- top[is.finite(.x)]
    x_lab <- "Fold enrichment (observed / expected)"
} else {
    top[, .x := .logp]
    x_lab <- expression(-log[10]~"(FDR)")
}
if (nrow(top) == 0L) {
    write_placeholder(args$output, plot_title,
                      "No result had both a usable p-value and effect size.")
    quit(save = "no", status = 0)
}
# With free y scales each panel shows only its own levels, so ordering the
# global factor by (facet, x) makes every panel read ascending internally.
top[, .wrapped := factor(as.character(.wrapped),
                         levels = as.character(.wrapped)[order(.facet, .x)])]

if (faceted) {
    per <- dt[, .N, by = .facet][order(.facet)]
    tested_desc <- paste(sprintf("%s %d", per$.facet, per$N), collapse = ", ")
    head_desc <- sprintf(
        "top %d per collection, each corrected within itself and on its own x axis (tested: %s)",
        args$top_n, tested_desc)
} else if (per_group) {
    head_desc <- sprintf(
        "strongest hit from each of %d knowledgebases (%d features tested); x shows effect size, which is comparable across knowledgebases",
        nrow(top), nrow(dt))
} else {
    head_desc <- sprintf("top %d of %d tested by FDR", nrow(top), nrow(dt))
}
subtitle <- sprintf(
    "%s; %d at FDR < %g%s%s",
    head_desc, n_sig, args$fdr_threshold,
    if (n_sig == 0L) " (no significant results -- all points hollow)" else "",
    if (n_capped > 0L) sprintf(
        "\n%d p-value%s underflowed to zero, so the x axis shows %s",
        n_capped, if (n_capped == 1L) "" else "s",
        if (use_effect) "fold enrichment" else
            sprintf("-log10(p) capped at %.0f", LOGP_CAP)) else "")

has_groups <- !all(is.na(top$.group)) && uniqueN(top$.group) > 1L

# Built up front rather than patched onto p$layers[[i]] afterwards: the layer
# index is not stable (the FDR reference line is also a layer), and indexing
# into it silently coloured the wrong thing.
point_layer <- if (has_groups) {
    geom_point(aes(shape = .signif, colour = .group), stroke = 0.9)
} else {
    geom_point(aes(shape = .signif), colour = spec$colour, stroke = 0.9)
}

p <- ggplot(top, aes(x = .x, y = .wrapped, size = .count, alpha = .count)) +
    # The threshold that decides solid vs hollow, on the axis itself.
    geom_vline(xintercept = if (use_effect) NA_real_
                            else -log10(args$fdr_threshold),
               linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    point_layer +
    scale_shape_manual(
        values = c(`TRUE` = 19, `FALSE` = 1),
        breaks = c(TRUE, FALSE),
        labels = c(sprintf("FDR < %g", args$fdr_threshold),
                   sprintf("FDR >= %g", args$fdr_threshold)),
        drop = FALSE, name = NULL) +
    scale_alpha(range = c(0.65, 1), guide = "none") +
    scale_size(range = c(1.5, 6), name = paste("Number of", spec$unit)) +
    labs(x = x_lab, y = spec$axis,
         title = plot_title,
         # Wrapped per line: the faceted subtitle is long enough to be
         # clipped at this figure width, but strwrap() on the whole string
         # would also swallow the explicit break before the capping note.
         subtitle = paste(vapply(strsplit(subtitle, "\n", fixed = TRUE)[[1]],
                                 function(ln) paste(strwrap(ln, width = 95),
                                                    collapse = "\n"),
                                 character(1), USE.NAMES = FALSE),
                          collapse = "\n")) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          legend.direction = "vertical",
          plot.subtitle = element_text(colour = "grey30", size = 10),
          axis.text.y = element_text(size = 9, lineheight = 1))

if (has_groups) {
    n_groups <- uniqueN(top$.group)
    p <- p + guides(colour = guide_legend(
        title = tools::toTitleCase(gsub("_", " ", spec$group)),
        override.aes = list(size = 2.5, alpha = 1, shape = 19)))
    # Dark2 carries eight colours; fall back to viridis beyond that.
    p <- p + if (n_groups <= 8L) {
        scale_colour_brewer(palette = "Dark2")
    } else {
        scale_colour_viridis_d(option = "turbo", end = 0.9)
    }
}

if (faceted && uniqueN(top$.facet) > 1L) {
    # facet_wrap, not facet_grid: facet_grid with rows only SHARES the x axis
    # across panels, and GO p-values are routinely orders of magnitude smaller
    # than KEGG's (they can underflow outright), which would crush the KEGG
    # panel onto the axis origin and lose its ranking entirely. facet_wrap
    # gives each panel its own x axis.
    #
    # Separate axes are also the honest choice here: each collection is
    # corrected within itself, so the two panels are NOT directly comparable
    # and should not share a scale that implies they are. Significance is read
    # from the point fill, which is per-collection, not from the x position.
    #
    # The cost is equal panel heights even when one collection has fewer terms;
    # proportional heights would need ggh4x::force_panelsizes, which is not
    # worth a dependency for the cosmetics.
    p <- p + facet_wrap(~ .facet, ncol = 1, scales = "free")
}

# Height grows with the number of rows so wrapped labels do not collide.
height <- max(3.2, 1.6 + 0.52 * nrow(top) +
                   0.16 * sum(lengths(regmatches(as.character(top$.wrapped),
                                                 gregexpr("\n", as.character(top$.wrapped))))))
ggsave(args$output, plot = p, width = 9.5, height = height, dpi = 300,
       bg = "white", limitsize = FALSE)
message(sprintf("Wrote %s (%d rows, %d significant, height %.1f in)",
                basename(args$output), nrow(top), n_sig, height))
