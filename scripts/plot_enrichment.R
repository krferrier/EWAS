#!/usr/bin/env Rscript
#
# One summary plot per enrichment analysis, in one of two geometries.
#
# KYCG features use the shape of knowYourCG::KYCG_plotEnrichAll: every tested
# knowledgebase gets a block along the x axis, its hits sit above it, y is
# -log10(FDR) and point size is the effect size. Because each knowledgebase
# occupies its own block, nothing in the figure implies a ranking ACROSS
# knowledgebases -- which is what per-knowledgebase FDR correction requires,
# and which a single pooled ranking would have implied. Sets with no hits keep
# their block, so "tested and found nothing" is visible rather than absent.
#
# GO/KEGG terms and EWAS Atlas traits use a dot plot of the top N, with
# -log10(FDR) on the x axis, point size and transparency carrying how many
# CpGs or genes drive the result, and solid versus hollow points separating
# what passes FDR from what does not.
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
    # Label placement for the feature plots; the dot-plot kinds do not need it.
    library(ggrepel)
})

parser <- argparse::ArgumentParser(
    description = "Plot the top enrichment results for one analysis")
parser$add_argument('--input-file', required = TRUE,
                    help = "Enrichment results table from one of the enrich_* rules")
parser$add_argument('--kind', required = TRUE,
                    choices = c("features", "features_qc", "pathways", "traits"),
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
# `geom` picks the geometry: "enrichall" for the KYCG features, "dot" for the
# pathway and trait tables.
#
# `facet` matters statistically, not just visually. enrich_pathways.R calls
# gometh once per collection and keeps each collection's own FDR, so GO and
# KEGG are SEPARATE testing families -- and GO carries some 22,000 terms
# against KEGG's ~350, so a pooled top N is GO-dominated and hides KEGG
# entirely. Those get one panel each, with top N taken within each. The
# feature table is corrected per knowledgebase for the same reason, but with
# up to 17 of them faceting is unreadable, so it uses the enrichall geometry
# instead: one block per knowledgebase along the x axis.
SPEC <- list(
    features = list(
        geom = "enrichall", role = "biological",
        label = "feature", group = "knowledgebase", facet = NULL,
        count = "n_overlap",
        unit = "CpGs", axis = "Feature",
        title = "KYCG feature enrichment", colour = "#1B9E77"),
    # The design and QC knowledgebases, kept out of the biological figure.
    # Enrichment here is not a finding about biology: it says the hit list
    # tracks array design, artefact-prone regions or CpG density. An empty
    # version of this plot is the good outcome.
    features_qc = list(
        geom = "enrichall", role = "qc",
        label = "feature", group = "knowledgebase", facet = NULL,
        count = "n_overlap",
        unit = "CpGs", axis = "Feature",
        title = "KYCG post-hoc QC checks", colour = "#999999",
        qc_note = paste("These sets are technical, not biological. Enrichment",
                        "here suggests the hits track probe design,",
                        "mappability or local CpG density rather than biology,",
                        "so it is a reason to be careful with the other",
                        "enrichment results -- not a finding in itself. No",
                        "points is the expected outcome.")),
    pathways = list(
        geom = "dot",
        label = "term", group = NULL, facet = "collection",
        count = "n_significant_genes",
        unit = "Genes", axis = "Term",
        title = "GO and KEGG enrichment", colour = "#D95F02"),
    traits = list(
        geom = "dot",
        label = "trait", group = NULL, facet = NULL, count = "n_overlap",
        unit = "CpGs", axis = "Trait",
        title = "EWAS Atlas trait enrichment", colour = "#7570B3")
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

#' One block per tested knowledgebase along x, -log10(FDR) on y.
#'
#' After knowYourCG::KYCG_plotEnrichAll. Block width grows with log(set size)
#' so a 1188-member set does not swallow the axis and a 2-member set stays
#' findable, and every tested knowledgebase keeps its block whether or not it
#' has hits.
plot_enrich_all <- function(dt, spec, args, plot_title, path) {
    CAP <- 40          # -log10(FDR) beyond this is capped; the reference uses 25
    N_LABEL <- 14

    dt <- dt[!is.na(fdr) & !is.na(fold_enrichment) & is.finite(fold_enrichment)]
    if (nrow(dt) == 0L) {
        write_placeholder(path, plot_title,
                          "No result had both a usable FDR and effect size.")
        return(invisible(NULL))
    }

    # Smaller sets first, so the wide ones do not crowd the left edge.
    gp_size <- sort(table(dt$knowledgebase))
    gp_width <- log(2 + gp_size)
    dt[, .kb := factor(knowledgebase, levels = names(gp_size))]
    setorder(dt, .kb, feature)
    dt[, .inc := as.numeric((gp_width / gp_size)[as.character(.kb)])]
    dt[, .step := c(0, ifelse(head(as.character(.kb), -1) !=
                              tail(as.character(.kb), -1), 1, 0))]
    dt[, .xpos := cumsum(.inc + .step)]

    # Only enriched, significant features are drawn. Depletion is not shown:
    # the test is one-sided (phyper lower.tail = FALSE), so a fold below 1
    # carries no evidence here.
    pts <- dt[fold_enrichment > 1 & fdr < args$fdr_threshold]
    n_capped <- pts[neg_log10_fdr > CAP, .N]
    if (nrow(pts)) pts[, .y := pmin(neg_log10_fdr, CAP * 1.08)]

    blocks <- dt[, .(beg = min(.xpos), middle = mean(.xpos), end = max(.xpos),
                     n = .N), by = .kb]
    blocks[, lab := sprintf("%s (%d)", .kb, n)]
    # A one-feature set has beg == end, which draws a zero-length segment and
    # so renders as a label with no block. Give every block a floor width.
    min_w <- 0.006 * max(1e-9, diff(range(dt$.xpos)))
    blocks[, half := pmax((end - beg) / 2, min_w)]
    blocks[, `:=`(beg = middle - half, end = middle + half)]

    # The axis follows the data rather than the cap: on the QC figure the
    # points sit near the threshold, and defaulting to the cap would leave the
    # panel almost empty.
    y_top <- if (nrow(pts)) max(6, max(pts$.y) * 1.14) else 6
    sub <- sprintf("%d feature%s enriched at FDR < %g, of %d tested across %d knowledgebases; FDR corrected within each",
                   nrow(pts), if (nrow(pts) == 1L) "" else "s",
                   args$fdr_threshold, nrow(dt), uniqueN(dt$.kb))
    if (!is.null(spec$qc_note)) sub <- paste0(sub, "\n", spec$qc_note)
    if (nrow(pts) == 0L) {
        sub <- paste0(sub, "\nNothing reached the threshold in any set.",
                      if (!is.null(spec$qc_note)) "" else
                          " Blocks below the axis show what was tested.")
    }
    if (n_capped > 0L) {
        sub <- paste0(sub, sprintf(
            "\n%d point%s capped at -log10(FDR) = %d (dotted line).",
            n_capped, if (n_capped == 1L) "" else "s", CAP))
    }

    p <- ggplot(pts, aes(.xpos, .y)) +
        geom_hline(yintercept = -log10(args$fdr_threshold), linetype = "dashed",
                   colour = "grey55", linewidth = 0.3)
    if (n_capped > 0L) {
        p <- p + geom_hline(yintercept = CAP, linetype = "dotted",
                            colour = "grey60", linewidth = 0.3)
    }
    if (nrow(pts)) {
        p <- p +
            geom_point(aes(size = log2_odds_ratio, colour = .kb), alpha = 0.65) +
            geom_text_repel(
                data = pts[order(-neg_log10_fdr)][seq_len(min(N_LABEL, .N))],
                aes(label = feature, colour = .kb), size = 2.9,
                direction = "y", nudge_y = 0.25, max.overlaps = 100,
                segment.colour = "grey70", segment.size = 0.2,
                show.legend = FALSE)
    }
    p <- p +
        # One segment per tested knowledgebase: a set with no hits is still
        # shown to have been tested.
        geom_segment(data = blocks,
                     aes(x = beg, xend = end, y = -0.055 * y_top,
                         yend = -0.055 * y_top, colour = .kb),
                     linewidth = 2.5) +
        geom_text(data = blocks,
                  aes(middle, -0.075 * y_top, label = lab, colour = .kb),
                  vjust = 1, hjust = 1, angle = 35, size = 2.9,
                  show.legend = FALSE) +
        scale_colour_discrete(guide = "none") +
        scale_size_continuous(range = c(1.5, 6),
                              name = expression(log[2] ~ "(odds ratio)")) +
        scale_x_continuous(expand = expansion(mult = 0.045)) +
        coord_cartesian(clip = "off", ylim = c(-0.06 * y_top, y_top)) +
        labs(x = NULL, y = expression(-log[10] ~ "(FDR)"), title = plot_title,
             subtitle = paste(vapply(strsplit(sub, "\n", fixed = TRUE)[[1]],
                                     function(l) paste(strwrap(l, 95),
                                                       collapse = "\n"),
                                     character(1), USE.NAMES = FALSE),
                              collapse = "\n")) +
        theme_bw(base_size = 13) +
        theme(legend.position = "right",
              plot.subtitle = element_text(colour = "grey30", size = 10),
              plot.margin = margin(6, 6, 62, 16),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_blank(), axis.ticks.x = element_blank())

    ggsave(path, plot = p, width = 9.5, height = 6.6, dpi = 300, bg = "white",
           limitsize = FALSE)
    message(sprintf("Wrote %s (%d points, %d knowledgebase blocks, %d capped)",
                    basename(path), nrow(pts), nrow(blocks), n_capped))
}

plot_title <- sprintf("%s: %s", args$assoc, spec$title)

res <- tryCatch(fread(args$input_file), error = function(e) NULL)
if (is.null(res) || nrow(res) == 0L) {
    write_placeholder(args$output, plot_title,
                      paste("No results to plot. The enrichment table is empty --",
                            "see the run log for why the test was skipped."))
    quit(save = "no", status = 0)
}
need <- c(spec$label, spec$count, "p_value", "fdr")
if (identical(spec$geom, "enrichall")) {
    need <- c(need, "knowledgebase", "fold_enrichment", "log2_odds_ratio",
              "neg_log10_fdr")
}
for (nm in need) {
    if (!nm %in% names(res)) {
        stop(sprintf("column '%s' missing from %s", nm, args$input_file))
    }
}

dt <- as.data.table(res)
dt <- dt[!is.na(p_value) & is.finite(p_value)]

# The feature table holds both the biological and the QC knowledgebases, so
# each figure takes its own slice. A table written before the role column
# existed is treated as all-biological.
if (!is.null(spec$role)) {
    if ("role" %in% names(dt)) {
        dt <- dt[role == spec$role]
    } else if (!identical(spec$role, "biological")) {
        write_placeholder(args$output, plot_title,
                          paste("No role column in the enrichment table, so",
                                "the QC knowledgebases cannot be separated",
                                "out. Re-run enrich_features."))
        quit(save = "no", status = 0)
    }
}

if (nrow(dt) == 0L) {
    write_placeholder(args$output, plot_title,
                      if (is.null(spec$role))
                          "No results carried a usable p-value."
                      else sprintf(
                          "No %s knowledgebases were tested. Check enrichment.%s in the config.",
                          spec$role,
                          if (identical(spec$role, "qc")) "qc_knowledgebases"
                          else "knowledgebases"))
    quit(save = "no", status = 0)
}

if (identical(spec$geom, "enrichall")) {
    plot_enrich_all(dt, spec, args, plot_title, args$output)
    quit(save = "no", status = 0)
}

# Break ties on effect size where there is one, so a saturated p-value column
# still produces a stable, meaningful top N.
if ("fold_enrichment" %in% names(dt)) {
    setorder(dt, p_value, -fold_enrichment)
} else {
    setorder(dt, p_value)
}
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

top[, .label := as.character(get(spec$label))]
top[, .count := suppressWarnings(as.numeric(get(spec$count)))]
top[, .count := fifelse(is.na(.count), 0, .count)]
# A factor carrying BOTH levels, not a bare logical: the shape guide should
# always show the solid/hollow meaning, including on a plot where every point
# happens to be significant. drop = FALSE on the scale cannot keep a level that
# a logical vector has no way to declare, so a one-sided plot would otherwise
# get a one-entry legend that explains nothing.
top[, .signif := factor(!is.na(fdr) & fdr < args$fdr_threshold,
                        levels = c(TRUE, FALSE))]

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

if (!is.null(spec$group) && spec$group %in% names(top)) {
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

n_sig <- sum(top$.signif == "TRUE")

# Fallback for a collapsed p-value axis: if the p-values available to this plot
# have underflowed, capping them puts every point on one vertical line and
# hides the ranking, so plot effect size instead. Both the axis title and the
# subtitle say which is shown, so the two cases are never confused.
#
# With neg_log10_p present this is unreachable; it remains for the pathways
# table, which carries only gometh's linear p-value.
use_effect <- n_capped > nrow(top) / 2 && "fold_enrichment" %in% names(top)
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

# A fully invisible layer carrying BOTH significance levels, so the shape guide
# always explains solid vs hollow -- including on a plot where every point
# passed the threshold, or none did. drop = FALSE on the scale is not enough:
# ggplot2 (4.0.3) will emit the key label for an unused level but has no data
# row to draw a glyph from, so the entry renders as a bare label. override.aes
# cannot fill that in either; the key needs to exist. Anchored on the first
# real point rather than at Inf (which bleeds a clipped mark into the panel
# corner), at size 0 and alpha 0, and carrying .facet so facet_wrap does not
# open an extra NA panel. Verified not to alter the panel ranges.
key_layer <- geom_point(
    data = data.table(.x = top$.x[1], .wrapped = top$.wrapped[1],
                      .facet = top$.facet[1],
                      .signif = factor(c(TRUE, FALSE), levels = c(TRUE, FALSE))),
    aes(x = .x, y = .wrapped, shape = .signif),
    size = 0, alpha = 0, inherit.aes = FALSE)

p <- ggplot(top, aes(x = .x, y = .wrapped, size = .count, alpha = .count)) +
    # The threshold that decides solid vs hollow, on the axis itself.
    geom_vline(xintercept = if (use_effect) NA_real_
                            else -log10(args$fdr_threshold),
               linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    point_layer +
    key_layer +
    scale_shape_manual(
        values = c(`TRUE` = 19, `FALSE` = 1),
        breaks = c("TRUE", "FALSE"),
        limits = c("TRUE", "FALSE"),
        labels = c(sprintf("FDR < %g", args$fdr_threshold),
                   sprintf("FDR >= %g", args$fdr_threshold)),
        drop = FALSE, name = NULL) +
    scale_alpha(range = c(0.65, 1), guide = "none") +
    # Point size and alpha are data-mapped, so the key for a level with no rows
    # in this plot has nothing to draw from and renders as a bare label. Fix
    # the key glyphs explicitly. Colour is neutral when colour already encodes
    # the knowledgebase, so the shape key cannot be read as a group.
    guides(shape = guide_legend(
        override.aes = list(size = 2.5, alpha = 1,
                            colour = if (has_groups) "grey25" else spec$colour))) +
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
