#!/usr/bin/env Rscript
#
# One summary plot per enrichment analysis, in one of two geometries.
#
# Similar plotting to knowYourCG::KYCG_plotEnrichAll: every tested
# knowledgebase has a block along the x axis and y is -log10(FDR).
# Point size is the effect size. Sets with no hits still have a
# block, so "tested and found nothing" is visible.
#
# GO/KEGG terms and EWAS Atlas traits use a dot plot of the top N, with
# -log10(FDR) on the x axis, point size and transparency showing how many
# CpGs or genes drive the result, and solid versus hollow points separating
# what passes FDR from what does not.
#
# The plot is ALWAYS written, including when the enrichment table is empty or
# nothing reached significance. This keeps the output trackable by Snakemake. 
# The same solid versus hollow dots to show significance is used.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(ggtext)
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

# Define specs aligned with the OUT_COLS in the enrich_* scripts.
SPEC <- list(
    features = list(
        geom = "enrichall", role = "biological",
        label = "feature", group = "knowledgebase", facet = NULL,
        count = "n_overlap",
        unit = "CpGs", axis = "Feature",
        title = "KYCG feature enrichment", colour = "#1B9E77"),
    # The design and QC knowledgebases are plotted separately from
    # biological knowledgebases.
    features_qc = list(
        geom = "enrichall", role = "qc",
        label = "feature", group = "knowledgebase", facet = NULL,
        count = "n_overlap",
        unit = "CpGs", axis = "Feature",
        title = "KYCG post-hoc QC checks", colour = "#999999"),
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

#' Wrap long labels
wrap_labels <- function(x, width) {
    vapply(as.character(x), function(s) {
        if (is.na(s) || !nzchar(s)) return("(unnamed)")
        paste(strwrap(s, width = width), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
}

#' Placeholder used when there is nothing to plot.
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

# Display order for the knowledgebase blocks, grouped so that sets showing
# related information are adjacent: the two chromatin-state models
# beside the histone marks, the two repeat resolutions together, and so on. 
# Alternate shading occurs between these groups.
#
# A set not listed here (MSA has 15 features that EPIC does not) is appended in
# an "other" group, alphabetically.
KB_GROUPS <- list(
    `chromatin state`   = c("ChromHMM", "REMCChromHMM", "ChromHMMfullStack", "HM"),
    `protein binding`   = c("TFBSrm", "TFBSconsensus", "CTCFbind"),
    `gene context`      = c("CGI", "MetagenePC", "GeneFeatures"),
    `large-scale domain` = c("ABCompartment", "PMD"),
    repeats             = c("rmsk1", "rmsk2"),
    `sequence context`  = c("Tetranuc2", "nFlankCG", "GCfrac"),
    imprinting          = c("ImprintingDMR"),
    `array design`      = c("ProbeType", "InfiniumChemistry", "Blacklist")
)

#' One block per tested knowledgebase along x, -log10(FDR) on y.
#'
#' After knowYourCG::KYCG_plotEnrichAll. Block width grows with log(set size)
#' so a large feature set does not dominate the axis.
plot_enrich_all <- function(dt, spec, args, plot_title, path) {
    N_LABEL <- 14

    dt <- dt[!is.na(fdr) & !is.na(fold_enrichment) & is.finite(fold_enrichment)]
    if (nrow(dt) == 0L) {
        write_placeholder(path, plot_title,
                          "No result had both a usable FDR and effect size.")
        return(invisible(NULL))
    }

    # Group the knowledgebases, then order blocks by group.
    grp <- data.table(
        kb = unlist(KB_GROUPS, use.names = FALSE),
        grp = rep(names(KB_GROUPS), lengths(KB_GROUPS)),
        gi = rep(seq_along(KB_GROUPS), lengths(KB_GROUPS)))
    grp[, ki := seq_len(.N)]
    present <- data.table(kb = sort(unique(dt$knowledgebase)))
    present <- grp[present, on = "kb"]
    present[is.na(gi), `:=`(grp = "other", gi = length(KB_GROUPS) + 1L,
                            ki = nrow(grp) + seq_len(.N))]
    setorder(present, gi, ki)

    gp_size <- table(dt$knowledgebase)[present$kb]
    gp_width <- log(2 + gp_size)
    dt[, .kb := factor(knowledgebase, levels = present$kb)]
    setorder(dt, .kb, feature)
    dt[, .inc := as.numeric((gp_width / gp_size)[as.character(.kb)])]
    dt[, .step := c(0, ifelse(head(as.character(.kb), -1) !=
                              tail(as.character(.kb), -1), 1, 0))]
    dt[, .xpos := cumsum(.inc + .step)]

    # Every enriched feature is shown: the solid/hollow encoding is what 
    # shows significant vs not.
    #
    # Features with fold enrichment <= 1 are left out. The hypergeometric here
    # is one-sided (phyper lower.tail = FALSE), so depletion is not tested and
    # such a feature carries no evidence either way.
    pts <- dt[fold_enrichment > 1]
    n_sig <- pts[fdr < args$fdr_threshold, .N]
    pts[, .y := neg_log10_fdr]
    # A factor carrying BOTH levels, so the legend explaining solid
    # vs hollow is always shown in the output plot.
    pts[, .signif := factor(!is.na(fdr) & fdr < args$fdr_threshold,
                            levels = c(TRUE, FALSE))]

    blocks <- dt[, .(beg = min(.xpos), middle = mean(.xpos), end = max(.xpos),
                     n = .N), by = .kb]
    blocks <- present[, .(.kb = kb, grp, gi)][blocks, on = ".kb"]
    setorder(blocks, middle)
    blocks[, lab := sprintf("%s (%d)", .kb, n)]
    # A one-feature set has beg == end, so its band would have zero width.
    # Give every block a base width, and split the gap between neighbours so
    # the shading tiles the axis without overlapping.
    min_w <- 0.006 * max(1e-9, diff(range(dt$.xpos)))
    blocks[, half := pmax((end - beg) / 2, min_w)]
    blocks[, `:=`(beg = middle - half, end = middle + half)]
    pad <- if (nrow(blocks) > 1L) min(diff(blocks$middle)) * 0.12 else min_w
    blocks[, `:=`(rmin = beg - pad, rmax = end + pad)]

    # Shade alternate GROUPS, so a band covers the related sets together.
    bands <- blocks[, .(rmin = min(rmin), rmax = max(rmax)), by = .(gi, grp)]
    setorder(bands, rmin)
    bands[, shade := seq_len(.N) %% 2L == 0L]

    # Colours are assigned explicitly rather than left to the default scale so
    # the same values can colour the axis labels: the tick label under a block
    # then matches the points above it. 
    kb_cols <- scales::hue_pal()(nrow(blocks))
    names(kb_cols) <- as.character(blocks$.kb)
    blocks[, lab_md := sprintf("<span style='color:%s'>%s</span>",
                               kb_cols[as.character(.kb)], lab)]

    # Headroom for the repelled labels, which sit above the highest points.
    y_top <- if (nrow(pts)) max(6, max(pts$.y) * 1.18) else 6

    p <- ggplot(pts, aes(.xpos, .y)) +
        geom_rect(data = bands[shade == TRUE],
                  aes(xmin = rmin, xmax = rmax, ymin = -Inf, ymax = Inf),
                  fill = "grey92", colour = NA, alpha = 0.55,
                  inherit.aes = FALSE)
    if (nrow(pts)) {
        # Label the strongest, significant hit in each knowledgebase.
        # Only significant hits are labelled.
        sig_pts <- pts[.signif == "TRUE"]
        labs_dt <- if (nrow(sig_pts)) {
            lab_kb <- sig_pts[, .(best = max(neg_log10_fdr)), by = .kb
                              ][order(-best)][seq_len(min(N_LABEL, .N)), .kb]
            sig_pts[.kb %in% lab_kb][order(-neg_log10_fdr),
                                     head(.SD, 1L), by = .kb]
        } else {
            sig_pts[0]
        }
        # layer so the hollow/filled legend is always created.
        key_layer <- geom_point(
            data = data.table(.xpos = pts$.xpos[1], .y = pts$.y[1],
                              .signif = factor(c(TRUE, FALSE),
                                               levels = c(TRUE, FALSE))),
            aes(x = .xpos, y = .y, shape = .signif),
            size = 0, alpha = 0, inherit.aes = FALSE)
        p <- p +
            geom_point(aes(size = log2_odds_ratio, colour = .kb,
                           shape = .signif), alpha = 0.65) +
            key_layer +
            # Repelled upward off its own block's points. point.padding has to
            # cover the largest points (6 mm) because ggrepel knows nothing of
            # point size; direction = "y" keeps a label over the block it
            # belongs to instead of drifting across a neighbour, and
            # min.segment.length = 0 keeps the connector drawn however short.
            geom_text_repel(
                data = labs_dt, aes(label = feature, colour = .kb), size = 2.9,
                direction = "y", nudge_y = 0.04 * y_top,
                point.padding = 0.6, box.padding = 0.35,
                min.segment.length = 0, force = 4, force_pull = 0.3,
                max.overlaps = Inf, seed = 1,
                segment.colour = "grey60", segment.size = 0.2,
                show.legend = FALSE)
    }
    y_breaks <- pretty(c(0, y_top))
    y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= y_top]
    p <- p +
        scale_colour_manual(values = kb_cols, guide = "none") +
        scale_shape_manual(
            values = c(`TRUE` = 19, `FALSE` = 1),
            breaks = c("TRUE", "FALSE"), limits = c("TRUE", "FALSE"),
            labels = c(sprintf("FDR < %g", args$fdr_threshold),
                       sprintf("FDR >= %g", args$fdr_threshold)),
            drop = FALSE, name = NULL) +
        scale_size_continuous(range = c(1.5, 6),
                              name = expression(log[2] ~ "(odds ratio)")) +
        guides(shape = guide_legend(
            override.aes = list(size = 2.5, alpha = 1, colour = "grey25"))) +
        # Knowledgebase axis labels are coloured to match their block's points.
        scale_x_continuous(breaks = blocks$middle, labels = blocks$lab_md,
                           expand = expansion(mult = 0.02)) +
        scale_y_continuous(breaks = y_breaks, expand = expansion(mult = 0.02)) +
        # x limits come from the blocks, so feature sets with no significant
        # points aren't dropped from the plot.
        coord_cartesian(xlim = c(min(bands$rmin), max(bands$rmax)),
                        ylim = c(0, y_top)) +
        labs(x = "Knowledgebase (features tested)",
             y = expression(-log[10] ~ "(FDR)"), title = plot_title) +
        theme_bw(base_size = 13) +
        theme(legend.position = "right",
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.text.x = element_markdown(angle = 40, hjust = 1, vjust = 1,
                                             size = 9),
              axis.ticks.x = element_line(colour = "grey70", linewidth = 0.3))

    ggsave(path, plot = p, width = 9.5, height = 6.6, dpi = 300, bg = "white",
           limitsize = FALSE)
    message(sprintf(paste("Wrote %s (%d blocks in %d groups; %d enriched",
                          "features drawn, %d of them at FDR < %g)"),
                    basename(path), nrow(blocks), uniqueN(blocks$grp),
                    nrow(pts), n_sig, args$fdr_threshold))
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
# was defined is treated as all-biological.
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
# A factor with both significant and non-significant levels so the legend
# always shows the solid/hollow definition.
top[, .signif := factor(!is.na(fdr) & fdr < args$fdr_threshold,
                        levels = c(TRUE, FALSE))]

# A strong enrichment over a large feature can produce a p-value below the
# smallest positive double, where p_value underflows to exactly 0 and
# -log10(0) is Inf -- which ggplot drops silently, giving an empty panel.
#
# enrich_features.R and enrich_traits.R therefore also emit neg_log10_p,
# computed on the log scale, which is exact; use it when present. gometh
# returns only a linear p-value, so the pathways table has no such column and
# still needs the cap below.
# The axis is -log10(FDR), following knowYourCG plots.
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
# A missing name falls back to the ID, and any label still repeated gets its
# ID appended: y positions are factor levels, and a repeated level is an error
# (and would otherwise put two results on one row).
id_col <- intersect(c("term_id", "feature", "trait"), names(top))[1]
if (!is.na(id_col)) {
    ids <- as.character(top[[id_col]])
    blank <- is.na(top$.label) | !nzchar(trimws(top$.label)) | top$.label == "NA"
    top[blank, .label := ids[blank]]
    dup <- duplicated(top$.label) | duplicated(top$.label, fromLast = TRUE)
    top[dup, .label := paste0(.label, " (", ids[dup], ")")]
}
top[, .label := make.unique(.label, sep = " ")]
top[, .wrapped := wrap_labels(.label, args$label_width)]
top[, .wrapped := make.unique(as.character(.wrapped), sep = " ")]
top[, .wrapped := factor(.wrapped, levels = .wrapped[order(.logp)])]

n_sig <- sum(top$.signif == "TRUE")

# Fallback for a collapsed p-value axis: if the p-values available to this plot
# have underflowed, capping them puts every point on one vertical line and
# hides the ranking, so plot effect size instead. The axis title says which
# quantity is shown, so the two cases are never confused.
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
    guides(shape = guide_legend(
        override.aes = list(size = 2.5, alpha = 1,
                            colour = if (has_groups) "grey25" else spec$colour))) +
    scale_size(range = c(1.5, 6), name = paste("Number of", spec$unit)) +
    labs(x = x_lab, y = spec$axis, title = plot_title) +
    theme_bw(base_size = 13) +
    theme(legend.position = "right",
          legend.direction = "vertical",
          axis.text.y = element_text(size = 9, lineheight = 1))

if (has_groups) {
    n_groups <- uniqueN(top$.group)
    p <- p + guides(colour = guide_legend(
        title = tools::toTitleCase(gsub("_", " ", spec$group)),
        override.aes = list(size = 2.5, alpha = 1, shape = 19)))
    p <- p + if (n_groups <= 8L) {
        scale_colour_brewer(palette = "Dark2")
    } else {
        scale_colour_viridis_d(option = "turbo", end = 0.9)
    }
}

if (faceted && uniqueN(top$.facet) > 1L) {
    p <- p + facet_wrap(~ .facet, ncol = 1, scales = "free")
}

# Height grows with the number of rows so wrapped labels do not collide.
height <- max(3.2, 1.6 + 0.52 * nrow(top) +
                   0.16 * sum(lengths(regmatches(as.character(top$.wrapped),
                                                 gregexpr("\n", as.character(top$.wrapped))))))
ggsave(args$output, plot = p, width = 9.5, height = height, dpi = 300,
       bg = "white", limitsize = FALSE)
message(sprintf("Wrote %s (%s; %d significant at FDR < %g; height %.1f in)",
                basename(args$output), head_desc, n_sig, args$fdr_threshold,
                height))
