#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(cowplot)
  library(jsonlite)
})

# ================================================================
# Arguments
# ================================================================

parser <- argparse::ArgumentParser(
  description = "Create DMR Manhattan and optional zoomed DMR plots with UCSC refGene gene tracks"
)
parser$add_argument("--slk-file",
                    required = TRUE,
                    help = "Path to comb-p/SLK CpG-level BED file, e.g. cotinine_ewas.slk.bed")
parser$add_argument("--regions-file",
                    required = TRUE,
                    help = "Path to significant DMR regions BED file, e.g. cotinine_ewas.regions-p.bed.gz")
parser$add_argument("--out-dir",
                    required = TRUE,
                    help = "Output directory")
parser$add_argument("--assoc",
                    required = TRUE,
                    help = "Association name used in output filenames, e.g. cotinine")
parser$add_argument("--min-probes",
                    type = "integer",
                    default = 2,
                    help = "Minimum number of CpGs required for a DMR to be highlighted")
parser$add_argument("--max-y",
                    type = "double",
                    default = -1,
                    help = "Optional y-axis maximum for genome-wide plot. If <= 0, calculated from data.")
parser$add_argument("--make-zoom",
                    choices = c("yes", "no"),
                    default = "yes",
                    help = "Whether to create zoomed-in plots for clusters of nearby multi-CpG DMRs")
parser$add_argument("--zoom-padding",
                    type = "integer",
                    default = 2000,
                    help = "Number of bp to extend on each side of a zoomed DMR cluster")
parser$add_argument("--cluster-gap",
                    type = "integer",
                    default = 3000,
                    help = "Maximum gap in bp between DMRs for grouping into the same zoom cluster")
parser$add_argument("--point-jitter-bp",
                    type = "integer",
                    default = 25,
                    help = "Horizontal jitter in bp for highlighted DMR CpGs in zoomed plots")
parser$add_argument("--zoom-midlines",
                    choices = c("yes", "no"),
                    default = "yes",
                    help = "Draw one faint vertical midpoint line per DMR in zoomed plots")
parser$add_argument("--genome-build",
                    choices = c("hg19", "hg38"),
                    default = "hg38",
                    help = "Genome build for UCSC refGene annotation")
parser$add_argument("--gene-label-size",
                    type = "double",
                    default = 3.0,
                    help = "Gene label size in zoomed gene track")
parser$add_argument("--make-combined",
                    choices = c("yes", "no"),
                    default = "yes",
                    help = "Whether to create a combined multi-panel figure from the Manhattan and zoomed plots")
parser$add_argument("--combined-formats",
                    default = "pdf",
                    help = "Comma-separated output formats for the combined figure, e.g. svg,pdf,png,jpg")
parser$add_argument("--combined-width-cm",
                    type = "double",
                    default = 18.3,
                    help = "Width of combined multi-panel figure in cm")
parser$add_argument("--combined-manhattan-height-cm",
                    type = "double",
                    default = 8.5,
                    help = "Target height in cm allocated to the Manhattan panel in the combined figure")
parser$add_argument("--combined-region-height-cm",
                    type = "double",
                    default = 10,
                    help = "Target height in cm allocated to each row of regional panels in the combined figure")
parser$add_argument("--panel-label-size",
                    type = "double",
                    default = 16,
                    help = "Font size for multi-panel labels A, B, C, ...")

args <- parser$parse_args()


# ================================================================
# Helpers
# ================================================================

chrom_to_num <- function(chrom) {
  x <- gsub("^chr", "", chrom, ignore.case = TRUE)

  out <- suppressWarnings(as.numeric(x))
  out[toupper(x) == "X"] <- 23
  out[toupper(x) == "Y"] <- 24
  out[toupper(x) %in% c("M", "MT")] <- 25

  out
}

get_gene_track_data_hgnc <- function(chr_string,
                                     region_start,
                                     region_end,
                                     genome_build = "hg38") {

  # UCSC uses 0-based starts internally.
  ucsc_start <- max(0, region_start - 1)
  ucsc_end <- region_end

  api_url <- paste0(
    "https://api.genome.ucsc.edu/getData/track?",
    "genome=", genome_build,
    ";track=hgnc",
    ";chrom=", chr_string,
    ";start=", ucsc_start,
    ";end=", ucsc_end
  )

  res <- jsonlite::fromJSON(api_url)

  if (!"hgnc" %in% names(res) || length(res$hgnc) == 0) {
    return(list(
      genes = data.frame(),
      exons = data.frame()
    ))
  }

  hgnc_df <- as.data.frame(res$hgnc)

  if (nrow(hgnc_df) == 0) {
    return(list(
      genes = data.frame(),
      exons = data.frame()
    ))
  }

  gene_df <- hgnc_df %>%
    transmute(
      gene_id = name,
      symbol = symbol,
      chr = chrom,
      # Convert UCSC 0-based chromStart to 1-based plotting coordinate.
      gene_start = as.numeric(chromStart) + 1,
      gene_end = as.numeric(chromEnd),
      strand = strand
    ) %>%
    filter(
      chr == chr_string,
      gene_start <= region_end,
      gene_end >= region_start
    ) %>%
    distinct(symbol, gene_start, gene_end, strand, .keep_all = TRUE) %>%
    arrange(gene_start, gene_end)

  # HGNC track is gene-span based, not transcript/exon based.
  # For plotting, draw one block per gene span.
  exon_df <- gene_df %>%
    transmute(
      gene_id = gene_id,
      exon_start = gene_start,
      exon_end = gene_end
    )

  list(
    genes = gene_df,
    exons = exon_df
  )
}

flag_cpgs_in_dmrs <- function(slk_df, dmr_df) {
  slk_df$in_multi_cpg_dmr <- FALSE

  if (nrow(dmr_df) == 0) {
    return(slk_df)
  }

  for (i in seq_len(nrow(dmr_df))) {
    hit <- slk_df$CHR == dmr_df$CHR[i] &
      slk_df$MAPINFO >= dmr_df$start[i] &
      slk_df$MAPINFO <= dmr_df$end[i]

    slk_df$in_multi_cpg_dmr[hit] <- TRUE
  }

  slk_df
}

cluster_dmrs <- function(dmr_df, max_gap = 3000) {
  if (nrow(dmr_df) == 0) {
    dmr_df$cluster_id <- integer(0)
    return(dmr_df)
  }

  dmr_df <- dmr_df %>%
    arrange(CHR, start, end)

  cluster_id <- integer(nrow(dmr_df))
  current_cluster <- 0L
  current_chr <- NA
  current_cluster_end <- NA_real_

  for (i in seq_len(nrow(dmr_df))) {
    start_new_cluster <- FALSE

    if (i == 1) {
      start_new_cluster <- TRUE
    } else if (dmr_df$CHR[i] != current_chr) {
      start_new_cluster <- TRUE
    } else if ((dmr_df$start[i] - current_cluster_end) > max_gap) {
      start_new_cluster <- TRUE
    }

    if (start_new_cluster) {
      current_cluster <- current_cluster + 1L
      current_cluster_end <- dmr_df$end[i]
    } else {
      current_cluster_end <- max(current_cluster_end, dmr_df$end[i])
    }

    cluster_id[i] <- current_cluster
    current_chr <- dmr_df$CHR[i]
  }

  dmr_df$cluster_id <- cluster_id
  dmr_df
}

format_list_arg <- function(x) {
  out <- unlist(strsplit(x, ",", fixed = TRUE))
  out <- trimws(tolower(out))
  out[nzchar(out)]
}

panel_label <- function(i) {
  # 1 -> A, 2 -> B, ..., 26 -> Z, 27 -> AA, etc.
  stopifnot(i >= 1)
  label <- ""
  while (i > 0) {
    i <- i - 1
    label <- paste0(LETTERS[(i %% 26) + 1], label)
    i <- i %/% 26
  }
  label
}

make_panel_labels <- function(n) {
  vapply(seq_len(n), panel_label, character(1))
}

save_plot_multi_format <- function(plot,
                                   file_prefix,
                                   formats,
                                   width,
                                   height,
                                   units = "cm",
                                   dpi = 300) {
  formats <- format_list_arg(formats)

  if (length(formats) == 0) {
    stop("--combined-formats must include at least one format, e.g. svg,pdf")
  }

  saved_files <- character(0)

  for (fmt in formats) {
    out_file <- paste0(file_prefix, ".", fmt)

    if (fmt == "svg") {
      if (!requireNamespace("svglite", quietly = TRUE)) {
        stop("Saving SVG output requires the svglite package. ",
             "Install it or remove svg from --combined-formats.")
      }

      ggsave(
        filename = out_file,
        plot = plot,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        device = svglite::svglite
      )
    } else if (fmt == "pdf") {
      ggsave(
        filename = out_file,
        plot = plot,
        width = width,
        height = height,
        units = units,
        dpi = dpi,
        device = grDevices::cairo_pdf
      )
    } else {
      ggsave(
        filename = out_file,
        plot = plot,
        width = width,
        height = height,
        units = units,
        dpi = dpi
      )
    }

    saved_files <- c(saved_files, out_file)
  }

  saved_files
}

make_combined_dmr_figure <- function(manhattan_plot,
                                     regional_plots = list(),
                                     width_cm = 18.3,
                                     manhattan_height_cm = 8.5,
                                     region_row_height_cm = 9.5,
                                     panel_label_size = 16) {
  regional_plots <- Filter(Negate(is.null), regional_plots)
  n_regions <- length(regional_plots)
  n_panels <- 1 + n_regions

  if (n_panels > 26) {
    warning(
      "Combined figure has more than 26 panels. ",
      "Labels after Z will continue as AA, AB, etc."
    )
  }

  labels <- make_panel_labels(n_panels)

  manhattan_panel <- cowplot::plot_grid(
    manhattan_plot,
    labels = labels[1],
    label_size = panel_label_size,
    label_fontface = "bold"
  )

  if (n_regions == 0) {
    return(list(
      plot = manhattan_panel,
      width_cm = width_cm,
      height_cm = manhattan_height_cm
    ))
  }

  region_ncol <- ifelse(n_regions == 1, 1, 2)
  region_nrow <- ceiling(n_regions / region_ncol)

  regional_grid <- cowplot::plot_grid(
    plotlist = regional_plots,
    ncol = region_ncol,
    labels = labels[-1],
    label_size = panel_label_size,
    label_fontface = "bold",
    align = "hv",
    axis = "tblr"
  )

  regional_height_cm <- region_nrow * region_row_height_cm
  total_height_cm <- manhattan_height_cm + regional_height_cm

  combined <- cowplot::plot_grid(
    manhattan_panel,
    regional_grid,
    ncol = 1,
    rel_heights = c(manhattan_height_cm, regional_height_cm),
    align = "v",
    axis = "lr"
  )

  list(
    plot = combined,
    width_cm = width_cm,
    height_cm = total_height_cm
  )
}


assign_gene_lanes <- function(gene_df,
                              plot_start,
                              plot_end,
                              label_bp_per_char = 0.018,
                              gap_frac = 0.004) {
  if (nrow(gene_df) == 0) return(gene_df)

  plot_span <- plot_end - plot_start
  gap_bp <- plot_span * gap_frac

  gene_df <- gene_df %>%
    mutate(
      plot_gene_start = pmax(gene_start, plot_start),
      plot_gene_end = pmin(gene_end, plot_end),
      gene_mid = (plot_gene_start + plot_gene_end) / 2,
      label_width_bp = pmax(
        plot_gene_end - plot_gene_start,
        nchar(symbol) * plot_span * label_bp_per_char
      ),
      label_start = pmax(plot_start, gene_mid - label_width_bp / 2),
      label_end = pmin(plot_end, gene_mid + label_width_bp / 2),
      interval_start = pmin(plot_gene_start, label_start),
      interval_end = pmax(plot_gene_end, label_end)
    ) %>%
    arrange(interval_start, interval_end)

  lane_end <- numeric(0)
  lane <- integer(nrow(gene_df))

  for (i in seq_len(nrow(gene_df))) {
    placed <- FALSE

    if (length(lane_end) > 0) {
      for (j in seq_along(lane_end)) {
        if (gene_df$interval_start[i] > lane_end[j] + gap_bp) {
          lane[i] <- j
          lane_end[j] <- gene_df$interval_end[i]
          placed <- TRUE
          break
        }
      }
    }

    if (!placed) {
      lane[i] <- length(lane_end) + 1
      lane_end <- c(lane_end, gene_df$interval_end[i])
    }
  }

  gene_df$lane <- lane
  gene_df$label_y <- gene_df$lane + 0.30

  gene_df
}

make_gene_track_plot <- function(gene_data,
                                 dmr_sub,
                                 dmr_cols,
                                 chr_string,
                                 plot_start,
                                 plot_end,
                                 gene_label_size = 3.0) {
  gene_df <- gene_data$genes
  exon_df <- gene_data$exons

  if (nrow(gene_df) == 0) {
    return(
      ggplot() +
        geom_rect(
          data = dmr_sub,
          aes(xmin = start, xmax = end,
              ymin = -Inf, ymax = Inf,
              fill = dmr_id),
          inherit.aes = FALSE,
          alpha = 0.4,
          color = NA
        ) +
        scale_fill_manual(values = dmr_cols, guide = "none") +
        scale_x_continuous(
          limits = c(plot_start, plot_end),
          labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
          expand = c(0, 0)
        ) +
        theme_bw(base_size = 12) +
        theme(
          panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        ) +
        labs(
          x = paste0(chr_string, " position"),
          y = "Genes\n(UCSC HGNC hg38)"
        )
    )
  }

  gene_df <- assign_gene_lanes(
    gene_df = gene_df,
    plot_start = plot_start,
    plot_end = plot_end
  )

  exon_df <- exon_df %>%
    left_join(
      gene_df %>%
        dplyr::select(gene_id, lane),
      by = "gene_id"
    ) %>%
    mutate(
      exon_start = pmax(exon_start, plot_start),
      exon_end = pmin(exon_end, plot_end)
    ) %>%
    filter(exon_start <= exon_end)

  max_lane <- max(gene_df$lane, na.rm = TRUE)

  ggplot() +
    geom_rect(
      data = dmr_sub,
      aes(xmin = start, xmax = end,
          ymin = -Inf, ymax = Inf,
          fill = dmr_id),
      inherit.aes = FALSE,
      alpha = 0.2,
      color = NA
    ) +
    geom_segment(
      data = gene_df,
      aes(x = plot_gene_start,
          xend = plot_gene_end,
          y = lane,
          yend = lane),
      color = "navy",
      linewidth = 0.45
    ) +
    geom_rect(
      data = exon_df,
      aes(xmin = exon_start,
          xmax = exon_end,
          ymin = lane - 0.12,
          ymax = lane + 0.12),
      fill = "navy",
      color = "navy",
      linewidth = 0.35
    ) +
    geom_text(
      data = gene_df,
      aes(x = gene_mid,
          y = label_y,
          label = symbol),
      size = gene_label_size,
      fontface = "bold",
      vjust = 0,
      color = "black"
    ) +
    scale_fill_manual(values = dmr_cols, guide = "none") +
    scale_x_continuous(
      limits = c(plot_start, plot_end),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.6, max_lane + 0.85),
      expand = c(0, 0)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 5.5, 5.5, 5.5)
    ) +
    labs(
      x = paste0(chr_string, " position"),
      y = "Genes\n(HGNC)"
    )
}

make_dmr_zoom_plot <- function(slk_df,
                               dmr_df,
                               chr,
                               chr_string,
                               region_start,
                               region_end,
                               out_file,
                               gene_table_file,
                               padding = 2000,
                               point_jitter_bp = 25,
                               show_dmr_midlines = TRUE,
                               title = NULL,
                               genome_build = "hg38",
                               gene_label_size = 3.0) {
  plot_start <- max(1, region_start - padding)
  plot_end <- region_end + padding

  slk_sub <- slk_df %>%
    filter(
      CHR == chr,
      MAPINFO >= plot_start,
      MAPINFO <= plot_end
    )

  dmr_sub <- dmr_df %>%
    filter(
      CHR == chr,
      end >= plot_start,
      start <= plot_end
    ) %>%
    arrange(start, end) %>%
    mutate(
      dmr_id = factor(seq_len(n())),
      mid = (start + end) / 2
    )

  if (nrow(slk_sub) == 0) {
    message("No CpGs found in zoom window: ", chr_string, ":",
            plot_start, "-", plot_end)
    return(NULL)
  }

  if (is.null(title)) {
    title <- paste0(
      "Zoomed DMR plot: ", chr_string, ":",
      scales::comma(plot_start), "-",
      scales::comma(plot_end)
    )
  }

  slk_sub$dmr_id <- NA_character_

  if (nrow(dmr_sub) > 0) {
    for (i in seq_len(nrow(dmr_sub))) {
      hit <- slk_sub$MAPINFO >= dmr_sub$start[i] &
        slk_sub$MAPINFO <= dmr_sub$end[i]

      slk_sub$dmr_id[hit & is.na(slk_sub$dmr_id)] <-
        as.character(dmr_sub$dmr_id[i])
    }
  }

  n_dmrs <- nrow(dmr_sub)
  dmr_cols <- scales::hue_pal()(n_dmrs)
  names(dmr_cols) <- levels(dmr_sub$dmr_id)

  gene_data <- get_gene_track_data_hgnc(
    chr_string = chr_string,
    region_start = plot_start,
    region_end = plot_end,
    genome_build = genome_build
  )

  if (nrow(gene_data$genes) > 0) {
    fwrite(
      gene_data$genes %>%
        arrange(gene_start, gene_end),
      gene_table_file,
      sep = "\t"
    )

    message(
      "Genes shown in ", basename(out_file), ": ",
      paste(unique(gene_data$genes$symbol), collapse = ", ")
    )
  } else {
    message("No refGene genes found in ", basename(out_file))
  }

  y_max <- max(1, max(slk_sub$neglogp, na.rm = TRUE) * 1.05)

  midline_layer <- if (show_dmr_midlines) {
    geom_vline(
      data = dmr_sub,
      aes(xintercept = mid, color = dmr_id),
      inherit.aes = FALSE,
      alpha = 0.30,
      linewidth = 0.85
    )
  } else {
    NULL
  }

  p_top <- ggplot(slk_sub, aes(x = MAPINFO, y = neglogp)) +
    geom_rect(
      data = dmr_sub,
      aes(xmin = start, xmax = end,
          ymin = -Inf, ymax = Inf,
          fill = dmr_id),
      inherit.aes = FALSE,
      alpha = 0.20,
      color = NA
    ) +
    midline_layer +
    geom_point(
      color = "grey45",
      alpha = 0.45,
      size = 1.4
    ) +
    geom_point(
      data = slk_sub %>% filter(!is.na(dmr_id)),
      aes(fill = dmr_id),
      shape = 21,
      color = "black",
      stroke = 0.7,
      size = 2.7,
      alpha = 0.75,
      position = position_jitter(
        width = point_jitter_bp,
        height = 0,
        seed = 1
      )
    ) +
    scale_fill_manual(values = dmr_cols, guide = "none") +
    scale_color_manual(values = dmr_cols, guide = "none") +
    scale_x_continuous(
      limits = c(plot_start, plot_end),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0.03, 0.05))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(5.5, 5.5, 3, 5.5)
    ) +
    labs(
      y = expression(-log[10]("regional P-value")),
      title = title
    )

  p_gene <- make_gene_track_plot(
    gene_data = gene_data,
    dmr_sub = dmr_sub,
    dmr_cols = dmr_cols,
    chr_string = chr_string,
    plot_start = plot_start,
    plot_end = plot_end,
    gene_label_size = gene_label_size
  )

  n_lanes <- if (nrow(gene_data$genes) > 0) {
    max(assign_gene_lanes(gene_data$genes, plot_start, plot_end)$lane)
  } else {
    1
  }

  gene_rel_height <- max(1.15, 0.75 + 0.35 * n_lanes)
  out_height <- max(10.5, 10 + 0.45 * n_lanes)

  combined_plot <- cowplot::plot_grid(
    p_top,
    p_gene,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(3.5, gene_rel_height)
  )

  ggsave(
    filename = out_file,
    plot = combined_plot,
    width = 18.3,
    height = out_height,
    units = "cm",
    dpi = 300
  )

  combined_plot
}

# ================================================================
# Read SLK CpG-level file
# ================================================================

slk <- data.table::fread(args$slk_file)

names(slk)[names(slk) == "#chrom"] <- "chrom"
names(slk)[names(slk) == "region-p"] <- "region_p"

if (!all(c("chrom", "start", "end", "p", "region_p") %in% names(slk)) &&
    ncol(slk) >= 5) {
  names(slk)[1:5] <- c("chrom", "start", "end", "p", "region_p")
}

stopifnot(all(c("chrom", "start", "end", "region_p") %in% names(slk)))

slk <- slk %>%
  mutate(
    CHR = chrom_to_num(chrom),
    CHR_label = gsub("^chr", "", chrom, ignore.case = TRUE),
    MAPINFO = as.numeric(start),
    Pvalue = as.numeric(region_p),
    neglogp = -log10(pmax(Pvalue, .Machine$double.xmin))
  ) %>%
  filter(
    !is.na(CHR),
    !is.na(MAPINFO),
    !is.na(Pvalue),
    Pvalue > 0
  )

# ================================================================
# Read significant DMR region file
# ================================================================

regions <- data.table::fread(args$regions_file)

names(regions)[names(regions) == "#chrom"] <- "chrom"

if (!all(c("chrom", "start", "end", "n_probes") %in% names(regions)) &&
    ncol(regions) >= 5) {
  names(regions)[1:5] <- c("chrom", "start", "end", "min_p", "n_probes")
}

stopifnot(all(c("chrom", "start", "end", "n_probes") %in% names(regions)))

dmrs_to_highlight <- regions %>%
  mutate(
    CHR = chrom_to_num(chrom),
    start = as.numeric(start),
    end = as.numeric(end),
    n_probes = as.integer(n_probes),
    mid = (start + end) / 2
  ) %>%
  filter(
    !is.na(CHR),
    n_probes >= args$min_probes
  ) %>%
  arrange(CHR, start, end)

dmrs_to_highlight <- cluster_dmrs(
  dmrs_to_highlight,
  max_gap = args$cluster_gap
)

slk <- flag_cpgs_in_dmrs(slk, dmrs_to_highlight)

message("DMRs plotted: ", nrow(dmrs_to_highlight))
message("CpGs highlighted in DMRs: ", sum(slk$in_multi_cpg_dmr))

# ================================================================
# Genome-wide DMR Manhattan plot
# ================================================================

chr.pos <- slk %>%
  group_by(CHR) %>%
  summarize(chr_len = max(MAPINFO, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  dplyr::select(CHR, tot)

chr.labels <- slk %>%
  group_by(CHR) %>%
  summarize(CHR_label = dplyr::first(CHR_label), .groups = "drop")

slk <- slk %>%
  left_join(chr.pos, by = "CHR") %>%
  arrange(CHR, MAPINFO) %>%
  mutate(POS = MAPINFO + tot)

dmrs_to_highlight <- dmrs_to_highlight %>%
  left_join(chr.pos, by = "CHR") %>%
  mutate(
    start_POS = start + tot,
    end_POS = end + tot,
    mid_POS = mid + tot
  )

x_axis <- slk %>%
  group_by(CHR) %>%
  summarize(center = mean(range(POS, na.rm = TRUE)), .groups = "drop") %>%
  left_join(chr.labels, by = "CHR") %>%
  arrange(CHR)

ymax <- if (is.null(args$max_y) || args$max_y <= 0) {
  ceiling(max(slk$neglogp, na.rm = TRUE))
} else {
  args$max_y
}

dmr_manhattan <- ggplot(slk, aes(x = POS, y = neglogp)) +
  geom_vline(
    data = dmrs_to_highlight,
    aes(xintercept = mid_POS),
    inherit.aes = FALSE,
    color = "red3",
    linewidth = 0.45,
    alpha = 0.45
  ) +
  geom_point(
    aes(color = as.factor(CHR)),
    alpha = 0.35,
    size = 0.55
  ) +
  geom_point(
    data = slk %>% filter(in_multi_cpg_dmr),
    color = "red3",
    alpha = 0.75,
    size = 1.05
  ) +
  scale_color_manual(
    values = rep(c("grey65", "grey35"), length.out = nrow(x_axis))
  ) +
  scale_x_continuous(
    labels = x_axis$CHR_label,
    breaks = x_axis$center,
    guide = guide_axis(check.overlap = TRUE),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, ymax),
    expand = expansion(mult = c(0, 0.03))
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10]("regional P-value")),
    #title = " "
    title = paste0(args$assoc, " DMR Manhattan plot")
  )

genome_out <- file.path(
  args$out_dir,
  paste0(args$assoc, "_dmr_manhattan.jpg")
)

ggsave(
  filename = genome_out,
  plot = dmr_manhattan,
  width = 18.3,
  height = 12,
  units = "cm",
  dpi = 300
)

message("Saved genome-wide DMR plot: ", genome_out)

# ================================================================
# Optional zoomed plots with UCSC refGene track
# ================================================================

zoom_plots <- list()

if (args$make_zoom == "yes" && nrow(dmrs_to_highlight) > 0) {

  zoom_windows <- dmrs_to_highlight %>%
    group_by(cluster_id, CHR, chrom) %>%
    summarize(
      cluster_start = min(start),
      cluster_end = max(end),
      n_dmrs = n(),
      .groups = "drop"
    ) %>%
    arrange(CHR, cluster_start)

  for (i in seq_len(nrow(zoom_windows))) {
    chr_i <- zoom_windows$CHR[i]
    chr_string <- zoom_windows$chrom[i]
    #title_i <- " "
    title_i <- paste0(
     "Zoomed DMR cluster ", i,
     " (", chr_string, "; ",
     zoom_windows$n_dmrs[i], " DMR",
     ifelse(zoom_windows$n_dmrs[i] > 1, "s", ""),
     ")"
    )

    out_i <- file.path(
      args$out_dir,
      paste0(args$assoc, "_dmr_zoom_cluster_", i, ".jpg")
    )

    gene_table_i <- file.path(
      args$out_dir,
      paste0(args$assoc, "_dmr_zoom_cluster_", i, ".refGene_genes.tsv")
    )

    zoom_plot_i <- make_dmr_zoom_plot(
      slk_df = slk,
      dmr_df = dmrs_to_highlight,
      chr = chr_i,
      chr_string = chr_string,
      region_start = zoom_windows$cluster_start[i],
      region_end = zoom_windows$cluster_end[i],
      out_file = out_i,
      gene_table_file = gene_table_i,
      padding = args$zoom_padding,
      point_jitter_bp = args$point_jitter_bp,
      show_dmr_midlines = args$zoom_midlines == "yes",
      title = title_i,
      genome_build = args$genome_build,
      gene_label_size = args$gene_label_size
    )

    if (!is.null(zoom_plot_i)) {
      zoom_plots[[length(zoom_plots) + 1]] <- zoom_plot_i
    }

    message("Saved zoomed DMR plot: ", out_i)
  }
}

# ================================================================
# Optional combined multi-panel figure
# ================================================================

if (args$make_combined == "yes") {
  if (args$make_zoom == "no" || length(zoom_plots) == 0) {
    message(
      "No regional zoom plots were available for the combined figure. ",
      "Creating a single-panel combined figure with the Manhattan plot only."
    )
  }

  combined_result <- make_combined_dmr_figure(
    manhattan_plot = dmr_manhattan,
    regional_plots = zoom_plots,
    width_cm = args$combined_width_cm,
    manhattan_height_cm = args$combined_manhattan_height_cm,
    region_row_height_cm = args$combined_region_height_cm,
    panel_label_size = args$panel_label_size
  )

  combined_prefix <- file.path(
    args$out_dir,
    paste0(args$assoc, "_dmr_combined")
  )

  combined_files <- save_plot_multi_format(
    plot = combined_result$plot,
    file_prefix = combined_prefix,
    formats = args$combined_formats,
    width = combined_result$width_cm,
    height = combined_result$height_cm,
    units = "cm",
    dpi = 300
  )

  message(
    "Saved combined DMR figure: ",
    paste(combined_files, collapse = ", ")
  )
}
