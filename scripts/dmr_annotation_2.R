suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(data.table)
})
# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for adding annotation data to DMR results")
parser$add_argument('--dmr-regions',
                    required=TRUE,
                    help="Path to DMR regions-p results data")
parser$add_argument('--ewas-bed',
                    required=TRUE,
                    help="Path to EWAS results bed used as input to dmr analysis")
parser$add_argument('--hgnc',
                    required=TRUE,
                    help="Path to HGNC annotation data")
parser$add_argument('--refgene',
                    required=TRUE,
                    help="Path to UCSC RefGene annotation data")
parser$add_argument('--cpgIslandExt',
                    required=TRUE,
                    help="Path to UCSC cpgIslandExt annotation data")
parser$add_argument('--out-dir',
                    required=TRUE,
                    help="Path to output directory")
parser$add_argument("--assoc",
                    required=TRUE,
                    type="character",
                    nargs=1,
                    help="Association variable EWAS was performed with.")


#parse arguments
args <- parser$parse_args()
dmr_regions <- args$dmr_regions
ewas_bed <- args$ewas_bed
hgnc <- args$hgnc
refgene <- args$refgene
cpgIslandExt <- args$cpgIslandExt
out_dir <- args$out_dir
assoc <- args$assoc


# dmr_regions <- "results_test/dmr/BMI_ewas.regions-p.bed.gz"
# ewas_bed <- "results_test/BMI_ewas_annotated_results.bed"
# hgnc <- "resources/ucsc/hg38/latest/hgnc.approved.bed.gz"
# refgene <- "resources/ucsc/hg38/latest/refGene.bed.gz"
# cpgIslandExt <- "resources/ucsc/hg38/latest/cpgIslandExt.bed.gz"
# out_dir <- "results_test/dmr"
# assoc <- "BMI"


# Define outfile name
out_file  <- paste0(out_dir, "/", assoc, "_dmr_annotated_results.tsv")

# -------------------------
# Read input files
# -------------------------

# Read DMR annotation bed file output from comb-p
dmr <- fread(dmr_regions, sep = "\t", header = TRUE)
setnames(dmr, 1, "chrom")
dmr[, dmr_row_id := .I]

if (!"region_id" %in% names(dmr)) {
  dmr[, region_id := paste0("DMR_", sprintf("%06d", dmr_row_id))]
}

# Read EWAS annotation bed file used as input to comb-p
cpg <- fread(ewas_bed, sep = "\t", header = TRUE)
setnames(cpg, 1, "chrom")

# Read HGNC annotation data
hgnc_anno <- fread(hgnc, sep = "\t", header = F)
colnames(hgnc_anno) <- c("chrom", "start", "end", "HGNC_gene", "HGNC_symbol", "HGNC_gene_name",
                         "HGNC_locus_group", "HGNC_locus_type")

# Read in refgene data
refgene_anno <- fread(refgene, sep = "\t", header = F)
colnames(refgene_anno) <- c("chrom", "start", "end", "refGene_gene", "refGene_ID", "refGene_strand")

# Read in cpgIslandExt data
cpgIslandExt_anno <- fread(cpgIslandExt, sep = "\t", header = F)
colnames(cpgIslandExt_anno) <- c("chrom", "start", "end", "cpgIsland_ID")

# -------------------------
# EWAS CpG summary
# -------------------------
# Some files use 'MarkerName', others use 'cpgid'; otherwise use the 5th column name.
if ("MarkerName" %in% names(cpg)) {
  id_col <- "MarkerName"
} else if ("cpgid" %in% names(cpg)) {
  id_col <- "cpgid"
} else {
  id_col <- names(cpg)[5]
}
setnames(cpg, id_col, "cpg_id")

# Set the threshold for EWAS significance
threshold <- 0.05 / nrow(cpg)  # Bonferroni correction for epigenome‑wide significance

# Perform overlap join of EWAS-CpG information
setkey(dmr, chrom, start, end)
setkey(cpg, chrom, start, end)
ov <- foverlaps(cpg, dmr, type = "within", nomatch = 0L)

if (nrow(ov) > 0L) {
  ewas_summary <- ov[, .(
    CpGs              = paste(unique(cpg_id), collapse = ";"),
    min_cpg_pvalue   = min(pvals, na.rm = TRUE),
    num_EWAS_sig_CpGs = sum(pvals < threshold, na.rm = TRUE),
    detected_by_EWAS = ifelse(sum(pvals < threshold, na.rm = TRUE) > 0L, "Yes", "No")
  ), by = dmr_row_id]
} else {
  ewas_summary <- data.table(
    dmr_row_id = integer(),
    CpGs = character(),
    min_cpg_pvalue = numeric(),
    num_EWAS_sig_CpGs = integer(),
    detected_by_EWAS = character()
  )
}

# -------------------------
# Helper functions
# -------------------------
signed_interval_distance <- function(region_start, region_end, ann_start, ann_end) {
  # 0 = overlap
  # negative = annotation is upstream/left of region
  # positive = annotation is downstream/right of region
  fifelse(
    ann_end > region_start & ann_start < region_end,
    0L,
    fifelse(
      ann_end <= region_start,
      ann_end - region_start,
      ann_start - region_end
    )
  )
}

cpg_feature_from_distance <- function(distance_string, shore_dist = 3000L) {
  if (is.na(distance_string) || distance_string %in% c("", "NA")) return("")

  tokens <- unlist(strsplit(as.character(distance_string), ";", fixed = TRUE))
  tokens <- trimws(tokens)
  tokens <- tokens[tokens != ""]

  if (length(tokens) == 0L) return("")

  d <- suppressWarnings(abs(as.integer(tokens)))
  d <- d[!is.na(d)]

  if (length(d) == 0L) return("")

  features <- ifelse(d == 0L, "island",
                     ifelse(d <= shore_dist, "shore", ""))

  features <- sort(unique(features[features != ""]))
  if (length(features) == 0L) return("")
  paste(features, collapse = ";")
}

make_nearest_summary <- function(dmr_dt, anno_dt, value_cols, distance_col) {
  anno_by_chr <- split(anno_dt, anno_dt$chrom, drop = TRUE)

  res <- vector("list", nrow(dmr_dt))

  for (i in seq_len(nrow(dmr_dt))) {
    chr <- dmr_dt$chrom[i]
    anns <- anno_by_chr[[chr]]

    row_out <- list(dmr_row_id = dmr_dt$dmr_row_id[i])

    if (is.null(anns) || nrow(anns) == 0L) {
      for (j in seq_along(value_cols)) row_out[[value_cols[j]]] <- "NA"
      row_out[[distance_col]] <- "NA"
      res[[i]] <- as.data.table(row_out)
      next
    }

    dist <- signed_interval_distance(
      region_start = dmr_dt$start[i],
      region_end   = dmr_dt$end[i],
      ann_start    = anns$start,
      ann_end      = anns$end
    )

    best_abs <- min(abs(dist), na.rm = TRUE)
    keep <- abs(dist) == best_abs
    hits <- anns[keep]
    dist_string <- paste(sort(unique(dist[keep])), collapse = ";")

    for (j in seq_along(value_cols)) {
      row_out[[value_cols[j]]] <- paste(hits[[value_cols[j]]], collapse = ";")
    }

    row_out[[distance_col]] <- dist_string
    res[[i]] <- as.data.table(row_out)
  }

  rbindlist(res, fill = TRUE)
}

# -------------------------
# HGNC annotation summary
# -------------------------

hgnc_nearest <- make_nearest_summary(
  dmr_dt = dmr,
  anno_dt = hgnc_anno,
  value_cols = c("HGNC_gene", "HGNC_symbol", "HGNC_gene_name",
                 "HGNC_locus_group", "HGNC_locus_type"),
  distance_col = "hgnc_distance_bp"
)

# -------------------------
# refGene annotation summary
# -------------------------

refgene_nearest <- make_nearest_summary(
  dmr_dt = dmr,
  anno_dt = refgene_anno,
  value_cols = c("refGene_gene"),
  distance_col = "refgene_distance_bp"
)

# -------------------------
# cpgIslandExt nearest annotation
# -------------------------

cpgIsland_nearest <- make_nearest_summary(
  dmr_dt = dmr,
  anno_dt = cpgIslandExt_anno,
  value_cols = c("cpgIsland_ID"),
  distance_col = "cpgIslandExt_distance"
)

cpgIsland_nearest[,
  cpgIslandExt_feature := vapply(
    cpgIslandExt_distance,
    cpg_feature_from_distance,
    character(1)
  )
]

# -------------------------
# Combine all annotations
# -------------------------

annotation_tables <- list(
  ewas_summary,
  hgnc_nearest,
  refgene_nearest,
  cpgIsland_nearest
)

dmr_summary <- Reduce(
  function(x, y) merge(x, y, by = "dmr_row_id", all.x = TRUE, sort = FALSE),
  c(list(dmr), annotation_tables)
)

# Remove internal helper columns before writing.
internal_cols <- c("dmr_row_id", "region_id", "HGNC_symbol", "start_i", "end_i")
dmr_summary[, (intersect(internal_cols, names(dmr_summary))) := NULL]
dmr_summary <- dmr_summary %>% filter(n_probes >= 2)

fwrite(dmr_summary, file = out_file, sep = "\t", quote = FALSE, na = "NA")

cat("Output written to:", out_file, "\n")
