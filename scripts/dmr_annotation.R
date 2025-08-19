library(data.table)

# Read input files
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript add_cpgs_and_ewas.R <dmr_file> <cpg_file> [output_file]")
}
dmr_file  <- args[1]
cpg_file  <- args[2]
out_file  <- sub("(\\.bed)?$", ".with_cpgs.bed", dmr_file)


# Read DMR annotation bed file output from comb-p
dmr <- fread(dmr_file, sep = "\t", header = TRUE)
setnames(dmr, 1, "chrom")
dmr[, start := as.integer(start)]
dmr[, end   := as.integer(end)]
dmr[, region_id := .I]

# Read EWAS annotation bed file used as input to comb-p
cpg <- fread(cpg_file, sep = "\t", header = TRUE)
setnames(cpg, 1, "chrom")
cpg[, start := as.integer(start)]
cpg[, end   := as.integer(end)]
cpg[, pvals := as.numeric(pvals)]

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

# Perform overlap join
setkey(dmr, chrom, start, end)
setkey(cpg, chrom, start, end)
ov <- foverlaps(cpg, dmr, type = "within", nomatch = 0L)

# Aggregate CpGs and EWAS summary by region
summary_dt <- ov[, .(
  CpGs            = paste(cpg_id, collapse = ","),
  min_cpg_pvalue     = min(pvals),
  num_EWAS_sig_CpGs  = sum(pvals < threshold),
  detected_by_EWAS   = ifelse(sum(pvals < threshold) > 0L, "Yes", "No")
), by = region_id]

# Merge back to DMR data and write output
dmr_summary <- merge(dmr, summary_dt, by = "region_id", all.x = TRUE)
setorder(dmr_summary, region_id)
dmr_summary[, region_id := NULL]
fwrite(dmr_summary, file = out_file, sep = "\t")

cat("Output written to:", out_file, "\n")