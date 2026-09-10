# Script for adding annotation data to EWAS results.
#
# Annotation sources
# ------------------
# Genes / promoters : Zhou lab Infinium manifest, GENCODE-annotated, fetched by
#                     rules/annotate.smk from zhou-lab/InfiniumAnnotationData
#                     and pinned by release tag. Columns: CpG_chrm, CpG_beg,
#                     CpG_end, probe_strand, probeID, genesUniq, geneNames,
#                     transcriptTypes, transcriptIDs, distToTSS.
# CpG islands       : called here from the UCSC cpgIslandExt BED cached by
#                     rules/dmr.smk. GENCODE v41 dropped the CGI and
#                     CGIposition columns that v36 supplied, and Zhou's CGI
#                     annotation was itself derived from this UCSC track, so
#                     recomputing keeps the EWAS and DMR outputs consistent.
# eQTM              : BIOS cis-eQTM table with HGNC-resolved GRCh38 symbols.
#
# SNP annotation is intentionally absent. Zhou's EPIC.hg38.commonsnp.tsv.gz has
# no successor in the v8.1 release, and SNP-affected probes are expected to be
# masked before the EWAS is run.
# Import libraries
suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(data.table)
})
# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for adding annotation data to ewas results")
parser$add_argument('--input-file', '-i',
                    required=TRUE,
                    help="Path to ewas results data")
parser$add_argument('--gene-anno',
                    required=TRUE,
                    help=paste("Path to the Zhou lab GENCODE-annotated Infinium manifest",
                               "(e.g. EPIC.hg38.manifest.gencode.v41.tsv.gz)"))
parser$add_argument('--cpg-islands',
                    required=TRUE,
                    help=paste("Path to the UCSC cpgIslandExt BED cache",
                               "(gzipped, headerless: chrom, start, end, name)"))
# Shores extend 2 kb from the island and shelves a further 2 kb (so 2-4 kb from
# the island). These are the standard definitions used across the methylation
# array literature and by minfi / sesame, so they are not exposed in config.yml;
# they remain arguments only so the boundaries can be varied for a sensitivity
# check without editing the script.
parser$add_argument('--shore-bp',
                    type="integer",
                    default=2000L,
                    help="Width in bp of the shore flanking each CpG island [default %(default)s]")
parser$add_argument('--shelf-bp',
                    type="integer",
                    default=2000L,
                    help="Width in bp of the shelf beyond each shore [default %(default)s]")
parser$add_argument('--eQTM-anno',
                    required=TRUE,
                    help="Path to eQTM annotation file")
parser$add_argument('--out-dir',
                    required=TRUE,
                    help="Path to output directory")
parser$add_argument('--stratified',
                    choices=c("yes", "no", "True", "False"),
                    default="no",
                    help="Results from a stratified analysis: yes or no")
parser$add_argument("--assoc", 
                    required=TRUE,
                    type="character", 
                    nargs=1, 
                    help="Association variable EWAS was performed with.")
parser$add_argument('--out-type', 
                    type="character",
                    choices=c(".csv", ".csv.gz"), 
                    nargs="?",                    
                    const=".csv",
                    default=".csv",  
                    help="Output file type: CSV or CSV.GZ")

# parse arguments
args <- parser$parse_args()
results <- args$input_file
out_dir <- args$out_dir
stratified <- args$stratified
assoc <- args$assoc
out_type <- args$out_type
gene_file <- args$gene_anno
island_file <- args$cpg_islands
shore_bp <- args$shore_bp
shelf_bp <- args$shelf_bp
eqtm_file <- args$eQTM_anno

#' Assign each probe to a CpG island, shore, shelf, or open sea.
#'
#' Reproduces the CGI and CGIposition columns that Zhou's GENCODE v36 manifest
#' carried and v41 dropped. Coordinates are 0-based half-open on both sides;
#' foverlaps treats intervals as closed, so ends are decremented by one before
#' the join. Where a probe falls in more than one zone, islands beat shores and
#' shores beat shelves, and ties within a tier go to the nearest island.
#'
#' @return data.table of cpgid, CGI, CGIposition for assigned probes only.
assign_cpg_islands <- function(island_bed, probes, shore_bp, shelf_bp) {
    isl <- fread(island_bed, header = FALSE, select = 1:3,
                 col.names = c("chrm", "istart", "iend"))
    isl <- isl[!is.na(chrm) & !is.na(istart) & !is.na(iend)]

    zone <- function(zstart, zend, label) {
        data.table(chrm = isl$chrm,
                   zstart = pmax(zstart, 0),
                   zend = pmax(zend, 0),
                   zone = label,
                   istart = isl$istart,
                   iend = isl$iend)
    }
    zones <- rbindlist(list(
        zone(isl$istart, isl$iend, "Island"),
        zone(isl$istart - shore_bp, isl$istart, "N_Shore"),
        zone(isl$iend, isl$iend + shore_bp, "S_Shore"),
        zone(isl$istart - shore_bp - shelf_bp, isl$istart - shore_bp, "N_Shelf"),
        zone(isl$iend + shore_bp, isl$iend + shore_bp + shelf_bp, "S_Shelf")
    ))
    # Drop zones collapsed to zero width by the clamp at the chromosome start.
    zones <- zones[zend > zstart]
    zones[, zend_c := zend - 1L]
    setkeyv(zones, c("chrm", "zstart", "zend_c"))

    q <- probes[!is.na(CpG_chrm) & !is.na(CpG_beg) & !is.na(CpG_end),
                .(cpgid, CpG_chrm, pstart = as.integer(CpG_beg),
                  pend = pmax(as.integer(CpG_end) - 1L, as.integer(CpG_beg)))]
    if (nrow(q) == 0L) {
        return(data.table(cpgid = character(), CGI = character(),
                          CGIposition = character()))
    }

    hits <- foverlaps(q, zones,
                      by.x = c("CpG_chrm", "pstart", "pend"),
                      type = "any", nomatch = NULL)
    if (nrow(hits) == 0L) {
        return(data.table(cpgid = character(), CGI = character(),
                          CGIposition = character()))
    }

    hits[, prio := fifelse(zone == "Island", 1L,
                    fifelse(zone %chin% c("N_Shore", "S_Shore"), 2L, 3L))]
    hits[, dist := fifelse(zone == "Island", 0L,
                    pmin(abs(istart - pstart), abs(pstart - iend)))]
    setorder(hits, cpgid, prio, dist)
    best <- unique(hits, by = "cpgid")   # setorder above makes this the best hit

    best[, .(cpgid,
             CGI = sprintf("CGI:%s:%d-%d", CpG_chrm, istart, iend),
             CGIposition = zone)]
}

# Read in EWAS summary statistics
ewas <- fread(results)

# Load gene / promoter annotation
gene_anno <- fread(gene_file)
if (!"probeID" %in% names(gene_anno)) {
    # Zhou's own column documentation calls this Probe_ID; the shipped v41
    # files use probeID. Accept either so a manifest refresh cannot break us.
    if ("Probe_ID" %in% names(gene_anno)) {
        setnames(gene_anno, "Probe_ID", "probeID")
    } else {
        stop("No probeID / Probe_ID column found in ", gene_file)
    }
}
setnames(gene_anno, "probeID", "cpgid")
setDT(gene_anno)

# Call islands, shores and shelves from the UCSC track
islands <- assign_cpg_islands(island_file, gene_anno, shore_bp, shelf_bp)
message(sprintf("CpG-island assignment: %d of %d probes in island/shore/shelf",
                nrow(islands), nrow(gene_anno)))

eqtm <- fread(eqtm_file) %>%
  dplyr::select(SNPName, HGNCName_GRCh38) %>%
  filter(!is.na(HGNCName_GRCh38)) %>%
  rename(cpgid = SNPName, BIOS_eQTM_gene = HGNCName_GRCh38) %>%
  group_by(cpgid) %>%
  summarize(BIOS_eQTM_genes = paste(unique(BIOS_eQTM_gene), collapse = ";")) %>%
  ungroup()

annotation <- left_join(gene_anno, islands, by = "cpgid") %>%
              left_join(eqtm, by = "cpgid")
rm(gene_anno, islands)

# Guard against a platform mismatch between the EWAS results and the manifest.
# Every Zhou platform manifest has the same columns, so a wrong
# annotation.array_platform does not fail loudly -- it just joins nothing.
# EPICv2 and MSA probe IDs also carry a design suffix (cg00000029_TC21), so a
# matrix built with those suffixes stripped will not match either.
result_ids <- if (stratified == "no" || stratified == "False") {
    ewas$cpgid
} else {
    ewas$MarkerName
}
matched <- sum(result_ids %in% annotation$cpgid)
match_rate <- if (length(result_ids) > 0) matched / length(result_ids) else 0
message(sprintf("Annotation match: %d of %d result rows (%.1f%%)",
                matched, length(result_ids), 100 * match_rate))
if (match_rate < 0.50) {
    stop(sprintf(paste0(
        "Only %.1f%% of result CpGs matched the annotation manifest.\n",
        "  Check that annotation.array_platform in config.yml matches the array\n",
        "  the M-value matrix came from, and that probe-ID suffixes agree\n",
        "  (EPICv2 and MSA use IDs such as cg00000029_TC21)."), 100 * match_rate))
} else if (match_rate < 0.95) {
    warning(sprintf(
        "Only %.1f%% of result CpGs matched the annotation manifest.", 100 * match_rate))
}

if(stratified=="no" | stratified == "False"){
    ewas <- left_join(ewas, annotation, by = "cpgid")
    ewas <- ewas[order(ewas$bacon.pval),]
} else{
        ewas <- left_join(ewas, annotation, by = c("MarkerName"="cpgid"))
        ewas <- ewas %>% dplyr::select(-Allele1, -Allele2)
        ewas$"P-value" <- as.numeric(ewas$"P-value")
        ewas <- ewas[order(ewas$"P-value"),]
}
file_name <- paste0(out_dir,"/", assoc, "_ewas_annotated_results", out_type)
fwrite(ewas, file = file_name)
