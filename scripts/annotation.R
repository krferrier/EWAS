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
parser$add_argument('--kycg-dir',
                    default=NULL,
                    help=paste("Directory of normalised KYCG .cm feature sets for this",
                               "platform. Omit to skip the feature columns."))
parser$add_argument('--probe-order',
                    default=NULL,
                    help=paste("Path to <platform>.ordering.tsv.gz. Required with",
                               "--kycg-dir: the .cm files carry no probe IDs and are",
                               "row-aligned to this file."))
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
kycg_dir <- args$kycg_dir
order_file <- args$probe_order

# KYCG set prefix -> output column name. Must stay in step with
# ConfigWizard.KYCG_FEATURE_SETS in helper_fxns.py, which decides what
# rules/annotate.smk downloads. A set listed here but absent from the cache is
# skipped, so the two lists going briefly out of step is not fatal.
KYCG_FEATURE_SETS <- list(
    ChromHMM       = "chromHMM_state",
    PMD            = "PMD",
    ABCompartment  = "AB_compartment",
    rmsk1          = "repeat_class",
    ImprintingDMR  = "imprinting_DMR",
    CTCFbind       = "CTCF_binding",
    Blacklist      = "ENCODE_blacklist"
)

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

#' Read the KYCG feature sets into one probe-keyed table.
#'
#' The .cm files carry no probe IDs -- they are row-aligned to
#' <platform>.ordering.tsv.gz -- and come in two packings, which `yame unpack -a`
#' renders uniformly:
#'   * a single column of state labels (ChromHMM, ImprintingDMR, Blacklist,
#'     CTCFbind), used as the column value directly; and
#'   * one 0/1 column per record (PMD, ABCompartment, rmsk1), collapsed to a
#'     single label here. These are near-exclusive but not strictly so -- 326
#'     EPIC probes carry more than one rmsk1 class -- so multiples are joined
#'     with ";" rather than silently dropped.
#'
#' Any set a platform does not publish is simply absent and is skipped.
#'
#' @return data.table keyed on cpgid, or NULL when nothing could be read.
read_kycg_features <- function(kycg_dir, order_file, feature_names) {
    if (is.null(kycg_dir) || is.null(order_file)) return(NULL)
    if (!dir.exists(kycg_dir) || !file.exists(order_file)) return(NULL)
    if (Sys.which("yame") == "") {
        warning("yame not found on PATH; skipping KYCG feature columns.")
        return(NULL)
    }

    probe_order <- fread(order_file, select = 1L)
    setnames(probe_order, 1L, "cpgid")
    out <- data.table(cpgid = probe_order$cpgid)

    for (prefix in names(feature_names)) {
        cm <- file.path(kycg_dir, paste0(prefix, ".cm"))
        if (!file.exists(cm)) next
        column <- feature_names[[prefix]]

        values <- tryCatch({
            tbl <- fread(cmd = paste("yame unpack -a", shQuote(cm)),
                         header = FALSE, sep = "\t", colClasses = "character")
            if (nrow(tbl) != nrow(out)) {
                stop(sprintf("%d rows, expected %d from the ordering file",
                             nrow(tbl), nrow(out)))
            }
            if (ncol(tbl) == 1L) {
                # already a state label per probe
                v <- tbl[[1]]
                v[v %in% c("NA", "", ".")] <- NA_character_
                v
            } else {
                # one 0/1 column per record; recover names from `yame info`
                info <- fread(cmd = paste("yame info", shQuote(cm)), sep = "\t")
                labels <- as.character(info[[2]])
                if (length(labels) != ncol(tbl)) {
                    stop(sprintf("%d record names for %d columns",
                                 length(labels), ncol(tbl)))
                }
                # Record names live in the .cm.idx sidecar. Without it yame
                # falls back to 1..N, which would put meaningless integers in
                # the results, so refuse the set rather than degrade quietly.
                if (all(grepl("^[0-9]+$", labels))) {
                    stop(sprintf(paste("record names are bare indices -- the %s.cm.idx",
                                       "sidecar is missing from the cache"), prefix))
                }
                m <- as.matrix(tbl) == "1"
                apply(m, 1L, function(hits) {
                    if (!any(hits)) NA_character_ else paste(labels[hits], collapse = ";")
                })
            }
        }, error = function(e) {
            warning(sprintf("KYCG set '%s' skipped: %s", prefix, conditionMessage(e)))
            NULL
        })

        if (!is.null(values)) {
            out[[column]] <- values
            message(sprintf("KYCG %-16s -> %-16s %d of %d probes annotated",
                            prefix, column, sum(!is.na(values)), length(values)))
        }
    }

    if (ncol(out) == 1L) NULL else out
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

# KYCG feature columns (chromatin state, PMD, A/B compartment, repeat class,
# imprinting DMR, CTCF binding, ENCODE blacklist). Skipped without error when
# the cache is absent or the platform does not publish the sets.
kycg <- read_kycg_features(kycg_dir, order_file, KYCG_FEATURE_SETS)
if (!is.null(kycg)) {
    annotation <- left_join(annotation, kycg, by = "cpgid")
    rm(kycg)
} else {
    message("No KYCG feature sets loaded; annotated results will omit those columns.")
}

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
