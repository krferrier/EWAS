# Make a bed file of the EWAS results

suppressPackageStartupMessages({
    library(R.utils)
    library(argparse)
    library(data.table)
    library(tibble)
    library(dplyr)
})


# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for making BED file from EWAS results")
parser$add_argument("--results",
                    required=TRUE,
                    help="Path or URL to EWAS results file.")
parser$add_argument('--out-dir', 
                    type="character",
                    nargs="?",                    
                    const="~/", 
                    default="~/",  
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
           

# parse arguments
args <- parser$parse_args()
results <- args$results
out_dir <- args$out_dir
stratified <- args$stratified
assoc <- args$assoc

# Read in EWAS results files
res <- fread(results)

# Wrangle results into BED format (chr, start, end, pvalue)
if(stratified=="no" | stratified == "False"){
    res <- res  %>% 
            dplyr::select(CpG_chrm,CpG_beg,CpG_end, bacon.pval, cpgid) %>%
            dplyr::rename("#chrom" = "CpG_chrm",
                          "start" = "CpG_beg",
                          "end" = "CpG_end",
                          "pvals" = "bacon.pval")  
} else{
    res <- res  %>% 
        dplyr::select(CpG_chrm,CpG_beg,CpG_end, "P-value", MarkerName)  %>%
        dplyr::rename("#chrom" = "CpG_chrm",
                      "start" = "CpG_beg",
                      "end" = "CpG_end",
                      "pvals" = "P-value")
                      
}

# comb-p needs every CpG to have a position, and the file sorted by chromosome
# in plain character (C locale) order -- it stops with "chromosomes must be
# sorted as characters" otherwise. dplyr::arrange() is not safe for this: before
# dplyr 1.1 it sorted in the system locale, which puts "chr1_KI270711v1_random"
# before "chr10", and the DMR environment carries such a dplyr. order(method =
# "radix") always sorts in C order, whatever the locale or package versions.
# Plain data.frame: data.table's [ ] treats i and j differently from base R.
res <- as.data.frame(res)
chrom_col <- names(res)[1]
no_pos <- is.na(res[[chrom_col]]) | res[[chrom_col]] %in% c("", "NA") | is.na(res$start)
if (any(no_pos)) {
    message(sprintf("make_bed: dropped %d CpG(s) with no genomic position", sum(no_pos)))
    res <- res[!no_pos, ]
}
res <- res[order(res[[chrom_col]], res$start, method = "radix"), ]

# Export BED file
file_name <- paste0(out_dir, "/", assoc, "_ewas_annotated_results.bed")
write.table(res, file = file_name, append = FALSE, sep = "\t",
             row.names = FALSE, col.names = TRUE, quote = FALSE)
