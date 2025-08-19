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
                    choices=c("yes", "no"), 
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
if(stratified=="no"){
    res <- res  %>% 
            dplyr::select(CpG_chrm,CpG_beg,CpG_end, bacon.pval, cpgid) %>%
            arrange(CpG_chrm, CpG_beg)  %>% 
            dplyr::rename("#chrom" = "CpG_chrm",
                          "start" = "CpG_beg",
                          "end" = "CpG_end",
                          "pvals" = "bacon.pval")  
} else{
    res <- res  %>% 
        dplyr::select(CpG_chrm,CpG_beg,CpG_end, "P-value", MarkerName)  %>%
        arrange(CpG_chrm, CpG_beg)  %>% 
        dplyr::rename("#chrom" = "CpG_chrm",
                      "start" = "CpG_beg",
                      "end" = "CpG_end",
                      "pvals" = "P-value")
                      
}

# Export BED file
file_name <- paste0(out_dir, assoc, "_ewas_annotated_results.bed")
write.table(res, file = file_name, append = FALSE, sep = "\t",
             row.names = FALSE, col.names = TRUE, quote = FALSE)
