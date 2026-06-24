suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(dplyr)
})

parser <- argparse::ArgumentParser(description = paste("Add current HGNC annotations to BIOS cis-eQTM results using",
        "the Ensembl gene IDs in the BIOS ProbeName column."))
parser$add_argument('--bios-eqtm',
                    required = TRUE,
                    help = "Raw BIOS cis-eQTM table")
parser$add_argument("--hgnc-complete-set",
                    required = TRUE,
                    help = "HGNC complete set file")
parser$add_argument("--output",
                    required = TRUE,
                    help = "Output annotated BIOS cis-eQTM table")

args <- parser$parse_args()

# Read raw source tables.
bios <- fread(args$bios_eqtm, colClasses = "character", quote = "")

hgnc <- fread(args$hgnc_complete_set, colClasses = "character", quote = "")

# Preserve the original BIOS symbol under an explicit name.
setnames(bios, old = "HGNCName", new = "HGNCName_BIOS_source")

# The BIOS ProbeName field contains Ensembl gene IDs.
# Remove an optional Ensembl version suffix before matching.
bios$ensembl_gene_id <- sub("\\..*$", "", bios$ProbeName)

# Retain approved HGNC records with Ensembl gene identifiers.
hgnc <- hgnc %>%
    filter(status == "Approved", !is.na(ensembl_gene_id), ensembl_gene_id != "") %>%
    transmute(
        ensembl_gene_id = sub("\\..*$", "", ensembl_gene_id),
        HGNC_ID = hgnc_id,
        HGNCName_GRCh38 = symbol
    )
hgnc <- unique(hgnc)

# Resolve rare one-to-many Ensembl ID -> approved HGNC mappings
# deterministically by retaining the lowest numeric HGNC ID.
hgnc <- hgnc %>%
    mutate(hgnc_id_num = as.integer(sub("^HGNC:", "", HGNC_ID))) %>%
    arrange(ensembl_gene_id, hgnc_id_num, HGNCName_GRCh38) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    select(-hgnc_id_num)

# Match preserves BIOS row order.
hgnc_index <- match(bios$ensembl_gene_id, hgnc$ensembl_gene_id)

bios$HGNC_ID <- hgnc$HGNC_ID[hgnc_index]
bios$HGNCName_GRCh38 <- hgnc$HGNCName_GRCh38[hgnc_index]

bios <- bios %>% select(-ensembl_gene_id)


fwrite(bios, args$output, sep = "\t", quote = FALSE, na = "")