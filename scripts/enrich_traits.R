# Trait enrichment of the significant CpG set against the EWAS Atlas.
#
# Replaces looking hits up in the Atlas web interface by hand: for every trait
# with enough probes on this array, test whether the significant set overlaps
# its reported CpGs more than expected.
#
# Two things to keep in mind when reading the output. The Atlas is a catalogue
# of published associations, so its coverage reflects what has been studied --
# smoking, ageing and sex dominate it, and an overlap with those traits partly
# reflects publication volume. And each trait's probe list is restricted here
# to probes present in this analysis's background, so counts will not match the
# Atlas website, which reports against all arrays at once.

suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
})
source(file.path(dirname(sub("--file=", "",
       grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])),
       "enrichment_common.R"))

parser <- argparse::ArgumentParser(
    description = "EWAS Atlas trait enrichment for the significant CpG set")
parser$add_argument('--input-file', required = TRUE)
parser$add_argument('--ewas-atlas', required = TRUE,
                    help = "Cached EWAS_Atlas_associations.tsv")
parser$add_argument('--stratified', required = TRUE)
parser$add_argument('--significance', default = "fdr")
parser$add_argument('--threshold', type = "double", default = 0.05)
parser$add_argument('--min-set-size', type = "integer", default = 20L)
parser$add_argument('--max-overlap-listed', type = "integer", default = 25L,
                    help = "Cap on probes listed per trait [default %(default)s]")
parser$add_argument('--output', required = TRUE)
args <- parser$parse_args()

OUT_COLS <- c("trait", "n_universe", "n_significant", "n_trait_probes",
              "n_overlap", "expected", "fold_enrichment", "neg_log10_p",
              "p_value", "fdr",
              "n_studies", "pmids", "overlapping_probes")

sets <- read_cpg_sets(args$input_file, args$stratified,
                      args$significance, args$threshold)
if (length(sets$foreground) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       "no CpGs passed the significance threshold")
    quit(save = "no", status = 0)
}

# The Atlas table is not valid UTF-8 (some trait names carry Latin-1 bytes),
# so read it as Latin-1 and convert rather than letting fread mangle it.
atlas <- fread(args$ewas_atlas, sep = "\t", quote = "",
               select = c("probe_ID", "trait", "study_ID", "PMID"),
               encoding = "Latin-1", showProgress = FALSE)
atlas[, trait := trimws(enc2utf8(as.character(trait)))]
atlas <- atlas[!is.na(probe_ID) & nzchar(trait)]
message(sprintf("EWAS Atlas: %d associations, %d probes, %d traits",
                nrow(atlas), uniqueN(atlas$probe_ID), uniqueN(atlas$trait)))

# Restrict the Atlas to probes this analysis actually tested, so the universe
# is the same one used for the other enrichment rules.
bg <- unique(sets$background)
fg <- unique(sets$foreground)
atlas <- atlas[probe_ID %in% bg]
if (nrow(atlas) == 0L) {
    write_empty_result(args$output, OUT_COLS, paste(
        "no Atlas probe is present in the tested background -- the Atlas is",
        "keyed on bare cg identifiers, so an EPICv2/MSA run with suffixed",
        "probe IDs will not match"))
    quit(save = "no", status = 0)
}

n_univ <- length(bg)
n_sig  <- length(fg)

per_trait <- atlas[, .(
    n_trait_probes = uniqueN(probe_ID),
    n_studies      = uniqueN(study_ID),
    pmids          = paste(unique(stats::na.omit(PMID))[1:min(10L, uniqueN(PMID))],
                           collapse = ";"),
    overlap_probes = list(intersect(unique(probe_ID), fg))
), by = trait]

per_trait[, n_overlap := lengths(overlap_probes)]
per_trait <- per_trait[n_trait_probes >= args$min_set_size]
if (nrow(per_trait) == 0L) {
    write_empty_result(args$output, OUT_COLS,
                       sprintf("no trait reached min_set_size=%d", args$min_set_size))
    quit(save = "no", status = 0)
}

per_trait[, n_universe := n_univ]
per_trait[, n_significant := n_sig]
per_trait[, expected := n_sig * (n_trait_probes / n_univ)]
per_trait[, fold_enrichment := fifelse(expected > 0, n_overlap / expected, NA_real_)]
.ht <- hyper_test(n_overlap = per_trait$n_overlap, n_query = per_trait$n_sig,
                  n_mask = per_trait$n_trait_probes, n_univ = per_trait$n_univ)
per_trait[, p_value := .ht$p_value]
# Exact magnitude even where p_value has underflowed to 0; see hyper_test().
per_trait[, neg_log10_p := .ht$neg_log10_p]
per_trait[, fdr := p.adjust(p_value, method = "BH")]
per_trait[, overlapping_probes := vapply(overlap_probes, function(p) {
    if (length(p) == 0L) return(NA_character_)
    shown <- p[seq_len(min(length(p), args$max_overlap_listed))]
    paste(c(shown, if (length(p) > length(shown)) sprintf("...(+%d)", length(p) - length(shown))),
          collapse = ";")
}, character(1))]
per_trait[, overlap_probes := NULL]

setorder(per_trait, p_value, -fold_enrichment)
fwrite(per_trait[, ..OUT_COLS], args$output, sep = "\t")
message(sprintf("Tested %d traits; %d at FDR < 0.05",
                nrow(per_trait), sum(per_trait$fdr < 0.05, na.rm = TRUE)))
