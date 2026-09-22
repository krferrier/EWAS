suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

# Genomic inflation factor (lambda): the median observed chi-square statistic
# over its expected value under the null. Identical to QCEWAS::P_lambda
# (QCEWAS 1.2-3, copied from its source), defined here so the environment does
# not depend on QCEWAS: conda-forge's r-qcewas has no build newer than R 4.3,
# while bioconductor-bacon >= 1.32 requires R >= 4.4, so the two cannot be
# installed together. This was the only QCEWAS function the workflow used.
p_lambda <- function(p) {
    p <- p[!is.na(p)]
    if (length(p) < 2L) {
        stop("'p' does not contain sufficient non-missing values to calculate lambda")
    }
    median(qchisq(p, df = 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}


# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for plotting ewas results")
parser$add_argument('--input-file', '-i',
                    required=TRUE,
                    help="Path to annotated ewas results data")
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


# parse arguments
args <- parser$parse_args()
results <- args$input_file
out_dir <- args$out_dir
stratified <- args$stratified
assoc <- args$assoc

# Read in annotated EWAS results
ewas <- fread(results)

# Chromosomes drawn on the Manhattan plot, in axis order. CHR is the position in
# this vector (X = 23, Y = 24), found by matching rather than as.numeric(), which
# turned X and Y into NA with a coercion warning and drew them as a block
# labelled "NA" after chr22. Anything else -- chrM and unplaced or alt contigs
# such as chrUn_KI270748v1 -- also gets NA and is left off the Manhattan plot
# only; those CpGs stay in the QQ plot and in lambda, since they were tested.
CHROMOSOMES <- c(as.character(1:22), "X", "Y")

# Wrangle data for plotting
if (stratified == "yes" | stratified == "True"){
  ewas <- ewas %>% 
    dplyr::select(MarkerName, CpG_chrm, CpG_beg, "P-value") %>% 
    rename("Pvalue" = "P-value") %>% 
    mutate("CHR" = match(sub("^chr", "", CpG_chrm), CHROMOSOMES),
           "MAPINFO" = CpG_beg)
} else{
  ewas <- ewas %>% 
    dplyr::select(cpgid, CpG_chrm, CpG_beg, bacon.pval)%>% 
    rename(Pvalue = bacon.pval) %>% 
    mutate(CHR = match(sub("^chr", "", CpG_chrm), CHROMOSOMES),
           MAPINFO = CpG_beg)
}

# CpGs on the Manhattan plot: a standard chromosome and a position.
manh <- ewas %>% filter(!is.na(CHR), !is.na(MAPINFO))
n_off <- nrow(ewas) - nrow(manh)
if (n_off > 0) {
  off <- ewas %>% filter(is.na(CHR) | is.na(MAPINFO)) %>%
    count(chrom = ifelse(is.na(CpG_chrm) | CpG_chrm == "", "<no position>", CpG_chrm))
  message(sprintf("Manhattan plot omits %d CpG(s) not on chr1-22, X or Y: %s",
                  n_off, paste(sprintf("%s %d", off$chrom, off$n), collapse = ", ")))
}

#################################################################
#                  MANHATTAN & QQ PLOT                          #
#################################################################


# Calculate the cumulative position of each chromosome
chr.pos <- manh %>% 
  group_by(CHR) %>% 
  summarize(chr_len=max(as.numeric(MAPINFO))) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>% 
  dplyr::select(-chr_len) 


# Add cumulative position data to results
manh <- left_join(manh, chr.pos, by = "CHR") %>% 
  arrange(CHR, as.numeric(MAPINFO)) %>% 
  mutate(POS = as.numeric(MAPINFO)+tot)


x_axis <- manh %>% group_by(CHR) %>% 
  summarise(center = (max(POS) + min(POS))/2) 

manh %>% 
  filter(-log10(Pvalue)>1) %>% 
  ggplot(aes(x=POS, y=-log10(Pvalue))) +
  geom_point(aes(color=as.factor(CHR)), alpha = 0.8, size = 2) +
  scale_color_manual(values= rep(c("steelblue1", "steelblue4"), 
                                 length.out = length(CHROMOSOMES))) + 
  scale_x_continuous(label = CHROMOSOMES[x_axis$CHR], breaks = x_axis$center,
                     guide = guide_axis(check.overlap = T), expand=c(0,0)) +
  # Lower limit only. A fixed upper limit (this was 28) silently dropped every
  # CpG with p < 1e-28 -- the strongest hits -- with nothing but ggplot's
  # "Removed N rows ... outside the scale range" warning to show for it. The
  # small top expansion keeps the highest point from being cut by the frame.
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                     limits = c(1, NA)) +
  geom_hline(yintercept = -log10(0.05/nrow(ewas)),
             linetype = 'solid',
             color = "red",
             linewidth = 0.5) +
  theme_bw(base_size = 16) +
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(y= expression(-log[10]("P-value")),
       x= "Chromosome") -> manh.plot


StatQQplot <- ggproto("StatQQplot", Stat,
                      default_aes = aes(y = stat(observed), x = stat(expected)),
                      
                      required_aes = c("observed"),
                      
                      compute_group = function(data, scales, dparams = list(),
                                               na.rm = FALSE) {
                        
                        observed <- data$observed#[!is.na(data$x)]
                        N <- length(observed)
                        
                        ## expected
                        expected <- sort(-log10((1:N)/N-1/(2*N)))
                        observed <- sort(-log10(observed))
                        data.frame(observed, expected)
                        
                      }
)

stat_qqplot <- function(mapping = NULL, data = NULL, geom = "point",
                        position = "identity", na.rm = FALSE, show.legend = NA, 
                        inherit.aes = TRUE, ...) {
  layer(
    stat = StatQQplot, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

lambda <- p_lambda(ewas$Pvalue)
lambda_label <- paste0("lambda==", round(lambda, digits = 2))

ggplot(ewas, aes(observed= Pvalue)) +
  stat_qqplot() +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw(base_size = 16) +
  annotate("text", x=1, y=25, label = lambda_label, parse=T, size = 6)+
  labs(y= expression(Observed ~ -log[10]("P-value")),
      x = expression(Expected ~ -log[10]("P-value"))) -> qq.plot

left.panel <- plot_grid(NULL, qq.plot, labels= c("B", ""), label_size = 22, 
                        ncol=1, rel_heights = c(1,5))
full.plot <- plot_grid(manh.plot, left.panel, labels = c("A", ""), ncol = 2,
          rel_widths = c(2.5,1), label_size = 22)

filename <- paste0(out_dir, "/", assoc, "_ewas_manhattan_qq_plots.jpg")
ggsave(filename,
      plot = full.plot,
      width = 32,
      height = 12,
      units = "cm")
