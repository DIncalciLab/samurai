#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(dplyr, quietly = TRUE)
}
)

parser <- arg_parser("Adjust logR.", hide.opts = TRUE)

parser <- add_argument(parser, "--seg",
                       help = "Segmentation file from which extract information.",
                       nargs = Inf)

parser <- add_argument(parser, "--ploidy",
                       help = "Ploidy Summary file from which extract information.",
                       nargs = Inf)

parser <- add_argument(parser, "--project",
                       help = "Samplename")

args <- parse_args(parser)

message("Adjusting log-Ratio...")
# TO DO: Conditional correction of logR
#Create df for Signature Extraction
############### adjusted logR (https://github.com/lima1/PureCN/issues/40)
#Solution in:
#https://github.com/maitnnguyen/Oseq_CNV_Paper/blob/main/src/pipeline/step4_CNVcall.R # nolint
df <- read_tsv(args$seg)
df_ploidy <- read_tsv(args$ploidy)

df <- df %>%  rename(all_of(c(sample = "ID", chromosome = "chrom")))
df_ploidy <- df_ploidy %>% 
            rename(all_of(c(purity = "Tumor Fraction", ploidy = "Ploidy")))

R <- 2^(df$seg.median.logR)

df$adj.seg <- (df_ploidy$purity*df_ploidy$ploidy*R + 2*(1-df_ploidy$purity)*(R - 1))/(df_ploidy$purity*df_ploidy$ploidy)

df_gistic <- df %>%
  dplyr::select(sample, chromosome, start, end, num.mark, adj.seg) %>%
   na.omit()

readr::write_tsv(df_gistic, file = paste0(args$project, "_df_gistic.seg"),
                 quote = "needed")