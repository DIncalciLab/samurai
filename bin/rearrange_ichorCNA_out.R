#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr, quietly = TRUE)
  library(argparser, quietly = TRUE)
  library(dplyr, quietly = TRUE)
}
)
parser <- arg_parser("Rearrange ichoRCNA Output.", hide.opts = TRUE)

parser <- add_argument(parser, "--seg_file",
                       help = "Segmentation File to rearrange.",
                       nargs = Inf)
args <- parse_args(parser)
df <- read_tsv(args$seg_file)

samplename <- gsub("\\..*", "", colnames(df)[[4]])
df <- df %>%
        select(chr, start, end, ends_with("logR_Copy_Number")) %>%
        mutate(sample = samplename) %>%
        rename("chromosome" = "chr") %>%
        rename(segVal = ends_with("logR_Copy_Number")) %>% na.omit()

readr::write_tsv(df, file = paste0(samplename, "_df_signatures.seg"),
                 quote = "needed")
