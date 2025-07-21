#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(dplyr, quietly = TRUE)
}
)


adjust_log_ratio <- function(copy_number, ploidy, min.ratio = 2^-8) {

    adjusted <- copy_number / ploidy
    adjusted <- pmax(min.ratio, adjusted)
    adjusted <- log2(adjusted)
    return(adjusted)
}

parser <- arg_parser("Adjust logR.", hide.opts = TRUE)

parser <- add_argument(parser, "--seg",
                       help = "Segmentation file from which extract information.", #nolint
                       nargs = Inf)

parser <- add_argument(parser, "--ploidy",
                       help = "Ploidy Summary file from which extract information.", #nolint
                       nargs = Inf)

parser <- add_argument(parser, "--project",
                       help = "Samplename")

args <- parse_args(parser)

message("Adjusting log-Ratio...")


# TO DO: Conditional correction of logR

# IMPORT TWO DATAFRAMES
df <- read_tsv(args$seg)
df_ploidy <- read_tsv(args$ploidy)

#RENAME COLUMNS

df <- df %>%  rename(all_of(c(chromosome = "chrom", sample = "ID")))

df_ploidy <- df_ploidy %>%
  rename(all_of(c(purity = "Tumor Fraction",
                  ploidy = "Ploidy",
                  sample = "samplename")))

# CREATE A TEMPORARY COMPLETE DF TO CORRECT LOGR
df_tmp <- df %>% select(sample) %>% left_join(df_ploidy, by = 'sample') %>% select(-sample)
df <- cbind(df, df_tmp)

# adjusted logR (https://github.com/lima1/PureCN/issues/40)
#Solution in:
#https://github.com/maitnnguyen/Oseq_CNV_Paper/blob/main/src/pipeline/step4_CNVcall.R # nolint

df <- df %>%
    dplyr::mutate(adj.seg = adjust_log_ratio(logR_Copy_Number,
        ploidy = df_ploidy$ploidy))

# df$adj.seg <- (df$purity*df$ploidy*R + 2*(1-df$purity)*(R - 1))/(df$purity*df$ploidy) #nolint

# FINAL DF WITH SELECTED COLUMNS
df_gistic <- df %>%
            dplyr::select(sample, chromosome, start, end, num.mark, adj.seg) %>%
            rowwise %>%
            filter(!is.infinite(adj.seg))

readr::write_tsv(df_gistic,
                file = "segments_logR_corrected_gistic.seg",
                quote = "needed")
