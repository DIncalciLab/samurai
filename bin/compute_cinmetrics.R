#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(CINmetrics, quietly = TRUE)
  library(argparser, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(dplyr, quietly = TRUE)
})

parser <- arg_parser("Compute CINmetrics", hide.opts = TRUE)

parser <- add_argument(parser, "--seg_file",
  help = "Input segmentation file",
  nargs = Inf)
parser <- add_argument(parser, "--output",
  help = "Output TSV file",
  default = "cinmetrics_output.tsv")
parser <- add_argument(parser, "--segmentMean_tai",
  help = "Threshold for TAI",
  type = "numeric", default = 0.2)
parser <- add_argument(parser, "--segmentMean_cna",
  help = "Threshold for CNA",
  type = "numeric", default = 0.2)
parser <- add_argument(parser, "--segmentMean_base_segments",
  help = "Threshold for base segments",
  type = "numeric", default = 0.2)
parser <- add_argument(parser, "--segmentMean_break_points",
  help = "Threshold for break points",
  type = "numeric", default = 0.2)
parser <- add_argument(parser, "--segmentMean_fga",
  help = "Threshold for FGA",
  type = "numeric", default = 0.2)
parser <- add_argument(parser, "--numProbes",
  help = "Minimum number of probes",
  type = "numeric", default = NA)
parser <- add_argument(parser, "--segmentDistance_cna",
  help = "Segment distance for CNA",
  type = "numeric", default = 0.2)
parser <- add_argument(parser, "--minSegSize_cna",
  help = "Minimum segment size for CNA",
  type = "numeric", default = 10)
parser <- add_argument(parser, "--genomeSize_fga",
  help = "Genome size for FGA",
  type = "numeric", default = 2873203431)

args <- parse_args(parser)

input_file <- args$seg_file
cat("Reading:", input_file, "\n")
df <- read_tsv(input_file, 
               show_col_types = FALSE,
               col_names = c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean"),
               col_types = cols(Sample        = col_character(),
                                Chromosome    = col_character(),
                                Start         = col_double(),
                                End           = col_double(),
                                Num_Probes    = col_double(),
                                Segment_Mean  = col_double())
) %>%
  dplyr::mutate(
    Chromosome = gsub("^chr", "", Chromosome)
  )

cat("CINmetrics calculation...\n")
metrics <- CINmetrics(
  cnvData                    = df,
  segmentMean_tai            = args$segmentMean_tai,
  segmentMean_cna            = args$segmentMean_cna,
  segmentMean_base_segments  = args$segmentMean_base_segments,
  segmentMean_break_points   = args$segmentMean_break_points,
  segmentMean_fga            = args$segmentMean_fga,
  numProbes                  = args$numProbes,
  segmentDistance_cna        = args$segmentDistance_cna,
  minSegSize_cna             = args$minSegSize_cna,
  genomeSize_fga             = args$genomeSize_fga
) %>%
  mutate(tai = coalesce(tai, 0))

#metrics$tai[is.na(metrics$tai)] <- 0

write_tsv(metrics, args$output)
cat("Output file:", args$output, "\n")