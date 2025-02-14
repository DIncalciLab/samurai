#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(HRDCNA, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(argparser, quietly = TRUE)
    library(readr, quietly = TRUE)
    library(sigminer, quietly = TRUE)
}
)

sigminercopy_hg38 <-function(data, genome_build = "hg38"){
  copy <- sigminer::read_copynumber(data,
                                    seg_cols = c("chromosome", "start", "end", "segVal"),
                                    genome_build="hg38", complement = FALSE, verbose = TRUE)
  
  
  cnall_tally_W <- sigminer::sig_tally(copy, method = "W")
  nmf <- cnall_tally_W$nmf_matrix
  nmf <- as.data.frame(nmf)
  nmf$sample <- rownames(nmf)
  rownames(nmf) <- NULL
  return(nmf=nmf)
}

parser <- arg_parser("Run HRDCNA", hide.opts = TRUE)

parser <- add_argument(parser, "--seg_file",
                       help = "Input segmentation file",
                       nargs = Inf)
parser <- add_argument(parser, "--genome", type = "string",
    default = "hg38",
    help = "Genome to use")
parser <- add_argument(parser, "--hrdcna_threshold", type = "float",
    default = 0.2,
    help = "Threshold to classify samples into 'HRD' or 'HRP'.")

args <- parse_args(parser)

message("Starting HRDCNA score analysis...")

segments <- read_tsv(args$seg_file, show_col_types = FALSE)
# Dataframe MUST have chr, start, end, segVal, sample
segments$segVal <- round(segments$segVal, 0)
# Perform HRD analysis
# Perform HRD analysis based on genome build
nmfcn_wgs <- if (args$genome == "hg38") {
  sigminercopy_hg38(data = segments, genome_build = "hg38")
} else if (args$genome == "hg19") {
  sigminercopy(data = segments, genome_build = "hg19")
} else {
  stop("Unsupported genome build. Please specify either 'hg38' or 'hg19'.")
}

message("Saving RDS with Sigminer features activity...")
saveRDS(nmfcn_wgs, "hrdcna_features_activity.rds")
message("Sigminer features activity RDS has been saved...")

readr::write_tsv(nmfcn_wgs, file = "hrdcna_features_activity.tsv", quote = "needed")

score_wgs <- HRDprediction(data = nmfcn_wgs)

# Assign HRD status based on threshold declared in the original paper
score_wgs$hrd_status_predicted <- ifelse(score_wgs$HRDCNAScore >= args$hrdcna_threshold, "HRD+", "HRD-")
# Order columns for multiqc
score_wgs <- score_wgs %>% dplyr::select(sample, HRDCNAScore, hrd_status_predicted)

message("Saving RDS with HRDCNA Scores...")
saveRDS(score_wgs, "hrdcna_scores.rds")
message("HRDCNA Scores RDS has been saved...")

readr::write_tsv(score_wgs, file = "hrdcna_summary_mqc.tsv", quote = "needed")

message("Complete.")