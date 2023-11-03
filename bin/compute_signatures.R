#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(CINSignatureQuantification, quietly = TRUE)
    library(argparser, quietly = TRUE)
    library(dplyr, quietly = TRUE)
}
)

create_pdf <- function(template, project) {
    pdf_file <- sprintf(template, project)
    pdf(paste(pdf_file, sep = "/"), width = 15, height = 10)
}

parser <- arg_parser("Compute Signatures", hide.opts = TRUE)

parser <- add_argument(parser, "--seg_file",
    help = "Input segmentation file",
    nargs = Inf)
parser <- add_argument(parser, "--cpus", type = "integer", default = 1,
    help = "Number of cores to use.")
parser <- add_argument(parser, "--projectname", type = "string",
    default = "Signatures",
    help = "Flag to include in output file(s).")
parser <- add_argument(parser, "--build", type = "string",
    default = "hg38",
    help = "Genome to use")
args <- parse_args(parser)

message("Starting Signature Extraction...")

cnobj <- quantifyCNSignatures(object = args$seg_file,
    experimentName = args$projectname,
    method = "drews",
    cores = args$cpus,
    build = args$build)

message("Saving RDS dataset...")
saveRDS(cnobj, paste0(args$project, "_signatures.rds"))

# Extract df with norm. expr. of signatures
df_activity <- as.data.frame(cnobj@activities[["thresholdAct2"]])
df_activity <- df_activity %>%
    mutate(sample = rownames(df_activity)) %>%
    relocate((sample))

readr::write_tsv(df_activity, file = paste0(args$project, "_activity.txt"),
    quote = "needed")


message("Starting Platinum Clinical Prediction...")
clin_pred <- clinPredictionPlatinum(object = cnobj)
saveRDS(clin_pred, paste0(args$project, "_platinum_prediction.rds"))

create_pdf("%s_plot_by_component.pdf", args$project)
plotSampleByComponent(object = cnobj)
dev.off()

create_pdf("%s_plot_activities.pdf", args$project)
plotActivities(object = cnobj, type = "threshold")
dev.off()

message("Complete.")
