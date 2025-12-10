#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(CINSignatureQuantification, quietly = TRUE)
    library(argparser, quietly = TRUE)
    library(dplyr, quietly = TRUE)
    library(readr, quietly=TRUE)
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
parser <- add_argument(parser, "--cpus",
                       default = 1,
                       help = "Number of cores to use.")
parser <- add_argument(parser, "--genome",
                       default = "hg38",
                       help = "Genome to use")
parser <- add_argument(parser, "--projectname",
                       default = "signatures",
                       help = "Prefix to include in output file(s).")
args <- parse_args(parser)

message("Starting Signature Extraction...")

segments <- read_tsv(args$seg_file, show_col_types = FALSE)
cnobj <- quantifyCNSignatures(object = segments,
    experimentName = args$projectname,
    method = "drews",
    cores = args$cpus,
    build = args$genome)

message("Saving RDS dataset...")
saveRDS(cnobj, paste0(args$project, "_signatures.rds"))

# Extract df with norm. expr. of signatures
df_activity <- as.data.frame(getActivities(cnobj, type = "threshold"))
df_activity <- df_activity %>%
    mutate(sample = rownames(df_activity)) %>%
    relocate((sample))

readr::write_tsv(df_activity, file = paste0(args$projectname, "_activity.txt"),
    quote = "needed")

# Only plot if samples are > 1
total_samples <- segments %>% pull(sample) %>% unique() %>% length()
png(filename = "ascat_sc_plot_by_component.png")
if (total_samples > 1) {
    plotSampleByComponent(object = cnobj)
}
dev.off()


png("signatures_summary_mqc.png")
plotActivities(object = cnobj, type = "threshold")
dev.off()

message("Complete.")
