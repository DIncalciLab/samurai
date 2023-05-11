#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(argparser, quietly = TRUE)
}
)

append_seg <- function(sample_list, file_name) {
        sample <- read.table(sample_list[1], header = TRUE)
        for (file in sample_list[-1]) {
                sample <- rbind(sample, read.table(file, header = TRUE))
                }
        write.table(sample,
                paste0(file_name, ".seg"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE)
}

parser <- arg_parser("Concatenate segments", hide.opts = TRUE)

parser <- add_argument(parser, "--segfiles",
                        help = "Input Segmentation files to be concatenated",
                        nargs = Inf)

parser <- add_argument(parser, "--output",
                        help = "Name of output segmentation file",
                        default = "segments")
args <- parse_args(parser)

message("Concatenating segments..")

options("R.cache.rootPath" = "~/.Rcache")
R.cache::getCacheRootPath()

append_seg(args$segfiles, args$output)

message("Complete!")