#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(maftools, quietly = TRUE)
    library(argparser, quietly = TRUE)
}
)

create_pdf <- function(template, project) {
    pdf_file <- sprintf(template, project)
    pdf(paste(pdf_file, sep = "/"), width = 15, height = 10)
}

parser <- arg_parser("Run Maftools", hide.opts = TRUE)

parser <- add_argument(parser, "--all_lesions",
                    help = "Paths to the GISTIC lesions file.",
                    nargs = Inf)
parser <- add_argument(parser, "--amp_genes",
                    help = "Paths to the GISTIC amplified genes file.",
                    nargs = Inf)
parser <- add_argument(parser, "--del_genes",
                    help = "Paths to the GISTIC deleted genes file.",
                    nargs = Inf)
parser <- add_argument(parser, "--gistic_scores",
                    help = "Paths to the GISTIC scores file.",
                    nargs = Inf)

args <- parse_args(parser)

gistic_obj <- readGistic(gisticAllLesionsFile = args$all_lesions,
                gisticAmpGenesFile = args$amp_genes,
                gisticDelGenesFile = args$del_genes,
                gisticScoresFile = args$gistic_scores)

create_pdf("%s_gistic_chrom.pdf", "plot")
gisticChromPlot(gistic = gistic_obj, markBands = "all")
dev.off()

create_pdf("%s_gistic_bubble.pdf", "plot")
gisticBubblePlot(gistic = gistic_obj)
dev.off()
