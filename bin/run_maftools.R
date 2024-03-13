#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(maftools, quietly = TRUE)
    library(argparser, quietly = TRUE)
}
)

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

png(filename="maftools_summary_mqc.png")
gisticChromPlot(gistic = gistic_obj, markBands = "all")
dev.off()

png(filename="maftools_bubble.png")
gisticBubblePlot(gistic = gistic_obj)
dev.off()
