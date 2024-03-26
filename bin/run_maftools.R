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
parser <- add_argument(parser, "--ref_build",
    help = "Reference genome build", default="hg38")

args <- parse_args(parser)

gistic_obj <- readGistic(gisticAllLesionsFile = args$all_lesions,
                        gisticAmpGenesFile = args$amp_genes,
                        gisticDelGenesFile = args$del_genes,
                        gisticScoresFile = args$gistic_scores)

png(filename="maftools_summary_mqc.png", width = 800, height = 480, res=90)
par(srt = 30, xpd = TRUE)
gisticChromPlot(gistic = gistic_obj, 
                markBands = "all", 
                ref.build = args$ref_build, cytobandOffset = 0.4, txtSize = 0.9, cytobandTxtSize = 0.5)
dev.off()

png(filename="maftools_bubble.png", width = 2000, height = 600, units = "px", res=150)
gisticBubblePlot(gistic = gistic_obj)
dev.off()
