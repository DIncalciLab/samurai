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
    help = "Reference genome build", default = "hg38")
parser <- add_argument(parser, "--cutoff",
    help = "q-value cutoff to display bands", default = 0.05)

args <- parse_args(parser)

gistic_obj <- readGistic(gisticAllLesionsFile = args$all_lesions,
                        gisticAmpGenesFile = args$amp_genes,
                        gisticDelGenesFile = args$del_genes,
                        gisticScoresFile = args$gistic_scores)

png(filename="maftools_summary_mqc.png", width = 1600, height = 800, res = 150)
par(srt = 30, xpd = TRUE)
gisticChromPlot(gistic = gistic_obj,
                markBands = "all",
                ref.build = args$ref_build,
                cytobandOffset = 0.4,
                txtSize = 0.9, cytobandTxtSize = 0.5,
                fdrCutoff = args$cutoff,
                y_lims = c(-10, 10))
dev.off()

png(filename="maftools_bubble.png", width = 2000, height = 600, units = "px", res=150)
gisticBubblePlot(gistic = gistic_obj)
dev.off()
