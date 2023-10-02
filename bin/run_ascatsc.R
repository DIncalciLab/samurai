#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ASCAT.sc, quietly = TRUE)
    library(argparser, quietly = TRUE)
    library(copynumber, quietly = TRUE)
    library(future, quietly = TRUE)
}
)

parser <- arg_parser("Run ASCAT.sc", hide.opts = TRUE)

parser <- add_argument(parser, "--tumour_bams",
                       help = "Paths to the bams (sorted and indexed).",
                       nargs = Inf)
parser <- add_argument(parser, "--binsize", type = "integer",
                       default = 30000,
                       help = "Bin size in bp (default: 30000).")
parser <- add_argument(parser, "--cpus", type = "integer", default = 1,
                       help = "Number of cores to use.")
parser <- add_argument(parser, "--projectname", type = "string",
                       default = "ASCATsc",
                       help = "Flag to include in output file(s).")
parser <- add_argument(parser, "--build", type = "string",
                       default = "hg38",
                       help = "Genome to use")
parser <- add_argument(parser, "--outdir",
                       default = "./",
                       help = "Destination directory to save data to")
parser <- add_argument(parser, "--allchr", type = "string",
                       default = paste0("chr", (1:22)),
                       help = "Chromosome names.")
parser <- add_argument(parser, "--chrstring_bam", type = "string",
                       default = "chr",
                       help = "Chromosome prefix present in BAMs.")
parser <- add_argument(parser, "--sex", type = "string",
                       default = "female",
                       help = "Patient gender (default: female)")
parser <- add_argument(parser, "--segmentation_alpha", type = "float",
                       default = 0.01,
                       help = "Segmentation parameter for CBS.")
parser <- add_argument(parser, "--predict_refit", type = "boolean",
                       default = TRUE,
                       help = "xgboost predictor to predict refitted ploidy \n
                                and refit the profiles automatically")
parser <- add_argument(parser, "--max_tumor_ploidy", type = "float",
                       default = 5,
                       help = "Maximum tumour ploidy above which solutions \n
                       will be masked")


args <- parse_args(parser)

message("Starting analysis...")

res <- run_sc_sequencing(tumour_bams = args$tumour_bams,
                         allchr = paste0("chr", (1:22)),
                         sex = args$sex,
                         binsize = args$binsize,
                         chrstring_bam = args$chrstring_bam,
                         purs = seq(0.05, 1, 0.01),
                         ploidies = seq(1.7, 5, 0.01),
                         maxtumourpsi = 5,
                         build = args$build,
                         MC.CORES = args$cpus,
                         projectname = args$projectname,
                         segmentation_alpha = args$segmentation_alpha,
                         predict_refit = args$predict_refit)

for (i in seq_along(length(args$tumour_bams))) {
  if (args$predict_refit == TRUE)
    {
      df_summary <- as.data.frame(res$summary$allSols.refitted)
      df_segment <- as.data.frame(res[["allProfiles.refitted.auto"]][[1]])
      df_seg <- as.data.frame(res[["allProfiles.refitted.auto"]][[i]])
      df_final <- rbind(df_segment, df_seg)
    } else {
        df_summary <- as.data.frame(res$summary$allSols)
        df_segment <- as.data.frame(res[["allProfiles"]][[1]])
        df_seg <- as.data.frame(res[["allProfiles"]][[i]])
        df_final <- rbind(df_segment, df_seg)
        }
}
rownames(df_summary) <- NULL
rownames(df_final) <- NULL

write.table(df_summary, file = paste0(args$project, "_summary.txt"),
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
write.table(df_final,
            file = paste0(args$project, "_segments.seg"),
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

df_gistic <- df_final[c("chromosome", "start", "end", "num.mark", "logr")]
df_gistic$logr <- round(as.numeric(df_gistic$logr), 5)
write.table(df_gistic, file = paste0(args$project, "_gistic.seg"),
            quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
message("Complete.")
