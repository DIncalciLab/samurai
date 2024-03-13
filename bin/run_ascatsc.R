#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ASCAT.sc, quietly = TRUE)
    library(argparser, quietly = TRUE)
    library(copynumber, quietly = TRUE)
    library(readr, quietly = TRUE)
    library(dplyr, quietly = TRUE)
}
)

parser <- arg_parser("Run ASCAT.sc", hide.opts = TRUE)

parser <- add_argument(parser, "--tumor_bams",
    help = "Paths to the tumor BAMs (sorted and indexed)",
    nargs = Inf)
parser <- add_argument(parser, "--binsize", type = "integer",
    default = 30000,
    help = "Bin size in bp (default: 30000).")
parser <- add_argument(parser, "--cpus", type = "integer", default = 1,
    help = "Number of cores to use.")
parser <- add_argument(parser, "--projectname", type = "string",
    default = "ASCATsc",
    help = "Flag to include in output file(s).")
parser <- add_argument(parser, "--genome", type = "string",
    default = "hg38",
    help = "Genome to use")
parser <- add_argument(parser, "--outdir",
    default = "./",
    help = "Destination directory to save data to")
parser <- add_argument(parser, "--chromosomes",
    help = "Chromosomes to evaluate",
    default = "autosomes_x")
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
parser <- add_argument(parser, "--min-ploidy", type = "float",
    help = "Minimum ploidy to consider",
    default = 1.7)
parser <- add_argument(parser, "--max-ploidy", type = "float",
    help = "Maximum ploidy to consider",
    default = 5)
parser <- add_argument(parser, "--min-purity", type = "float",
    help = "Minimum purity to consider (fractional scale)",
    default = 0.05)
parser <- add_argument(parser, "--max-purity", type = "float",
    help = "Maximum purity to consider (fractional scale)",
    default = 1)
parser <- add_argument(parser, "--max-tumor-ploidy", type = "float",
    help = "Maximum tumour ploidy above which solutions will be masked
        (distance set to infinity in the grid search). Can reduce runtime.",
    default = 5)

args <- parse_args(parser)

allchr <- switch(
    args$chromosomes,
    "autosomes" = paste0(args$chrstring_bam, (1:22)),
    "autosomes_x" = paste0(args$chrstring_bam, c(1:22, "X")),
    "all" = paste0(args$chrstring_bam, c(1:22, "X", "Y")),
    stop("Invalid 'chromsomes' value. Use 'autosomes', 'autosomes_x' or 'all'.")
)

if (args$min_ploidy > args$max_ploidy) {
    stop("Minimum ploidy cannot be higher than maximum ploidy")
} else if (args$min_purity > args$max_purity) {
    stop("Minimum purity cannot be higher than maximum purity")
}

message("Starting analysis...")

multipcf <- FALSE

if (length(args$tumour_bams) > 1) {
    multipcf <- TRUE
}

purities <- seq(args$min_purity, args$max_purity, 0.01)
ploidies <- seq(args$min_ploidy, args$max_ploidy, 0.01)

res <- run_sc_sequencing(
    tumour_bams = args$tumor_bams,
    allchr = allchr,
    sex = args$sex,
    binsize = args$binsize,
    chrstring_bam = args$chrstring_bam,
    purs = purities,
    ploidies = ploidies,
    maxtumourpsi = args$max_tumor_ploidy,
    build = args$genome,
    MC.CORES = args$cpus,
    projectname = args$projectname,
    segmentation_alpha = args$segmentation_alpha,
    predict_refit = args$predict_refit,
    multipcf = multipcf,
    outdir = args$outdir
)

message("Creation of final segmentation dataframe...")
# create segmentation dataframe
df_final <- as.data.frame(ifelse(args$predict_refit,
        res[["allProfiles.refitted.auto"]],
        res[["allProfiles"]])) %>%
    mutate(samplename = res$summary$allSols$samplename,
           logr = as.numeric(logr))

message("Creation of summary df with purity and ploidy")
df_summary <- as.data.frame(res$summary$allSols) # Purity and ploidy are the same in refitted and non-refitted

# save output files
message("Saving output into rds...")
saveRDS(res, paste0(args$project, "_ASCAT.rds"))

readr::write_tsv(df_summary, file = paste0(args$project, "_summary.txt"),
    quote = "needed")
readr::write_tsv(df_final, file = paste0(args$project, "_segments.seg"),
    quote = "needed")

# Create df for Signature Extraction
message("Creating dataframe for signature extraction...")
df_sig <- df_final %>%
    dplyr::select(chromosome, start, end, total_copy_number_logr, samplename) %>%
    dplyr::rename(segVal = total_copy_number_logr,
                  sample = samplename) %>%
    na.omit()

readr::write_tsv(df_sig, file = paste0(args$project, "_df_signatures.seg"),
    quote = "needed")

# TODO: Conditional correction of logR
# Create df for GISTIC
############### adjusted logR (https://github.com/lima1/PureCN/issues/40)
# Solution in:
#https://github.com/maitnnguyen/Oseq_CNV_Paper/blob/main/src/pipeline/step4_CNVcall.R # nolint

message("Creating df for GISTIC Analysis...")
message("Computing R value...")
R <- 2^(df_final$logr)

message("Adding column with adjusted logR values...")
df_final$adj.seg <- (df_summary$purity * df_summary$ploidy * R + 2 * (1 - df_summary$purity) * (R - 1)) /
                     (df_summary$purity * df_summary$ploidy)
message("Creating final GISTIC df...")
df_gistic <- df_final %>%
    dplyr::select(samplename, chromosome, start, end, num.mark, adj.seg) %>%
    na.omit()
readr::write_tsv(df_gistic, file = paste0(args$projectname, "_gistic.seg"),
    quote = "needed")

message("Complete.")
