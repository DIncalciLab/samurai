#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Biobase, quietly = TRUE)
    library(parallel, quietly = TRUE)
    library(methods, quietly = TRUE)
    library(QDNAseq, quietly = TRUE)
    library(argparser, quietly = TRUE)
}
)

binsToUse <- function(object) {

    if ("use" %in% colnames(fData(object))) {
        return(fData(object)$use)
    } else if ("filter" %in% colnames(fData(object))) {
        fData(object)$use <- fData(object)$filter
        fData(object)$filter <- NULL
        return(fData(object)$use)
    }
    rep(TRUE, times = nrow(object))
}


log2offset <- function(offset = .Machine$double.xmin) {
    # Get offset?
    if (missing(offset)) {
        offset <- getOption("QDNAseq::log2offset", .Machine$double.xmin)
        offset <- as.double(offset)
        stopifnot(is.finite(offset))
        return(offset)
    }
    # Reset offset?
    if (is.null(offset)) offset <- eval(formals(log2offset)$offset)
    # Set offset
    stopifnot(length(offset) == 1L)
    offset <- as.double(offset)
    stopifnot(is.finite(offset))
    options("QDNAseq::log2offset" = offset)
    offset
}

log2adhoc <- function(x, offset = log2offset(), inv = FALSE) {
    if (!inv) {
        x[x < 0] <- 0
        x <- x + offset
        log2(x)
    } else {
        x <- 2^x
        x - offset
    }
}

factor2int <- function(x) {
    return(as.integer(as.character(x)))
}

factor2double <- function(x) {
    return(as.double(as.character(x)))
}

exportSegments <- function(object, filter = TRUE, include_calls = FALSE) {
    require(QDNAseq)
    require(Biobase)
    if (filter) {
        object <- object[binsToUse(object), ]
    }
    calls <- assayDataElement(object, "calls")
    segments <- log2adhoc(assayDataElement(object, "segmented"))
    fnames <- pData(object)$name
    fd <- fData(object)
    datalist <- list()
    segment_table <- NULL
    for (i in 1:ncol(calls)) {
        d <- cbind(fd[, 1:3], calls[, i], segments[, i])
        sel <- !is.na(d[, 4])
        dsel <- d[sel, ]
        rleD <- rle(paste(d[sel, 1], d[sel, 4], sep = ":"))
        endI <- cumsum(rleD$lengths)
        posI <- c(1, endI[-length(endI)] + 1)
        chr <- dsel[posI, 1]
        pos <- dsel[posI, 2]
        end <- dsel[endI, 3]
        score <- dsel[posI, 4]
        segVal <- round(dsel[posI, 5], 2)
        bins <- rleD$lengths
        options(scipen = 100)
        if (include_calls) {
            out <- cbind(fnames[i], paste0("chr", chr), pos, end,
                bins, segVal, score)
            colnames(out) <- c("ID", "chrom", "start", "end", "num.mark",
                "seg.mean", "call")
        } else {
            out <- cbind(fnames[i], paste0("chr", chr), pos, end,
            bins, segVal)
            colnames(out) <- c("ID", "chrom", "start", "end", "num.mark",
                "seg.mean")
        }
        out <- as.data.frame(out)
        datalist[[fnames[i]]] <- out
    }
    segment_table <- do.call(rbind, datalist)
    segment_table$start <- factor2int(segment_table$start)
    segment_table$end <- factor2int(segment_table$end)
    segment_table[["num.mark"]] <- factor2int(segment_table[["num.mark"]])
    segment_table[["seg.mean"]] <- factor2double(segment_table[["seg.mean"]])
    if (include_calls) {
        segment_table[["call"]] <- factor2int(segment_table[["call"]])
    }
    return(segment_table)
}

create_pdf <- function(output, template, project) {
    pdf_file <- sprintf(template, project)
    pdf(paste(output, pdf_file, sep = "/"), width = 15, height = 10)
}

parser <- arg_parser("Run QDNAseq", hide.opts = TRUE)

parser <- add_argument(parser, "--cpus", type = "integer", default = 10,
                       help = "Number of CPU cores to use")
parser <- add_argument(parser, "--bin-size", type = "integer", default = 30,
                        help = "Bin size in kbp (default: 30)")
parser <- add_argument(parser, "--bin-data",
                        help = "Path to RDS containing bin data",
                        default = NULL)
parser <- add_argument(parser, "--project-name", default = "sWGS",
                        help = "Set the name for the collected data")
parser <- add_argument(parser, "--min-mapq", type = "integer",
                        default = 1,
                        help = "Minimum mapping quality (default: 1)")
parser <- add_argument(parser, "--purity", type = "integer",
                        default = 1,
                        help = "Set purity for the calling algorithm (default: 1)")
parser <- add_argument(parser, "--genome", default = "hg38",
                       help = "Genome to use")
parser <- add_argument(parser, "--paired-end",
                        help = "Whether reads are paired or not",
                        flag = TRUE)
parser <- add_argument(parser, "--source",
                        help = "Source BAM files (one or more)",
                        nargs = Inf)
parser <- add_argument(parser, "destination",
                        help = "Destination directory to save data to")

args <- parse_args(parser)

is_paired <- args$paired_end
bamfiles <- args$source

if (length(bamfiles) == 0) {
    stop("No BAM files found")
}

message(sprintf("Found %d BAM files.", length(bamfiles)))

output <- args$destination
message(sprintf("Destination directory: %s", output))

message("Initializing...")
# Initialize multisession
future::plan("multisession", workers = args$cpus)

options("R.cache.rootPath" = "~/.Rcache")
R.cache::getCacheRootPath()


if (!is.null(args$bin_data)) {
    bins <- readRDS(args$bin_data)
} else {
    bins <- getBinAnnotations(binSize = args$bin_size)
}

message("Reading data...")
read_counts <- binReadCounts(bins, bamfiles = bamfiles,
                             cache = TRUE,
                             minMapq = args$min_mapq,
                             pairedEnds = args$paired_end,
                             isPaired = TRUE,
                             isProperPair = TRUE,
                             chunkSize = TRUE)

message("Filtering read counts...")
reads_filtered <- applyFilters(read_counts, residual = TRUE,
                                blacklist = TRUE,
                                chromosome = c("X", "Y", "MT"))
reads_filtered <- estimateCorrection(reads_filtered)
dest_rds <- paste(output, sprintf("%s_reads_filtered.rds", args$project),
                  sep = "/")
saveRDS(reads_filtered, dest_rds)

message("Correcting and smoothing copy numbers...")
corrected_bins <- correctBins(reads_filtered)
corrected_bins <- normalizeBins(corrected_bins)#, method = 'median' as defaut)
corrected_bins <- smoothOutlierBins(corrected_bins)

create_pdf(output, template = "%s_bin_plot.pdf", args$project)
plot(corrected_bins)
dev.off()

message("Saving bins...")
exportBins(corrected_bins, file = "%s_bins.bed", format = "bed", filter = TRUE)
dest_rds <- paste(output, sprintf("%s_bins.rds", args$project),
                  sep = "/")
saveRDS(corrected_bins, dest_rds)

message("Segmenting data...")
segmented <- segmentBins(corrected_bins, transformFun = "log2")

create_pdf(output, "%s_segment_plot.pdf", args$project)
plot(segmented)
dev.off()

# Create Summary table
col_names <- c("sample", "nr_reads", "binsize", "expected_sd")
summary_table <- data.frame(reads_filtered@phenoData@data, row.names = NULL)[ c("name", "total.reads", "expected.variance")]
summary_table$expected.variance <- round(summary_table$expected.variance, 5)
summary_table$expected_sd <- round(sqrt(summary_table$expected.variance), 5)
summary_table$binsize <- bins@data$end[1]
write.table(summary_table, file = paste0(args$project, "_summary.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# This is useless, but QDNAseq won't export data if there isn't call information
called <- callBins(segmented, method = "CGHcall", build = args$genome,
                   cellularity = args$purity, nclass = 3)

message("Saving RDS data set...")
dest_rds <- paste(output, sprintf("%s.rds", args$project),
                  sep = "/")
saveRDS(called, dest_rds)

message("Exporting segments to SEG format...")

segments <- exportSegments(called, include_calls = TRUE)
segments_nocall <- exportSegments(called, include_calls = FALSE)

savedfiles <- exportBins(called, type = "segments",
                         file = "%s_filt.seg",
                         format = "seg", filter = TRUE)

write.table(segments, paste0(args$project, ".calls.seg"), sep = "\t",
            row.names = FALSE)
write.table(segments_nocall, paste0(args$project, "_.seg"), sep = "\t",
            row.names = FALSE)

message("Complete.")
