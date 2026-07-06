#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

## arguments

p <- arg_parser("Plot WisecondorX output (bin-level log2 ratio + segments)")
p <- add_argument(p, "--id", help = "Sample ID")
p <- add_argument(p, "--seg_file", help = "Segmented file with calls from WisecondorX", nargs = Inf)
p <- add_argument(p, "--binfile", help = "Bin-level file from WisecondorX")
p <- add_argument(p, "--outdir", help = "Output directory", default = ".")
p <- add_argument(p, "--ratio_limit", help = "Y-axis limit (+-) in log2(ratio) scale; values beyond this are visually clipped", default = 1)
argv <- parse_args(p)

sample_id      <- argv$id
bin_path       <- argv$binfile
output_dir     <- argv$outdir
ratio_limit    <- argv$ratio_limit

# --seg_file puo' ricevere piu' di un path se, a monte, un glob troppo
# permissivo (es. "*.seg") matcha anche il file "*_gistic.seg" generato
# nello stesso step (che NON ha la colonna 'call', e' pensato per GISTIC).
# Scartiamo esplicitamente quel file e teniamo quello con i call reali.
seg_candidates <- argv$seg_file
if (length(seg_candidates) > 1) {
  message("Ricevuti piu' file per --seg_file (", paste(seg_candidates, collapse = ", "),
          "): scarto quelli che terminano in '_gistic.seg'.")
  seg_candidates <- seg_candidates[!grepl("_gistic\\.seg$", seg_candidates)]
  if (length(seg_candidates) == 0) {
    stop("Nessun file valido rimasto per --seg_file dopo aver escluso '_gistic.seg'.")
  }
}
seg_path <- seg_candidates[1]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper functions ------------------------------------------------------

# Fix column names
read_and_normalize <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE) %>%
    rename_with(tolower)

  if ("chrom" %in% names(df)) df <- rename(df, chr = chrom)
  if ("seg.mean.adj" %in% names(df)) df <- rename(df, ratio = `seg.mean.adj`)

  df
}

# Rimappa i call gain/loss/neut di WisecondorX in etichette leggibili per il plot
classify_segments <- function(df) {
  call_map <- c(gain = "GAIN", loss = "LOSS", neut = "NEUTRAL")
  df %>% mutate(call = recode(tolower(call), !!!call_map, .default = "NEUTRAL"))
}

# Plot style formatting, genomic coordinates
compute_genomic_coords <- function(bins, segs) {
  chr_order <- c(as.character(1:22), "X", "Y")
  present_chr <- intersect(chr_order, unique(c(bins$chr, segs$chr)))

  chr_lengths <- bins %>%
    filter(chr %in% present_chr) %>%
    group_by(chr) %>%
    summarise(len = max(end), .groups = "drop") %>%
    mutate(chr = factor(chr, levels = present_chr)) %>%
    arrange(chr) %>%
    mutate(offset = lag(cumsum(len), default = 0))

  add_offset <- function(df) {
    df %>%
      filter(chr %in% present_chr) %>%
      mutate(chr = factor(chr, levels = present_chr)) %>%
      left_join(chr_lengths %>% select(chr, offset), by = "chr") %>%
      mutate(start_g = start + offset, end_g = end + offset)
  }

  list(
    bins = add_offset(bins) %>% mutate(pos_g = (start_g + end_g) / 2),
    segs = add_offset(segs),
    chr_lengths = chr_lengths %>% mutate(xmax = offset + len)
  )
}

# Plot wisecondorx plot
plot_wisecondorx_cnv <- function(bins, segs, chr_lengths, sample_id, ratio_limit, output_dir) {
  chr_ticks <- chr_lengths %>% mutate(mid = offset + len / 2)
  chr_boundaries <- chr_lengths$xmax[-nrow(chr_lengths)]

  color_mapping <- c(
    "NEUTRAL" = "#377eb8",
    "GAIN" = "#e41a1c",
    "LOSS" = "#4daf4a"
  )

  p <- ggplot() +
    geom_vline(xintercept = chr_boundaries, color = "grey85", linewidth = 0.4, linetype = "dotted") +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.4) +
    geom_point(
      data = bins,
      aes(x = pos_g, y = pmin(pmax(ratio, -ratio_limit), ratio_limit)),
      color = "grey75", size = 0.25, alpha = 0.35
    ) +
    geom_segment(
      data = segs,
      aes(x = start_g, xend = end_g, y = ratio, yend = ratio, color = call),
      linewidth = 1.6, lineend = "round"
    ) +
    scale_color_manual(
      values = color_mapping,
      name = "Copy Number Call",
      guide = guide_legend(override.aes = list(linewidth = 4))
    ) +
    scale_x_continuous(breaks = chr_ticks$mid, labels = chr_ticks$chr, expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(-ratio_limit, ratio_limit)) +
    labs(
      x = "Chromosome", y = expression(log[2](ratio)),
      title = paste0("Copy Number Profile", if (!is.null(sample_id)) paste0(" - ", sample_id) else "")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text = element_text(color = "grey30"),
      axis.title = element_text(color = "grey20"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40", size = 10),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )

  out_png <- file.path(output_dir, paste0(sample_id, ".copy_number.png"))
  out_svg <- file.path(output_dir, paste0(sample_id, ".copy_number.svg"))
  ggsave(out_png, plot = p, width = 14, height = 6, dpi = 300)
  ggsave(out_svg, plot = p, width = 14, height = 6)
  message("Plot (PNG) saved to: ", out_png)
  message("Plot (SVG) saved to: ", out_svg)
}

# ---- main -------------------------------------------------------------

message("Processing sample: ", sample_id)

bins <- read_and_normalize(bin_path) %>%
  filter(!is.na(ratio))
segs <- read_and_normalize(seg_path) %>%
  classify_segments()

message("Data loaded - Bins: ", nrow(bins), " Segments: ", nrow(segs))

coords <- compute_genomic_coords(bins, segs)

plot_wisecondorx_cnv(coords$bins,
                     coords$segs,
                     coords$chr_lengths,
                     sample_id,
                     ratio_limit,
                     output_dir)

message("Plot generation completed!")