#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparser)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# Silence R CMD check / lintr warnings
utils::globalVariables(c(
  "midpoint", "chrom_num", "max_pos", "cumulative_start", "cumulative_end",
  "chrom_center", "plot_position", "Corrected_Copy_Number",
  "call", "start", "end"
))

# Argument parsing
p <- arg_parser("Plot IchorCNA output")
p <- add_argument(p, "--id", help = "Sample ID")
p <- add_argument(p, "--seg_file", help = "Segmented CN file")
p <- add_argument(p, "--binfile", help = "Binned CN file")
p <- add_argument(p, "--summary", help = "IchorCNA summary file")
p <- add_argument(p, "--outdir", help = "Output directory", default = ".")
argv <- parse_args(p)

sample_id <- argv$id
seg_file <- argv$seg_file
bin_file <- argv$binfile
summary_file <- argv$summary
output_dir <- argv$outdir

# Helper functions

# Read summary info for specific sample
get_summary <- function(file, sample_id) {
  lines <- readLines(file)
  header <- unlist(strsplit(lines[1], "\t"))

  # Find the row for our sample
  sample_row <- NULL
  for (i in 2:length(lines)) {
    data <- unlist(strsplit(lines[i], "\t"))
    if (length(data) > 0 && data[1] == sample_id) {
      sample_row <- data
      break
    }
  }

  if (is.null(sample_row)) {
    stop(paste("Sample", sample_id, "not found in summary file"))
  }

  df <- as.data.frame(t(sample_row), stringsAsFactors = FALSE)
  colnames(df) <- header

  # TODO: Consider outputting TF/ploidy as a machine-readable format
  df$`Tumor Fraction` <- as.numeric(df$`Tumor Fraction`)
  df$Ploidy <- as.numeric(df$Ploidy)

  # Look for MAD value in the file
  mad_line <- grep("GC-Map correction MAD:", lines, value = TRUE)
  if (length(mad_line) > 0) {
    mad_val <- as.numeric(sub("GC-Map correction MAD:\\s*", "", mad_line[1]))
    df$MAD <- mad_val
  } else {
    df$MAD <- NA
  }

  return(df)
}

# Clean column names in bins file - remove sample prefix
clean_bins_colnames <- function(df, sample_id) {
  cols <- colnames(df)

  # Keep chr, start, end as-is, remove sample prefix from others
  new_cols <- cols
  for (i in 1:length(cols)) {
    if (startsWith(cols[i], paste0(sample_id, "."))) {
      new_cols[i] <- gsub(paste0("^", sample_id, "\\."), "", cols[i])
    }
  }

  # Rename chr to chrom for consistency with segments
  new_cols[new_cols == "chr"] <- "chrom"

  colnames(df) <- new_cols
  return(df)
}

prepare_plot_data <- function(df) {
  # Filter first to keep only autosomes
  df <- df %>%
    filter(!grepl("X|Y|M|MT", chrom)) %>%
    mutate(
      # Extract numeric chromosome identifier
      chrom_num = as.numeric(gsub("chr", "", chrom)),
      # Calculate midpoint for plotting
      midpoint = (start + end) / 2
    ) %>%
    # Remove rows with NA in any *_Copy_Number column
    filter(!if_any(matches("_Copy_Number$"), is.na)) %>%
    # Sort by chromosome and position
    arrange(chrom_num, midpoint)

  # Compute the maximum position per chromosome to determine its length
  chrom_lengths <- df %>%
    group_by(chrom_num, chrom) %>%
    summarise(
      max_pos = max(end),  # end coordinate is treated as chromosome length
      .groups = "drop"
    ) %>%
    arrange(chrom_num)

  # Define the full set of chromosomes 1 to 22 (ensures consistency)
  full_chrom_df <- data.frame(
    chrom_num = 1:22,
    chrom = paste0("chr", 1:22)
  )

  # Join with actual chromosome lengths to ensure completeness
  chrom_lengths_full <- full_chrom_df %>%
    left_join(chrom_lengths, by = c("chrom_num", "chrom")) %>%
    mutate(
      # Fill missing chromosome lengths with 0 (in case some are absent in data)
      max_pos = ifelse(is.na(max_pos), 0, max_pos),

      # Calculate cum. start pos. of each chrom for concatenated plotting
      cumulative_start = lag(cumsum(as.numeric(max_pos)), default = 0),

      # Cumulative end = start + length
      cumulative_end = cumulative_start + max_pos,

      # Chromosome center (used for axis labels)
      chrom_center = cumulative_start + max_pos / 2
    )

  # Add the cumulative_start position back to the main dataframe
  df <- df %>%
    left_join(
      chrom_lengths_full %>%
        select(chrom_num, cumulative_start),
      by = "chrom_num"
    ) %>%

    # Compute the final plotting position as midpoint within concatenated genome
    mutate(
      plot_position = cumulative_start + midpoint
    )

  # Return both processed data and chromosome metadata
  return(list(df = df, chrom_lengths = chrom_lengths_full))
}

plot_combined_copy_number <- function(df_bins, df_segments, chrom_lengths,
                                      tf, ploidy, mad, sample_id, output_dir) {
  # Color map for calls
  color_mapping <- c(
    "NEUTRAL" = "#4b83d8",
    "AMPLIFICATION" = "#F24236",
    "DELETION" = "#08a50c"
  )

  # Define Y-axis range
  copy_numbers <- df_bins$logR_Copy_Number

  if (length(copy_numbers) > 0) {
    # Use actual data range with sensible bounds
    data_min <- min(copy_numbers)
    data_max <- max(copy_numbers)

    # Set reasonable limits
    y_lower <- max(data_min - 0.5, -1)  # Don't go below -1
    y_upper <- min(data_max + 0.5, 12)  # Don't go above 12
  } else {
    # Fallback for empty data
    y_lower <- -1
    y_upper <- 6
  }

  # Base plot
  p <- ggplot() +
    # Grey background for alternate chromosomes
    geom_rect(
      data = chrom_lengths %>% filter(chrom_num %% 2 == 0),
      aes(xmin = cumulative_start, xmax = cumulative_end,
          ymin = -Inf, ymax = Inf),
      fill = "gray95", alpha = 1
    ) +

    # Raw logR points
    geom_point(
      data = df_bins,
      aes(x = plot_position, y = logR_Copy_Number),
      color = "grey39", size = 0.4, alpha = 0.3, stroke = 0
    ) +

    # Segmented corrected copy number calls
    geom_segment(
      data = df_segments,
      aes(x = cumulative_start + start,
          xend = cumulative_start + end,
          y = Corrected_Copy_Number,
          yend = Corrected_Copy_Number,
          color = call),
      linewidth = 1.7, alpha = 0.9, lineend = "round"
    ) +

    # Axis and color configuration
    scale_color_manual(
      values = color_mapping,
      name = "Copy Number Call",
      guide = guide_legend(override.aes = list(linewidth = 4, alpha = 1))
    ) +
    scale_x_continuous(
      breaks = chrom_lengths$chrom_center,
      labels = gsub("chr", "", chrom_lengths$chrom),
      expand = c(0.005, 0)
    ) +
    scale_y_continuous(
      limits = c(y_lower, y_upper),
      breaks = pretty_breaks(n = 10),
      expand = expansion(mult = c(0.02, 0.02))
    ) +

    # Labels
    labs(
      title = paste("Sample:", sample_id),
      subtitle = paste0("Tumor Fraction: ", round(tf, 3),
                        " - Ploidy: ", round(ploidy, 2),
                        " - MAD: ", round(mad, 4)),
      x = "Chromosome",
      y = "Copy Number"
    ) +

    # Styling
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.5, linetype = "dotted"),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "grey25"),
      axis.text.y = element_text(size = 9, color = "grey25"),
      axis.title = element_text(size = 10, color = "grey15", face = "bold"),
      axis.line = element_line(color = "grey35", linewidth = 0.5),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold", color = "grey10", margin = margin(b = 6)),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey35", margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_text(size = 10, face = "bold", color = "grey15"),
      legend.text = element_text(size = 9, color = "grey25"),
      legend.key.size = unit(0.7, "cm"),
      legend.margin = margin(t = 10),
      legend.box.spacing = unit(0.2, "cm"),
      plot.margin = margin(t = 15, r = 15, b = 12, l = 15)
    )

  # Save the figure
  out_file <- file.path(output_dir, paste0(sample_id, ".copy_number.png"))
  out_file_svg <- file.path(output_dir, paste0(sample_id, ".copy_number.svg"))
  ggsave(out_file, plot = p, width = 14, height = 6, dpi = 300)
  ggsave(out_file_svg, plot = p, width = 14, height = 6)
  message("Plot (PNG) saved to:", out_file, "\n")
  message("Plot (SVG) saved to:", out_file_svg, "\n")
}

# Main execution
message("Processing sample:", sample_id, "\n")

# Read and process summary
summary_df <- get_summary(summary_file, sample_id)
ploidy <- summary_df$Ploidy
tf <- summary_df$`Tumor Fraction`
mad <- ifelse(is.na(summary_df$MAD), 0, summary_df$MAD)

message("Sample info - TF:", tf, "Ploidy:", ploidy, "MAD:", mad, "\n")

# Read bins file and clean column names
df_bins <- read_tsv(bin_file, show_col_types = FALSE)
df_bins <- clean_bins_colnames(df_bins, sample_id)

# Create dynamic calls for bins based on ploidy
df_bins <- df_bins %>%
  dplyr::mutate(call = case_when(
    Corrected_Copy_Number > round(ploidy) ~ "AMPLIFICATION",
    Corrected_Copy_Number < round(ploidy) ~ "DELETION",
    TRUE ~ "NEUTRAL"
  ))

# Read segments file (already has standard column names)
df_segments <- read_tsv(seg_file, show_col_types = FALSE)

# Create dynamic calls for segments based on ploidy
df_segments <- df_segments %>%
  dplyr::mutate(call = case_when(
    Corrected_Copy_Number > round(ploidy) ~ "AMPLIFICATION",
    Corrected_Copy_Number < round(ploidy) ~ "DELETION",
    TRUE ~ "NEUTRAL"
  ))

message("Data loaded - Bins:", nrow(df_bins), "Segments:", nrow(df_segments), "\n")

# Prepare plot data
prep_bins <- prepare_plot_data(df_bins)
prep_segments <- prepare_plot_data(df_segments)

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Generate plot
plot_combined_copy_number(prep_bins$df, prep_segments$df, prep_bins$chrom_lengths, tf, ploidy, mad, sample_id, output_dir) # nolint: line_length_linter.

message("Plot generation completed!\n")
