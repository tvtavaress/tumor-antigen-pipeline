#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(tools)
})

# =========================
# ARGUMENTS
# =========================
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input CSV file with specificity summary"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "specificity_groups_heatmap.svg",
    help = "Output heatmap SVG file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("Input CSV file does not exist: ", opt$input)
}

# =========================
# LOAD AND VALIDATE DATA
# =========================
df <- read.csv(
  opt$input,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "celltypist_labels",
  "specificity_profile",
  "n_genes"
)

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

df <- df %>%
  filter(
    !is.na(specificity_profile),
    specificity_profile != "Not detected"
  )

if (nrow(df) == 0) {
  stop("No specificity profiles available after filtering.")
}

# =========================
# FORMAT CELL TYPE NAMES
# =========================
df <- df %>%
  mutate(
    celltypist_labels = case_when(
      celltypist_labels == "epithelial cell of thymus" ~ "TEC",
      celltypist_labels == "plasmacytoid dendritic cell, human" ~ "pDCs",
      celltypist_labels == "hematopoietic multipotent progenitor cell" ~ "MPPs",
      str_detect(celltypist_labels, regex("cd4", ignore_case = TRUE)) ~
        toupper(celltypist_labels),
      str_detect(celltypist_labels, regex("cd8", ignore_case = TRUE)) ~
        toupper(celltypist_labels),
      TRUE ~ tools::toTitleCase(celltypist_labels)
    ),
    specificity_profile = case_when(
      specificity_profile == "Diferentiation" ~ "Differentiation",
      TRUE ~ specificity_profile
    ),
    n_genes = as.numeric(n_genes)
  ) %>%
  filter(
    !is.na(celltypist_labels),
    !is.na(n_genes)
  )

# =========================
# MERGE CATEGORIES
# =========================
df <- df %>%
  group_by(
    celltypist_labels,
    specificity_profile
  ) %>%
  summarise(
    n_genes = sum(n_genes),
    .groups = "drop"
  )

# =========================
# BUILD MATRIX
# =========================
mat_raw <- df %>%
  pivot_wider(
    names_from = specificity_profile,
    values_from = n_genes,
    values_fill = 0
  ) %>%
  column_to_rownames("celltypist_labels") %>%
  as.matrix()

mat_raw <- t(mat_raw)
storage.mode(mat_raw) <- "numeric"

# Ensure consistent order of specificity profiles
preferred_order <- c(
  "Enriched",
  "Group Enriched",
  "Enhanced",
  "Low specificity"
)

ordered_rows <- c(
  intersect(preferred_order, rownames(mat_raw)),
  setdiff(rownames(mat_raw), preferred_order)
)

mat_raw <- mat_raw[ordered_rows, , drop = FALSE]

# =========================
# LOG TRANSFORM
# =========================
mat <- log10(mat_raw + 1)

# =========================
# COLOR SCALE
# =========================
count_breaks <- c(0, 10, 100, 1000, 3000, 6000, 10000)

breaks <- log10(count_breaks + 1)

base_colors <- c(
  "#F2F2F2",
  "#D9D9D9",
  "#BFBFBF",
  "#A3A3A3",
  "#858585",
  "#666666",
  "#454545"
)

col_fun <- circlize::colorRamp2(
  breaks,
  base_colors
)

# =========================
# DISPLAY NUMBERS
# =========================
display_mat <- matrix(
  "",
  nrow = nrow(mat),
  ncol = ncol(mat),
  dimnames = dimnames(mat)
)

target_row <- "Enriched"
target_col <- "Malignant Cell"

if (
  target_row %in% rownames(mat_raw) &&
  target_col %in% colnames(mat_raw)
) {
  display_mat[target_row, target_col] <- format(
    mat_raw[target_row, target_col],
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

# =========================
# CELL FUNCTION
# =========================
cell_fun <- function(j, i, x, y, width, height, fill) {
  
  is_target <-
    rownames(mat)[i] == target_row &&
    colnames(mat)[j] == target_col
  
  if (is_target) {
    grid.rect(
      x = x,
      y = y,
      width = width,
      height = height,
      gp = gpar(
        fill = "#9F3B4A",
        col = "#9F3B4A",
        lwd = 1
      )
    )
  }
  
  if (display_mat[i, j] != "") {
    grid.text(
      display_mat[i, j],
      x = x,
      y = y,
      gp = gpar(
        fontsize = if (is_target) 14 else 10,
        fontface = if (is_target) "bold" else "plain",
        fontfamily = "Arial",
        col = if (is_target) "white" else "black"
      )
    )
  }
}

# =========================
# HEATMAP
# =========================
ht <- Heatmap(
  mat,
  name = "Gene count",
  col = col_fun,
  
  cluster_rows = nrow(mat) > 1,
  cluster_columns = ncol(mat) > 1,
  
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  
  rect_gp = gpar(col = NA),
  
  width = unit(ncol(mat) * 12, "mm"),
  height = unit(65, "mm"),
  
  column_names_side = "bottom",
  row_names_side = "left",
  row_dend_side = "left",
  
  row_names_gp = gpar(
    fontsize = 11,
    fontfamily = "Arial"
  ),
  
  column_names_gp = gpar(
    fontsize = 11,
    fontfamily = "Arial"
  ),
  
  row_names_max_width = unit(100, "mm"),
  row_dend_width = unit(10, "mm"),
  column_dend_height = unit(5, "mm"),
  
  heatmap_legend_param = list(
    at = breaks,
    labels = c("0", "10", "100", "1k", "3k", "6k", "10k"),
    legend_width = unit(8, "mm"),
    legend_height = unit(5, "cm"),
    
    title_gp = gpar(
      fontsize = 10,
      fontface = "plain",
      fontfamily = "Arial"
    ),
    
    labels_gp = gpar(
      fontsize = 9,
      fontfamily = "Arial"
    ),
    
    title_position = "topcenter"
  ),
  
  cell_fun = cell_fun
)

# =========================
# EXPORT SVG
# =========================
output_dir <- dirname(opt$output)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

svg(
  filename = opt$output,
  width = 15,
  height = 5,
  bg = "transparent"
)

draw(
  ht,
  heatmap_legend_side = "right",
  background = "transparent",
  padding = unit(c(5, 20, 5, 20), "mm")
)

dev.off()

message("Specificity heatmap saved to: ", opt$output)
