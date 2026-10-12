#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(ComplexUpset)
  library(ggplot2)
  library(dplyr)
})

# ==========================================================
# ARGUMENTS
# ==========================================================

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  
  while (i <= length(args)) {
    key <- args[i]
    
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    
    parsed[[substring(key, 3)]] <- args[i + 1]
    i <- i + 2
  }
  
  parsed
}

params <- parse_args(args)

if (is.null(params$input)) {
  stop("Required argument: --input")
}

if (is.null(params$output_dir)) {
  stop("Required argument: --output_dir")
}

input_file <- params$input
output_dir <- params$output_dir

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================================
# LOAD
# ==========================================================

message("Reading binary candidate matrix: ", input_file)

df <- read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c("Gene_ID", "Gene_Name")

missing_columns <- setdiff(required_columns, colnames(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# ==========================================================
# DEFINE SETS
# ==========================================================

sets <- setdiff(
  colnames(df),
  c("Gene_ID", "Gene_Name")
)

if (length(sets) == 0) {
  stop("No tumor columns were found in the input matrix.")
}

# Standardize tumor names for display.
sets_title <- tools::toTitleCase(tolower(sets))

# Avoid duplicated column names after title capitalization.
if (anyDuplicated(sets_title)) {
  stop(
    "Tumor names become duplicated after capitalization: ",
    paste(sets_title, collapse = ", ")
  )
}

colnames(df)[match(sets, colnames(df))] <- sets_title
sets <- sets_title

# Ensure that the tumor columns contain binary values.
for (tumor in sets) {
  df[[tumor]] <- as.numeric(df[[tumor]])
  
  if (any(
    !is.na(df[[tumor]]) &
    !df[[tumor]] %in% c(0, 1)
  )) {
    stop(
      "Column '", tumor,
      "' contains values other than 0 and 1."
    )
  }
  
  df[[tumor]][is.na(df[[tumor]])] <- 0
}

# ==========================================================
# CREATE GENE LIST BY TUMOR COMBINATION
# ==========================================================

df$Combination <- apply(
  as.matrix(df[, sets, drop = FALSE]),
  1,
  function(x) {
    paste(
      sets[x == 1],
      collapse = "_"
    )
  }
)

# Exclude genes that are not detected in any tumor.
df_combination <- df %>%
  filter(Combination != "")

gene_combination_table <- df_combination %>%
  group_by(Combination) %>%
  summarise(
    N_genes = n(),
    N_tumors = n_distinct(
      unlist(strsplit(Combination, "_", fixed = TRUE))
    ),
    Genes = paste(
      ifelse(
        is.na(Gene_Name) | Gene_Name == "",
        Gene_ID,
        Gene_Name
      ),
      collapse = "; "
    ),
    .groups = "drop"
  ) %>%
  arrange(desc(N_genes), desc(N_tumors))

summary_file <- file.path(
  output_dir,
  "genes_by_tumor_combination_summary.csv"
)

write_csv(
  gene_combination_table,
  summary_file
)

message("Saved: ", summary_file)

# ==========================================================
# OUTPUT PATH
# ==========================================================

plot_file <- file.path(
  output_dir,
  "upset_candidate_genes.svg"
)

# ==========================================================
# CREATE UPSET PLOT
# ==========================================================

# Convert the tumor columns to logical values for ComplexUpset.
df_plot <- df

for (tumor in sets) {
  df_plot[[tumor]] <- df_plot[[tumor]] == 1
}

p <- upset(
  df_plot,
  intersect = sets,
  min_size = 3,
  sort_intersections_by = "cardinality",
  width_ratio = 0.22,
  
  base_annotations = list(
    "Intersection size" =
      intersection_size(
        counts = TRUE,
        bar_number_threshold = 0.85,
        mapping = aes(fill = "Intersection")
      ) +
      scale_fill_manual(
        values = c("Intersection" = "#454545")
      ) +
      theme(
        legend.position = "none"
      )
  ),
  
  set_sizes =
    upset_set_size(
      mapping = aes(fill = "Set size")
    ) +
    scale_fill_manual(
      values = c("Set size" = "#454545")
    ) +
    theme(
      legend.position = "none"
    )
) +
  theme_minimal(
    base_size = 14,
    base_family = "Arial"
  ) +
  theme(
    panel.grid.major = element_line(
      colour = "#F5F6F7",
      linewidth = 0.3
    ),
    
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    
    text = element_text(
      family = "Arial",
      colour = "black",
      size = 11
    ),
    
    axis.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    axis.title = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    strip.text = element_text(
      family = "Arial",
      size = 11,
      colour = "black"
    ),
    
    plot.title = element_text(
      family = "Arial",
      size = 11,
      hjust = 0.5
    ),
    
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    
    panel.background = element_rect(
      fill = "#F5F6F7",
      colour = NA
    )
  )

# Save SVG with transparent outer background.
ggsave(
  filename = plot_file,
  plot = p,
  device = "svg",
  width = 10.7,
  height = 4.7,
  bg = "transparent"
)

message("Saved: ", plot_file)
message("Analysis completed successfully.")