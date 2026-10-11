
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(optparse)
  library(stringr)
  library(forcats)
  library(tools)
})

# =========================
# ARGUMENTS
# =========================
option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    help = "Input CSV file"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "specificity_groups_by_type_heatmap.svg",
    help = "Output SVG file"
  ),
  make_option(
    c("--merged_output"),
    type = "character",
    default = "merged_summed_genes.csv",
    help = "Output CSV with merged rows"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

# =========================
# VALIDATE INPUT
# =========================
if (is.null(opt$input) || !file.exists(opt$input)) {
  stop("Input CSV does not exist. Provide a valid --input path.")
}

# =========================
# LOAD DATA
# =========================
df <- read_csv(
  opt$input,
  show_col_types = FALSE,
  na = c("", "NA", "na", "Na", "N/A", "n/a")
)

required_columns <- c(
  "celltypist_labels",
  "type",
  "specificity_profile",
  "gene_name"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# =========================
# CLEAN LABELS
# =========================
df <- df %>%
  mutate(
    celltypist_labels = as.character(celltypist_labels),
    celltypist_labels = str_trim(celltypist_labels),
    celltypist_labels = replace_na(celltypist_labels, "Unknown"),
    celltypist_labels = if_else(
      celltypist_labels == "",
      "Unknown",
      celltypist_labels
    ),
    type = as.character(type),
    specificity_profile = as.character(specificity_profile),
    gene_name = as.character(gene_name)
  )

# =========================
# FORMAT CELL TYPE NAMES
# =========================
df <- df %>%
  mutate(
    celltypist_labels = case_when(
      celltypist_labels == "epithelial cell of thymus" ~ "TEC",
      celltypist_labels == "plasmacytoid dendritic cell, human" ~ "pDCs",
      celltypist_labels == "hematopoietic multipotent progenitor cell" ~ "MPPs",
      str_detect(str_to_lower(celltypist_labels), "cd4") ~
        str_to_upper(celltypist_labels),
      str_detect(str_to_lower(celltypist_labels), "cd8") ~
        str_to_upper(celltypist_labels),
      TRUE ~ tools::toTitleCase(celltypist_labels)
    ),
    type = case_when(
      str_to_lower(type) == "diferentiation" ~ "Differentiation",
      TRUE ~ tools::toTitleCase(type)
    )
  )

# =========================
# MERGE HETEROGENEOUS GROUPS
# =========================
df <- df %>%
  mutate(
    specificity_profile = case_when(
      specificity_profile %in% c(
        "Enriched",
        "Heterogeneous Enriched"
      ) ~ "Enriched",
      
      specificity_profile %in% c(
        "Enhanced",
        "Heterogeneous Enhanced"
      ) ~ "Enhanced",
      
      specificity_profile %in% c(
        "Group enriched",
        "Group Enriched",
        "Heterogeneous Group Enriched"
      ) ~ "Group Enriched",
      
      TRUE ~ specificity_profile
    )
  )

# =========================
# SAVE MERGED RAW DATA
# =========================
df_merged <- df %>%
  filter(
    specificity_profile %in% c(
      "Enriched",
      "Enhanced",
      "Group Enriched"
    )
  ) %>%
  filter(
    !str_detect(str_to_lower(type), "hla")
  )

if (nrow(df_merged) == 0) {
  stop(
    "No rows remain after filtering specificity profiles and HLA categories."
  )
}

# =========================
# OUTPUT DIRECTORIES
# =========================
dir.create(
  dirname(opt$output),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dirname(opt$merged_output),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  df_merged,
  opt$merged_output
)

# =========================
# HEATMAP SUMMARY
# =========================
heatmap_df <- df_merged %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  group_by(type, celltypist_labels) %>%
  summarise(
    total = n_distinct(gene_name),
    .groups = "drop"
  ) %>%
  complete(
    type,
    celltypist_labels,
    fill = list(total = 0)
  ) %>%
  mutate(
    type = fct_rev(factor(type))
  )

if (nrow(heatmap_df) == 0) {
  stop("No genes available to construct the heatmap.")
}

# =========================
# DEFINE TEXT COLOR
# =========================
dark_threshold <- quantile(
  heatmap_df$total,
  probs = 0.80,
  na.rm = TRUE
)

heatmap_df <- heatmap_df %>%
  mutate(
    text_color = if_else(
      total >= dark_threshold & total > 0,
      "white",
      "black"
    )
  )

# =========================
# PLOT
# =========================
p <- ggplot(
  heatmap_df,
  aes(
    x = celltypist_labels,
    y = type,
    fill = total
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(
      label = ifelse(total > 0, total, ""),
      color = text_color
    ),
    size = 3.5
  ) +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = c(
      "#F2F2F2",
      "#D9D9D9",
      "#BFBFBF",
      "#A3A3A3",
      "#858585",
      "#666666",
      "#454545"
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 12,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      colour = "black"
    ),
    axis.title.x = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      colour = "black"
    ),
    legend.title = element_text(
      size = 12,
      colour = "black"
    ),
    legend.text = element_text(
      colour = "black"
    ),
    panel.grid = element_blank(),
    plot.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "transparent",
      colour = NA
    )
  ) +
  labs(
    x = NULL,
    y = "Antigen category",
    fill = "Enriched gene categories*"
  )

# =========================
# EXPORT SVG
# =========================
ggsave(
  filename = opt$output,
  plot = p,
  device = "svg",
  width = 13,
  height = 4,
  units = "in",
  bg = "transparent"
)

message("Heatmap saved to: ", opt$output)
message("Merged CSV saved to: ", opt$merged_output)
