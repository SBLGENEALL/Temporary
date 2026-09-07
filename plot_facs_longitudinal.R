# CHO 5'UTR transposase FACS longitudinal analysis
# Preferred input: Plate_Map.csv + FACS_Data.csv
# Secondary input: Plate_Map.tsv + FACS_Data.tsv
# Fallback input: CHO_5UTR_FACS_D1-D11_Input_Template.xlsx
# Required packages: dplyr, tidyr, ggplot2, stringr, ggrepel
# readxl is required only when the XLSX fallback is used.

CSV_PLATE_MAP_FILE <- "Plate_Map.csv"
CSV_FACS_DATA_FILE <- "FACS_Data.csv"
TSV_PLATE_MAP_FILE <- "Plate_Map.tsv"
TSV_FACS_DATA_FILE <- "FACS_Data.tsv"
XLSX_INPUT_FILE <- "CHO_5UTR_FACS_D1-D11_Input_Template.xlsx"

csv_available <- file.exists(CSV_PLATE_MAP_FILE) && file.exists(CSV_FACS_DATA_FILE)
tsv_available <- file.exists(TSV_PLATE_MAP_FILE) && file.exists(TSV_FACS_DATA_FILE)
xlsx_available <- file.exists(XLSX_INPUT_FILE)

if (csv_available) {
  INPUT_MODE <- "csv"
} else if (tsv_available) {
  INPUT_MODE <- "tsv"
} else if (xlsx_available) {
  INPUT_MODE <- "xlsx"
} else {
  stop(
    "Input files were not found. Put either:\n",
    "  1) Plate_Map.csv and FACS_Data.csv,\n",
    "  2) Plate_Map.tsv and FACS_Data.tsv, or\n",
    "  3) CHO_5UTR_FACS_D1-D11_Input_Template.xlsx\n",
    "in the same directory as this R script."
  )
}

required_packages <- c("dplyr", "tidyr", "ggplot2", "stringr", "ggrepel")
if (INPUT_MODE == "xlsx") required_packages <- c("readxl", required_packages)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ggrepel)
  if (INPUT_MODE == "xlsx") library(readxl)
})

# ----------------------------- User settings -----------------------------
OUTPUT_DIR <- "FACS_plot_results_v3_prism"
ORIGINAL_NAME <- "Original"
DAY_ORDER <- c("D1", "D3", "D5", "D7", "D9", "D11")
# Usually leave blank. If a legacy Korean Windows TSV is garbled, set "CP949".
TSV_ENCODING_OVERRIDE <- ""

# Leave empty for automatic selection of Original + top candidates at latest day.
# Example: HIGHLIGHT_CONSTRUCTS <- c("Original", "TOP3", "TOP11", "TOP24")
HIGHLIGHT_CONSTRUCTS <- character(0)
N_AUTO_HIGHLIGHT <- 7

# TRUE: use only GFP_GeoMean among GFP-positive cells, if that is how FlowJo exported it.
# Record the exact parent population in the workbook Notes column.
GFP_GEOMEAN_IS_WITHIN_POSITIVE_GATE <- TRUE
# ------------------------------------------------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

clean_numeric <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, ",", "")
  x <- str_replace_all(x, "%", "")
  suppressWarnings(as.numeric(x))
}

clean_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, c("POSIXct", "POSIXlt"))) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  suppressWarnings(as.Date(as.character(x)))
}

normalize_well <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  str_remove(x, regex("\\.fcs$", ignore_case = TRUE))
}

detect_tsv_encoding <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con))
  bom <- as.integer(readBin(con, what = "raw", n = 3))
  if (length(bom) >= 2 && identical(bom[1:2], c(255L, 254L))) return("UTF-16LE")
  if (length(bom) >= 2 && identical(bom[1:2], c(254L, 255L))) return("UTF-16BE")
  if (length(bom) >= 3 && identical(bom[1:3], c(239L, 187L, 191L))) return("UTF-8")
  "UTF-8"
}

read_tsv_flexible <- function(path) {
  encoding <- if (nzchar(TSV_ENCODING_OVERRIDE)) TSV_ENCODING_OVERRIDE else detect_tsv_encoding(path)
  x <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    na.strings = c("", "NA"),
    quote = "\"",
    comment.char = "",
    fileEncoding = encoding,
    stringsAsFactors = FALSE
  )
  names(x) <- sub(paste0("^", intToUtf8(65279)), "", names(x))
  x
}

read_csv_utf8 <- function(path) {
  x <- read.csv(
    path,
    header = TRUE,
    check.names = FALSE,
    na.strings = c("", "NA"),
    quote = "\"",
    comment.char = "",
    fileEncoding = "UTF-8",
    stringsAsFactors = FALSE
  )
  names(x) <- sub(paste0("^", intToUtf8(65279)), "", names(x))
  x
}

if (INPUT_MODE == "csv") {
  plate_map_input <- read_csv_utf8(CSV_PLATE_MAP_FILE)
  facs_input <- read_csv_utf8(CSV_FACS_DATA_FILE)
  message("Input mode: CSV")
} else if (INPUT_MODE == "tsv") {
  plate_map_input <- read_tsv_flexible(TSV_PLATE_MAP_FILE)
  facs_input <- read_tsv_flexible(TSV_FACS_DATA_FILE)
  message("Input mode: TSV")
} else {
  plate_map_input <- read_excel(XLSX_INPUT_FILE, sheet = "Plate_Map")
  facs_input <- read_excel(XLSX_INPUT_FILE, sheet = "FACS_Data")
  message("Input mode: XLSX")
}

require_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

require_columns(
  plate_map_input,
  c("Well", "Sample_File", "Construct", "Replicate", "DNA_Prep", "DNA_Batch", "Notes"),
  "Plate map input"
)
require_columns(
  facs_input,
  c("Day", "Condition", "Acquisition_Date", "Well", "Sample_File",
    "Live_Cell_Pct", "GFP_Positive_Pct", "GFP_GeoMean", "Notes"),
  "FACS data input"
)

plate_map <- plate_map_input %>%
  transmute(
    Well = normalize_well(Well),
    Sample_File_Map = as.character(Sample_File),
    Construct = str_trim(as.character(Construct)),
    Replicate = suppressWarnings(as.integer(Replicate)),
    DNA_Prep = na_if(str_trim(as.character(DNA_Prep)), ""),
    DNA_Batch = na_if(str_trim(as.character(DNA_Batch)), ""),
    Plate_Notes = na_if(as.character(Notes), "")
  ) %>%
  filter(!is.na(Well), Well != "")

facs_raw <- facs_input %>%
  transmute(
    Day = toupper(str_trim(as.character(Day))),
    Condition = na_if(str_trim(as.character(Condition)), ""),
    Acquisition_Date = clean_date(Acquisition_Date),
    Well = normalize_well(if_else(
      is.na(Well) | str_trim(as.character(Well)) == "",
      as.character(Sample_File),
      as.character(Well)
    )),
    Sample_File = as.character(Sample_File),
    Live_Cell_Pct = clean_numeric(Live_Cell_Pct),
    GFP_Positive_Pct = clean_numeric(GFP_Positive_Pct),
    GFP_GeoMean = clean_numeric(GFP_GeoMean),
    FACS_Notes = na_if(as.character(Notes), "")
  ) %>%
  mutate(Condition = replace_na(Condition, "Unspecified")) %>%
  filter(if_any(c(Live_Cell_Pct, GFP_Positive_Pct, GFP_GeoMean), ~ !is.na(.x)))

if (nrow(facs_raw) == 0) {
  stop("No numeric FACS values were found in the FACS_Data sheet.")
}

for (pct_col in c("Live_Cell_Pct", "GFP_Positive_Pct")) {
  observed <- facs_raw[[pct_col]][!is.na(facs_raw[[pct_col]])]
  if (length(observed) > 0 && max(observed) <= 1) {
    warning(pct_col, " values are all <= 1. Enter 50% as 50, not 0.50.")
  }
}

unknown_days <- setdiff(unique(facs_raw$Day), DAY_ORDER)
if (length(unknown_days) > 0) {
  stop("Unexpected Day values: ", paste(unknown_days, collapse = ", "))
}

dat <- facs_raw %>%
  left_join(plate_map, by = "Well") %>%
  mutate(
    Day = factor(Day, levels = DAY_ORDER, ordered = TRUE),
    Day_Number = as.integer(str_remove(as.character(Day), "D")),
    Construct = na_if(Construct, ""),
    Construct = if_else(is.na(Construct), paste0("UNMAPPED_", Well), Construct)
  )

# Preserve the first-appearance order from Plate_Map exactly as entered by the user.
# Repeated wells/replicates do not create duplicate axis entries.
input_construct_order <- plate_map %>%
  filter(!is.na(Construct), Construct != "") %>%
  pull(Construct) %>%
  unique()
input_construct_order <- c(
  input_construct_order,
  setdiff(unique(as.character(dat$Construct)), input_construct_order)
)

if (any(str_starts(dat$Construct, "UNMAPPED_"))) {
  warning("Some wells have no Construct in Plate_Map. They are retained as UNMAPPED_<well>.")
}

duplicate_rows <- dat %>%
  count(Day, Condition, Well, name = "n_rows") %>%
  filter(n_rows > 1)
if (nrow(duplicate_rows) > 0) {
  warning("Duplicate Day/Condition/Well rows detected. Check qc_duplicate_rows.csv.")
  write.csv(duplicate_rows, file.path(OUTPUT_DIR, "qc_duplicate_rows.csv"), row.names = FALSE)
}

summary_dat <- dat %>%
  group_by(Day, Day_Number, Condition, Construct, DNA_Prep) %>%
  summarise(
    n = sum(!is.na(GFP_GeoMean)),
    Live_Cell_mean = mean(Live_Cell_Pct, na.rm = TRUE),
    Live_Cell_sd = sd(Live_Cell_Pct, na.rm = TRUE),
    Live_Cell_n = sum(!is.na(Live_Cell_Pct)),
    Live_Cell_sem = Live_Cell_sd / sqrt(Live_Cell_n),
    GFP_Positive_mean = mean(GFP_Positive_Pct, na.rm = TRUE),
    GFP_Positive_sd = sd(GFP_Positive_Pct, na.rm = TRUE),
    GFP_Positive_n = sum(!is.na(GFP_Positive_Pct)),
    GFP_Positive_sem = GFP_Positive_sd / sqrt(GFP_Positive_n),
    GFP_GeoMean_mean = mean(GFP_GeoMean, na.rm = TRUE),
    GFP_GeoMean_sd = sd(GFP_GeoMean, na.rm = TRUE),
    GFP_GeoMean_sem = GFP_GeoMean_sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_mean"), ~ ifelse(is.nan(.x), NA_real_, .x)))

original_summary <- summary_dat %>%
  filter(str_to_lower(Construct) == str_to_lower(ORIGINAL_NAME)) %>%
  select(
    Day, Condition,
    Original_Live_Cell = Live_Cell_mean,
    Original_GFP_Positive = GFP_Positive_mean,
    Original_GFP_GeoMean = GFP_GeoMean_mean
  )

if (nrow(original_summary) == 0) {
  stop("Original control was not found. Check Plate_Map Construct and ORIGINAL_NAME.")
}

summary_vs_original <- summary_dat %>%
  left_join(original_summary, by = c("Day", "Condition")) %>%
  mutate(
    MFI_Fold_vs_Original = GFP_GeoMean_mean / Original_GFP_GeoMean,
    MFI_Pct_vs_Original = 100 * (MFI_Fold_vs_Original - 1),
    MFI_SEM_Fold = GFP_GeoMean_sem / Original_GFP_GeoMean,
    GFP_Positive_pp_vs_Original = GFP_Positive_mean - Original_GFP_Positive,
    Live_Cell_pp_vs_Original = Live_Cell_mean - Original_Live_Cell
  )

dat_vs_original <- dat %>%
  left_join(original_summary, by = c("Day", "Condition")) %>%
  mutate(MFI_Fold_vs_Original = GFP_GeoMean / Original_GFP_GeoMean)

write.csv(dat, file.path(OUTPUT_DIR, "facs_cleaned_long_data.csv"), row.names = FALSE)
write.csv(summary_dat, file.path(OUTPUT_DIR, "facs_summary_mean_sd_sem.csv"), row.names = FALSE)
write.csv(summary_vs_original, file.path(OUTPUT_DIR, "facs_summary_vs_original.csv"), row.names = FALSE)

latest_day <- max(summary_vs_original$Day_Number[!is.na(summary_vs_original$GFP_GeoMean_mean)], na.rm = TRUE)
if (length(HIGHLIGHT_CONSTRUCTS) == 0) {
  auto_hits <- summary_vs_original %>%
    filter(Day_Number == latest_day, !is.na(MFI_Fold_vs_Original)) %>%
    arrange(desc(MFI_Fold_vs_Original)) %>%
    distinct(Construct) %>%
    slice_head(n = N_AUTO_HIGHLIGHT) %>%
    pull(Construct)
  selected_highlights <- unique(c(ORIGINAL_NAME, auto_hits))
  HIGHLIGHT_CONSTRUCTS <- input_construct_order[input_construct_order %in% selected_highlights]
}

prism_theme <- theme_classic(base_size = 13) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.7),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, color = "#444444", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.margin = margin(12, 14, 12, 12)
  )

# Keep the exact first-appearance order from Plate_Map; do not sort TOP numbers.
construct_order <- input_construct_order
candidate_order <- construct_order[
  str_to_lower(construct_order) != str_to_lower(ORIGINAL_NAME)
]

summary_plot <- summary_vs_original %>%
  mutate(
    Construct = factor(Construct, levels = construct_order),
    Bar_group = ifelse(str_to_lower(as.character(Construct)) == str_to_lower(ORIGINAL_NAME),
                       "Original", "Candidate")
  )

raw_plot <- dat_vs_original %>%
  mutate(
    Construct = factor(Construct, levels = construct_order),
    Bar_group = ifelse(str_to_lower(as.character(Construct)) == str_to_lower(ORIGINAL_NAME),
                       "Original", "Candidate")
  )

bar_dir_mfi <- file.path(OUTPUT_DIR, "01_MFI_barplots_by_day")
bar_dir_gfp <- file.path(OUTPUT_DIR, "02_GFP_positive_barplots_by_day")
dir.create(bar_dir_mfi, showWarnings = FALSE, recursive = TRUE)
dir.create(bar_dir_gfp, showWarnings = FALSE, recursive = TRUE)

safe_filename <- function(x) str_replace_all(as.character(x), "[^A-Za-z0-9_-]+", "_")

# Prism-style bar plots: mean +/- SD plus all individual replicate points.
available_panels <- summary_plot %>%
  filter(!is.na(GFP_GeoMean_mean) | !is.na(GFP_Positive_mean)) %>%
  distinct(Day, Day_Number, Condition) %>%
  arrange(Day_Number, Condition)

for (i in seq_len(nrow(available_panels))) {
  one_day <- available_panels$Day[i]
  one_condition <- available_panels$Condition[i]

  sum_sub <- summary_plot %>% filter(Day == one_day, Condition == one_condition)
  raw_sub <- raw_plot %>% filter(Day == one_day, Condition == one_condition)
  file_tag <- paste0(safe_filename(one_day), "_", safe_filename(one_condition))

  p_mfi_bar <- ggplot(sum_sub, aes(Construct, GFP_GeoMean_mean, fill = Bar_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.45) +
    geom_errorbar(
      aes(ymin = pmax(0, GFP_GeoMean_mean - GFP_GeoMean_sd),
          ymax = GFP_GeoMean_mean + GFP_GeoMean_sd),
      width = 0.22, linewidth = 0.65
    ) +
    geom_point(
      data = raw_sub,
      aes(Construct, GFP_GeoMean),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.11, height = 0),
      shape = 21, size = 2.1, stroke = 0.45, fill = "white", color = "black"
    ) +
    scale_fill_manual(values = c("Candidate" = "#78A9DC", "Original" = "#333333"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
    labs(
      title = paste0(one_day, " | ", one_condition, " | GFP GeoMean"),
      subtitle = "Bars: mean | error bars: SD | points: individual replicates",
      x = NULL, y = "GFP GeoMean"
    ) +
    prism_theme +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9))
  ggsave(file.path(bar_dir_mfi, paste0(file_tag, "_MFI.png")), p_mfi_bar,
         width = 15, height = 7, dpi = 320, bg = "white")

  p_gfp_bar <- ggplot(sum_sub, aes(Construct, GFP_Positive_mean, fill = Bar_group)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.45) +
    geom_errorbar(
      aes(ymin = pmax(0, GFP_Positive_mean - GFP_Positive_sd),
          ymax = pmin(100, GFP_Positive_mean + GFP_Positive_sd)),
      width = 0.22, linewidth = 0.65
    ) +
    geom_point(
      data = raw_sub,
      aes(Construct, GFP_Positive_Pct),
      inherit.aes = FALSE,
      position = position_jitter(width = 0.11, height = 0),
      shape = 21, size = 2.1, stroke = 0.45, fill = "white", color = "black"
    ) +
    scale_fill_manual(values = c("Candidate" = "#78A9DC", "Original" = "#333333"), guide = "none") +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste0(one_day, " | ", one_condition, " | GFP-positive population"),
      subtitle = "Bars: mean | error bars: SD | points: individual replicates",
      x = NULL, y = "GFP-positive cells (%)"
    ) +
    prism_theme +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9))
  ggsave(file.path(bar_dir_gfp, paste0(file_tag, "_GFP_positive.png")), p_gfp_bar,
         width = 15, height = 7, dpi = 320, bg = "white")
}

# Clear longitudinal plots: Original plus automatically selected top candidates.
highlight_levels <- HIGHLIGHT_CONSTRUCTS
highlight_summary <- summary_vs_original %>%
  filter(Construct %in% highlight_levels) %>%
  mutate(Construct = factor(Construct, levels = highlight_levels))
highlight_raw <- dat_vs_original %>%
  filter(Construct %in% highlight_levels) %>%
  mutate(Construct = factor(Construct, levels = highlight_levels))

candidate_levels <- setdiff(highlight_levels, ORIGINAL_NAME)
candidate_colors <- if (length(candidate_levels) > 0) {
  setNames(grDevices::hcl.colors(length(candidate_levels), palette = "Dark 3"), candidate_levels)
} else {
  character(0)
}
line_colors <- c(setNames("#111111", ORIGINAL_NAME), candidate_colors)

p_mfi_line <- ggplot(
  highlight_summary,
  aes(Day_Number, MFI_Fold_vs_Original, color = Construct, group = Construct)
) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.7, color = "#555555") +
  geom_errorbar(
    aes(ymin = MFI_Fold_vs_Original - MFI_SEM_Fold,
        ymax = MFI_Fold_vs_Original + MFI_SEM_Fold),
    width = 0.15, linewidth = 0.55
  ) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3) +
  geom_point(
    data = highlight_raw,
    aes(Day_Number, MFI_Fold_vs_Original, color = Construct),
    inherit.aes = FALSE, alpha = 0.28, size = 1.7,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  facet_wrap(~Condition) +
  scale_color_manual(values = line_colors, drop = FALSE) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = paste0("D", c(1, 3, 5, 7, 9, 11))) +
  labs(
    title = "GFP GeoMean over time",
    subtitle = "Mean +/- SEM; values are normalized to Original within each day and condition",
    x = "Day", y = "GFP GeoMean fold vs Original", color = NULL
  ) +
  prism_theme
ggsave(file.path(OUTPUT_DIR, "03_MFI_lineplot_top_candidates.png"), p_mfi_line,
       width = 12, height = 7, dpi = 320, bg = "white")

p_gfp_line <- ggplot(
  highlight_summary,
  aes(Day_Number, GFP_Positive_mean, color = Construct, group = Construct)
) +
  geom_errorbar(
    aes(ymin = pmax(0, GFP_Positive_mean - GFP_Positive_sem),
        ymax = pmin(100, GFP_Positive_mean + GFP_Positive_sem)),
    width = 0.15, linewidth = 0.55
  ) +
  geom_line(linewidth = 1.15) +
  geom_point(size = 3) +
  geom_point(
    data = dat %>% filter(Construct %in% highlight_levels),
    aes(Day_Number, GFP_Positive_Pct, color = Construct),
    inherit.aes = FALSE, alpha = 0.28, size = 1.7,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  facet_wrap(~Condition) +
  scale_color_manual(values = line_colors, drop = FALSE) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = paste0("D", c(1, 3, 5, 7, 9, 11))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "GFP-positive population over time",
    subtitle = "Mean +/- SEM; faint points are individual replicates",
    x = "Day", y = "GFP-positive cells (%)", color = NULL
  ) +
  prism_theme
ggsave(file.path(OUTPUT_DIR, "04_GFP_positive_lineplot_top_candidates.png"), p_gfp_line,
       width = 12, height = 7, dpi = 320, bg = "white")

# One small panel per construct allows every candidate trajectory to be inspected
# without assigning 32 competing colors to one graph.
condition_levels <- unique(as.character(summary_vs_original$Condition))
condition_colors <- setNames(
  grDevices::hcl.colors(length(condition_levels), palette = "Dark 3"),
  condition_levels
)

p_all_small <- summary_vs_original %>%
  filter(str_to_lower(Construct) != str_to_lower(ORIGINAL_NAME)) %>%
  mutate(Construct = factor(Construct, levels = candidate_order)) %>%
  ggplot(aes(Day_Number, MFI_Fold_vs_Original, color = Condition, group = Condition)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.45, color = "#555555") +
  geom_errorbar(
    aes(ymin = MFI_Fold_vs_Original - MFI_SEM_Fold,
        ymax = MFI_Fold_vs_Original + MFI_SEM_Fold),
    width = 0.18, linewidth = 0.35
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~Construct, ncol = 4) +
  scale_color_manual(values = condition_colors) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = paste0("D", c(1, 3, 5, 7, 9, 11))) +
  labs(
    title = "MFI trajectory for every candidate",
    subtitle = "Dashed line = Original; mean +/- SEM",
    x = "Day", y = "Fold vs Original"
  ) +
  prism_theme +
  theme(
    strip.text = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.spacing = grid::unit(0.8, "lines")
  )
ggsave(file.path(OUTPUT_DIR, "05_MFI_lineplot_all_candidates_small_multiples.png"), p_all_small,
       width = 13, height = 20, dpi = 320, bg = "white")

message("Analysis complete. Results written to: ", normalizePath(OUTPUT_DIR))
message("Highlighted constructs: ", paste(HIGHLIGHT_CONSTRUCTS, collapse = ", "))
