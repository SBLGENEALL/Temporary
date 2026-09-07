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
OUTPUT_DIR <- "FACS_plot_results_v2"
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
  HIGHLIGHT_CONSTRUCTS <- unique(c(ORIGINAL_NAME, auto_hits))
}

plot_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 5)),
    plot.subtitle = element_text(color = "#4D4D4D", size = 11, margin = margin(b = 10)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "#222222"),
    strip.background = element_rect(fill = "#EAF0F6", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    panel.spacing = grid::unit(1.1, "lines"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(12, 16, 12, 12)
  )

# Global construct order emphasizes sustained MFI performance across measured days.
construct_order <- summary_vs_original %>%
  filter(str_to_lower(Construct) != str_to_lower(ORIGINAL_NAME), !is.na(MFI_Pct_vs_Original)) %>%
  group_by(Construct) %>%
  summarise(Order_score = mean(MFI_Pct_vs_Original, na.rm = TRUE), .groups = "drop") %>%
  arrange(Order_score) %>%
  pull(Construct)
construct_order <- c(construct_order, ORIGINAL_NAME)

heatmap_dat <- summary_vs_original %>%
  mutate(Construct_plot = factor(Construct, levels = construct_order))

highlight_dat <- summary_vs_original %>%
  filter(Construct %in% HIGHLIGHT_CONSTRUCTS) %>%
  mutate(Construct = factor(Construct, levels = HIGHLIGHT_CONSTRUCTS))

replicate_highlight <- dat_vs_original %>%
  filter(Construct %in% HIGHLIGHT_CONSTRUCTS) %>%
  mutate(Construct = factor(Construct, levels = HIGHLIGHT_CONSTRUCTS))

candidate_levels <- setdiff(HIGHLIGHT_CONSTRUCTS, ORIGINAL_NAME)
candidate_colors <- if (length(candidate_levels) > 0) {
  setNames(grDevices::hcl.colors(length(candidate_levels), palette = "Dark 3"), candidate_levels)
} else {
  character(0)
}
highlight_colors <- c(setNames("#111111", ORIGINAL_NAME), candidate_colors)

# 1. Main overview: all constructs and all time points. Labels appear only for
# differences of at least 5%, keeping the heatmap readable.
p1 <- ggplot(heatmap_dat, aes(Day, Construct_plot, fill = MFI_Pct_vs_Original)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = ifelse(!is.na(MFI_Pct_vs_Original) & abs(MFI_Pct_vs_Original) >= 5,
                       sprintf("%+.0f", MFI_Pct_vs_Original), "")),
    size = 2.5, color = "#202020"
  ) +
  facet_wrap(~Condition, scales = "free_x") +
  scale_fill_gradient2(
    low = "#C44E52", mid = "#F7F7F7", high = "#2C7FB8", midpoint = 0,
    na.value = "#ECECEC", name = "% vs Original"
  ) +
  labs(
    title = "GFP GeoMean relative to Original",
    subtitle = "Blue: higher than Original | Red: lower than Original | cell labels show differences >=5%",
    x = NULL, y = NULL
  ) +
  plot_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.key.width = grid::unit(2.2, "cm")
  )
ggsave(file.path(OUTPUT_DIR, "01_MFI_heatmap_all_constructs.png"), p1, width = 12, height = 11, dpi = 320, bg = "white")

# Latest measured day is determined separately for each condition.
latest_by_condition <- summary_vs_original %>%
  filter(!is.na(MFI_Pct_vs_Original)) %>%
  group_by(Condition) %>%
  filter(Day_Number == max(Day_Number, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Rank_key = paste(Construct, Condition, sep = "___"),
    Rank_key = reorder(Rank_key, MFI_Pct_vs_Original)
  )

# 2. Clean ranked plot for the latest measured day in each condition.
p2 <- ggplot(latest_by_condition, aes(MFI_Pct_vs_Original, Rank_key)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6, color = "#555555") +
  geom_segment(aes(x = 0, xend = MFI_Pct_vs_Original, yend = Rank_key), color = "#D0D0D0", linewidth = 0.55) +
  geom_point(aes(color = MFI_Pct_vs_Original >= 0), size = 3) +
  facet_wrap(~Condition, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
  scale_color_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "#B8B8B8"), guide = "none") +
  labs(
    title = "Latest-day ranking",
    subtitle = "Horizontal distance from zero is the mean GFP GeoMean difference versus Original",
    x = "GFP GeoMean difference vs Original (%)", y = NULL
  ) +
  plot_theme +
  theme(panel.grid.major.x = element_line(color = "#ECECEC"), axis.text.y = element_text(size = 9))
ggsave(file.path(OUTPUT_DIR, "02_latest_day_MFI_ranking.png"), p2, width = 12, height = 10.5, dpi = 320, bg = "white")

# 3. Longitudinal plot is restricted to the strongest candidates plus Original.
p3 <- ggplot(highlight_dat, aes(Day_Number, MFI_Fold_vs_Original, color = Construct, group = Construct)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.7, color = "#555555") +
  geom_point(
    data = replicate_highlight,
    aes(Day_Number, MFI_Fold_vs_Original, color = Construct),
    inherit.aes = FALSE, alpha = 0.25, size = 1.8,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  geom_errorbar(
    aes(ymin = MFI_Fold_vs_Original - MFI_SEM_Fold,
        ymax = MFI_Fold_vs_Original + MFI_SEM_Fold),
    width = 0.13, linewidth = 0.45, alpha = 0.75
  ) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.8) +
  geom_line(
    data = highlight_dat %>% filter(as.character(Construct) == ORIGINAL_NAME),
    linewidth = 1.6, color = "#111111"
  ) +
  facet_wrap(~Condition) +
  scale_color_manual(values = highlight_colors, drop = FALSE) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11), labels = paste0("D", c(1, 3, 5, 7, 9, 11))) +
  labs(
    title = "Top-candidate MFI trajectories",
    subtitle = "Large points: triplicate mean | error bars: SEM | faint points: individual replicates",
    x = "Day", y = "GFP GeoMean fold vs Original", color = "Construct"
  ) +
  plot_theme +
  theme(panel.grid.major.y = element_line(color = "#ECECEC"))
ggsave(file.path(OUTPUT_DIR, "03_top_candidate_MFI_trajectories.png"), p3, width = 12, height = 7, dpi = 320, bg = "white")

# 4. Latest-day decision map separates intensity gain from GFP-positive fraction.
decision_dat <- latest_by_condition %>%
  mutate(
    Decision = case_when(
      MFI_Pct_vs_Original >= 0 & GFP_Positive_pp_vs_Original >= 0 ~ "Both higher",
      MFI_Pct_vs_Original >= 0 & GFP_Positive_pp_vs_Original < 0 ~ "MFI higher / GFP+ lower",
      MFI_Pct_vs_Original < 0 & GFP_Positive_pp_vs_Original >= 0 ~ "MFI lower / GFP+ higher",
      TRUE ~ "Both lower"
    ),
    Point_label = ifelse(Construct %in% HIGHLIGHT_CONSTRUCTS, Construct, "")
  )

p4 <- ggplot(decision_dat, aes(GFP_Positive_pp_vs_Original, MFI_Pct_vs_Original)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_hline(yintercept = 0, linetype = 2, color = "#777777") +
  geom_point(aes(color = Decision), size = 3.2, alpha = 0.9) +
  geom_text_repel(aes(label = Point_label), size = 3.2, max.overlaps = Inf, box.padding = 0.35, show.legend = FALSE) +
  facet_wrap(~Condition, scales = "free") +
  scale_color_manual(values = c(
    "Both higher" = "#167D5A",
    "MFI higher / GFP+ lower" = "#2C7FB8",
    "MFI lower / GFP+ higher" = "#D99A2B",
    "Both lower" = "#B8B8B8"
  )) +
  labs(
    title = "Latest-day expression trade-off",
    subtitle = "Upper-right quadrant identifies candidates exceeding Original in both intensity and GFP-positive fraction",
    x = "GFP-positive difference vs Original (percentage points)",
    y = "GFP GeoMean difference vs Original (%)", color = NULL
  ) +
  plot_theme +
  theme(panel.grid.major = element_line(color = "#EEEEEE"))
ggsave(file.path(OUTPUT_DIR, "04_latest_day_expression_tradeoff.png"), p4, width = 12, height = 7.5, dpi = 320, bg = "white")

# 5. GFP-positive fraction is shown as a difference from Original, not as 32 overlapping lines.
p5 <- ggplot(heatmap_dat, aes(Day, Construct_plot, fill = GFP_Positive_pp_vs_Original)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(
    aes(label = ifelse(!is.na(GFP_Positive_pp_vs_Original) & abs(GFP_Positive_pp_vs_Original) >= 3,
                       sprintf("%+.1f", GFP_Positive_pp_vs_Original), "")),
    size = 2.4, color = "#202020"
  ) +
  facet_wrap(~Condition, scales = "free_x") +
  scale_fill_gradient2(
    low = "#C44E52", mid = "#F7F7F7", high = "#2C7FB8", midpoint = 0,
    na.value = "#ECECEC", name = "pp vs Original"
  ) +
  labs(
    title = "GFP-positive fraction relative to Original",
    subtitle = "Values are percentage-point differences; labels are shown for absolute differences >=3 pp",
    x = NULL, y = NULL
  ) +
  plot_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    legend.key.width = grid::unit(2.2, "cm")
  )
ggsave(file.path(OUTPUT_DIR, "05_GFP_positive_heatmap_all_constructs.png"), p5, width = 12, height = 11, dpi = 320, bg = "white")

message("Analysis complete. Results written to: ", normalizePath(OUTPUT_DIR))
message("Highlighted constructs: ", paste(HIGHLIGHT_CONSTRUCTS, collapse = ", "))
