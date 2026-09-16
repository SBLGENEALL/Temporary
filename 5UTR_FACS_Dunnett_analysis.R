#!/usr/bin/env Rscript

# 5′UTR FACS analysis
# Primary inference: two-sided Welch t-tests for each variant versus Original,
# followed by Holm family-wise error correction within each Day × Selection × Metric condition.
# Dunnett-adjusted results are retained as a secondary comparison.

required_packages <- c("dplyr", "tidyr", "ggplot2", "multcomp")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ", paste(missing_packages, collapse = ", "), "\n",
      "Install them before running the script."
    ),
    call. = FALSE
  )
}

get_argument <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop(paste("Missing value after", flag), call. = FALSE)
  args[[index + 1]]
}

has_flag <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)

input_file <- get_argument("--input", "5UTR_FACS_Data_Input.csv")
input_sheet <- get_argument("--sheet", "Data_Input")
output_dir <- get_argument("--output-dir", "5UTR_FACS_analysis_results")
control_name <- get_argument("--control", "Original")
transform_setting <- get_argument("--transform", "auto")
show_ns <- has_flag("--show-ns")

if (!transform_setting %in% c("auto", "none", "log2")) {
  stop("--transform must be auto, none, or log2", call. = FALSE)
}

replicate_columns <- c("Rep1", "Rep2", "Rep3")
key_columns <- c("Day", "Selection", "Metric", "Construct")

significance_symbol <- function(p_value) {
  dplyr::case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

choose_transform <- function(metric, setting) {
  if (setting != "auto") return(setting)
  if (grepl("MFI|intensity|fluorescence", metric, ignore.case = TRUE)) "log2" else "none"
}

natural_order <- function(values) {
  numeric_part <- suppressWarnings(as.numeric(gsub("[^0-9.-]", "", values)))
  order(is.na(numeric_part), numeric_part, values)
}

safe_filename <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", trimws(as.character(value)))
  value <- gsub("^_+|_+$", "", value)
  ifelse(nchar(value) == 0, "unnamed", value)
}

write_csv_bom <- function(data, path) {
  connection <- file(path, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw("\xEF\xBB\xBF"), connection)
  write.table(
    data,
    file = connection,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = TRUE,
    na = ""
  )
}

input_extension <- tolower(tools::file_ext(input_file))
raw_wide <- switch(
  input_extension,
  "tsv" = utils::read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8"
  ),
  "txt" = utils::read.delim(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8"
  ),
  "csv" = utils::read.csv(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA"),
    fileEncoding = "UTF-8"
  ),
  "xlsx" = {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("The readxl package is required only for XLSX input. Use the supplied CSV file instead.", call. = FALSE)
    }
    readxl::read_excel(input_file, sheet = input_sheet)
  },
  "xls" = {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("The readxl package is required only for XLS input. Use the supplied CSV file instead.", call. = FALSE)
    }
    readxl::read_excel(input_file, sheet = input_sheet)
  },
  stop("Unsupported input format. Use .tsv, .txt, .csv, .xlsx, or .xls", call. = FALSE)
)
names(raw_wide) <- sub("^\ufeff", "", trimws(names(raw_wide)))
missing_columns <- setdiff(c(key_columns, replicate_columns), names(raw_wide))
if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")), call. = FALSE)
}

raw_wide <- raw_wide |>
  dplyr::mutate(
    dplyr::across(dplyr::all_of(key_columns), ~ trimws(as.character(.x))),
    dplyr::across(dplyr::all_of(replicate_columns), ~ suppressWarnings(as.numeric(.x))),
    Input_order = dplyr::row_number()
  ) |>
  dplyr::filter(dplyr::if_any(dplyr::all_of(replicate_columns), ~ !is.na(.x)))

if (nrow(raw_wide) == 0) stop("No values were found in Rep1–Rep3", call. = FALSE)
if (any(is.na(raw_wide[key_columns])) || any(raw_wide[key_columns] == "")) {
  stop("Every row containing measurements must have Day, Selection, Metric, and Construct", call. = FALSE)
}

duplicates <- raw_wide |>
  dplyr::group_by(dplyr::across(dplyr::all_of(key_columns))) |>
  dplyr::summarise(Rows = dplyr::n(), .groups = "drop") |>
  dplyr::filter(Rows > 1)
if (nrow(duplicates) > 0) {
  print(duplicates)
  stop("Duplicate Day × Selection × Metric × Construct rows were detected", call. = FALSE)
}

raw_long <- raw_wide |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(replicate_columns),
    names_to = "Replicate",
    values_to = "Value",
    values_drop_na = TRUE
  )

condition_keys <- raw_long |>
  dplyr::distinct(Day, Selection, Metric)

all_results <- list()
all_diagnostics <- list()
warning_messages <- character()

for (condition_index in seq_len(nrow(condition_keys))) {
  current_key <- condition_keys[condition_index, ]
  condition <- raw_long |>
    dplyr::filter(
      Day == current_key$Day,
      Selection == current_key$Selection,
      Metric == current_key$Metric
    )

  condition_label <- paste(current_key$Day, current_key$Selection, current_key$Metric, sep = " | ")
  available_constructs <- unique(condition$Construct[order(condition$Input_order)])
  if (!control_name %in% available_constructs) {
    warning_messages <- c(
      warning_messages,
      paste0("Skipped ", condition_label, ": control '", control_name, "' is missing")
    )
    next
  }

  replicate_counts <- condition |>
    dplyr::count(Construct, name = "n")
  insufficient <- replicate_counts$Construct[replicate_counts$n < 2]
  if (length(insufficient) > 0) {
    warning_messages <- c(
      warning_messages,
      paste0(
        "Excluded from ", condition_label, " because n<2: ",
        paste(insufficient, collapse = ", ")
      )
    )
    condition <- condition |>
      dplyr::filter(!Construct %in% insufficient)
    available_constructs <- setdiff(available_constructs, insufficient)
  }

  variants <- setdiff(available_constructs, control_name)
  if (!control_name %in% available_constructs || length(variants) == 0) {
    warning_messages <- c(warning_messages, paste0("Skipped ", condition_label, ": insufficient groups"))
    next
  }

  transform_method <- choose_transform(current_key$Metric, transform_setting)
  if (transform_method == "log2" && any(condition$Value <= 0)) {
    warning_messages <- c(
      warning_messages,
      paste0("Skipped ", condition_label, ": log2 transformation requires values > 0")
    )
    next
  }
  condition <- condition |>
    dplyr::mutate(
      Test_value = if (transform_method == "log2") log2(Value) else Value,
      Construct = factor(Construct, levels = c(control_name, variants))
    )

  model <- stats::aov(Test_value ~ Construct, data = condition)
  dunnett_model <- multcomp::glht(
    model,
    linfct = multcomp::mcp(Construct = "Dunnett"),
    alternative = "two.sided"
  )
  dunnett_summary <- summary(dunnett_model)
  adjusted_p <- as.numeric(dunnett_summary$test$pvalues)
  test_statistics <- as.numeric(dunnett_summary$test$tstat)

  group_summary <- condition |>
    dplyr::group_by(Construct) |>
    dplyr::summarise(
      Input_order = min(Input_order),
      n = dplyr::n(),
      Mean = mean(Value),
      SD = stats::sd(Value),
      .groups = "drop"
    ) |>
    dplyr::mutate(Construct = as.character(Construct))

  control_mean <- group_summary$Mean[group_summary$Construct == control_name]
  variant_results <- data.frame(
    Construct = variants,
    Dunnett_statistic = test_statistics,
    Dunnett_p_adjusted = adjusted_p,
    stringsAsFactors = FALSE
  )

  control_test_values <- condition$Test_value[condition$Construct == control_name]
  welch_raw <- vapply(
    variants,
    function(variant) {
      stats::t.test(
        condition$Test_value[condition$Construct == variant],
        control_test_values,
        var.equal = FALSE,
        alternative = "two.sided"
      )$p.value
    },
    numeric(1)
  )
  variant_results$Welch_p_raw <- welch_raw
  variant_results$Welch_p_Holm <- stats::p.adjust(welch_raw, method = "holm")

  condition_results <- group_summary |>
    dplyr::left_join(variant_results, by = "Construct") |>
    dplyr::mutate(
      Day = current_key$Day,
      Selection = current_key$Selection,
      Metric = current_key$Metric,
      Transform = transform_method,
      Fold_vs_Original = Mean / control_mean,
      Dunnett_significance = ifelse(
        Construct == control_name,
        "control",
        significance_symbol(Dunnett_p_adjusted)
      ),
      Welch_Holm_significance = ifelse(
        Construct == control_name,
        "control",
        significance_symbol(Welch_p_Holm)
      ),
      Primary_method = "Welch t-test + Holm correction",
      Primary_p_adjusted = Welch_p_Holm,
      Primary_significance = Welch_Holm_significance
    ) |>
    dplyr::select(
      Day, Selection, Metric, Construct, Input_order, n, Mean, SD,
      Fold_vs_Original, Transform, Dunnett_statistic, Dunnett_p_adjusted,
      Dunnett_significance, Welch_p_raw, Welch_p_Holm, Welch_Holm_significance,
      Primary_method, Primary_p_adjusted, Primary_significance
    )

  anova_table <- summary(model)[[1]]
  anova_p <- anova_table[["Pr(>F)"]][1]
  fligner_p <- stats::fligner.test(Test_value ~ Construct, data = condition)$p.value
  diagnostic <- data.frame(
    Day = current_key$Day,
    Selection = current_key$Selection,
    Metric = current_key$Metric,
    Transform = transform_method,
    Groups_analyzed = length(available_constructs),
    ANOVA_F = anova_table[["F value"]][1],
    ANOVA_p = anova_p,
    Fligner_Killeen_p = fligner_p,
    Variance_warning = ifelse(
      fligner_p < 0.05,
      "Unequal variance signal: Welch-Holm is the primary analysis; interpret n=3 cautiously",
      ""
    ),
    stringsAsFactors = FALSE
  )

  all_results[[length(all_results) + 1]] <- condition_results
  all_diagnostics[[length(all_diagnostics) + 1]] <- diagnostic
}

if (length(all_results) == 0) {
  stop("No condition could be analyzed. Check the data and Original label.", call. = FALSE)
}

results <- dplyr::bind_rows(all_results) |>
  dplyr::arrange(Metric, Selection, Input_order)
diagnostics <- dplyr::bind_rows(all_diagnostics)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write_csv_bom(results, file.path(output_dir, "Welch_Holm_results.csv"))
write_csv_bom(results, file.path(output_dir, "Dunnett_results.csv"))
write_csv_bom(diagnostics, file.path(output_dir, "Condition_diagnostics.csv"))

plot_data <- raw_long |>
  dplyr::inner_join(
    results |>
      dplyr::select(
        Day, Selection, Metric, Construct, Input_order, Mean, SD,
        Significance = Primary_significance
      ),
    by = c("Day", "Selection", "Metric", "Construct", "Input_order")
  )

make_condition_plot <- function(condition_data, condition_results, title_text) {
  construct_levels <- condition_results$Construct[order(condition_results$Input_order)]
  condition_data <- condition_data |>
    dplyr::mutate(Construct = factor(Construct, levels = construct_levels))
  condition_results <- condition_results |>
    dplyr::mutate(
      Construct = factor(Construct, levels = construct_levels),
      Label = dplyr::case_when(
        Construct == control_name ~ "",
        Primary_significance == "ns" & !show_ns ~ "",
        TRUE ~ Primary_significance
      )
    )

  maximum_value <- max(condition_data$Value, condition_results$Mean + condition_results$SD, na.rm = TRUE)
  minimum_value <- min(condition_data$Value, na.rm = TRUE)
  data_span <- max(maximum_value - min(0, minimum_value), maximum_value * 0.08, 1e-9)
  label_offset <- data_span * 0.045
  condition_results <- condition_results |>
    dplyr::mutate(Label_y = Mean + SD + label_offset)
  upper_limit <- max(maximum_value, condition_results$Label_y[condition_results$Label != ""], na.rm = TRUE) + data_span * 0.10

  ggplot2::ggplot(condition_results, ggplot2::aes(x = Construct, y = Mean, fill = Construct == control_name)) +
    ggplot2::geom_col(width = 0.72, alpha = 0.92) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
      width = 0.22,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      data = condition_data,
      ggplot2::aes(x = Construct, y = Value),
      inherit.aes = FALSE,
      position = ggplot2::position_jitter(width = 0.13, height = 0, seed = 260915),
      size = 1.45,
      color = "#202020",
      alpha = 0.78
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = Label_y, label = Label),
      size = 3.0,
      vjust = 0,
      color = "#202020"
    ) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#7F7F7F", `FALSE` = "#4472C4"), guide = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, upper_limit),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = "Mean ± SD with individual transfection wells; Welch t-tests with Holm-adjusted p-values vs Original",
      x = NULL,
      y = as.character(condition_results$Metric[[1]])
    ) +
    ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#555555"),
      axis.text.x = ggplot2::element_text(angle = 67, hjust = 1, vjust = 1, size = 7.5),
      axis.line = ggplot2::element_line(linewidth = 0.5, color = "#333333"),
      panel.grid.major.y = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

condition_combinations <- results |>
  dplyr::distinct(Day, Selection, Metric)

for (index in seq_len(nrow(condition_combinations))) {
  current <- condition_combinations[index, ]
  condition_results <- results |>
    dplyr::filter(
      Day == current$Day,
      Selection == current$Selection,
      Metric == current$Metric
    ) |>
    dplyr::arrange(Input_order)
  condition_data <- plot_data |>
    dplyr::filter(
      Day == current$Day,
      Selection == current$Selection,
      Metric == current$Metric
    )
  plot <- make_condition_plot(
    condition_data,
    condition_results,
    paste0(current$Metric, ": ", current$Day, ", ", current$Selection)
  )
  filename <- paste(
    safe_filename(current$Metric),
    safe_filename(current$Day),
    safe_filename(current$Selection),
    sep = "__"
  )
  ggplot2::ggsave(
    filename = file.path(figure_dir, paste0(filename, ".png")),
    plot = plot,
    width = 18,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
}

summary_combinations <- results |>
  dplyr::distinct(Metric, Selection)

for (index in seq_len(nrow(summary_combinations))) {
  current <- summary_combinations[index, ]
  combined_results <- results |>
    dplyr::filter(Metric == current$Metric, Selection == current$Selection)
  combined_data <- plot_data |>
    dplyr::filter(Metric == current$Metric, Selection == current$Selection)

  day_levels <- unique(combined_results$Day)
  day_levels <- day_levels[natural_order(day_levels)]
  construct_levels <- combined_results |>
    dplyr::arrange(Input_order) |>
    dplyr::distinct(Construct) |>
    dplyr::pull(Construct)
  combined_results <- combined_results |>
    dplyr::mutate(
      Day = factor(Day, levels = day_levels),
      Construct = factor(Construct, levels = construct_levels),
      Label = dplyr::case_when(
        as.character(Construct) == control_name ~ "",
        Primary_significance == "ns" & !show_ns ~ "",
        TRUE ~ Primary_significance
      )
    ) |>
    dplyr::group_by(Day) |>
    dplyr::mutate(
      Panel_max = max(Mean + SD, na.rm = TRUE),
      Label_y = Mean + SD + max(Panel_max * 0.045, 1e-9)
    ) |>
    dplyr::ungroup()
  combined_data <- combined_data |>
    dplyr::mutate(
      Day = factor(Day, levels = day_levels),
      Construct = factor(Construct, levels = construct_levels)
    )

  combined_plot <- ggplot2::ggplot(
    combined_results,
    ggplot2::aes(x = Construct, y = Mean, fill = Construct == control_name)
  ) +
    ggplot2::geom_col(width = 0.72, alpha = 0.92) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD),
      width = 0.20,
      linewidth = 0.35
    ) +
    ggplot2::geom_point(
      data = combined_data,
      ggplot2::aes(x = Construct, y = Value),
      inherit.aes = FALSE,
      position = ggplot2::position_jitter(width = 0.12, height = 0, seed = 260915),
      size = 0.95,
      color = "#202020",
      alpha = 0.72
    ) +
    ggplot2::geom_text(
      ggplot2::aes(y = Label_y, label = Label),
      size = 2.25,
      vjust = 0,
      color = "#202020"
    ) +
    ggplot2::facet_wrap(~Day, ncol = 2, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#7F7F7F", `FALSE` = "#4472C4"), guide = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, NA),
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = paste0(current$Metric, ": ", current$Selection),
      subtitle = "Mean ± SD with individual transfection wells; Welch t-tests with Holm-adjusted p-values vs Original",
      x = NULL,
      y = current$Metric
    ) +
    ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 9, color = "#555555"),
      strip.background = ggplot2::element_rect(fill = "#D9E2F3", color = NA),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(angle = 67, hjust = 1, vjust = 1, size = 6.2),
      panel.grid.major.y = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing = grid::unit(1.0, "lines")
    )

  ggplot2::ggsave(
    filename = file.path(
      figure_dir,
      paste0("COMBINED__", safe_filename(current$Metric), "__", safe_filename(current$Selection), ".png")
    ),
    plot = combined_plot,
    width = 24,
    height = max(7, 5.8 * ceiling(length(day_levels) / 2)),
    dpi = 300,
    bg = "white"
  )
}

analysis_notes <- c(
  "Primary inference: two-sided Welch t-tests versus Original with Holm family-wise error correction.",
  "The model is fitted separately within every Day × Selection × Metric condition.",
  "MFI-like metrics are log2-transformed for testing in the default auto mode; figures remain on the original scale.",
  "Bars show mean ± SD and points show independently transfected wells.",
  "FACS events are not replicates; each independently transfected well contributes one value.",
  "Dunnett-adjusted p-values are retained in the result table as a secondary comparison.",
  "Welch_Holm_results.csv is the primary result file; Dunnett_results.csv is retained as a compatibility copy.",
  "",
  "Warnings:",
  if (length(warning_messages) == 0) "None" else warning_messages
)
writeLines(analysis_notes, file.path(output_dir, "Analysis_notes.txt"), useBytes = TRUE)

message("Analysis complete: ", normalizePath(output_dir))
message("Conditions analyzed: ", nrow(diagnostics))
message("Variant comparisons: ", sum(results$Construct != control_name))
if (length(warning_messages) > 0) {
  message("Warnings: ", length(warning_messages), " (see Analysis_notes.txt)")
}
