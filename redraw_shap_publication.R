# Regenerate publication-ready SHAP figures without retraining models.
required_packages <- c(
  "ggplot2", "dplyr", "forcats", "shapviz", "ragg", "svglite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

output_dir <- Sys.getenv("MOLAN_OUTPUT_DIR", unset = "results")
input_root <- Sys.getenv(
  "MOLAN_SHAP_INPUT",
  unset = file.path(output_dir, "explainability", "SHAP")
)
output_root <- Sys.getenv(
  "MOLAN_SHAP_PUBLICATION_OUTPUT",
  unset = file.path(output_dir, "explainability", "SHAP_publication")
)
if (!dir.exists(input_root)) {
  stop("SHAP input directory not found: ", input_root)
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

models <- c("CoxPH", "GBM", "XGBoost")
model_colours <- c(CoxPH = "#3B6FB6", GBM = "#159D89", XGBoost = "#D67A2C")
value_colours <- c("#4B146E", "#C43C63", "#FDB515")

pretty_label <- function(x) {
  trimws(gsub("\\s+", " ", gsub("[._]+", " ", x)))
}

canonical_label <- function(x) {
  sub(" (Yes|No|Low|High)$", "", pretty_label(x))
}

theme_shap <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.key.height = grid::unit(14, "mm"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(colour = "#555555"),
      plot.margin = margin(7, 9, 7, 7)
    )
}

save_publication <- function(p, base, width, height) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(base, ".pdf"), p, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
  ggsave(paste0(base, ".svg"), p, width = width, height = height,
         units = "in", device = svglite::svglite, bg = "white")
  ragg::agg_tiff(paste0(base, ".tiff"), width = width, height = height,
                 units = "in", res = 600, compression = "lzw")
  print(p)
  dev.off()
  ragg::agg_png(paste0(base, ".png"), width = width, height = height,
                units = "in", res = 300)
  print(p)
  dev.off()
}

importance_list <- list()

for (model_name in models) {
  model_input <- file.path(input_root, model_name)
  model_output <- file.path(output_root, model_name)
  required_files <- file.path(
    model_input,
    c("shapviz_object.rds", "SHAP_mean_absolute_importance.csv")
  )
  if (any(!file.exists(required_files))) {
    stop(
      model_name, " SHAP files are missing: ",
      paste(required_files[!file.exists(required_files)], collapse = ", ")
    )
  }
  dir.create(model_output, recursive = TRUE, showWarnings = FALSE)

  sv <- readRDS(file.path(model_input, "shapviz_object.rds"))
  importance <- read.csv(
    file.path(model_input, "SHAP_mean_absolute_importance.csv"),
    check.names = FALSE
  ) %>% arrange(desc(mean_abs_shap))
  importance_list[[model_name]] <- importance
  order_top <- importance$variable
  order_bottom <- rev(order_top)

  bar_data <- importance %>%
    mutate(label = factor(pretty_label(variable), levels = pretty_label(order_bottom)))
  p_bar <- ggplot(bar_data, aes(mean_abs_shap, label)) +
    geom_col(width = 0.68, fill = model_colours[[model_name]]) +
    geom_text(aes(label = formatC(mean_abs_shap, digits = 3, format = "f")),
              hjust = 1.08, colour = "white", size = 2.8, family = "Arial") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Mean |SHAP value|", y = NULL,
         title = paste0(model_name, " — global feature importance"),
         subtitle = "External validation cohort") +
    theme_shap()
  save_publication(p_bar, file.path(model_output, "SHAP_importance"),
                   6.2, max(4, 0.27 * nrow(importance) + 1.8))

  set.seed(26)
  bee <- bind_rows(lapply(seq_along(order_bottom), function(i) {
    feature <- order_bottom[[i]]
    raw <- sv$X[[feature]]
    numeric_raw <- if (is.numeric(raw)) as.numeric(raw) else as.numeric(factor(raw))
    limits <- range(numeric_raw, na.rm = TRUE)
    scaled <- if (!all(is.finite(limits)) || diff(limits) == 0) {
      rep(0.5, length(numeric_raw))
    } else {
      (numeric_raw - limits[1]) / diff(limits)
    }
    data.frame(
      feature = feature,
      shap = as.numeric(sv$S[, feature]),
      value = scaled,
      y = i + runif(length(scaled), -0.22, 0.22)
    )
  }))
  p_bee <- ggplot(bee, aes(shap, y, colour = value)) +
    geom_vline(xintercept = 0, colour = "#A9A9A9", linewidth = 0.55) +
    geom_point(size = 1.55, alpha = 0.82, stroke = 0) +
    scale_y_continuous(breaks = seq_along(order_bottom),
                       labels = pretty_label(order_bottom),
                       expand = expansion(add = 0.6)) +
    scale_colour_gradientn(colours = value_colours, limits = c(0, 1),
                           breaks = c(0, 1), labels = c("Low", "High"),
                           name = "Feature value") +
    labs(x = "SHAP value (impact on predicted risk)", y = NULL,
         title = paste0(model_name, " — SHAP summary"),
         subtitle = "Positive values increase predicted risk; negative values decrease it") +
    theme_shap() + theme(legend.position = "right")
  save_publication(p_bee, file.path(model_output, "SHAP_summary"),
                   7.2, max(4.2, 0.28 * nrow(importance) + 2))

  dependence_dir <- file.path(model_output, "dependence")
  dir.create(dependence_dir, recursive = TRUE, showWarnings = FALSE)
  for (feature in colnames(sv$S)) {
    p_dep <- sv_dependence(sv, v = feature, color_var = NULL) +
      geom_hline(yintercept = 0, colour = "#A9A9A9", linewidth = 0.45) +
      labs(title = paste0(model_name, " — ", pretty_label(feature)),
           subtitle = "SHAP dependence in the external validation cohort",
           x = pretty_label(feature), y = "SHAP value") +
      theme_shap()
    safe_name <- substr(gsub("[^[:alnum:]_.-]+", "_", feature), 1, 120)
    save_publication(p_dep, file.path(dependence_dir, paste0(safe_name, "_dependence")),
                     5.4, 4.2)
  }

  local_index_file <- file.path(model_input, "local", "local_explanation_rows.csv")
  if (file.exists(local_index_file)) {
    local_index <- read.csv(local_index_file, check.names = FALSE)
    for (i in seq_len(nrow(local_index))) {
      label <- as.character(local_index$explanation[i])
      row_id <- as.integer(local_index$shap_row[i])
      p_local <- sv_waterfall(sv, row_id = row_id, max_display = ncol(sv$S)) +
        labs(title = paste0(model_name, " — ", gsub("_", " ", label), " patient"),
             subtitle = "Local SHAP explanation of predicted risk") +
        theme_shap()
      save_publication(p_local,
                       file.path(model_output, "local", paste0(label, "_waterfall")),
                       7.2, max(4.2, 0.28 * ncol(sv$S) + 2))
    }
  }
}

comparison <- bind_rows(importance_list) %>%
  mutate(
    model = factor(model, levels = models),
    variable_label = canonical_label(variable),
    variable_label = forcats::fct_reorder(variable_label, importance_share)
  )
p_comparison <- ggplot(comparison, aes(importance_share, variable_label, colour = model)) +
  geom_point(position = position_dodge(width = 0.55), size = 2.3) +
  scale_colour_manual(values = model_colours) +
  labs(x = "Within-model share of mean absolute SHAP", y = NULL, colour = "Model",
       title = "SHAP importance comparison",
       subtitle = "Importance is normalized within each model") +
  theme_shap() + theme(legend.position = "top")
save_publication(p_comparison, file.path(output_root, "SHAP_model_comparison"),
                 7.2, max(4.2, 0.32 * length(unique(comparison$variable_label)) + 2))

write.csv(bind_rows(importance_list),
          file.path(output_root, "SHAP_importance_source_data.csv"), row.names = FALSE)
cat("Publication SHAP figures written to: ", normalizePath(output_root), "\n")

