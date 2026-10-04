# ============================================================
# Integrated CNPS wildfire analysis and figure pipeline
# ============================================================

options(stringsAsFactors = FALSE)

# ---- User settings ----
# For test data, set input_dir to: "C:/Users/79840/Desktop/Data2"
input_dir <- "C:/Users/79840/Desktop/Data2"
results_dir <- file.path(input_dir, "Integrated_Results")
figures_dir <- file.path(results_dir, "Figures")
tables_dir <- file.path(results_dir, "Tables")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Set TRUE to rebuild process-level summary tables from raw TPM matrices.
RUN_PREPROCESSING <- TRUE

# Network analysis settings are kept at full inference by default.
QUICK_NETWORK_TEST <- FALSE

# ---- Packages ----
required_packages <- c(
  "tidyverse", "patchwork", "randomForest", "nlme", "lme4", "lmerTest",
  "vegan", "permute", "sf", "rnaturalearth", "ggrepel", "cowplot", "igraph",
  "scales", "RColorBrewer", "gridExtra", "ggalluvial", "ape"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(randomForest)
  library(nlme)
  library(lme4)
  library(lmerTest)
  library(vegan)
  library(sf)
  library(rnaturalearth)
  library(ggrepel)
  library(cowplot)
  library(igraph)
  library(scales)
  library(RColorBrewer)
  library(gridExtra)
  library(ggalluvial)
  library(ape)
  library(grid)
})

if (!dir.exists(input_dir)) {
  stop("input_dir does not exist: ", input_dir)
}
setwd(input_dir)

required_input_files <- c(
  "Basic_ALL.txt",
  "group_ALL.txt", "group_SX.txt",
  "C-classification.txt", "N-classification.txt",
  "P-classification.txt", "S-classification.txt",
  "TPM_C_ALL.tsv", "TPM_N_ALL.tsv",
  "TPM_P_ALL.tsv", "TPM_S_ALL.tsv",
  "TPM_C_SX.tsv", "TPM_N_SX.tsv",
  "TPM_P_SX.tsv", "TPM_S_SX.tsv"
)

missing_input_files <- required_input_files[
  !file.exists(file.path(input_dir, required_input_files))
]

if (length(missing_input_files)) {
  stop(
    "Missing required input files:\n",
    paste0("  - ", missing_input_files, collapse = "\n"),
    "\n\nCopy these files into input_dir before running the full Figure 1-24 pipeline."
  )
}


figure_file <- function(number, label) {
  file.path(
    figures_dir,
    sprintf("Figure_%02d_%s.pdf", as.integer(number), label)
  )
}

subfigure_file <- function(tag, label) {
  file.path(
    figures_dir,
    sprintf("Figure_%s_%s.pdf", tag, label)
  )
}

save_figure <- function(number, label, plot, width, height, family = NULL) {
  args <- list(
    filename = figure_file(number, label),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
  if (!is.null(family)) args$family <- family
  do.call(ggsave, args)
}

format_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  out <- rep("p = NA", length(p))
  out[is.finite(p) & p < 0.001] <- "p < 0.001"
  idx <- is.finite(p) & p >= 0.001
  out[idx] <- paste0("p = ", formatC(p[idx], format = "f", digits = 3))
  out
}

sig_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

clean_text <- function(data) {
  data %>%
    mutate(across(where(is.character), ~ trimws(gsub("\r", "", .x))))
}

# Robust nitrogen-process recoding.
recode_n_process <- function(x) {
  dplyr::case_when(
    x == "ANR" ~ "Assimilatory nitrate reduction",
    x == "DNF" ~ "Denitrification",
    x == "DNRA" ~ "Dissimilatory nitrate reduction",
    x == "N-fix" ~ "N-fixation",
    x == "NIT" ~ "Nitrification",
    x == "Org" ~ "Organic synthesis & degradation",
    TRUE ~ x
  )
}

# ---- Shared process-level TPM aggregation ----
first_existing <- function(candidates) {
  hits <- candidates[file.exists(file.path(input_dir, candidates))]
  if (!length(hits)) return(NA_character_)
  hits[[1]]
}

# Shared process-level TPM aggregation.
# Includes deduplication, metadata matching and sample-set checks.
aggregate_process_tpm <- function(
    tpm_file, classification_file, group_file, output_file,
    place_filter = NULL
) {
  if (any(is.na(c(tpm_file, classification_file, group_file)))) return(invisible(NULL))
  
  tpm <- read.table(
    file.path(input_dir, tpm_file),
    header = TRUE, row.names = 1, sep = "\t",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  classification <- read.table(
    file.path(input_dir, classification_file),
    header = FALSE, sep = "\t", fill = TRUE, strip.white = TRUE,
    col.names = c("gene", "process"), stringsAsFactors = FALSE,
    quote = "", comment.char = ""
  ) %>%
    mutate(
      gene = trimws(gene),
      process = trimws(gsub('"', "", process))
    ) %>%
    filter(
      !is.na(gene), gene != "",
      !is.na(process), process != "",
      !tolower(process) %in% c("deleted", "others", "other", "unclassified")
    ) %>%
    distinct(gene, process)

  if (grepl("^N(-|cycle_)", basename(classification_file), ignore.case = TRUE)) {
    classification <- classification %>%
      mutate(process = recode_n_process(process))
  }

  group_info <- read.table(
    file.path(input_dir, group_file),
    header = TRUE, sep = "\t",
    check.names = FALSE, stringsAsFactors = FALSE
  ) %>% 
    clean_text() %>%
    distinct(sample, .keep_all = TRUE)
  
  if (!setequal(colnames(tpm), group_info$sample)) {
    warning(sprintf("TPM samples and metadata samples are not identical for %s.", tpm_file))
  }
  
  summary_df <- tpm %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "TPM") %>%
    mutate(TPM = suppressWarnings(as.numeric(TPM))) %>%
    filter(sample %in% group_info$sample) %>%
    inner_join(classification, by = "gene") %>%
    inner_join(group_info, by = "sample")
  
  if (!is.null(place_filter)) {
    summary_df <- summary_df %>% filter(place %in% place_filter)
  }
  
  grouping_vars <- intersect(
    c("sample", "place", "time", "process", "treatment"),
    names(summary_df)
  )

  summary_df <- summary_df %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    select(any_of(c("sample", "place", "time", "process", "treatment", "TPM")))
  
  write.table(
    summary_df,
    file.path(input_dir, output_file),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  invisible(summary_df)
}

if (RUN_PREPROCESSING) {
  elements <- c("C", "N", "P", "S")
  
  for (element in elements) {
    regional_tpm <- first_existing(
      if (element == "N") c("TPM_N_ALL.tsv", "TPM_by_category.tsv")
      else paste0("TPM_", element, "_ALL.tsv")
    )
    regional_class <- first_existing(
      if (element == "N") c("N-classification.txt", "Ncycle_classification.txt")
      else paste0(element, "-classification.txt")
    )
    
    if (!is.na(regional_tpm) && !is.na(regional_class) &&
        file.exists(file.path(input_dir, "group_ALL.txt"))) {
      aggregate_process_tpm(
        regional_tpm, regional_class, "group_ALL.txt",
        paste0("TPM_summary_", element, "_ALL.txt")
      )
    }
    
    sx_tpm <- first_existing(paste0("TPM_", element, "_SX.tsv"))
    sx_class <- first_existing(paste0(element, "-classification.txt"))
    
    if (!is.na(sx_tpm) && !is.na(sx_class) &&
        file.exists(file.path(input_dir, "group_SX.txt"))) {
      aggregate_process_tpm(
        sx_tpm, sx_class, "group_SX.txt",
        paste0("TPM_summary_", element, "_SX.txt"),
        place_filter = "XC"
      )
    }
  }
}


# ============================================================
# Figure 1. National sampling sites
# ============================================================

line_w <- 0.72 
tick_len <- unit(0.25, "cm")

base_map_theme <- theme_bw(base_size = 18) +
  theme(
    text = element_text(size = 18),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", face = "bold", size = 18),
    axis.ticks = element_line(color = "black", linewidth = line_w),
    axis.ticks.length = tick_len,
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18),
    panel.background = element_rect(fill = "#f0f7fa"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    legend.background = element_rect(fill = alpha("white", 0.7), color = "grey80", linewidth = 0.3)
  )

site_data <- data.frame(
  Abbr = c("GH", "HH", "OQ", "BD", "JZ", "HZ", "YT", "SM", "DL", "CZ", "ZH", "SY", "XC"),
  Lon = c(122.02, 127.29, 123.20, 114.62, 113.00, 119.66, 117.07, 118.76, 100.34, 113.17, 113.54, 123.78, 102.19),
  Lat = c(51.72, 50.22, 49.54, 39.25, 37.39, 30.04, 28.24, 25.93, 25.63, 25.52, 22.33, 42.00, 28.24),
  MAT = c(-5.3, -3.0, -0.5, 13.4, 9.5, 16.7, 18.0, 17.4, 14.9, 18.0, 22.5, 9.7, 18.3),
  MAP = c(437, 425, 414, 499, 458, 1477, 1750, 1700, 1051, 1493, 2062, 896, 965)
)

site_data <- site_data %>%
  mutate(Group = ifelse(Abbr %in% c("SY", "XC"), "Special", "General"))

sites_sf <- st_as_sf(site_data, coords = c("Lon", "Lat"), crs = 4326)

world_map <- ne_countries(scale = "medium", returnclass = "sf")
china_full <- world_map %>% filter(admin %in% c("China", "Taiwan"))

main_plot <- ggplot() +
  geom_sf(data = world_map, fill = "white", color = "grey85", linewidth = 0.2) +
  geom_sf(data = china_full, fill = "white", color = "black", linewidth = line_w) +
  geom_sf(data = sites_sf, aes(size = MAP, fill = MAT, shape = Group), 
          color = "black", stroke = 1.0) +
  geom_text_repel(data = site_data, aes(x = Lon, y = Lat, label = Abbr),
                  size = 5.5, fontface = "bold", box.padding = 0.5, 
                  point.padding = 0.4, segment.color = "black", segment.linewidth = line_w) +
  scale_shape_manual(values = c("General" = 21, "Special" = 24), guide = "none") +
  scale_fill_gradientn(colors = c("#2166AC", "#D1E5F0", "#FFFFBF", "#F4A582", "#B2182B"),
                       name = "MAT (°C)") +
  scale_size_continuous(range = c(3, 8), name = "MAP (mm)") +
  coord_sf(xlim = c(73, 135), ylim = c(18, 54), expand = FALSE) +
  base_map_theme +
  labs(x = "Longitude (°E)", y = "Latitude (°N)")

inset_plot <- ggplot() +
  geom_sf(data = world_map, fill = "white", color = "grey85", linewidth = 0.1) +
  geom_sf(data = china_full, fill = "white", color = "black", linewidth = line_w) +
  coord_sf(xlim = c(106, 122), ylim = c(3, 25)) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    panel.background = element_rect(fill = "white", color = NA)
  )

final_output <- ggdraw(main_plot) +
  draw_plot(inset_plot, x = 0.16, y = 0.12, width = 0.14, height = 0.22)

ggsave(subfigure_file("1a", "Sampling_Sites_Map"), final_output, width = 10.5, height = 11, device = cairo_pdf, bg = "white")


# ============================================================
# Figure 2. Environmental responses to wildfire
# ============================================================
out_dir <- file.path(tables_dir, "Figure_02_Environmental_Responses")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

line_w <- 0.8

content_vars <- c("SOC", "TN", "TP", "TS", "DOC", "DON", "AN", "AP", "AS")
ratio_vars <- c(
  "SOCTN", "SOCTP", "TNTP", "SOCTS", "TNTS", "TPTS",
  "DOCAN", "DOCAP", "ANAP", "DOCAS", "ANAS", "APAS"
)
microbial_vars <- c("CO2", "MBC", "qCO2")
target_vars <- c(content_vars, ratio_vars, microbial_vars)

label_map <- c(
  SOC = "SOC", TN = "TN", TP = "TP", TS = "TS",
  DOC = "DOC", DON = "DON", AN = "Available N", AP = "Available P", AS = "Available S",
  SOCTN = "SOC:TN", SOCTP = "SOC:TP", TNTP = "TN:TP",
  SOCTS = "SOC:TS", TNTS = "TN:TS", TPTS = "TP:TS",
  DOCAN = "DOC:AN", DOCAP = "DOC:AP", ANAP = "AN:AP",
  DOCAS = "DOC:AS", ANAS = "AN:AS", APAS = "AP:AS",
  CO2 = "CO2", MBC = "MBC", qCO2 = "qCO2", pH = "pH"
)

group_map <- c(
  setNames(rep("Nutrient contents", length(content_vars)), content_vars),
  setNames(rep("Nutrient ratios", length(ratio_vars)), ratio_vars),
  setNames(rep("Microbial properties", length(microbial_vars)), microbial_vars),
  pH = "Soil pH"
)

basic <- read_tsv(
  file.path(input_dir, "Basic_ALL.txt"),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
) %>%
  mutate(sample = as.character(sample))

group <- read_tsv(
  file.path(input_dir, "group_ALL.txt"),
  show_col_types = FALSE,
  na = c("", "NA", "NaN")
) %>%
  transmute(
    sample = as.character(sample),
    place = factor(place),
    Fire = case_when(
      tolower(trimws(as.character(treatment))) %in%
        c("1ub", "ub", "unburned", "unburnt", "1_unburntsoil") ~ "Unburned",
      tolower(trimws(as.character(treatment))) %in%
        c("b", "burned", "burnt", "2_highburntsoil") ~ "Burned",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Fire)) %>%
  mutate(Fire = factor(Fire, levels = c("Unburned", "Burned")))

dat <- group %>%
  inner_join(
    basic %>% select(-any_of(c("place", "time", "treatment", "Fire"))),
    by = "sample"
  )

fit_effect <- function(data, variable, is_log1p = TRUE) {
  d <- data %>%
    transmute(
      place,
      Fire,
      response = suppressWarnings(as.numeric(.data[[variable]]))
    ) %>%
    drop_na()
  
  if (is_log1p) {
    if (any(d$response < 0)) return(NULL)
    d$response <- log1p(d$response)
  }
  
  fit <- tryCatch(
    lmer(response ~ Fire + (1 | place), data = d, REML = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(NULL)
  
  tab <- coef(summary(fit))
  if (!"FireBurned" %in% rownames(tab)) return(NULL)
  
  estimate <- tab["FireBurned", "Estimate"]
  se <- tab["FireBurned", "Std. Error"]
  
  # [Bug Fix: lmer subscript] Prevent crash if lmerTest fails Satterthwaite
  p_val <- if ("Pr(>|t|)" %in% colnames(tab)) tab["FireBurned", "Pr(>|t|)"] else NA_real_
  t_val <- if ("t value" %in% colnames(tab)) tab["FireBurned", "t value"] else NA_real_
  
  tibble(
    Variable = variable,
    Transform = if_else(is_log1p, "log1p", "raw"),
    Estimate = estimate,
    SE = se,
    t_value = t_val,
    p_value = p_val,
    CI_low = estimate - 1.96 * se,
    CI_high = estimate + 1.96 * se,
    n = nrow(d)
  )
}

main_results <- map_dfr(
  intersect(target_vars, names(dat)),
  ~ fit_effect(dat, .x, is_log1p = TRUE)
)

ph_result <- if ("pH" %in% names(dat)) {
  fit_effect(dat, "pH", is_log1p = FALSE)
} else {
  NULL
}

results <- bind_rows(main_results, ph_result) %>%
  mutate(
    Star = case_when(
      is.na(p_value) ~ "ns",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    Impact = case_when(
      !is.na(p_value) & p_value < 0.05 & Estimate > 0 ~ "Significant increase",
      !is.na(p_value) & p_value < 0.05 & Estimate < 0 ~ "Significant decrease",
      TRUE ~ "Non-significant"
    ),
    Label = unname(label_map[Variable]),
    Group = unname(group_map[Variable])
  )

write.csv(
  results,
  file.path(out_dir, "Environmental_LMM_Forest_results_rawP.csv"),
  row.names = FALSE
)

impact_colors <- c(
  "Significant increase" = "#D53E4F",
  "Significant decrease" = "#3288BD",
  "Non-significant" = "grey65"
)

base_theme_fig2 <- theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", face = "bold", size = 18),
    axis.ticks = element_line(color = "black", linewidth = line_w),
    axis.ticks.length = grid::unit(0.25, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dotted"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = line_w),
    strip.text = element_text(face = "bold", size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 16)
  )

plot_forest <- function(data, x_title) {
  data <- data %>%
    mutate(Label = forcats::fct_reorder(Label, Estimate))
  
  effect_range <- range(c(data$CI_low, data$CI_high), na.rm = TRUE)
  star_offset <- max(diff(effect_range) * 0.04, 0.02)
  
  ggplot(data, aes(x = Estimate, y = Label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 1) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high, color = Impact), height = 0.20, linewidth = 1) +
    geom_point(aes(fill = Impact), shape = 21, size = 4.5, color = "black", stroke = line_w) +
    geom_text(aes(x = CI_high + star_offset, label = Star), hjust = 0, size = 5.5, fontface = "bold") +
    scale_fill_manual(values = impact_colors) +
    scale_color_manual(values = impact_colors) +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.18))) +
    labs(x = x_title, y = NULL) +
    base_theme_fig2
}

p_main <- results %>%
  filter(Transform == "log1p") %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Nutrient contents", "Nutrient ratios", "Microbial properties")
    )
  ) %>%
  plot_forest("LMM Estimate (Burned vs Unburned)") +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Response of soil environmental factors to wildfire") +
  theme(legend.position = "bottom")

if (any(results$Transform == "raw")) {
  p_ph <- results %>%
    filter(Transform == "raw") %>%
    plot_forest("LMM Estimate for pH") +
    theme(legend.position = "none")
  
  final_plot <- p_main / p_ph +
    plot_layout(
      heights = c(max(nrow(main_results), 10), 2.5),
      guides = "collect"
    ) &
    theme(legend.position = "bottom")
} else {
  final_plot <- p_main
}

ggsave(subfigure_file("1b", "Environmental_Responses"), final_plot, width = 9, height = max(11, 0.40 * nrow(main_results) + 4), device = cairo_pdf, bg = "white")


# ============================================================
# Figure 2 and Supplementary Figures S1-S4, S12-S13. Multifunctionality and environmental controls
# ============================================================
out_dir <- file.path(tables_dir, "Figure_03_to_07_Multifunctionality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_n_samples <- 108L
expected_n_sites <- 11L
rf_ntree <- 1000L
rf_perm_n <- 199L
random_seed <- 123L
line_w <- 0.8

treat_cols <- c(Unburned = "#3288BD", Burned = "#D95A68")

files <- list(
  C = file.path(input_dir, "TPM_summary_C_ALL.txt"),
  N = file.path(input_dir, "TPM_summary_N_ALL.txt"),
  P = file.path(input_dir, "TPM_summary_P_ALL.txt"),
  S = file.path(input_dir, "TPM_summary_S_ALL.txt"),
  soil = file.path(input_dir, "Basic_ALL.txt"),
  metadata = file.path(input_dir, "group_ALL.txt")
)

missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files)) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))

soil_pc1_vars <- c("SOC", "TN", "TP", "TS", "DOC", "DON", "AN", "AP", "AS", "NH4", "NO3")
soil_predictor_vars <- c(
  "pH", "SOC", "TN", "TP", "TS", "DOC", "DON", "AN", "AP", "AS",
  "SOCTN", "SOCTP", "TNTP", "SOCTS", "TNTS", "TPTS",
  "DOCAN", "DOCAP", "ANAP", "DOCAS", "ANAS", "APAS",
  "CO2", "MBC", "qCO2"
)

soil_labels <- c(
  pH = "pH", SOC = "SOC", TN = "TN", TP = "TP", TS = "TS",
  DOC = "DOC", DON = "DON", AN = "AN", AP = "AP", AS = "AS",
  SOCTN = "SOC:TN", SOCTP = "SOC:TP", TNTP = "TN:TP",
  SOCTS = "SOC:TS", TNTS = "TN:TS", TPTS = "TP:TS",
  DOCAN = "DOC:AN", DOCAP = "DOC:AP", ANAP = "AN:AP",
  DOCAS = "DOC:AS", ANAS = "AN:AS", APAS = "AP:AS",
  CO2 = "CO2", MBC = "MBC", qCO2 = "qCO2",
  NH4 = "NH4", NO3 = "NO3"
)

base_theme_fig3 <- theme_classic(base_size = 18) +
  theme(
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", face = "bold", size = 18),
    axis.ticks = element_line(color = "black", linewidth = line_w),
    axis.ticks.length = grid::unit(0.25, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 18)
  )

clean_names <- function(x) {
  nm <- trimws(names(x))
  nm[is.na(nm) | nm == ""] <- paste0("empty_", seq_len(sum(is.na(nm) | nm == "")))
  names(x) <- make.unique(nm, sep = "_")
  x
}

read_cnps <- function(path, element) {
  readr::read_tsv(path, show_col_types = FALSE, na = c("", "NA", "NaN", "nan")) %>%
    clean_names() %>%
    transmute(
      sample = as.character(sample),
      Element = element,
      process = trimws(gsub('"', "", process)),
      TPM = suppressWarnings(as.numeric(TPM))
    ) %>%
    filter(
      !is.na(process), process != "", is.finite(TPM),
      !tolower(process) %in% c("na", "n/a", "nan", "null", "none", "others", "other", "deleted", "unclassified")
    ) %>%
    mutate(
      process = if_else(Element == "N", recode_n_process(process), process),
      Module = paste(Element, process, sep = ": ")
    ) %>%
    group_by(sample, Element, process, Module) %>%
    summarise(TPM = sum(TPM), .groups = "drop")
}

read_metadata <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, na = c("", "NA", "NaN", "nan")) %>%
    clean_names() %>%
    transmute(
      sample = as.character(sample),
      place = trimws(as.character(place)),
      treatment = trimws(as.character(treatment))
    ) %>%
    distinct(sample, .keep_all = TRUE)
}

read_soil <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, na = c("", "NA", "NaN", "nan")) %>%
    clean_names() %>%
    mutate(sample = as.character(sample)) %>%
    distinct(sample, .keep_all = TRUE) %>%
    select(-any_of(c("place", "time", "treatment"))) %>%
    mutate(across(-sample, ~ suppressWarnings(as.numeric(.x))))
}

calc_emf <- function(df, vars) {
  x <- df %>% select(all_of(vars)) %>% as.data.frame()
  x[] <- lapply(x, log1p)
  x <- x[, vapply(x, sd, numeric(1), na.rm = TRUE) > 0, drop = FALSE]
  rowMeans(scale(x), na.rm = TRUE)
}

calc_pc1 <- function(df, vars) {
  x <- df %>%
    select(all_of(vars)) %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(across(everything(), ~ replace_na(.x, median(.x, na.rm = TRUE)))) %>%
    as.data.frame()
  
  keep <- vapply(x, sd, numeric(1), na.rm = TRUE) > 0
  x <- x[, keep, drop = FALSE]
  z <- scale(x)
  pca <- prcomp(z, center = FALSE, scale. = FALSE)
  pc1 <- pca$x[, 1]
  
  anchor <- rowMeans(z)
  if (cor(pc1, anchor, use = "complete.obs") < 0) {
    pc1 <- -pc1
    pca$x[, 1] <- -pca$x[, 1]
    pca$rotation[, 1] <- -pca$rotation[, 1]
  }
  
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  write.csv(
    data.frame(Variable = rownames(pca$rotation), PC1_loading = pca$rotation[, 1]),
    file.path(out_dir, "Soil_nutrient_11var_PCA_loadings.csv"), row.names = FALSE
  )
  write.csv(
    data.frame(Axis = paste0("PC", seq_along(variance)), Explained_variance_percent = variance),
    file.path(out_dir, "Soil_nutrient_11var_PCA_variance.csv"), row.names = FALSE
  )
  
  list(score = pc1, variance = variance)
}

cnps_long <- bind_rows(
  read_cnps(files$C, "C"), read_cnps(files$N, "N"),
  read_cnps(files$P, "P"), read_cnps(files$S, "S")
)

process_check <- cnps_long %>% distinct(Element, process, Module)
if (nrow(process_check) != 33L) stop("Expected 33 CNPS processes, found ", nrow(process_check), ".")

cnps_wide <- cnps_long %>%
  select(sample, Module, TPM) %>%
  pivot_wider(names_from = Module, values_from = TPM, values_fill = 0, values_fn = sum)

metadata <- read_metadata(files$metadata)
soil <- read_soil(files$soil)

dat <- metadata %>%
  inner_join(cnps_wide, by = "sample") %>%
  inner_join(soil, by = "sample") %>%
  mutate(
    place = factor(place),
    Fire_factor = factor(treatment, levels = c("1UB", "B"), labels = c("Unburned", "Burned"))
  ) %>%
  filter(!is.na(Fire_factor))

if (nrow(dat) != expected_n_samples || nlevels(dat$place) != expected_n_sites) {
  stop(
    "Unexpected regional design: found ", nrow(dat),
    " samples and ", nlevels(dat$place), " sites; expected ",
    expected_n_samples, " samples and ", expected_n_sites, " sites."
  )
}

required_soil <- unique(c(soil_pc1_vars, soil_predictor_vars))
missing_soil <- setdiff(required_soil, names(dat))
if (length(missing_soil)) stop("Missing soil variables: ", paste(missing_soil, collapse = ", "))

cnps_process_vars <- process_check$Module
dat$CNPS_EMF <- calc_emf(dat, cnps_process_vars)
soil_pc <- calc_pc1(dat, soil_pc1_vars)
dat$Soil_nutrient_PC1 <- soil_pc$score

write.csv(
  data.frame(PC1 = "Soil nutrient PC1 (11 variables)", Explained_variance_percent = soil_pc$variance[1]),
  file.path(out_dir, "Fig2_PC1_explained_variance_summary.csv"), row.names = FALSE
)

# Fig. 2a
fit_a <- nlme::lme(CNPS_EMF ~ Fire_factor, random = ~ 1 | place, data = dat, method = "REML")
p_a_value <- summary(fit_a)$tTable["Fire_factorBurned", "p-value"]
write.csv(as.data.frame(summary(fit_a)$tTable), file.path(out_dir, "Fig2a_CNPS_EMF_LMM_results.csv"))

dat$Panel_Title <- "CNPS Overall"
y_range <- range(dat$CNPS_EMF, na.rm = TRUE)
y_pos <- y_range[2] + diff(y_range) * 0.12

p_a <- ggplot(dat, aes(Fire_factor, CNPS_EMF, fill = Fire_factor)) +
  geom_violin(color = "black", linewidth = line_w, trim = FALSE) +
  geom_boxplot(width = 0.18, fill = "white", outlier.shape = NA, linewidth = line_w) +
  geom_jitter(width = 0.12, size = 2.2, shape = 21, fill = "white", color = "black", alpha = 0.75) +
  annotate("text", x = 1.5, y = y_pos, label = format_p(p_a_value), fontface = "italic", size = 7) +
  facet_grid(. ~ Panel_Title) +
  scale_fill_manual(values = treat_cols) +
  labs(x = NULL, y = "Multifunctionality Index (Z-score EMF)") +
  coord_cartesian(clip = "off") +
  base_theme_fig3 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 22),
    strip.background = element_rect(fill = "#F5F5F5", color = NA),
    strip.text = element_text(face = "bold", size = 24)
  )

ggsave(subfigure_file("2a", "CNPS_Multifunctionality"), p_a, width = 5.5, height = 5.5, device = cairo_pdf, bg = "white")

# Fig. 2b
rf_data <- dat %>% select(CNPS_EMF, Fire_factor, all_of(soil_predictor_vars)) %>% drop_na()

run_rf <- function(data, fire_level) {
  d <- data %>% filter(Fire_factor == fire_level) %>% select(CNPS_EMF, all_of(soil_predictor_vars))
  
  constant <- names(d)[-1][vapply(d[-1], sd, numeric(1), na.rm = TRUE) == 0]
  # [Bug Fix: RF zero-variance] Remove constant predictors to prevent RF crash
  if (length(constant)) {
    warning(fire_level, ": zero-variance predictors removed: ", paste(constant, collapse = ", "))
    d <- d %>% select(-all_of(constant))
  }
  
  set.seed(random_seed)
  model <- randomForest(CNPS_EMF ~ ., data = d, importance = TRUE, ntree = rf_ntree)
  observed <- importance(model, type = 1, scale = FALSE) %>%
    as.data.frame() %>% rownames_to_column("Variable")
  names(observed)[names(observed) == "%IncMSE"] <- "IncMSE"
  
  set.seed(random_seed)
  null <- map_dfr(seq_len(rf_perm_n), function(i) {
    perm <- d
    perm$CNPS_EMF <- sample(perm$CNPS_EMF)
    fit <- randomForest(CNPS_EMF ~ ., data = perm, importance = TRUE, ntree = rf_ntree)
    x <- importance(fit, type = 1, scale = FALSE) %>%
      as.data.frame() %>% rownames_to_column("Variable")
    names(x)[names(x) == "%IncMSE"] <- "IncMSE_null"
    x
  })
  
  importance_table <- observed %>%
    mutate(
      Fire_factor = fire_level,
      p_perm = map2_dbl(Variable, IncMSE, ~ {
        z <- null %>% filter(Variable == .x) %>% pull(IncMSE_null)
        (sum(z >= .y, na.rm = TRUE) + 1) / (sum(is.finite(z)) + 1)
      }),
      star = sig_star(p_perm),
      Significant = p_perm < 0.05,
      Plot_label = unname(soil_labels[Variable])
    )
  
  summary_table <- tibble(
    Fire_factor = fire_level,
    n_complete = nrow(d),
    n_predictors = ncol(d) - 1,
    ntree = rf_ntree,
    response_permutations = rf_perm_n,
    OOB_explained_percent = tail(model$rsq, 1) * 100
  )
  
  list(importance = importance_table, summary = summary_table)
}

rf_unburned <- run_rf(rf_data, "Unburned")
rf_burned <- run_rf(rf_data, "Burned")
rf_importance <- bind_rows(rf_unburned$importance, rf_burned$importance)
rf_summary <- bind_rows(rf_unburned$summary, rf_burned$summary)

write.csv(rf_importance, file.path(out_dir, "Fig2b_random_forest_importance_by_fire.csv"), row.names = FALSE)
write.csv(rf_summary, file.path(out_dir, "Fig2b_random_forest_summary_by_fire.csv"), row.names = FALSE)

rf_order <- rf_importance %>%
  group_by(Variable) %>%
  summarise(mean_importance = mean(IncMSE, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_importance)) %>%
  pull(Variable)

rf_plot_data <- rf_importance %>%
  mutate(
    Variable = factor(Variable, levels = rev(rf_order)),
    Panel = paste0(Fire_factor, "\nExplained = ", round(rf_summary$OOB_explained_percent[match(Fire_factor, rf_summary$Fire_factor)]), "%")
  )

p_b <- ggplot(rf_plot_data, aes(Variable, IncMSE, fill = Significant)) +
  geom_col(width = 0.72, color = "black", linewidth = line_w) +
  geom_text(aes(label = star), hjust = -0.15, size = 5.2, fontface = "bold") +
  facet_wrap(~ Panel, nrow = 1) +
  coord_flip(clip = "off") +
  scale_x_discrete(labels = soil_labels) +
  scale_fill_manual(values = c(`TRUE` = "#FD8D3C", `FALSE` = "#FFF2DD")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = NULL, y = "Increase in MSE%") +
  base_theme_fig3 +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 15),
    strip.background = element_rect(fill = "#F5F5F5", color = "black", linewidth = line_w),
    strip.text = element_text(face = "bold", size = 19)
  )

ggsave(subfigure_file("2b", "Random_Forest_Importance"), p_b, width = 13.5, height = 10.5, device = cairo_pdf, bg = "white")

# Fig. 2c
regression_stats <- list()
interaction_stats <- list()

make_regression_plot <- function(data, variable, show_legend = FALSE) {
  d <- data %>%
    select(CNPS_EMF, Fire_factor, place, all_of(variable)) %>%
    drop_na() %>%
    mutate(Fire_factor = factor(Fire_factor, levels = c("Unburned", "Burned")), place = factor(place))
  
  group_stats <- map_dfr(levels(d$Fire_factor), function(g) {
    dg <- filter(d, Fire_factor == g)
    fit <- lm(reformulate(variable, response = "CNPS_EMF"), data = dg)
    s <- summary(fit)
    tibble(
      Predictor = variable,
      Fire_factor = g,
      n = nrow(dg),
      estimate = coef(fit)[2],
      SE = coef(s)[2, "Std. Error"],
      R2 = s$r.squared,
      adjusted_R2 = s$adj.r.squared,
      p_value = coef(s)[2, "Pr(>|t|)"]
    )
  })
  regression_stats[[variable]] <<- group_stats
  
  fit_int <- nlme::lme(
    fixed = as.formula(paste0("CNPS_EMF ~ `", variable, "` * Fire_factor")),
    random = ~ 1 | place,
    data = d,
    method = "REML",
    control = nlme::lmeControl(opt = "optim")
  )
  
  tab <- as.data.frame(summary(fit_int)$tTable) %>% rownames_to_column("Term")
  int_term <- grep(":Fire_factorBurned$|Fire_factorBurned:", tab$Term, value = TRUE)
  if (length(int_term) != 1) stop("Interaction term not found for ", variable, ".")
  
  int_row <- tab %>%
    filter(Term == int_term) %>%
    transmute(
      Predictor = variable,
      estimate = Value,
      SE = `Std.Error`,
      DF = DF,
      t_value = `t-value`,
      p_interaction = `p-value`
    )
  interaction_stats[[variable]] <<- int_row
  
  ann <- group_stats %>%
    mutate(
      Label = paste0(Fire_factor, ": R2 = ", formatC(R2, format = "f", digits = 3), ", ", format_p(p_value)),
      vjust_value = if_else(Fire_factor == "Unburned", 1.2, 2.7)
    )
  
  ggplot(d, aes(.data[[variable]], CNPS_EMF)) +
    geom_point(aes(fill = Fire_factor), shape = 21, size = 3, color = "black", alpha = 0.82) +
    geom_smooth(aes(color = Fire_factor, fill = Fire_factor), method = "lm", formula = y ~ x,
                se = TRUE, linewidth = 1, alpha = 0.16) +
    geom_text(
      data = ann,
      aes(x = -Inf, y = Inf, label = Label, color = Fire_factor, vjust = vjust_value),
      inherit.aes = FALSE, hjust = -0.05, size = 4.5, fontface = "bold"
    ) +
    annotate(
      "text", x = -Inf, y = -Inf,
      label = paste0("Interaction: ", format_p(int_row$p_interaction)),
      hjust = -0.05, vjust = -0.7, size = 4.5, fontface = "bold.italic"
    ) +
    scale_fill_manual(values = treat_cols, name = NULL) +
    scale_color_manual(values = treat_cols, name = NULL) +
    labs(x = unname(soil_labels[variable]), y = "Multifunctionality Index (Z-score EMF)") +
    coord_cartesian(clip = "off") +
    base_theme_fig3 +
    theme(
      legend.position = if (show_legend) "top" else "none",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold")
    )
}

p_c_list <- map(
  soil_predictor_vars,
  ~ make_regression_plot(dat, .x, TRUE)
)
names(p_c_list) <- soil_predictor_vars

regression_results <- bind_rows(regression_stats)
interaction_results <- bind_rows(interaction_stats)
write.csv(
  regression_results,
  file.path(out_dir, "Fig2c_all25_group_specific_regression_results.csv"),
  row.names = FALSE
)
write.csv(
  interaction_results,
  file.path(out_dir, "Fig2c_all25_Fire_x_environment_LMM_interaction_results.csv"),
  row.names = FALSE
)

# The 25 predictor panels are retained in their original order and split
# into five five-panel figures.
regression_panel_sets <- list(
  "2c" = c("pH", "SOC", "TN", "AS", "SOCTP","AN", "TS"),
  "S1" = c("pH", "TP", "DOC", "DON", "AP"),
  "S2" = c("SOCTN", "TNTP", "SOCTS", "TPTS", "TNTS"),
  "S3" = c("DOCAN", "DOCAP", "DOCAS", "ANAS", "APAS", "ANAP"),
  "S4" = c("CO2", "MBC", "qCO2")
)

regression_panel_labels <- c(
  "2c" = "Environmental_Regressions_Core_Soil_Properties",
  "S1" = "Multifunctionality_Correlations_pH_TP_DOC_DON_AP",
  "S2" = "Multifunctionality_Correlations_Bulk_Stoichiometry",
  "S3" = "Multifunctionality_Correlations_Available_Stoichiometry",
  "S4" = "Multifunctionality_Correlations_Microbial_Properties"
)

for (panel_tag in names(regression_panel_sets)) {
  vars_now <- regression_panel_sets[[panel_tag]]
  panel_ncol <- length(vars_now)
  panel_width <- 5 * panel_ncol

  panel_now <- wrap_plots(
    p_c_list[vars_now],
    ncol = panel_ncol,
    guides = "collect"
  ) & theme(legend.position = "top")

  ggsave(
    subfigure_file(panel_tag, regression_panel_labels[[panel_tag]]),
    panel_now,
    width = panel_width,
    height = 5.5,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )
}

# Fig. 2d
# [Bug Fix: cor.test zero variance] Prevents crash when features have 0 variance.
calc_correlations <- function(data, target, items, panel, log_items = FALSE) {
  map_dfr(items, function(item) {
    x <- data[[target]]
    y <- data[[item]]
    if (log_items) y <- log1p(y)
    keep <- is.finite(x) & is.finite(y)
    
    if (sum(keep) < 4 || sd(x[keep], na.rm = TRUE) == 0 || sd(y[keep], na.rm = TRUE) == 0) {
      return(tibble(Panel = panel, Item = item, r = NA_real_, p_value = NA_real_, n = sum(keep)))
    }
    
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    tibble(Panel = panel, Item = item, r = unname(test$estimate), p_value = test$p.value, n = sum(keep))
  }) %>%
    mutate(star = sig_star(p_value))
}

cor_cnps <- calc_correlations(dat, "CNPS_EMF", cnps_process_vars, "CNPS Multifunctionality (EMF)", TRUE)
cor_soil <- calc_correlations(
  dat,
  "Soil_nutrient_PC1",
  soil_pc1_vars,
  paste0("Soil nutrient PC1 (", round(soil_pc$variance[1]), "%)")
)

write.csv(
  bind_rows(cor_cnps, cor_soil),
  file.path(out_dir, "Fig2d_EMF33_and_11varPC1_Spearman_rawP_results.csv"),
  row.names = FALSE
)

heat_theme <- theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    strip.background = element_rect(fill = "white", color = "black", linewidth = line_w),
    strip.text.y.right = element_text(face = "bold", size = 16, angle = 270),
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 15)
  )

make_heatmap <- function(data, order, labels, legend = FALSE) {
  plot_data <- data %>% mutate(Item = factor(Item, levels = rev(order)), X = "")
  ggplot(plot_data, aes(X, Item)) +
    geom_tile(aes(fill = r), color = "black", linewidth = 0.55) +
    geom_text(aes(label = star), size = 5.3, fontface = "bold") +
    facet_grid(rows = vars(Panel), scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = labels) +
    scale_fill_gradient2(
      low = "#2C7BB6", mid = "white", high = "#D7191C",
      midpoint = 0, limits = c(-1, 1), oob = scales::squish, name = "Spearman r"
    ) +
    heat_theme +
    theme(legend.position = if (legend) "bottom" else "none")
}

cnps_order <- cor_cnps %>% arrange(desc(r)) %>% pull(Item)
cnps_labels <- setNames(str_replace(cnps_process_vars, "^([CNPS]):\\s*", "\\1: "), cnps_process_vars)

p_heat_cnps <- make_heatmap(cor_cnps, cnps_order, cnps_labels, TRUE)
p_heat_soil <- make_heatmap(cor_soil, soil_pc1_vars, soil_labels[soil_pc1_vars], FALSE)
p_d <- (p_heat_cnps | p_heat_soil) + plot_layout(widths = c(2.4, 1))

ggsave(subfigure_file("S12", "Multifunctionality_Correlation_Heatmaps"), p_d, width = 14, height = 16, device = cairo_pdf, bg = "white")


# Supplementary Fig. S5
s5_process_order <- sort(unique(cnps_process_vars))
s5_env_order <- soil_predictor_vars
s5_cycle_cols <- c(C = "#ADD8E8", N = "#EED1F9", P = "#ACD189", S = "#FCA894")

s5_sig_label <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# [Bug Fix: cor.test zero variance] Prevents crash on rare processes
s5_correlations <- map_dfr(s5_process_order, function(process_name) {
  map_dfr(s5_env_order, function(env_name) {
    x <- suppressWarnings(as.numeric(dat[[process_name]]))
    y <- suppressWarnings(as.numeric(dat[[env_name]]))
    keep <- is.finite(x) & is.finite(y)
    
    if (sum(keep) < 4 || sd(x[keep], na.rm = TRUE) == 0 || sd(y[keep], na.rm = TRUE) == 0) {
      return(tibble(Process = process_name, Environment = env_name, r = NA_real_, p_value = NA_real_, n = sum(keep)))
    }
    
    test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    tibble(Process = process_name, Environment = env_name, r = unname(test$estimate), p_value = test$p.value, n = sum(keep))
  })
}) %>%
  mutate(
    Element = stringr::str_extract(Process, "^[CNPS]"),
    Process_label = stringr::str_remove(Process, "^[CNPS]:\\s*"),
    significance = s5_sig_label(p_value)
  )

write.csv(
  s5_correlations,
  file.path(out_dir, "Supplementary_Fig_S5_CNPS33_environment25_Spearman_rawP_results.csv"),
  row.names = FALSE
)

s5_plot_data <- s5_correlations %>%
  mutate(
    Process_factor = factor(Process, levels = rev(s5_process_order)),
    Environment_factor = factor(Environment, levels = s5_env_order)
  )

s5_process_labels <- setNames(stringr::str_remove(s5_process_order, "^[CNPS]:\\s*"), s5_process_order)

p_s5_heatmap <- ggplot(s5_plot_data, aes(x = Environment_factor, y = Process_factor, fill = r)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = significance), size = 2.55, color = "black") +
  geom_hline(yintercept = c(8.5, 18.5, 25.5), color = "white", linewidth = 2.0) +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C",
    midpoint = 0, limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5),
    oob = scales::squish, name = "Spearman r",
    guide = guide_colorbar(title.position = "top", title.hjust = 0, barheight = grid::unit(30, "mm"), barwidth = grid::unit(4.5, "mm"))
  ) +
  scale_x_discrete(labels = soil_labels[s5_env_order], expand = expansion(add = 0)) +
  scale_y_discrete(labels = s5_process_labels, expand = expansion(add = 0)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(), panel.border = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(color = "black", size = 11.5, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "black", size = 11.5),
    legend.position = "left", legend.justification = c(0, 1),
    legend.title = element_text(face = "bold", size = 12.5), legend.text = element_text(size = 11),
    plot.margin = ggplot2::margin(t = 5, r = 2, b = 5, l = 2)
  )

s5_cycle_bounds <- tibble(
  Process = s5_process_order,
  Element = stringr::str_extract(s5_process_order, "^[CNPS]"),
  top_index = seq_along(s5_process_order),
  y_index = length(s5_process_order) - seq_along(s5_process_order) + 1
) %>%
  group_by(Element) %>%
  summarise(
    ymin = min(y_index) - 0.5, ymax = max(y_index) + 0.5, ymid = mean(c(min(y_index), max(y_index))),
    .groups = "drop"
  ) %>%
  mutate(Element = factor(Element, levels = c("C", "N", "P", "S")))

p_s5_cycle_strip <- ggplot(s5_cycle_bounds) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = Element), color = "grey35", linewidth = 0.55) +
  geom_text(aes(x = 0.5, y = ymid, label = paste0(Element, " cycle")), angle = 270, size = 4.0, fontface = "bold") +
  scale_fill_manual(values = s5_cycle_cols, drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(add = 0)) +
  scale_y_continuous(limits = c(0.5, length(s5_process_order) + 0.5), expand = expansion(add = 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = ggplot2::margin(t = 5, r = 2, b = 5, l = 0))

p_s5 <- (p_s5_heatmap | p_s5_cycle_strip) + plot_layout(widths = c(30, 1.25))

ggsave(subfigure_file("S13", "CNPS_Environment_Correlation_Heatmap"), p_s5, width = 17.5, height = 10.5, device = cairo_pdf, bg = "white")

# SEM input
sem_vars <- unique(c("sample", "place", "treatment", "Fire_factor", "CNPS_EMF", "Soil_nutrient_PC1", "pH", "qCO2", soil_pc1_vars))
write.csv(dat %>% select(all_of(sem_vars)), file.path(out_dir, "Fig2_SEM_input_table.csv"), row.names = FALSE)

# ============================================================
# Figure 3. CNPS functional composition and process responses
# ============================================================
out_dir <- file.path(tables_dir, "Figure_08_to_10_CNPS_Functions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_permutations <- 999L
random_seed <- 123L
line_w <- 0.8

treat_colors <- c("1UB" = "#3288BD", "B" = "#D53E4F")
treat_labels <- c("1UB" = "Unburned", "B" = "Burned")

element_order <- c("Carbon", "Nitrogen", "Phosphorus", "Sulfur")
group_order <- paste(element_order, "Cycle")

status_colors <- c(
  "Significant Increase" = "#D53E4F",
  "Significant Decrease" = "#3288BD",
  "Non-significant" = "#FFFFFF"
)

files <- list(
  metadata = file.path(input_dir, "group_ALL.txt"),
  C_matrix = file.path(input_dir, "TPM_C_ALL.tsv"),
  N_matrix = file.path(input_dir, "TPM_N_ALL.tsv"),
  P_matrix = file.path(input_dir, "TPM_P_ALL.tsv"),
  S_matrix = file.path(input_dir, "TPM_S_ALL.tsv"),
  C_class = file.path(input_dir, "C-classification.txt"),
  N_class = file.path(input_dir, "N-classification.txt"),
  P_class = file.path(input_dir, "P-classification.txt"),
  S_class = file.path(input_dir, "S-classification.txt")
)

base_theme_fig8 <- theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    axis.ticks = element_line(color = "black", linewidth = line_w),
    axis.ticks.length = grid::unit(0.25, "cm"),
    axis.text = element_text(color = "black", size = 18),
    axis.title = element_text(color = "black", face = "bold", size = 18),
    strip.background = element_rect(fill = "#F5F5F5", color = NA),
    strip.text = element_text(face = "bold", size = 18),
    plot.background = element_rect(fill = "white", color = NA)
  )

read_tpm_matrix <- function(path) {
  if(!file.exists(path)) return(NULL)
  x <- read.table(
    path, header = TRUE, sep = "\t", row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
  )
  x[] <- lapply(x, function(v) suppressWarnings(as.numeric(v)))
  as.data.frame(x, check.names = FALSE)
}

read_classification <- function(path) {
  if(!file.exists(path)) return(NULL)
  read.table(
    path, header = FALSE, sep = "\t",
    col.names = c("gene", "process"),
    stringsAsFactors = FALSE, quote = "", comment.char = ""
  )
}

clean_class <- function(x) {
  if(is.null(x)) return(NULL)
  x %>%
    mutate(
      gene = trimws(gene),
      process = trimws(gsub('"', "", process))
    ) %>%
    filter(
      !is.na(gene), gene != "",
      !is.na(process),
      !tolower(process) %in% c("", "others", "deleted", "unclassified")
    ) %>%
    distinct(gene, process)
}

aggregate_process <- function(mat, class_df, cycle_name, samples) {
  if(is.null(mat) || is.null(class_df)) return(tibble())
  mat[, samples, drop = FALSE] %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "TPM") %>%
    mutate(
      gene = trimws(gene),
      TPM = suppressWarnings(as.numeric(TPM))
    ) %>%
    inner_join(class_df, by = "gene") %>%
    group_by(sample, process) %>%
    summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    mutate(element_group = cycle_name)
}

meta <- read.table(
  files$metadata, header = TRUE, sep = "\t",
  check.names = FALSE, stringsAsFactors = FALSE
) %>%
  mutate(across(where(is.character), ~ trimws(gsub("\r", "", .)))) %>%
  filter(!is.na(sample), sample != "") %>%
  distinct(sample, .keep_all = TRUE)

matrices <- list(
  Carbon = read_tpm_matrix(files$C_matrix),
  Nitrogen = read_tpm_matrix(files$N_matrix),
  Phosphorus = read_tpm_matrix(files$P_matrix),
  Sulfur = read_tpm_matrix(files$S_matrix)
)

classes <- list(
  Carbon = clean_class(read_classification(files$C_class)),
  Nitrogen = clean_class(read_classification(files$N_class)),
  Phosphorus = clean_class(read_classification(files$P_class)),
  Sulfur = clean_class(read_classification(files$S_class))
)
if (!is.null(classes$Nitrogen)) {
  classes$Nitrogen <- classes$Nitrogen %>% mutate(process = recode_n_process(process))
}

common_samples <- Reduce(
  intersect,
  c(list(meta$sample), lapply(matrices, colnames))
)

meta_clean <- meta %>%
  filter(sample %in% common_samples) %>%
  mutate(
    treatment = factor(treatment, levels = c("1UB", "B")),
    place = factor(place)
  ) %>%
  filter(!is.na(treatment), !is.na(place)) %>%
  arrange(match(sample, common_samples))

common_samples <- meta_clean$sample
rownames(meta_clean) <- meta_clean$sample

df_all_cnps <- bind_rows(
  aggregate_process(matrices$Carbon, classes$Carbon, "Carbon Cycle", common_samples),
  aggregate_process(matrices$Nitrogen, classes$Nitrogen, "Nitrogen Cycle", common_samples),
  aggregate_process(matrices$Phosphorus, classes$Phosphorus, "Phosphorus Cycle", common_samples),
  aggregate_process(matrices$Sulfur, classes$Sulfur, "Sulfur Cycle", common_samples)
) %>%
  inner_join(meta_clean %>% select(sample, treatment, place), by = "sample") %>%
  filter(is.finite(TPM)) %>%
  mutate(element_group = factor(element_group, levels = group_order))

write.table(
  df_all_cnps,
  file.path(out_dir, "CNPS_process_level_TPM_ALL.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

calc_cycle_index <- function(data, cycle_name, samples) {
  wide <- data %>%
    filter(as.character(element_group) == cycle_name) %>%
    select(sample, process, TPM) %>%
    pivot_wider(names_from = process, values_from = TPM, values_fill = 0) %>%
    right_join(tibble(sample = samples), by = "sample") %>%
    arrange(match(sample, samples)) %>%
    mutate(across(-sample, ~ replace_na(as.numeric(.x), 0)))
  
  x <- log1p(as.matrix(wide %>% select(-sample)))
  
  keep <- apply(x, 2, sd, na.rm = TRUE) > 0
  
  if (sum(keep) < 2) stop("Fewer than two variable processes remain for ", cycle_name, ".")
  
  score <- rowMeans(scale(x[, keep, drop = FALSE]), na.rm = TRUE)
  names(score) <- wide$sample
  return(score)
}
# Fig. 3b: cycle functional indices
cycle_scores <- setNames(
  lapply(
    group_order,
    function(cycle_name) calc_cycle_index(df_all_cnps, cycle_name, common_samples)
  ),
  element_order
)

cycle_df <- tibble(sample = common_samples)
for (x in element_order) {
  cycle_df[[x]] <- as.numeric(cycle_scores[[x]][common_samples])
}

cycle_long <- cycle_df %>%
  inner_join(meta_clean %>% select(sample, treatment, place), by = "sample") %>%
  pivot_longer(
    cols = all_of(element_order),
    names_to = "Element_Function",
    values_to = "Score"
  ) %>%
  mutate(Element_Function = factor(Element_Function, levels = element_order))

cycle_lmm <- lapply(element_order, function(x) {
  d <- cycle_long %>% filter(Element_Function == x) %>% droplevels()
  fit <- lmer(Score ~ treatment + (1 | place), data = d, REML = TRUE)
  tab <- coef(summary(fit))
  
  # [Bug Fix: lmer subscript]
  p_val <- if ("Pr(>|t|)" %in% colnames(tab)) tab["treatmentB", "Pr(>|t|)"] else NA_real_
  df_val <- if ("df" %in% colnames(tab)) tab["treatmentB", "df"] else NA_real_
  t_val <- if ("t value" %in% colnames(tab)) tab["treatmentB", "t value"] else NA_real_
  
  tibble(
    Element_Function = x,
    estimate = tab["treatmentB", "Estimate"],
    SE = tab["treatmentB", "Std. Error"],
    df = df_val,
    t_value = t_val,
    p_value = p_val
  )
}) %>%
  bind_rows() %>%
  mutate(Element_Function = factor(Element_Function, levels = element_order))

write.csv(
  cycle_lmm,
  file.path(out_dir, "Fig3b_cycle_index_LMM_rawP.csv"),
  row.names = FALSE
)

cycle_labels <- cycle_lmm %>%
  left_join(
    cycle_long %>%
      group_by(Element_Function) %>%
      summarise(
        ymax = max(Score, na.rm = TRUE),
        ymin = min(Score, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Element_Function"
  ) %>%
  mutate(
    treatment = factor("B", levels = c("1UB", "B")),
    y = ymax + 0.12 * (ymax - ymin),
    label = vapply(p_value, format_p, character(1))
  )

p_cycle <- ggplot(cycle_long, aes(treatment, Score, fill = treatment)) +
  geom_violin(color = "black", linewidth = line_w, trim = FALSE) +
  geom_boxplot(
    width = 0.18, fill = "white", outlier.shape = NA,
    color = "black", linewidth = line_w
  ) +
  geom_jitter(
    width = 0.12, shape = 21, size = 1.6,
    fill = "white", color = "black", stroke = line_w, alpha = 0.75
  ) +
  geom_text(
    data = cycle_labels,
    aes(treatment, y, label = label),
    inherit.aes = FALSE,
    size = 5.8, fontface = "italic"
  ) +
  facet_wrap(~ Element_Function, ncol = 4) +
  scale_fill_manual(values = treat_colors) +
  scale_x_discrete(labels = treat_labels) +
  labs(x = NULL, y = "Cycle functional index (process-based Z score)") +
  coord_cartesian(clip = "off") +
  base_theme_fig8 +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

ggsave(
  subfigure_file("3b", "Cycle_Functional_Indices"),
  p_cycle, width = 13, height = 5, device = cairo_pdf, bg = "white"
)

# Fig. 3a: CPCoA
ordination_stats <- list()
ordination_plots <- list()

for (element_name in element_order) {
  mat <- t(matrices[[element_name]][, common_samples, drop = FALSE])
  mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  bray <- vegdist(mat, method = "bray")
  
  perm <- permute::how(nperm = n_permutations)
  permute::setBlocks(perm) <- meta_clean$place
  
  set.seed(random_seed)
  perm_res <- adonis2(bray ~ treatment, data = meta_clean, permutations = perm)
  
  disp <- betadisper(bray, group = meta_clean$treatment, type = "median", bias.adjust = TRUE)
  
  set.seed(random_seed)
  disp_res <- vegan::permutest(disp, permutations = perm)
  
  cap <- capscale(bray ~ treatment + Condition(place), data = meta_clean)
  scores_df <- as.data.frame(scores(cap, display = "sites")) %>%
    rownames_to_column("sample") %>%
    left_join(meta_clean %>% select(sample, treatment, place), by = "sample")
  
  cap1_percent <- 100 * cap$CCA$eig[1] / cap$tot.chi
  mds1_percent <- 100 * cap$CA$eig[1] / cap$tot.chi
  permanova_r2 <- as.numeric(perm_res$R2[1])
  permanova_p <- as.numeric(perm_res$`Pr(>F)`[1])
  
  ordination_stats[[element_name]] <- tibble(
    Element = element_name,
    CAP1_percent = cap1_percent,
    MDS1_percent = mds1_percent,
    PERMANOVA_R2 = permanova_r2,
    PERMANOVA_p = permanova_p,
    dispersion_F = as.numeric(disp_res$tab[1, "F"]),
    dispersion_p = as.numeric(disp_res$tab[1, "Pr(>F)"])
  )
  
  ordination_plots[[element_name]] <- ggplot(
    scores_df, aes(CAP1, MDS1, fill = treatment)
  ) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey65") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey65") +
    geom_point(shape = 21, size = 4.2, alpha = 0.88, color = "black", stroke = line_w) +
    annotate(
      "text", x = Inf, y = -Inf,
      label = paste0("PERMANOVA: R² = ", formatC(100 * permanova_r2, format = "f", digits = 2), "%, ", format_p(permanova_p)),
      hjust = 1.05, vjust = -0.4, size = 5.0, fontface = "italic"
    ) +
    scale_fill_manual(values = treat_colors, labels = treat_labels, name = "Treatment") +
    labs(
      title = paste(element_name, "Cycle"),
      x = paste0("CAP1 (", formatC(cap1_percent, format = "f", digits = 2), "%)"),
      y = paste0("MDS1 (", formatC(mds1_percent, format = "f", digits = 2), "%)")
    ) +
    base_theme_fig8 +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
}

write.csv(
  bind_rows(ordination_stats),
  file.path(out_dir, "Fig3a_CPCoA_PERMANOVA_rawP.csv"),
  row.names = FALSE
)

p_ordination <- wrap_plots(ordination_plots, ncol = 2) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(
  subfigure_file("3a", "CNPS_Ordination"),
  p_ordination, width = 14, height = 12, device = cairo_pdf, bg = "white"
)

# Fig. 3c: LMM Coefficient Forest Plot
lmm_results <- df_all_cnps %>%
  group_by(element_group, process) %>%
  group_modify(~{
    d <- droplevels(.x)
    if (n_distinct(d$treatment) < 2 || n_distinct(d$place) < 2) {
      return(tibble(estimate = NA_real_, SE = NA_real_, t_value = NA_real_, p_value = NA_real_))
    }
    
    fit <- tryCatch(
      lmer(log1p(TPM) ~ treatment + (1 | place), data = d, REML = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(tibble(estimate = NA_real_, SE = NA_real_, t_value = NA_real_, p_value = NA_real_))
    
    tab <- coef(summary(fit))
    # [Bug Fix: lmer subscript]
    p_val <- if ("Pr(>|t|)" %in% colnames(tab)) tab["treatmentB", "Pr(>|t|)"] else NA_real_
    t_val <- if ("t value" %in% colnames(tab)) tab["treatmentB", "t value"] else NA_real_
    
    tibble(
      estimate = tab["treatmentB", "Estimate"],
      SE       = tab["treatmentB", "Std. Error"],
      t_value  = t_val,
      p_value  = p_val
    )
  }) %>%
  ungroup() %>%
  mutate(
    ci_low = estimate - 1.96 * SE,
    ci_high = estimate + 1.96 * SE,
    sig_star = sig_star(p_value),
    sig_status = case_when(
      !is.na(p_value) & p_value < 0.05 & estimate > 0 ~ "Significant Increase",
      !is.na(p_value) & p_value < 0.05 & estimate < 0 ~ "Significant Decrease",
      TRUE ~ "Non-significant"
    ),
    element_group = factor(element_group, levels = group_order),
    sig_status = factor(sig_status, levels = c("Significant Increase", "Significant Decrease", "Non-significant"))
  ) %>%
  arrange(element_group, estimate)

write.csv(
  lmm_results,
  file.path(out_dir, "Fig3c_LMM_coefficients_rawP.csv"),
  row.names = FALSE
)

plot_lmm <- lmm_results %>% filter(is.finite(estimate), is.finite(ci_low), is.finite(ci_high))

p_forest <- ggplot(plot_lmm, aes(estimate, reorder(process, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_segment(aes(x = ci_low, xend = ci_high, y = reorder(process, estimate), yend = reorder(process, estimate)), linewidth = line_w, color = "black") +
  geom_point(aes(fill = sig_status), shape = 21, size = 3.8, color = "black", stroke = line_w) +
  geom_text(aes(x = ci_high, label = sig_star), hjust = -0.3, size = 6, fontface = "bold") +
  facet_grid(element_group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = status_colors, name = "Wildfire Response (LMM)", drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.15, 0.20))) +
  coord_cartesian(clip = "off") +
  labs(x = "LMM Estimate (Burned vs Unburned)", y = NULL) +
  base_theme_fig8 +
  theme(
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4),
    axis.text.y = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(5.5, 40, 5.5, 5.5)
  )

ggsave(
  subfigure_file("3c", "Process_Level_LMM_Effects"),
  p_forest, width = 12.4, height = 17, device = cairo_pdf, bg = "white"
)


# ============================================================
# Figure 5 and Supplementary Figures S9 and S14. CNPS network organization and rewiring
# ============================================================
out_dir <- file.path(tables_dir, "Figure_11_to_14_CNPS_Networks")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

safe_write_csv <- function(x, path, ...) {
  parent_dir <- dirname(path)
  dir.create(parent_dir, showWarnings = FALSE, recursive = TRUE)
  first_error <- NULL
  success <- tryCatch({readr::write_csv(x, path, ...); TRUE}, error = function(e) {first_error <<- conditionMessage(e); FALSE})
  if (success) return(invisible(path))
  
  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(basename(path))
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  alternative <- file.path(parent_dir, paste0(stem, "_", timestamp, ifelse(extension == "", "", paste0(".", extension))))
  tryCatch({
    readr::write_csv(x, alternative, ...)
    warning("Could not overwrite ", path, ".\nResults written to: ", alternative, call. = FALSE)
    invisible(alternative)
  }, error = function(e) stop("CSV output failed."))
}

cor_method           <- "spearman"
r_cutoff             <- 0.7            
q_cutoff_network     <- 0.001    
prevalence_threshold <- 0.20
trt_levels           <- c("Unburned", "Burned")

quick_test_mode      <- QUICK_NETWORK_TEST
n_mixing_null       <- if(quick_test_mode) 19 else 199
n_module_perm       <- if(quick_test_mode) 49 else 499
analysis_seed        <- 2026
min_module_plot_size <- 5
z_preservation_moderate <- 2
z_preservation_strong   <- 10

font_family  <- "Arial"
text_size_pt <- 18
text_size_mm <- text_size_pt / ggplot2::.pt
line_w       <- 1.0
tick_length  <- grid::unit(2.5, "mm")

el_cols  <- c(C = "#7DB9DE", N = "#D1B1F4", P = "#85C494", S = "#F4A482")
trt_cols <- c("Unburned" = "#4A90E2", "Burned" = "#E05353")

base_theme_net <- theme_bw(base_size = text_size_pt, base_family = font_family, base_line_size = line_w, base_rect_size = line_w) + 
  theme(
    text              = element_text(family = font_family, size = text_size_pt, color = "black"),
    axis.text         = element_text(color = "black", size = text_size_pt),
    axis.title        = element_text(color = "black", face = "bold", size = text_size_pt),
    axis.ticks        = element_line(color = "black", linewidth = line_w),
    axis.ticks.length = tick_length,
    plot.title        = element_text(face = "bold", size = text_size_pt, hjust = 0.5),
    plot.subtitle     = element_text(size = text_size_pt, color = "grey30", hjust = 0.5),
    plot.tag          = element_text(face = "bold", size = text_size_pt, color = "black"),
    legend.title      = element_text(face = "bold", size = text_size_pt),
    legend.text       = element_text(size = text_size_pt),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    strip.background  = element_rect(fill = "#F2F2F2", color = "black", linewidth = line_w), 
    strip.text        = element_text(face = "bold", size = text_size_pt, color = "black"),
    panel.spacing     = grid::unit(0.8, "lines"), 
    panel.grid.major  = element_blank(), 
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = line_w)
  )

find_input_file <- function(input_dir, pattern, label) {
  hits <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) return(NA)
  hits[1]
}

files <- list(
  C_abund_ALL = find_input_file(input_dir, "TPM_C_ALL.*\\.tsv$", "C abundance ALL"),
  N_abund_ALL = find_input_file(input_dir, "TPM_N_ALL.*\\.tsv$", "N abundance ALL"),
  P_abund_ALL = find_input_file(input_dir, "TPM_P_ALL.*\\.tsv$", "P abundance ALL"),
  S_abund_ALL = find_input_file(input_dir, "TPM_S_ALL.*\\.tsv$", "S abundance ALL"),
  meta_ALL    = find_input_file(input_dir, "group_ALL.*\\.txt$", "Metadata ALL"),
  C_map       = find_input_file(input_dir, "C-classification.*\\.txt$", "C mapping"),
  N_map       = find_input_file(input_dir, "N-classification.*\\.txt$", "N mapping"),
  P_map       = find_input_file(input_dir, "P-classification.*\\.txt$", "P mapping"),
  S_map       = find_input_file(input_dir, "S-classification.*\\.txt$", "S mapping")
)

meta_merged <- readr::read_tsv(files$meta_ALL, show_col_types = FALSE) %>%
  mutate(
    Treatment = case_when(
      str_detect(treatment, "1UB") ~ "Unburned",
      treatment == "B"              ~ "Burned",
      TRUE                          ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment)) %>%
  distinct(sample, Treatment) %>%
  rename(SampleID = sample) %>%
  mutate(Treatment = factor(Treatment, levels = trt_levels))

read_national_tpm <- function(file_all, valid_samples, prev_threshold = prevalence_threshold) {
  if (is.na(file_all)) return(tibble())
  df_all <- readr::read_tsv(file_all, show_col_types = FALSE)
  names(df_all)[1] <- "FeatureID"
  df_all <- df_all %>% mutate(FeatureID = as.character(FeatureID))
  
  sample_cols <- intersect(names(df_all), valid_samples)
  df_filtered <- df_all %>%
    select(FeatureID, all_of(sample_cols)) %>%
    mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.x)) %>% replace_na(0)))
  
  min_samples <- ceiling(prev_threshold * length(sample_cols))
  df_filtered <- df_filtered %>%
    filter(rowSums(select(., all_of(sample_cols)) > 0) >= min_samples)
  df_filtered
}

read_mapping <- function(path, element, remove_others = TRUE) {
  if (is.na(path)) return(tibble())
  mp <- readr::read_tsv(path, col_names = c("FeatureID", "Process"), show_col_types = FALSE) %>%
    mutate(
      FeatureID = as.character(FeatureID),
      Process   = as.character(Process),
      Process   = trimws(gsub('"', "", Process)),
      Element   = element
    ) %>%
    filter(!is.na(FeatureID), FeatureID != "", !is.na(Process), Process != "", !tolower(Process) %in% c("deleted", "unclassified"))
  
  if (remove_others) mp <- mp %>% filter(!tolower(Process) %in% c("na", "nan", "others", "other"))
  
  if (element == "N") {
    mp <- mp %>% mutate(Process = recode_n_process(Process))
  }
  mp %>% distinct(FeatureID, Element, Process)
}

clr_transform <- function(mat, pseudocount = NULL) {
  mat <- as.matrix(mat)
  if (is.null(pseudocount)) {
    min_pos <- suppressWarnings(min(mat[mat > 0], na.rm = TRUE))
    pseudocount = if (is.finite(min_pos)) min_pos / 2 else 1e-6
  }
  logx <- log(mat + pseudocount)
  sweep(logx, 1, rowMeans(logx, na.rm = TRUE), FUN = "-")
}

calc_cor_p <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  n <- nrow(mat)
  rmat <- suppressWarnings(cor(mat, method = method, use = "pairwise.complete.obs"))
  pmat <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat), dimnames = list(colnames(mat), colnames(mat)))
  
  # [Bug Fix: cor matrix precision] Cap rmat to prevent NaN in sqrt
  rmat_capped <- pmin(pmax(rmat, -1), 1)
  t_mat <- rmat_capped * sqrt((n - 2) / pmax((1 - rmat_capped^2), .Machine$double.eps))
  
  finite_r <- is.finite(rmat_capped) & abs(rmat_capped) < 1
  pmat[finite_r] <- 2 * pt(-abs(t_mat[finite_r]), df = n - 2)
  diag(rmat) <- 1; diag(pmat) <- 0
  list(r = rmat, p = pmat)
}

upper_pair_table <- function(rmat, pmat, gene_info, treatment) {
  mods <- colnames(rmat)
  if(length(mods) < 2) return(tibble())
  comb <- t(combn(mods, 2)) %>% as.data.frame(stringsAsFactors = FALSE)
  names(comb) <- c("Module1", "Module2")
  
  comb %>%
    mutate(
      Treatment = treatment,
      r = map2_dbl(Module1, Module2, ~ rmat[.x, .y]),
      p = map2_dbl(Module1, Module2, ~ pmat[.x, .y])
    ) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    left_join(gene_info %>% rename(Module1 = FeatureID, Element1 = Element, Process1 = Process), by = "Module1") %>%
    left_join(gene_info %>% rename(Module2 = FeatureID, Element2 = Element, Process2 = Process), by = "Module2") %>%
    mutate(
      Element_pair = map2_chr(as.character(Element1), as.character(Element2), ~ paste(sort(c(.x, .y)), collapse = "-")),
      Edge_type    = if_else(Element1 == Element2, "Within-element", "Cross-element"),
      Sign         = if_else(r >= 0, "Positive", "Negative")
    )
}

valid_nat_samples <- meta_merged$SampleID

raw_C <- read_national_tpm(files$C_abund_ALL, valid_nat_samples)
map_C <- read_mapping(files$C_map, "C", TRUE)
raw_N <- read_national_tpm(files$N_abund_ALL, valid_nat_samples)
map_N <- read_mapping(files$N_map, "N", TRUE)
raw_P <- read_national_tpm(files$P_abund_ALL, valid_nat_samples)
map_P <- read_mapping(files$P_map, "P", TRUE)
raw_S <- read_national_tpm(files$S_abund_ALL, valid_nat_samples)
map_S <- read_mapping(files$S_map, "S", TRUE)

gene_info <- bind_rows(map_C, map_N, map_P, map_S) %>% distinct(FeatureID, .keep_all = TRUE)

all_raw_wide <- list(raw_C, raw_N, raw_P, raw_S) %>%
  reduce(bind_rows) %>%
  filter(FeatureID %in% gene_info$FeatureID) %>%
  distinct(FeatureID, .keep_all = TRUE)

gene_cols <- all_raw_wide$FeatureID
mat_raw   <- all_raw_wide %>% column_to_rownames("FeatureID") %>% as.matrix() %>% t()
mat_raw   <- mat_raw[intersect(rownames(mat_raw), valid_nat_samples), ]
mat_clr   <- clr_transform(mat_raw)

clr_df <- as.data.frame(mat_clr, check.names = FALSE) %>%
  rownames_to_column("SampleID") %>%
  inner_join(meta_merged %>% select(SampleID, Treatment), by = "SampleID")

gene_info <- gene_info %>% filter(FeatureID %in% colnames(mat_clr))

# Network construction
net_resolution <- 1.0
net_weights    <- TRUE
net_seed       <- 42

all_longs <- list(); all_edges_raw <- list(); all_net_plots <- list()
all_cor_matrices <- list(); all_metrics <- list(); all_hubs <- list(); all_focal <- list()

for (trt in trt_levels) {
  mat_subset <- clr_df %>% filter(Treatment == trt) %>% select(all_of(gene_info$FeatureID))
  if (nrow(mat_subset) < 4) next
  
  c_res <- calc_cor_p(mat_subset, method = cor_method)
  all_cor_matrices[[trt]] <- c_res$r
  
  edge_raw <- upper_pair_table(c_res$r, c_res$p, gene_info, trt)
  all_longs[[trt]] <- edge_raw
  all_edges_raw[[trt]] <- edge_raw
  
  x <- edge_raw %>% filter(!is.na(r), abs(r) >= r_cutoff, q < q_cutoff_network)
  all_net_plots[[trt]] <- x
  
  if (nrow(x) == 0) {
    all_metrics[[trt]] <- tibble(Treatment = trt, Nodes = length(gene_info$FeatureID), Edges = 0, Average_degree = 0, Density = 0, Modularity = NA_real_, Mean_abs_r = NA_real_, Cross_element_edge_ratio = 0)
    all_hubs[[trt]]    <- tibble(Treatment = trt, Module = character(), Element = character(), Process = character(), Degree = numeric(), Betweenness = numeric())
    all_focal[[trt]]   <- gene_info %>% mutate(Treatment = trt, Focal_coupling = 0)
  } else {
    sub_vertices <- gene_info %>% filter(FeatureID %in% unique(c(x$Module1, x$Module2))) %>% rename(name = FeatureID)
    g <- graph_from_data_frame(x %>% transmute(from = Module1, to = Module2, weight = abs(r)), directed = FALSE, vertices = sub_vertices)
    
    if (!is.null(net_seed)) set.seed(net_seed)
    edge_w_vector <- if (net_weights) E(g)$weight else NULL
    cl <- cluster_leiden(g, objective_function = "modularity", weights = edge_w_vector, resolution_parameter = net_resolution)
    V(g)$module_id <- cl$membership
    
    all_metrics[[trt]] <- tibble(
      Treatment = trt, Nodes = vcount(g), Edges = ecount(g), Average_degree = mean(degree(g)), Density = edge_density(g),
      Modularity = if (ecount(g) > 0) modularity(g, cl$membership, weights = edge_w_vector) else NA_real_,
      Mean_abs_r = mean(abs(x$r)), Cross_element_edge_ratio = sum(x$Edge_type == "Cross-element") / ecount(g)
    )
    all_hubs[[trt]] <- tibble(Treatment = trt, Module = V(g)$name, Element = V(g)$Element, Process = V(g)$Process, Degree = degree(g), Betweenness = betweenness(g, normalized = TRUE))
    
    focal <- bind_rows(x %>% transmute(Module = Module1, abs_r = abs(r)), x %>% transmute(Module = Module2, abs_r = abs(r))) %>%
      group_by(Module) %>% summarise(Focal_coupling = mean(abs_r), .groups = "drop")
    all_focal[[trt]] <- gene_info %>% rename(Module = FeatureID) %>% mutate(Treatment = trt) %>% left_join(focal, by = "Module") %>% mutate(Focal_coupling = replace_na(Focal_coupling, 0))
  }
}

merged_longs     <- bind_rows(all_longs)
merged_net_plots <- bind_rows(all_net_plots) %>% mutate(Treatment = factor(Treatment, levels = trt_levels))
metrics_all      <- bind_rows(all_metrics) %>% mutate(Treatment = factor(Treatment, levels = trt_levels))
hubs_all         <- bind_rows(all_hubs) %>% mutate(Treatment = factor(Treatment, levels = trt_levels))
focal_all        <- bind_rows(all_focal) %>% mutate(Treatment = factor(Treatment, levels = trt_levels))

# Zi-Pi Analysis
z_threshold      <- 2.5
p_threshold      <- 0.62
all_zipi_results <- list()

for (trt in trt_levels) {
  trt_edges <- merged_net_plots %>% filter(Treatment == trt) %>% transmute(from = Module1, to = Module2, weight = abs(r))
  if (nrow(trt_edges) == 0) next
  
  trt_nodes <- gene_info %>% filter(FeatureID %in% unique(c(trt_edges$from, trt_edges$to))) %>% rename(name = FeatureID)
  g <- graph_from_data_frame(d = trt_edges, directed = FALSE, vertices = trt_nodes)
  
  if (!is.null(net_seed)) set.seed(net_seed)
  edge_w_vector <- if (net_weights) E(g)$weight else NULL
  cl <- cluster_leiden(g, objective_function = "modularity", weights = edge_w_vector, resolution_parameter = net_resolution)
  V(g)$module_id <- cl$membership
  
  node_list          <- V(g)$name
  module_assignments <- V(g)$module_id
  total_degree       <- degree(g)
  
  Zi <- rep(0, length(node_list)); Pi <- rep(0, length(node_list))
  
  for (i in seq_along(node_list)) {
    node_name         <- node_list[i]
    node_mod          <- module_assignments[i]
    ki                <- total_degree[i] 
    if (ki == 0) { Zi[i] <- 0; Pi[i] <- 0; next }
    
    neighbors_nodes   <- neighbors(g, node_name)$name
    neighbors_modules <- module_assignments[match(neighbors_nodes, node_list)]
    ki_raw            <- sum(neighbors_modules == node_mod)
    same_mod_nodes    <- node_list[which(module_assignments == node_mod)]
    
    k_inner_all <- sapply(same_mod_nodes, function(v) sum(module_assignments[match(neighbors(g, v)$name, node_list)] == node_mod))
    mean_k_inner <- mean(k_inner_all); sd_k_inner <- sd(k_inner_all)
    
    Zi[i] <- if (is.na(sd_k_inner) || sd_k_inner == 0) 0 else (ki_raw - mean_k_inner) / sd_k_inner
    
    mod_counts <- table(neighbors_modules)
    Pi[i]      <- 1 - sum((as.vector(mod_counts) / ki) ^ 2)
  }
  
  all_zipi_results[[trt]] <- tibble(
    GeneID = node_list, Treatment = trt, Element = V(g)$Element, Process = V(g)$Process,
    Module_Cluster = paste0("Cluster_", V(g)$module_id), Degree = as.numeric(total_degree[node_list]), Zi = Zi, Pi = Pi
  ) %>%
    left_join(focal_all %>% filter(Treatment == trt) %>% select(Module, Focal_coupling), by = c("GeneID" = "Module")) %>%
    mutate(
      Role = case_when(
        Zi <= z_threshold & Pi <= p_threshold ~ "Peripherals",
        Zi <= z_threshold & Pi >  p_threshold ~ "Inter-cluster Connectors",
        Zi >  z_threshold & Pi <= p_threshold ~ "Intra-cluster Hubs",
        Zi >  z_threshold & Pi >  p_threshold ~ "Network Hubs",
        TRUE                                  ~ "Peripherals"
      )
    )
}

merged_zipi <- bind_rows(all_zipi_results) %>%
  mutate(Treatment = factor(Treatment, levels = trt_levels), Role = factor(Role, levels = c("Peripherals", "Inter-cluster Connectors", "Intra-cluster Hubs", "Network Hubs")))

zipi_dir <- file.path(out_dir, "ZiPi_Analysis")
dir.create(zipi_dir, showWarnings = FALSE, recursive = TRUE)
safe_write_csv(merged_zipi, file.path(zipi_dir, "GCN_Zi_Pi_Role_Classification.csv"))

plot_network_multi <- function(out_file) {
  all_g_edges <- merged_net_plots %>% transmute(from = Module1, to = Module2, weight = abs(r))
  g_global    <- graph_from_data_frame(all_g_edges, directed = FALSE, vertices = gene_info %>% rename(name = FeatureID))
  set.seed(42)
  lay_mat    <- layout_with_fr(g_global, weights = E(g_global)$weight)
  layout_df  <- tibble(Module = V(g_global)$name, x = lay_mat[,1], y = lay_mat[,2])
  
  node_df <- tidyr::crossing(Treatment = factor(trt_levels, levels = trt_levels), gene_info %>% rename(Module = FeatureID)) %>%
    inner_join(layout_df, by = "Module") %>%
    left_join(focal_all %>% select(Treatment, Module, Focal_coupling), by = c("Treatment", "Module")) %>%
    left_join(merged_zipi %>% select(Module = GeneID, Treatment, Module_Cluster), by = c("Treatment", "Module")) %>%
    mutate(
      Module_Cluster = if_else(is.na(Module_Cluster) | Module_Cluster == "Loose Node", "Unclustered", as.character(Module_Cluster)),
      Module_Cluster = factor(Module_Cluster)
    )
  
  unique_mods <- sort(unique(as.character(node_df$Module_Cluster)))
  core_mods   <- unique_mods[unique_mods != "Unclustered"]
  
  sub_journal_palette <- c("#3A8FB7", "#E03C8A", "#439A86", "#E9A93A", "#7B5CB7", "#66A135", "#D65A31", "#3E517A", "#A13D63", "#4A7C59", "#C84B31", "#2D4059", "#8D93AB", "#E79E8D", "#516079")
  
  # [Bug Fix: Color Palette] Ensure enough colors for all clusters
  palette_extended <- if (length(core_mods) == 0) {
    character(0)
  } else if (length(core_mods) <= length(sub_journal_palette)) {
    sub_journal_palette[seq_along(core_mods)]
  } else {
    colorRampPalette(sub_journal_palette)(length(core_mods))
  }
  mod_colors <- setNames(palette_extended, core_mods)
  mod_colors["Unclustered"] <- "#DCDCDC"
  
  module_centers <- node_df %>%
    filter(Module_Cluster != "Unclustered") %>%
    group_by(Treatment, Module_Cluster) %>%
    summarise(cx = mean(x), cy = mean(y), .groups = "drop") %>%
    mutate(Short_Label = gsub("Cluster_", "C", Module_Cluster))
  
  eplot <- merged_net_plots %>%
    left_join(layout_df %>% rename(Module1 = Module, x1 = x, y1 = y), by = "Module1") %>%
    left_join(layout_df %>% rename(Module2 = Module, x2 = x, y2 = y), by = "Module2")
  
  p_net <- ggplot() +
    geom_curve(data = eplot, aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = abs(r)), color = "grey84", alpha = 0.20, curvature = 0.05) +
    geom_point(data = node_df, aes(x = x, y = y, fill = Module_Cluster, size = Focal_coupling), shape = 21, color = "black", stroke = line_w, alpha = 0.95) +
    geom_text_repel(data = module_centers, aes(x = cx, y = cy, label = Short_Label), family = font_family, size = text_size_mm, fontface = "bold", color = "black", bg.color = "white", bg.r = 0.15, seed = 42, segment.alpha = 0) +
    facet_wrap(~ Treatment, ncol = 2) +
    scale_fill_manual(values = mod_colors, guide = "none") + 
    scale_linewidth(range = c(0.12, 0.55), guide = "none") + 
    scale_size_continuous(range = c(1.5, 6.0), name = "Mean |r| (Focal Coupling)") +
    theme_void(base_size = text_size_pt, base_family = font_family) + 
    theme(text = element_text(family = font_family, size = text_size_pt), strip.text = element_text(face = "bold", size = text_size_pt, margin = ggplot2::margin(b = 10)), legend.position = "bottom", legend.title = element_text(size = text_size_pt, face = "bold"), legend.text = element_text(size = text_size_pt))
  
  top_hub_tbl   <- hubs_all %>% group_by(Treatment) %>% slice_max(Betweenness, n = 1, with_ties = FALSE) %>% transmute(Treatment, Top_hub = Module)
  metrics_table <- metrics_all %>% left_join(top_hub_tbl, by = "Treatment") %>%
    mutate(across(c(Mean_abs_r, Average_degree, Modularity, Cross_element_edge_ratio), ~ round(replace_na(.x, 0), 2))) %>%
    select(Treatment, Nodes, Edges, Mean_abs_r, Average_degree, Modularity, Top_hub)
  
  table_grob <- gridExtra::tableGrob(metrics_table, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = text_size_pt, base_family = font_family))
  p_final    <- p_net / patchwork::wrap_elements(full = table_grob) + patchwork::plot_layout(heights = c(8, 1.2))
  ggsave(out_file, p_final, width = 12, height = 8.5, device = cairo_pdf, family = font_family)
}

# Figure 11 export cancelled as requested.
# plot_network_multi(figure_file(11, "CNPS_Gene_Networks"))

role_colors <- c("Peripherals" = "#999999", "Inter-cluster Connectors" = "#E69F00", "Intra-cluster Hubs" = "#56B4E9", "Network Hubs" = "#D55E00")
max_zi <- max(c(z_threshold + 1, max(merged_zipi$Zi, na.rm = TRUE)))
min_zi <- min(c(-2, min(merged_zipi$Zi, na.rm = TRUE)))

p_scatter <- ggplot(merged_zipi, aes(x = Pi, y = Zi)) +
  annotate("rect", xmin = -Inf, xmax = p_threshold, ymin = z_threshold, ymax = Inf, fill = "#F7F9FC", alpha = 0.5) +
  annotate("rect", xmin = p_threshold, xmax = Inf, ymin = z_threshold, ymax = Inf, fill = "#FFF5EE", alpha = 0.5) +
  geom_vline(xintercept = p_threshold, linetype = "dashed", color = "grey30", linewidth = line_w) +
  geom_hline(yintercept = z_threshold, linetype = "dashed", color = "grey20", linewidth = line_w) +
  geom_point(aes(color = Role, size = Focal_coupling), alpha = 0.75, stroke = line_w) +
  facet_wrap(~ Treatment, ncol = 2) +
  scale_color_manual(values = role_colors, name = "Topological Role") +
  scale_size_continuous(range = c(0.8, 3.5), name = "Mean |r|") +
  scale_x_continuous(limits = c(-0.02, 0.95), breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(min_zi - 0.2, max_zi + 0.2)) +
  base_theme_net +
  labs(x = "Participation Coefficient (Pi)", y = "Within-module Degree Similarity (Zi)") +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = ggplot2::margin(t = -5, b = 0), panel.border = element_rect(color = "black", fill = NA, linewidth = line_w))

role_counts <- merged_zipi %>% count(Treatment, Role, .drop = FALSE) %>% pivot_wider(names_from = Treatment, values_from = n, values_fill = 0)
safe_write_csv(role_counts, file.path(zipi_dir, "A3_ZiPi_Role_Counts.csv"))

role_counts_long <- role_counts %>% pivot_longer(cols = all_of(trt_levels), names_to = "Treatment", values_to = "N_genes") %>% mutate(Treatment = factor(Treatment, levels = trt_levels), Role = factor(Role, levels = names(role_colors)))

p_role_counts <- ggplot(role_counts_long, aes(x = Role, y = N_genes, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.74), width = 0.66, color = "black", linewidth = line_w) +
  geom_text(aes(label = N_genes), position = position_dodge(width = 0.74), hjust = -0.12, family = font_family, size = text_size_mm) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = trt_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  base_theme_net + labs(x = NULL, y = "Number of functional genes", fill = NULL, title = "Redistribution of topological roles") + theme(legend.position = "top")

# Global Network Properties
prop_dir <- file.path(out_dir, "Network_Properties")
dir.create(prop_dir, showWarnings = FALSE, recursive = TRUE)
all_properties_list <- list()

for (trt in trt_levels) {
  x_sub <- merged_net_plots %>% filter(Treatment == trt)
  if(nrow(x_sub) == 0) next
  v_sub <- gene_info %>% filter(FeatureID %in% unique(c(x_sub$Module1, x_sub$Module2))) %>% rename(name = FeatureID)
  g     <- graph_from_data_frame(d = x_sub %>% transmute(from = Module1, to = Module2, weight = abs(r)), directed = FALSE, vertices = v_sub)
  N     <- vcount(g)
  E_num <- ecount(g)
  if (N < 10) next
  
  node_num         <- N
  edge_num         <- E_num
  edge_density_val <- edge_density(g)
  avg_degree       <- mean(degree(g))
  
  edge_w_vector <- E(g)$weight
  edge_distance <- 1 / pmax(edge_w_vector, .Machine$double.eps)
  
  local_weighted_clustering <- transitivity(g, type = "barrat", weights = edge_w_vector)
  avg_clustering <- mean(local_weighted_clustering, na.rm = TRUE)
  avg_path_length <- mean_distance(g, directed = FALSE, unconnected = TRUE, weights = edge_distance)
  net_diameter <- diameter(g, directed = FALSE, unconnected = TRUE, weights = edge_distance)
  
  d_dist  <- degree_distribution(g)
  degrees <- 0:(length(d_dist) - 1)
  nonzero <- d_dist > 0 & degrees > 0
  scale_free_r2 <- if (sum(nonzero) >= 3) summary(lm(log(d_dist[nonzero]) ~ log(degrees[nonzero])))$r.squared else 0
  
  if (!is.null(net_seed)) set.seed(net_seed)
  cl <- cluster_leiden(g, objective_function = "modularity", weights = edge_w_vector, resolution_parameter = net_resolution)
  modularity_real <- modularity(g, cl$membership, weights = edge_w_vector)
  
  rand_mods <- rep(0, 100); rand_cc <- rep(0, 100); rand_apl <- rep(0, 100)
  for(r_step in 1:100) {
    g_rand <- sample_gnm(n = N, m = E_num, directed = FALSE, loops = FALSE)
    if(vcount(g_rand) > 0 && ecount(g_rand) > 0) {
      
      # [Bug Fix: Sample scalar] Fixed incorrect sampling behaviour in R when length == 1
      w_vec <- E(g)$weight
      E(g_rand)$weight <- w_vec[sample(length(w_vec), size = ecount(g_rand), replace = TRUE)]
      
      rand_mods[r_step] <- tryCatch({
        cl_r <- cluster_leiden(g_rand, objective_function = "modularity", weights = E(g_rand)$weight, resolution_parameter = net_resolution)
        modularity(g_rand, cl_r$membership, weights = E(g_rand)$weight)
      }, error = function(e) 0)
      
      rand_cc[r_step] <- mean(transitivity(g_rand, type = "barrat", weights = E(g_rand)$weight), na.rm = TRUE)
      rand_apl[r_step] <- mean_distance(g_rand, directed = FALSE, unconnected = TRUE, weights = 1 / pmax(E(g_rand)$weight, .Machine$double.eps))
    }
  }
  
  mean_rand_mod <- mean(rand_mods[rand_mods > 0], na.rm = TRUE)
  mean_rand_cc  <- mean(rand_cc[!is.na(rand_cc)], na.rm = TRUE)
  mean_rand_apl <- mean(rand_apl[rand_apl > 0], na.rm = TRUE)
  
  relative_modularity <- if(!is.na(mean_rand_mod) && mean_rand_mod > 0) (modularity_real - mean_rand_mod) / mean_rand_mod else 0
  small_world_index <- if(replace_na(mean_rand_cc, 0) > 0 && replace_na(mean_rand_apl, 0) > 0 && !is.na(avg_clustering) && !is.na(avg_path_length)) (avg_clustering / mean_rand_cc) / (avg_path_length / mean_rand_apl) else 1
  
  all_properties_list[[trt]] <- tibble(
    Treatment = trt, `Node Number` = node_num, `Edge Number` = edge_num, `Edge Density` = edge_density_val,
    `Average Degree` = avg_degree, `Clustering Coefficient` = avg_clustering, `Average Path Length` = avg_path_length,
    `Network Diameter` = net_diameter, `Scale-free R2` = scale_free_r2, `Small-world Index` = small_world_index,
    `Modularity` = modularity_real, `Connectance` = E_num / (N * (N - 1) / 2), `Relative Modularity (RM)` = relative_modularity,
    `Complexity Index (Avg K)` = avg_degree
  )
}
df_properties <- bind_rows(all_properties_list) %>% mutate(Treatment = factor(Treatment, levels = trt_levels))
safe_write_csv(df_properties, file.path(prop_dir, "GCN_Macro_Global_Properties.csv"))

df_prop_long <- df_properties %>% pivot_longer(cols = -Treatment, names_to = "Property", values_to = "Value") %>% filter(Property %in% c("Edge Number", "Edge Density", "Clustering Coefficient", "Average Path Length", "Scale-free R2", "Small-world Index", "Relative Modularity (RM)", "Complexity Index (Avg K)")) %>% mutate(Category = case_when(Property %in% c("Edge Number", "Edge Density") ~ "1. Size & Density", Property %in% c("Clustering Coefficient", "Average Path Length") ~ "2. Connection Traits", Property %in% c("Scale-free R2", "Small-world Index") ~ "3. Structural Type", TRUE ~ "4. Complexity & Modularity"), Property = factor(Property, levels = c("Edge Number", "Edge Density", "Clustering Coefficient", "Average Path Length", "Scale-free R2", "Small-world Index", "Relative Modularity (RM)", "Complexity Index (Avg K)")), Text_Label = case_when(Property == "Edge Number" ~ sprintf("%.0f", Value), Property == "Edge Density" ~ sprintf("%.4f", Value), TRUE ~ sprintf("%.3f", Value)))

p_prop <- ggplot(df_prop_long, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_col(width = 0.52, color = "black", linewidth = line_w, alpha = 0.9) +
  geom_text(aes(label = Text_Label), position = position_stack(vjust = 0.5), family = font_family, color = "white", size = text_size_mm, fontface = "bold") +
  facet_wrap(Category ~ Property, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = trt_cols) +
  base_theme_net + labs(x = NULL, y = "Calculated Topological Value") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1, size = text_size_pt), strip.text = element_text(size = text_size_pt, face = "bold"), panel.spacing.y = grid::unit(1.5, "lines"), panel.spacing.x = grid::unit(1.0, "lines"), panel.border = element_rect(color = "grey85", fill = NA, linewidth = line_w))

# ============================================================
# Figure 5 and Supplementary Figures S9 and S14. Module composition, topology and rewiring
# ============================================================
## =============== 16. Module Functional Composition Analysis & Plot ======
## =========================================================================
message("Analyzing functional process composition within major network modules...")

func_dir <- file.path(out_dir, "Module_Functional_Composition")
dir.create(func_dir, showWarnings = FALSE, recursive = TRUE)

df_mod_composition <- merged_zipi %>%
  filter(!is.na(Module_Cluster), !Module_Cluster %in% c("Cluster_NA", "Loose Node")) %>%
  group_by(Treatment, Module_Cluster, Element, Process) %>%
  summarise(Node_Count = n(), .groups = "drop") %>%
  group_by(Treatment, Module_Cluster) %>%
  mutate(
    Module_Total_Nodes = sum(Node_Count),
    Proportion         = Node_Count / Module_Total_Nodes,
    Percent_Label      = sprintf("%.1f", Proportion * 100) 
  ) %>%
  ungroup()

safe_write_csv(df_mod_composition, file.path(func_dir, "GCN_Module_Functional_Process_Composition_All.csv"))

## =========================================================================
## =============== 16.2 Major-module functional heatmap ===================
## =========================================================================
message("Filtering major clusters and plotting publication-ready functional blueprints...")

min_module_size <- 10 

df_major_modules <- df_mod_composition %>%
  filter(Module_Total_Nodes >= min_module_size) %>% 
  mutate(
    Treatment      = factor(Treatment, levels = trt_levels),
    Module_Cluster = factor(Module_Cluster, levels = paste0("Cluster_", 1:100)),
    Process_Full   = paste0("[", Element, "] ", Process)
  ) %>%
  filter(!is.na(Module_Cluster))

module_process_levels <- df_major_modules %>%
  distinct(Element, Process, Process_Full) %>%
  arrange(factor(Element, levels = c("C", "N", "P", "S")), Process) %>%
  pull(Process_Full)

df_major_modules <- df_major_modules %>%
  mutate(
    Process_Full = factor(Process_Full, levels = rev(module_process_levels))
  )

p_func_heatmap <- ggplot(df_major_modules, aes(x = Module_Cluster, y = Process_Full)) +
  geom_tile(aes(fill = Proportion), color = "white", linewidth = line_w) +
  geom_text(
    aes(label = Percent_Label, color = Proportion > 0.45),
    family = font_family,
    size = text_size_mm,
    fontface = "bold",
    show.legend = FALSE
  ) +
  facet_wrap(~ Treatment, scales = "free_x", ncol = 2) +
  scale_fill_gradient(
    low = "#F7F9FC", high = "#1F77B4",
    labels = scales::percent,
    name = "Within-module\nnode proportion"
  ) +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey25")) +
  base_theme_net +
  labs(
    title = "Functional Process Blueprint of Major Network Modules",
    x     = "Inferred Major Network Modules (Original Cluster ID)",
    y     = "Soil Biogeochemical Processes (C-N-P-S)"
  ) +
  theme(
    axis.text.x  = element_text(
      angle = 45, hjust = 1, face = "bold", size = text_size_pt
    ),
    axis.text.y  = element_text(
      size = text_size_pt, family = font_family
    ),
    strip.text   = element_text(size = text_size_pt, face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w)
  )


ggsave(
  subfigure_file("S9", "Module_Functional_Composition"),
  p_func_heatmap,
  width = 20,
  height = 10,
  device = cairo_pdf,
  family = font_family,
  bg = "white",
  limitsize = FALSE
)

## =========================================================================

## =============== 17. Analysis 1: CNPS edge rewiring ======================
## =========================================================================
message("[Selected A1] Quantifying conserved, lost and gained CNPS edges...")

extended_dir <- file.path(out_dir, "Selected_CNPS_Coupling_Visual")
dir.create(extended_dir, showWarnings = FALSE, recursive = TRUE)

element_pair_levels <- c(
  "C-C", "C-N", "C-P", "C-S", "N-N",
  "N-P", "N-S", "P-P", "P-S", "S-S"
)

make_edge_key <- function(x, y) {
  paste(pmin(as.character(x), as.character(y)),
        pmax(as.character(x), as.character(y)), sep = "|||")
}

edge_u <- merged_net_plots %>%
  filter(Treatment == "Unburned") %>%
  transmute(
    PairID = make_edge_key(Module1, Module2),
    Gene1 = pmin(Module1, Module2),
    Gene2 = pmax(Module1, Module2),
    Element_pair,
    r_Unburned = r,
    w_Unburned = abs(r)
  ) %>%
  distinct(PairID, .keep_all = TRUE)

edge_b <- merged_net_plots %>%
  filter(Treatment == "Burned") %>%
  transmute(
    PairID = make_edge_key(Module1, Module2),
    Gene1 = pmin(Module1, Module2),
    Gene2 = pmax(Module1, Module2),
    Element_pair,
    r_Burned = r,
    w_Burned = abs(r)
  ) %>%
  distinct(PairID, .keep_all = TRUE)

edge_turnover <- full_join(
  edge_u, edge_b,
  by = c("PairID", "Gene1", "Gene2"),
  suffix = c("_U", "_B")
) %>%
  mutate(
    Element_pair = coalesce(Element_pair_U, Element_pair_B),
    Status = case_when(
      !is.na(w_Unburned) & !is.na(w_Burned) ~ "Conserved",
      !is.na(w_Unburned) &  is.na(w_Burned) ~ "Lost",
       is.na(w_Unburned) & !is.na(w_Burned) ~ "Gained",
      TRUE                                  ~ "Absent"
    ),
    w_Unburned = replace_na(w_Unburned, 0),
    w_Burned   = replace_na(w_Burned, 0),
    Delta_weight = w_Burned - w_Unburned
  ) %>%
  select(
    PairID, Gene1, Gene2, Element_pair, Status,
    r_Unburned, r_Burned, w_Unburned, w_Burned, Delta_weight
  )

element_sizes <- gene_info %>%
  count(Element, name = "N_features") %>%
  complete(Element = c("C", "N", "P", "S"), fill = list(N_features = 0))

possible_edges <- tibble(Element_pair = element_pair_levels) %>%
  separate(Element_pair, into = c("E1", "E2"), sep = "-", remove = FALSE) %>%
  left_join(element_sizes, by = c("E1" = "Element")) %>%
  rename(N1 = N_features) %>%
  left_join(element_sizes, by = c("E2" = "Element")) %>%
  rename(N2 = N_features) %>%
  mutate(
    Possible_edges = if_else(
      E1 == E2,
      N1 * (N1 - 1) / 2,
      N1 * N2
    )
  ) %>%
  select(Element_pair, Possible_edges)

rewiring_by_pair <- edge_turnover %>%
  filter(Status != "Absent") %>%
  count(Element_pair, Status, name = "N_edges") %>%
  complete(
    Element_pair = element_pair_levels,
    Status = c("Conserved", "Lost", "Gained"),
    fill = list(N_edges = 0)
  ) %>%
  left_join(possible_edges, by = "Element_pair") %>%
  group_by(Element_pair) %>%
  mutate(
    Union_edges = sum(N_edges),
    Fraction_of_union = if_else(Union_edges > 0, N_edges / Union_edges, 0),
    Fraction_of_possible = if_else(
      Possible_edges > 0, N_edges / Possible_edges, 0
    )
  ) %>%
  ungroup()

rewiring_indices <- edge_turnover %>%
  filter(Status != "Absent") %>%
  group_by(Element_pair) %>%
  summarise(
    Conserved = sum(Status == "Conserved"),
    Lost      = sum(Status == "Lost"),
    Gained    = sum(Status == "Gained"),
    Union     = n(),
    Jaccard   = if_else(Union > 0, Conserved / Union, NA_real_),
    Rewiring_index = 1 - Jaccard,
    Weighted_turnover = sum(abs(w_Burned - w_Unburned)) /
      sum(w_Burned + w_Unburned),
    .groups = "drop"
  ) %>%
  right_join(tibble(Element_pair = element_pair_levels), by = "Element_pair")

global_rewiring <- edge_turnover %>%
  filter(Status != "Absent") %>%
  summarise(
    Element_pair = "CNPS overall",
    Conserved = sum(Status == "Conserved"),
    Lost      = sum(Status == "Lost"),
    Gained    = sum(Status == "Gained"),
    Union     = n(),
    Jaccard   = Conserved / Union,
    Rewiring_index = 1 - Jaccard,
    Weighted_turnover = sum(abs(w_Burned - w_Unburned)) /
      sum(w_Burned + w_Unburned)
  )

status_cols <- c(
  "Conserved" = "#8A8A8A",
  "Lost"      = "#4A90E2",
  "Gained"    = "#E05353"
)

overall_turnover_plot_df <- global_rewiring %>%
  select(Element_pair, Conserved, Lost, Gained, Rewiring_index,
         Weighted_turnover) %>%
  pivot_longer(
    cols = c(Conserved, Lost, Gained),
    names_to = "Status",
    values_to = "N_edges"
  ) %>%
  mutate(
    Status = factor(Status, levels = c("Conserved", "Lost", "Gained")),
    Fraction = N_edges / sum(N_edges),
    Label = paste0(
      scales::comma(N_edges), "\n",
      scales::percent(Fraction, accuracy = 0.1)
    )
  )

overall_rewiring_subtitle <- paste0(
  "Rewiring index = ",
  sprintf("%.3f", global_rewiring$Rewiring_index),
  "; weighted turnover = ",
  sprintf("%.3f", global_rewiring$Weighted_turnover)
)

p_overall_turnover <- ggplot(
  overall_turnover_plot_df,
  aes(x = "CNPS overall", y = Fraction, fill = Status)
) +
  geom_col(width = 0.72, color = "black", linewidth = line_w) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    family = font_family,
    color = "white",
    fontface = "bold",
    size = text_size_mm,
    lineheight = 0.9
  ) +
  scale_fill_manual(values = status_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  coord_flip() +
  base_theme_net +
  labs(
    x = NULL,
    y = "Fraction of union-network edges",
    fill = "Edge fate"
  ) +
  ggtitle("Overall CNPS edge turnover", subtitle = overall_rewiring_subtitle) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top"
  )

rewiring_index_plot_df <- rewiring_indices %>%
  mutate(
    Net_edge_change = Gained - Lost,
    Net_direction = case_when(
      Net_edge_change > 0 ~ "Net gain",
      Net_edge_change < 0 ~ "Net loss",
      TRUE                ~ "Balanced"
    ),
    Element_pair = forcats::fct_reorder(
      Element_pair,
      Rewiring_index,
      .desc = FALSE
    ),
    Label = paste0(
      sprintf("%.2f", Rewiring_index),
      "  (", if_else(Net_edge_change > 0, "+", ""),
      Net_edge_change, ")"
    )
  )

net_change_cols <- c(
  "Net loss" = "#4A90E2",
  "Net gain" = "#E05353",
  "Balanced" = "#8A8A8A"
)

p_rewiring_index <- ggplot(
  rewiring_index_plot_df,
  aes(x = Element_pair, y = Rewiring_index, fill = Net_direction)
) +
  geom_col(width = 0.70, color = "black", linewidth = line_w) +
  geom_text(
    aes(label = Label),
    family = font_family,
    hjust = -0.08,
    size = text_size_mm
  ) +
  coord_flip() +
  scale_fill_manual(values = net_change_cols) +
  scale_y_continuous(
    limits = c(0, max(rewiring_index_plot_df$Rewiring_index, na.rm = TRUE) * 1.28),
    expand = expansion(mult = c(0, 0.02))
  ) +
  base_theme_net +
  labs(
    x = "CNPS element pair",
    y = "Rewiring index",
    fill = "Net edge balance"
  ) +
  ggtitle("Rewiring across CNPS element pairs") +
  theme(legend.position = "top")

## =========================================================================

## =============== 18. Analysis 2: CNPS element mixing ======================
## =========================================================================
message("[Selected A2] Calculating CNPS mixing matrices...")

build_treatment_graph <- function(trt) {
  x <- merged_net_plots %>% filter(Treatment == trt)
  v <- gene_info %>% rename(name = FeatureID)
  graph_from_data_frame(
    x %>% transmute(from = Module1, to = Module2, weight = abs(r)),
    directed = FALSE,
    vertices = v
  )
}

summarise_mixing <- function(g, treatment = NA_character_, iteration = NA_integer_) {
  if (ecount(g) == 0) {
    return(tibble(
      Treatment = treatment,
      Iteration = iteration,
      Element_pair = element_pair_levels,
      Edge_count = 0,
      Weight_sum = 0
    ))
  }
  
  edge_ends <- ends(g, E(g), names = FALSE)
  element_vec <- as.character(V(g)$Element)
  pair_type <- map2_chr(
    element_vec[edge_ends[, 1]],
    element_vec[edge_ends[, 2]],
    ~ paste(sort(c(.x, .y)), collapse = "-")
  )
  
  tibble(
    Treatment = treatment,
    Iteration = iteration,
    Element_pair = pair_type,
    Weight = E(g)$weight
  ) %>%
    group_by(Treatment, Iteration, Element_pair) %>%
    summarise(
      Edge_count = n(),
      Weight_sum = sum(Weight),
      .groups = "drop"
    ) %>%
    complete(
      Treatment,
      Iteration,
      Element_pair = element_pair_levels,
      fill = list(Edge_count = 0, Weight_sum = 0)
    )
}


mixing_observed <- list()
mixing_null <- list()

set.seed(analysis_seed)
for (trt in trt_levels) {
  g_obs <- build_treatment_graph(trt)
  
  obs_mix <- summarise_mixing(g_obs, trt, 0L) %>%
    left_join(possible_edges, by = "Element_pair") %>%
    mutate(
      Edge_density_normalized = Edge_count / Possible_edges,
      Weight_density_normalized = Weight_sum / Possible_edges
    )
  mixing_observed[[trt]] <- obs_mix
  
  null_list <- vector("list", n_mixing_null)
  
  for (b in seq_len(n_mixing_null)) {
    g_null <- igraph::rewire(
      g_obs,
      with = igraph::keeping_degseq(
        niter = max(10 * ecount(g_obs), 1000),
        loops = FALSE
      )
    )
    E(g_null)$weight <- sample(E(g_obs)$weight, ecount(g_null), replace = FALSE)
    
    null_list[[b]] <- summarise_mixing(g_null, trt, b)
    
    if (b %% 25 == 0 || b == n_mixing_null) {
      message("Mixing nulls: ", trt, " ", b, "/", n_mixing_null)
    }
  }
  
  null_mix <- bind_rows(null_list)
  mixing_null[[trt]] <- null_mix
  
}

mixing_observed_df <- bind_rows(mixing_observed)
mixing_null_df <- bind_rows(mixing_null)

mixing_null_summary <- mixing_null_df %>%
  group_by(Treatment, Element_pair) %>%
  summarise(
    Null_edge_mean = mean(Edge_count),
    Null_edge_sd = sd(Edge_count),
    Null_weight_mean = mean(Weight_sum),
    Null_weight_sd = sd(Weight_sum),
    .groups = "drop"
  )

mixing_enrichment <- mixing_observed_df %>%
  left_join(mixing_null_summary, by = c("Treatment", "Element_pair")) %>%
  mutate(
    Z_edge_enrichment = if_else(
      Null_edge_sd > 0,
      (Edge_count - Null_edge_mean) / Null_edge_sd,
      NA_real_
    ),
    Z_weight_enrichment = if_else(
      Null_weight_sd > 0,
      (Weight_sum - Null_weight_mean) / Null_weight_sd,
      NA_real_
    )
  )

mixing_delta <- mixing_enrichment %>%
  select(
    Treatment, Element_pair,
    Edge_density_normalized,
    Weight_density_normalized,
    Z_edge_enrichment,
    Z_weight_enrichment
  ) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(
      Edge_density_normalized,
      Weight_density_normalized,
      Z_edge_enrichment,
      Z_weight_enrichment
    )
  ) %>%
  mutate(
    Delta_edge_density =
      Edge_density_normalized_Burned - Edge_density_normalized_Unburned,
    Delta_weight_density =
      Weight_density_normalized_Burned - Weight_density_normalized_Unburned,
    Delta_Z_edge = Z_edge_enrichment_Burned - Z_edge_enrichment_Unburned,
    Delta_Z_weight = Z_weight_enrichment_Burned -
      Z_weight_enrichment_Unburned
  )

a1a2_summary <- bind_rows(rewiring_indices, global_rewiring) %>%
  mutate(Net_edge_change = Gained - Lost) %>%
  left_join(mixing_delta, by = "Element_pair")

# A single concise summary replaces the former edge-level and intermediate
# A1/A2 tables. The large edge-turnover object remains available in memory.
safe_write_csv(
  a1a2_summary,
  file.path(extended_dir, "A1A2_CNPS_Rewiring_Mixing_Summary.csv")
)

mixing_heatmap_df <- mixing_delta %>%
  separate(Element_pair, into = c("Element1", "Element2"), sep = "-") %>%
  bind_rows(
    mixing_delta %>%
      separate(Element_pair, into = c("Element2", "Element1"), sep = "-") %>%
      filter(Element1 != Element2)
  ) %>%
  mutate(
    Element1 = factor(Element1, levels = c("C", "N", "P", "S")),
    Element2 = factor(Element2, levels = c("C", "N", "P", "S"))
  )

mix_lim <- max(abs(mixing_heatmap_df$Delta_Z_weight), na.rm = TRUE)
if (!is.finite(mix_lim) || mix_lim == 0) mix_lim <- 1

p_mixing_delta <- ggplot(
  mixing_heatmap_df,
  aes(x = Element1, y = Element2, fill = Delta_Z_weight)
) +
  geom_tile(color = "white", linewidth = line_w) +
  geom_text(
    aes(label = sprintf("%.2f", Delta_Z_weight)),
    family = font_family,
    size = text_size_mm
  ) +
  scale_fill_gradient2(
    low = "#3B7FB6", mid = "white", high = "#D9544D",
    midpoint = 0, limits = c(-mix_lim, mix_lim),
    name = expression(Delta * " weighted enrichment")
  ) +
  coord_equal() +
  base_theme_net +
  labs(
    x = NULL, y = NULL,
    title = "Change in CNPS mixing enrichment"
  ) +
  theme(
    axis.text = element_text(face = "bold", size = text_size_pt)
  )

## =========================================================================

## =============== 22. Analysis 6: module preservation =====================
## =========================================================================
message("[Selected A6] Testing module overlap, preservation, split and fusion...")

module_wide <- merged_zipi %>%
  transmute(
    GeneID,
    Treatment = as.character(Treatment),
    Module = as.character(Module_Cluster)
  ) %>%
  pivot_wider(names_from = Treatment, values_from = Module) %>%
  filter(!is.na(Unburned), !is.na(Burned)) %>%
  mutate(
    Module_Unburned = paste0("UB_", Unburned),
    Module_Burned = paste0("B_", Burned)
  )

module_overlap_counts <- module_wide %>%
  count(Module_Unburned, Module_Burned, name = "Shared_nodes")

module_size_u <- module_wide %>%
  count(Module_Unburned, name = "Size_Unburned")
module_size_b <- module_wide %>%
  count(Module_Burned, name = "Size_Burned")
module_N <- nrow(module_wide)

# Aggregate very small modules only for the alluvial plot. All calculations
# and exported overlap statistics below retain the original module labels.
small_modules_u <- module_size_u %>%
  filter(Size_Unburned < min_module_plot_size) %>%
  pull(Module_Unburned)
small_modules_b <- module_size_b %>%
  filter(Size_Burned < min_module_plot_size) %>%
  pull(Module_Burned)

module_alluvial_counts <- module_wide %>%
  mutate(
    Module_Unburned_plot = if_else(
      Module_Unburned %in% small_modules_u,
      "UB_Small_modules",
      Module_Unburned
    ),
    Module_Burned_plot = if_else(
      Module_Burned %in% small_modules_b,
      "B_Small_modules",
      Module_Burned
    )
  ) %>%
  count(
    Module_Unburned_plot,
    Module_Burned_plot,
    name = "Shared_nodes"
  )

module_overlap_stats <- module_overlap_counts %>%
  left_join(module_size_u, by = "Module_Unburned") %>%
  left_join(module_size_b, by = "Module_Burned") %>%
  mutate(
    Jaccard = Shared_nodes /
      (Size_Unburned + Size_Burned - Shared_nodes),
    Overlap_fraction_Unburned = Shared_nodes / Size_Unburned,
    Overlap_fraction_Burned = Shared_nodes / Size_Burned,
    P_hypergeometric = phyper(
      Shared_nodes - 1,
      Size_Burned,
      module_N - Size_Burned,
      Size_Unburned,
      lower.tail = FALSE
    ),
    Q_BH = p.adjust(P_hypergeometric, method = "BH")
  )

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(sum(tab))
  
  if (total_pairs == 0) return(NA_real_)
  expected <- sum_rows * sum_cols / total_pairs
  max_index <- 0.5 * (sum_rows + sum_cols)
  if (max_index == expected) return(NA_real_)
  (sum_cells - expected) / (max_index - expected)
}

normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(
    pxy[nz] *
      log(pxy[nz] / (px[nz[, 1]] * py[nz[, 2]]))
  )
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

module_partition_similarity <- tibble(
  N_common_nodes = module_N,
  Adjusted_Rand_Index = adjusted_rand_index(
    module_wide$Module_Unburned,
    module_wide$Module_Burned
  ),
  Normalized_Mutual_Information = normalized_mutual_information(
    module_wide$Module_Unburned,
    module_wide$Module_Burned
  )
)

upper_abs_vector <- function(rmat, ids) {
  sub <- abs(rmat[ids, ids, drop = FALSE])
  sub[upper.tri(sub)]
}

rmat_u <- all_cor_matrices[["Unburned"]]
rmat_b <- all_cor_matrices[["Burned"]]
common_matrix_features <- intersect(colnames(rmat_u), colnames(rmat_b))

reference_modules <- split(
  module_wide$GeneID,
  module_wide$Module_Unburned
)

set.seed(analysis_seed + 20)
module_preservation <- map_dfr(names(reference_modules), function(mod_name) {
  ids <- intersect(reference_modules[[mod_name]], common_matrix_features)
  module_size <- length(ids)
  
  if (module_size < 4) {
    return(tibble(
      Reference_module = mod_name,
      Module_size = module_size,
      Adjacency_cor_observed = NA_real_,
      Mean_abs_r_Unburned = NA_real_,
      Mean_abs_r_Burned = NA_real_,
      Z_preservation = NA_real_,
      P_preservation = NA_real_,
      Preservation_class = "Too small"
    ))
  }
  
  vec_u <- upper_abs_vector(rmat_u, ids)
  vec_b <- upper_abs_vector(rmat_b, ids)
  observed_cor <- suppressWarnings(
    cor(vec_u, vec_b, method = "spearman", use = "pairwise.complete.obs")
  )
  
  null_cor <- replicate(n_module_perm, {
    random_b_ids <- sample(
      common_matrix_features,
      size = module_size,
      replace = FALSE
    )
    random_b_vec <- upper_abs_vector(rmat_b, random_b_ids)
    suppressWarnings(
      cor(
        vec_u,
        random_b_vec,
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )
  })
  
  null_mean <- mean(null_cor, na.rm = TRUE)
  null_sd <- sd(null_cor, na.rm = TRUE)
  z_pres <- if_else(
    is.finite(observed_cor) & is.finite(null_sd) & null_sd > 0,
    (observed_cor - null_mean) / null_sd,
    NA_real_
  )
  p_pres <- if (is.finite(observed_cor)) {
    (
      1 + sum(null_cor >= observed_cor, na.rm = TRUE)
    ) / (sum(is.finite(null_cor)) + 1)
  } else {
    NA_real_
  }
  
  tibble(
    Reference_module = mod_name,
    Module_size = module_size,
    Adjacency_cor_observed = observed_cor,
    Mean_abs_r_Unburned = mean(vec_u, na.rm = TRUE),
    Mean_abs_r_Burned = mean(vec_b, na.rm = TRUE),
    Z_preservation = z_pres,
    P_preservation = p_pres,
    Preservation_class = case_when(
      is.na(z_pres) ~ "Unresolved",
      z_pres >= z_preservation_strong   ~ "Strong",
      z_pres >= z_preservation_moderate ~ "Moderate",
      TRUE          ~ "Not preserved"
    )
  )
})

# Remove obsolete A6 outputs from earlier script versions so that the result
# folder contains only the integrated figure and one concise summary table.
legacy_a6_outputs <- file.path(
  extended_dir,
  c(
    "A6_Module_Overlap_Jaccard.csv",
    "A6_Module_Partition_Similarity.csv",
    "A6_Module_Preservation_Permutation.csv",
    "A6_Module_Split_Fusion_Alluvial.pdf",
    "A6_Module_Preservation.pdf"
  )
)
invisible(file.remove(legacy_a6_outputs[file.exists(legacy_a6_outputs)]))

# One module-level table replaces the former three A6 tables. It records the
# strongest cross-treatment match, split/fusion extent and preservation test.
best_module_overlap <- module_overlap_stats %>%
  group_by(Module_Unburned) %>%
  slice_max(order_by = Jaccard, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    Reference_module = Module_Unburned,
    Best_Burned_module = Module_Burned,
    Shared_nodes_best = Shared_nodes,
    Best_Jaccard = Jaccard,
    Best_overlap_fraction_of_reference = Overlap_fraction_Unburned,
    Best_overlap_fraction_of_burned = Overlap_fraction_Burned,
    Best_match_Q_BH = Q_BH
  )

module_split_counts <- module_overlap_counts %>%
  count(Module_Unburned, name = "N_Burned_destinations") %>%
  rename(Reference_module = Module_Unburned)

module_fusion_counts <- module_overlap_counts %>%
  count(Module_Burned, name = "N_Unburned_sources")

a6_module_summary <- module_preservation %>%
  left_join(best_module_overlap, by = "Reference_module") %>%
  left_join(module_split_counts, by = "Reference_module") %>%
  left_join(
    module_fusion_counts,
    by = c("Best_Burned_module" = "Module_Burned")
  ) %>%
  mutate(
    Adjusted_Rand_Index =
      module_partition_similarity$Adjusted_Rand_Index,
    Normalized_Mutual_Information =
      module_partition_similarity$Normalized_Mutual_Information,
    N_common_nodes = module_partition_similarity$N_common_nodes
  ) %>%
  select(
    Reference_module, Module_size, N_Burned_destinations,
    Best_Burned_module, N_Unburned_sources,
    Shared_nodes_best, Best_Jaccard,
    Best_overlap_fraction_of_reference,
    Best_overlap_fraction_of_burned, Best_match_Q_BH,
    Adjacency_cor_observed,
    Mean_abs_r_Unburned, Mean_abs_r_Burned,
    Z_preservation, P_preservation, Preservation_class,
    Adjusted_Rand_Index, Normalized_Mutual_Information,
    N_common_nodes
  )

safe_write_csv(
  a6_module_summary,
  file.path(extended_dir, "A6_Module_Reorganization_Summary.csv")
)

p_module_alluvial <- ggplot(
  module_alluvial_counts,
  aes(
    axis1 = Module_Unburned_plot,
    axis2 = Module_Burned_plot,
    y = Shared_nodes
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = Module_Unburned_plot),
    width = 0.17, alpha = 0.65
  ) +
  ggalluvial::geom_stratum(
    width = 0.17, fill = "grey92",
    color = "black", linewidth = line_w
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    aes(
      label = stringr::str_remove(
        after_stat(stratum),
        "^(UB_|B_)"
      )
    ),
    family = font_family,
    size = text_size_mm
  ) +
  scale_x_discrete(
    limits = c("Unburned modules", "Burned modules"),
    expand = c(0.08, 0.08)
  ) +
  scale_fill_brewer(palette = "Set3") +
  base_theme_net +
  labs(
    x = NULL,
    y = "Number of shared functional features",
    title = "Module splitting and fusion after wildfire"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

major_u_order <- module_size_u %>%
  filter(Size_Unburned >= min_module_plot_size) %>%
  arrange(readr::parse_number(Module_Unburned), Module_Unburned) %>%
  pull(Module_Unburned)

major_b_order <- module_size_b %>%
  filter(Size_Burned >= min_module_plot_size) %>%
  arrange(readr::parse_number(Module_Burned), Module_Burned) %>%
  pull(Module_Burned)

module_overlap_heatmap <- tidyr::expand_grid(
  Module_Unburned = major_u_order,
  Module_Burned = major_b_order
) %>%
  left_join(
    module_overlap_stats %>%
      select(
        Module_Unburned, Module_Burned,
        Shared_nodes, Jaccard, Q_BH
      ),
    by = c("Module_Unburned", "Module_Burned")
  ) %>%
  mutate(
    Shared_nodes = replace_na(Shared_nodes, 0L),
    Jaccard = replace_na(Jaccard, 0),
    Module_Unburned_label = factor(
      stringr::str_remove(Module_Unburned, "^UB_"),
      levels = stringr::str_remove(major_u_order, "^UB_")
    ),
    Module_Burned_label = factor(
      stringr::str_remove(Module_Burned, "^B_"),
      levels = rev(stringr::str_remove(major_b_order, "^B_"))
    ),
    Cell_label = if_else(
      Shared_nodes > 0,
      paste0(sprintf("%.2f", Jaccard), "\n(n=", Shared_nodes, ")"),
      ""
    )
  )

overlap_fill_max <- max(
  module_overlap_heatmap$Jaccard,
  na.rm = TRUE
)
if (!is.finite(overlap_fill_max) || overlap_fill_max <= 0) {
  overlap_fill_max <- 1
}

p_module_overlap <- ggplot(
  module_overlap_heatmap,
  aes(
    x = Module_Unburned_label,
    y = Module_Burned_label,
    fill = Jaccard
  )
) +
  geom_tile(color = "white", linewidth = line_w) +
  geom_text(
    aes(label = Cell_label),
    family = font_family,
    size = text_size_mm,
    lineheight = 0.9
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#5B3A8E",
    limits = c(0, overlap_fill_max),
    name = "Jaccard"
  ) +
  coord_equal() +
  base_theme_net +
  labs(
    x = "Unburned reference module",
    y = "Burned module",
    title = "Cross-treatment module overlap",
    subtitle = "Cell labels: Jaccard (shared features)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

preservation_plot_df <- module_preservation %>%
  filter(is.finite(Z_preservation)) %>%
  mutate(
    Reference_module_label =
      stringr::str_remove(Reference_module, "^UB_"),
    Z_label = sprintf("%.2f", Z_preservation)
  )

p_module_preservation <- ggplot(
  preservation_plot_df,
  aes(
    x = reorder(Reference_module_label, Z_preservation),
    y = Z_preservation,
    fill = Preservation_class
  )
) +
  geom_hline(
    yintercept = c(z_preservation_moderate, z_preservation_strong),
    linetype = "dashed",
    color = "grey45",
    linewidth = line_w
  ) +
  geom_col(width = 0.68, color = "black", linewidth = line_w) +
  geom_text(
    aes(
      label = Z_label,
      hjust = if_else(Z_preservation >= 0, -0.10, 1.10)
    ),
    family = font_family,
    size = text_size_mm,
    show.legend = FALSE
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Strong" = "#2A9D8F",
      "Moderate" = "#E9C46A",
      "Not preserved" = "#E76F51",
      "Unresolved" = "grey70",
      "Too small" = "grey85"
    )
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.18))) +
  base_theme_net +
  labs(
    x = "Unburned reference module",
    y = "Permutation-based preservation Z",
    fill = "Preservation",
    title = "Preservation of within-module coupling",
    subtitle = paste0(
      "Dashed thresholds: Z = ",
      z_preservation_moderate, " and ", z_preservation_strong
    )
  ) +
  theme(legend.position = "top")

partition_plot_df <- module_partition_similarity %>%
  select(Adjusted_Rand_Index, Normalized_Mutual_Information) %>%
  pivot_longer(
    everything(),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      Adjusted_Rand_Index = "Adjusted Rand index",
      Normalized_Mutual_Information =
        "Normalized mutual information"
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Adjusted Rand index",
        "Normalized mutual information"
      )
    ),
    Value_label = sprintf("%.3f", Value)
  )

partition_y_min <- min(
  0,
  partition_plot_df$Value,
  na.rm = TRUE
)
partition_y_max <- max(
  1,
  partition_plot_df$Value,
  na.rm = TRUE
)
partition_y_padding <- max(
  0.08,
  0.08 * (partition_y_max - partition_y_min)
)

p_partition_similarity <- ggplot(
  partition_plot_df,
  aes(x = Metric, y = Value, fill = Metric)
) +
  geom_hline(yintercept = 0, color = "grey45", linewidth = line_w) +
  geom_col(width = 0.58, color = "black", linewidth = line_w) +
  geom_text(
    aes(
      label = Value_label,
      vjust = if_else(Value >= 0, -0.45, 1.35)
    ),
    family = font_family,
    fontface = "bold",
    size = text_size_mm
  ) +
  scale_fill_manual(
    values = c(
      "Adjusted Rand index" = "#5B8FF9",
      "Normalized mutual information" = "#61B15A"
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(
      partition_y_min - partition_y_padding,
      partition_y_max + partition_y_padding
    ),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  base_theme_net +
  labs(
    x = NULL,
    y = "Partition similarity",
    title = "Overall module-partition similarity",
    subtitle = paste0(
      "Common features = ", module_N,
      "; modules: ", nrow(module_size_u),
      " unburned to ", nrow(module_size_b), " burned"
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.subtitle = element_text(size = text_size_pt)
  )

## =========================================================================


# Final Figure 5 and supplementary network panels

# Figure 5a: global network architecture
p_fig5a <- p_prop +
  labs(tag = "a", title = "Wildfire-induced changes in global network architecture") +
  theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(8, 8, 12, 8)
  )

ggsave(
  subfigure_file("5a", "Global_Network_Architecture"),
  p_fig5a,
  width = 16,
  height = 9,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Figure 5b: overall edge turnover
p_fig5b <- p_overall_turnover +
  labs(tag = "b", title = "Overall edge turnover") +
  theme(legend.position = "bottom")

ggsave(
  subfigure_file("5b", "Overall_Edge_Turnover"),
  p_fig5b,
  width = 9,
  height = 6,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Figure 5c: Zi-Pi distribution and role counts
p_fig5c_left <- p_scatter +
  labs(title = "Topological roles of CNPS functional genes") +
  theme(
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 12, 8, 8)
  )

p_fig5c_right <- p_role_counts +
  theme(
    legend.position = "top",
    plot.margin = ggplot2::margin(8, 8, 8, 12)
  )

p_fig5c <- (p_fig5c_left | p_fig5c_right) +
  patchwork::plot_layout(widths = c(2.15, 1)) +
  patchwork::plot_annotation(tag_levels = NULL)

ggsave(
  subfigure_file("5c", "ZiPi_Roles"),
  p_fig5c,
  width = 20,
  height = 9,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Figure 5d: element-pair rewiring
p_fig5d <- p_rewiring_index +
  labs(tag = "d", title = "Element-pair rewiring") +
  theme(legend.position = "bottom")

ggsave(
  subfigure_file("5d", "Element_Pair_Rewiring"),
  p_fig5d,
  width = 11,
  height = 6,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Figure 5e: change in element mixing
p_fig5e <- p_mixing_delta +
  labs(tag = "e", title = "Change in element mixing") +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = grid::unit(4.2, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  )

ggsave(
  subfigure_file("5e", "Mixing_DeltaZWE"),
  p_fig5e,
  width = 8,
  height = 6,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Figure 5f: cross-treatment module overlap
p_fig5f <- p_module_overlap +
  labs(tag = "f", title = "Cross-treatment module overlap") +
  theme(legend.position = "bottom") +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = grid::unit(3.8, "cm"),
      barheight = grid::unit(0.45, "cm")
    )
  )

ggsave(
  subfigure_file("5f", "Module_Overlap"),
  p_fig5f,
  width = 9,
  height = 7,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Figure 5g: within-module coupling preservation
p_fig5g <- p_module_preservation +
  labs(tag = "g", title = "Within-module coupling preservation") +
  theme(legend.position = "bottom")

ggsave(
  subfigure_file("5g", "Module_Preservation"),
  p_fig5g,
  width = 10,
  height = 7,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)

# Supplementary Figure S11 cancelled.

# Supplementary Figure S14: module splitting and fusion
p_figS14 <- p_module_alluvial +
  labs(title = "Module splitting and fusion")

ggsave(
  subfigure_file("S14", "Module_Alluvial"),
  p_figS14,
  width = 13,
  height = 8,
  device = cairo_pdf,
  family = font_family,
  limitsize = FALSE
)


# ============================================================
# Figure 4 and Supplementary Figure S11. XC temporal CNPS dynamics
# ============================================================
out_dir <- file.path(tables_dir, "Figure_15_Temporal_CNPS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260729)

line_w <- 0.8
n_permutations <- 9999

time_levels <- c(
  "0week", "1week", "6months",
  "712months", "824months", "936months"
)
postfire_times <- time_levels[-1]

time_labels <- c(
  "0week" = "Prefire",
  "1week" = "1 week",
  "6months" = "6 months",
  "712months" = "1 year",
  "824months" = "2 years",
  "936months" = "3 years"
)

metric_order <- c("C cycle", "N cycle", "P cycle", "S cycle")

facet_colors <- c(
  "C cycle" = "#ADD8E8",
  "N cycle" = "#EED1F9",
  "P cycle" = "#ACD189",
  "S cycle" = "#FCA894"
)

point_colors <- c(
  "Increase ***" = "#B03040",
  "Increase **" = "#D53E4F",
  "Increase *" = "#E88B96",
  "Decrease ***" = "#224D70",
  "Decrease **" = "#377EB8",
  "Decrease *" = "#7EB1D6",
  "Non-significant" = "#FFFFFF",
  "Prefire" = "#FFFFFF"
)

sig_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

clean_text <- function(data) {
  data %>%
    mutate(across(where(is.character), ~ trimws(gsub("\r", "", .x))))
}

base_theme_fig15 <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
    axis.ticks = element_line(color = "black", linewidth = line_w),
    axis.ticks.length = unit(0.25, "cm"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 11),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(face = "bold", size = 13)
  )

apply_strip_colors <- function(grob, colors, metrics, right_strip = FALSE) {
  pattern <- if (right_strip) "strip-r|strip-right" else "strip-t|strip-top|strip"
  strip_rows <- which(grepl(pattern, grob$layout$name))

  edit_fill <- function(node, fill_color) {
    if (inherits(node, "rect")) {
      node$gp$fill <- fill_color
      node$gp$col <- "black"
      node$gp$lwd <- 1
      return(node)
    }
    if (!is.null(node$children)) {
      node$children <- lapply(node$children, edit_fill, fill_color = fill_color)
    }
    if (!is.null(node$grobs)) {
      node$grobs <- lapply(node$grobs, edit_fill, fill_color = fill_color)
    }
    node
  }

  for (i in seq_along(strip_rows)) {
    if (i > length(metrics)) break
    fill_color <- unname(colors[metrics[i]])
    if (is.na(fill_color)) fill_color <- "#FFFFFF"
    idx <- strip_rows[i]
    grob$grobs[[idx]] <- edit_fill(grob$grobs[[idx]], fill_color)
  }

  grob
}

# Input data -----------------------------------------------------------------

group <- read.delim(
  "group_SX.txt",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  clean_text() %>%
  drop_na(sample, place, time, treatment) %>%
  filter(place == "XC") %>%
  mutate(
    time = factor(time, levels = time_levels),
    treatment = factor(
      treatment,
      levels = c("1_Unburntsoil", "2_Highburntsoil")
    ),
    burn_status = factor(
      if_else(treatment == "1_Unburntsoil", "Unburnt", "Burnt"),
      levels = c("Unburnt", "Burnt")
    )
  ) %>%
  filter(
    (time == "0week" & treatment == "1_Unburntsoil") |
      (time != "0week" & treatment == "2_Highburntsoil")
  ) %>%
  distinct(sample, .keep_all = TRUE)

if (anyNA(group$time) || anyNA(group$treatment)) {
  stop("Unknown time or treatment labels were detected.")
}

summary_files <- c(
  "C cycle" = "TPM_summary_C_SX.txt",
  "N cycle" = "TPM_summary_N_SX.txt",
  "P cycle" = "TPM_summary_P_SX.txt",
  "S cycle" = "TPM_summary_S_SX.txt"
)

read_summary <- function(file, metric) {
  read.delim(
    file,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    clean_text() %>%
    filter(place == "XC", !grepl("others", process, ignore.case = TRUE)) %>%
    mutate(TPM = as.numeric(TPM), Metric = metric)
}

summary_list <- Map(read_summary, unname(summary_files), names(summary_files))
names(summary_list) <- names(summary_files)

common_samples <- Reduce(
  intersect,
  c(list(group$sample), lapply(summary_list, function(x) unique(x$sample)))
)

group <- group %>%
  filter(sample %in% common_samples) %>%
  arrange(time, sample)

analysis_samples <- group$sample
rownames(group) <- group$sample

summary_list <- lapply(summary_list, function(x) {
  x %>%
    filter(sample %in% analysis_samples) %>%
    select(-any_of(c("time", "treatment"))) %>%
    left_join(group %>% select(sample, time, treatment), by = "sample")
})

df_all <- bind_rows(summary_list) %>%
  mutate(Metric = factor(Metric, levels = metric_order))

if (nrow(df_all %>% distinct(Metric, process)) != 33L) {
  stop("The four input tables must contain 33 non-Others CNPS processes in total.")
}

sample_counts <- group %>% count(time, treatment, name = "n")
if (any(sample_counts$n != 10)) {
  warning("At least one sampling stage does not contain 10 aligned samples.")
}

# Panel a: PCoA and pairwise PERMANOVA ---------------------------------------

p_palette <- RColorBrewer::brewer.pal(9, "Blues")[3:8]
names(p_palette) <- time_levels

pcoa_plots <- list()
pcoa_stats <- list()
shared_legend <- NULL

for (metric in metric_order) {
  process_matrix <- summary_list[[metric]] %>%
    group_by(sample, process) %>%
    summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = process, values_from = TPM, values_fill = 0) %>%
    column_to_rownames("sample")

  process_matrix <- process_matrix[analysis_samples, , drop = FALSE]
  process_matrix <- process_matrix[, colSums(process_matrix) > 0, drop = FALSE]

  if (ncol(process_matrix) < 2) {
    stop("Fewer than two non-zero processes remain for ", metric, ".")
  }

  bray <- vegan::vegdist(process_matrix, method = "bray")
  pcoa <- ape::pcoa(bray, correction = "cailliez")
  coordinates <- if (!is.null(pcoa$vectors.cor)) pcoa$vectors.cor else pcoa$vectors
  eigenvalues <- pcoa$values$Rel_corr_eig
  if (is.null(eigenvalues) || all(!is.finite(eigenvalues))) {
    eigenvalues <- pcoa$values$Relative_eig
  }

  pcoa_data <- data.frame(
    sample = rownames(coordinates),
    PCoA1 = coordinates[, 1],
    PCoA2 = coordinates[, 2]
  ) %>%
    left_join(group %>% select(sample, time, treatment, burn_status), by = "sample") %>%
    mutate(Metric = metric)

  bray_matrix <- as.matrix(bray)

  pcoa_stats[[metric]] <- bind_rows(lapply(postfire_times, function(stage) {
    keep <- group %>%
      filter(time %in% c("0week", stage)) %>%
      pull(sample)

    meta_pair <- group[keep, , drop = FALSE] %>%
      mutate(
        comparison = factor(
          if_else(time == "0week", "Prefire", "Postfire"),
          levels = c("Prefire", "Postfire")
        )
      )

    pair_distance <- as.dist(bray_matrix[keep, keep, drop = FALSE])

    permanova <- vegan::adonis2(
      pair_distance ~ comparison,
      data = meta_pair,
      permutations = n_permutations
    )

    dispersion <- vegan::permutest(
      vegan::betadisper(pair_distance, meta_pair$comparison),
      permutations = n_permutations
    )

    data.frame(
      Metric = metric,
      Comparison = paste0(time_labels[[stage]], " vs Prefire"),
      Time = stage,
      n_prefire = sum(meta_pair$comparison == "Prefire"),
      n_postfire = sum(meta_pair$comparison == "Postfire"),
      pseudo_F = permanova$F[1],
      R2 = permanova$R2[1],
      p_PERMANOVA = permanova$`Pr(>F)`[1],
      p_dispersion = dispersion$tab$`Pr(>F)`[1]
    )
  }))

  p <- ggplot(pcoa_data, aes(PCoA1, PCoA2)) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60", linewidth = 0.5) +
    geom_point(
      aes(fill = time, shape = burn_status),
      size = 4.2,
      alpha = 0.88,
      color = "black",
      stroke = line_w
    ) +
    facet_wrap(~ Metric) +
    scale_fill_manual(
      values = p_palette,
      breaks = time_levels,
      labels = time_labels,
      name = "Sampling time",
      guide = guide_legend(
        override.aes = list(shape = 21, size = 4.2, stroke = line_w, color = "black")
      )
    ) +
    scale_shape_manual(
      values = c("Unburnt" = 21, "Burnt" = 24),
      labels = c("Unburnt" = "UB", "Burnt" = "B"),
      name = "Treatment"
    ) +
    labs(
      x = paste0("PCoA1 (", round(eigenvalues[1] * 100, 2), "%)"),
      y = paste0("PCoA2 (", round(eigenvalues[2] * 100, 2), "%)")
    ) +
    base_theme_fig15

  if (is.null(shared_legend)) {
    shared_legend <- cowplot::get_legend(
      p + theme(
        legend.position = "center",
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 13)
      )
    )
  }

  pcoa_plots[[metric]] <- apply_strip_colors(
    ggplotGrob(p + theme(legend.position = "none")),
    facet_colors,
    metric
  )
}

pcoa_stats <- bind_rows(pcoa_stats) %>%
  mutate(
    significance = sig_stars(p_PERMANOVA),
    dispersion = if_else(
      p_dispersion < 0.05,
      "Different dispersion",
      "No detected dispersion difference"
    )
  )

write.csv(pcoa_stats, file.path(out_dir, "Fig5_PCoA_pairwise_PERMANOVA.csv"), row.names = FALSE)

pcoa_grid <- plot_grid(plotlist = pcoa_plots, ncol = 2, align = "hv")

# Panel b: cycle functional indices ------------------------------------------

calc_cycle_index <- function(data, metric, samples) {
  wide <- data %>%
    filter(as.character(Metric) == metric) %>%
    group_by(sample, process) %>%
    summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = process, values_from = TPM, values_fill = 0) %>%
    right_join(tibble(sample = samples), by = "sample") %>%
    arrange(match(sample, samples)) %>%
    mutate(across(-sample, ~ replace_na(as.numeric(.x), 0)))

  x <- as.matrix(wide %>% select(-sample))
  keep <- apply(log1p(x), 2, sd, na.rm = TRUE) > 0

  if (sum(keep) < 2) {
    stop("Fewer than two variable processes remain for ", metric, ".")
  }

  score <- rowMeans(scale(log1p(x[, keep, drop = FALSE])), na.rm = TRUE)
  names(score) <- wide$sample
  score
}

cycle_scores <- data.frame(sample = analysis_samples)
for (metric in metric_order) {
  cycle_scores[[metric]] <- calc_cycle_index(df_all, metric, analysis_samples)[analysis_samples]
}

write.csv(
  cycle_scores %>% left_join(group %>% select(sample, time, treatment), by = "sample"),
  file.path(out_dir, "Fig5_process_based_cycle_indices.csv"),
  row.names = FALSE
)

cycle_data <- cycle_scores %>%
  left_join(group %>% select(sample, time, treatment), by = "sample") %>%
  pivot_longer(all_of(metric_order), names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = metric_order),
    time = factor(time, levels = time_levels)
  )

cycle_tests <- bind_rows(lapply(metric_order, function(metric) {
  prefire <- cycle_data %>%
    filter(Metric == metric, time == "0week") %>%
    pull(Value)
  prefire <- prefire[is.finite(prefire)]

  bind_rows(lapply(postfire_times, function(stage) {
    postfire <- cycle_data %>%
      filter(Metric == metric, time == stage) %>%
      pull(Value)
    postfire <- postfire[is.finite(postfire)]

    p <- suppressWarnings(
      wilcox.test(postfire, prefire, paired = FALSE, exact = FALSE)$p.value
    )

    mean_prefire <- mean(prefire)
    mean_postfire <- mean(postfire)
    direction <- ifelse(mean_postfire > mean_prefire, "Increase", "Decrease")
    star <- sig_stars(p)

    data.frame(
      Metric = metric,
      time = stage,
      n_prefire = length(prefire),
      n_postfire = length(postfire),
      mean_prefire = mean_prefire,
      mean_postfire = mean_postfire,
      p_value = p,
      sig_star = star,
      fill_key = ifelse(star == "", "Non-significant", paste(direction, star))
    )
  }))
})) %>%
  mutate(time = factor(time, levels = time_levels))

write.csv(cycle_tests, file.path(out_dir, "Fig5_cycle_level_tests.csv"), row.names = FALSE)

cycle_summary <- cycle_data %>%
  group_by(treatment, Metric, time) %>%
  summarise(
    n = sum(is.finite(Value)),
    mean = mean(Value, na.rm = TRUE),
    SE = if_else(n > 1, sd(Value, na.rm = TRUE) / sqrt(n), 0),
    .groups = "drop"
  )

plot_cycle <- bind_rows(
  cycle_summary %>% filter(time == "0week"),
  cycle_summary %>% filter(time != "0week", treatment == "2_Highburntsoil")
) %>%
  left_join(cycle_tests %>% select(Metric, time, fill_key), by = c("Metric", "time")) %>%
  mutate(fill_key = if_else(time == "0week", "Prefire", fill_key))

p_cycle <- ggplot(plot_cycle, aes(time, mean, group = Metric)) +
  geom_errorbar(aes(ymin = mean - SE, ymax = mean + SE), width = 0.15, linewidth = line_w) +
  geom_line(linewidth = line_w, alpha = 0.8) +
  geom_point(
    aes(fill = fill_key),
    shape = 24,
    color = "black",
    size = 4.2,
    stroke = line_w
  ) +
  scale_fill_manual(values = point_colors, guide = "none", drop = FALSE) +
  scale_x_discrete(labels = time_labels) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = "Cycle functional index\n(process-based Z score; Mean +/- SE)"
  ) +
  base_theme_fig15 +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1))

cycle_grob <- apply_strip_colors(
  ggplotGrob(p_cycle),
  facet_colors,
  metric_order
)

# Panel c: process-level lnRR heatmap ----------------------------------------

process_tests <- list()
combinations <- df_all %>% distinct(Metric, process)

for (i in seq_len(nrow(combinations))) {
  metric <- combinations$Metric[i]
  process <- combinations$process[i]

  data_now <- df_all %>% filter(Metric == metric, process == !!process)
  prefire <- data_now %>% filter(time == "0week") %>% pull(TPM)
  prefire <- prefire[is.finite(prefire)]
  prefire_mean <- mean(prefire)

  if (length(prefire) < 2 || !is.finite(prefire_mean) || prefire_mean <= 0) next

  for (stage in postfire_times) {
    postfire <- data_now %>% filter(time == stage) %>% pull(TPM)
    postfire <- postfire[is.finite(postfire)]
    postfire_mean <- mean(postfire)

    if (length(postfire) < 2 || !is.finite(postfire_mean) || postfire_mean <= 0) next

    p <- suppressWarnings(
      wilcox.test(postfire, prefire, paired = FALSE, exact = FALSE)$p.value
    )

    process_tests[[length(process_tests) + 1]] <- data.frame(
      Metric = metric,
      process = process,
      time = stage,
      n_prefire = length(prefire),
      n_postfire = length(postfire),
      mean_prefire = prefire_mean,
      mean_postfire = postfire_mean,
      lnRR = log(postfire_mean / prefire_mean),
      p_value = p,
      sig_star = sig_stars(p)
    )
  }
}

heatmap_data <- bind_rows(process_tests) %>%
  mutate(
    time = factor(time, levels = postfire_times),
    Metric = factor(Metric, levels = metric_order)
  )

write.csv(heatmap_data, file.path(out_dir, "Fig5_process_lnRR_and_tests.csv"), row.names = FALSE)

lnrr_limit <- quantile(abs(heatmap_data$lnRR), 0.95, na.rm = TRUE, names = FALSE)
if (!is.finite(lnrr_limit) || lnrr_limit <= 0) {
  stop("Unable to calculate the lnRR color scale.")
}

p_heat <- ggplot(heatmap_data, aes(time, process)) +
  geom_tile(aes(fill = lnRR), color = "grey88", linewidth = 0.5) +
  geom_text(aes(label = sig_star), size = 5.2, fontface = "bold", vjust = 0.75) +
  facet_grid(Metric ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = time_labels[postfire_times]) +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "#F7F7F7",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-lnrr_limit, lnrr_limit),
    oob = scales::squish,
    name = "lnRR"
  ) +
  labs(x = NULL, y = NULL) +
  base_theme_fig15 +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.position = "left",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11)
  )

heat_legend <- cowplot::get_legend(p_heat)
heat_grob <- apply_strip_colors(
  ggplotGrob(p_heat + theme(legend.position = "none")),
  facet_colors,
  metric_order,
  right_strip = TRUE
)

# Final Figure 4 and Supplementary Figure S11 -------------------------------

# Supplementary Figure S11: XC CNPS PCoA with shared legend
p_figS11 <- plot_grid(
  pcoa_grid,
  shared_legend,
  ncol = 2,
  rel_widths = c(1, 0.18)
)

ggsave(
  subfigure_file("S11", "XC_CNPS_PCoA"),
  p_figS11,
  width = 15.5,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

# Figure 4: cycle-level temporal responses above process-level lnRR heatmap
figure4_bottom <- plot_grid(
  heat_legend,
  heat_grob,
  ncol = 2,
  rel_widths = c(0.18, 1)
)

figure4 <- plot_grid(
  cycle_grob,
  figure4_bottom,
  ncol = 1,
  rel_heights = c(0.82, 1.18),
  labels = c("(a)", "(b)"),
  label_size = 20,
  label_fontface = "bold"
)

ggsave(
  figure_file(4, "XC_CNPS_Temporal_Functions"),
  figure4,
  width = 15,
  height = 18,
  device = cairo_pdf,
  bg = "white"
)


# ============================================================
# Supplementary Figures S5-S8. Site-level CNPS process responses
# ============================================================
out_dir <- file.path(tables_dir, "Figure_17_to_20_Site_Level_CNPS")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

line_w <- 0.8
sites <- c("GH","HH","OQ","JZ","BD","DL","HZ","SM","YT","CZ","ZH")
treat_cols <- c(Unburned="#3288BD", Burned="#D53E4F")
cycle_cols <- c(C="#ADD8E8", N="#EED1F9", P="#ACD189", S="#FCA894")

process_order <- list(
  C=c("Labile C degradation","Cellulose degradation","Hemicellulose degradation",
      "Pectin degradation","Necromass degradation","Oxidative C degradation",
      "C substrate binding","Glycan biosynthesis"),
  N=c("Anammox","Assimilatory nitrate reduction","Denitrification",
      "Dissimilatory nitrate reduction","N-fixation","Nitrification",
      "Organic synthesis & degradation"),
  P=c("Transporters","Phosphonate and phosphinate metabolism","Purine metabolism",
      "Pyrimidine metabolism","Organic phosphoesters hydrolysis","Two component system",
      "Pentose phosphate pathway","Pyruvate metabolism","Oxidative phosphorylation",
      "Phosphotransferase system"),
  S=c("Organic sulfur transformation",
      "Link between inorganic and organic sulfur transformation",
      "Sulfur reduction","Dissimilatory sulfur reduction and oxidation",
      "Sulfur oxidation","Assimilatory sulfate reduction",
      "SOX systems","Sulfur disproportionation")
)

recode_treatment <- function(x) case_when(
  as.character(x) %in% c("1UB","1_Unburntsoil","Unburned","Unburnt") ~ "Unburned",
  as.character(x) %in% c("B","2_Highburntsoil","Burned","Burnt") ~ "Burned",
  TRUE ~ NA_character_
)

standardize_cycle <- function(d, cycle) {
  if (cycle == "N") {
    d <- d %>% mutate(process = recode_n_process(process))
  }

  d %>%
    mutate(
      TPM = suppressWarnings(as.numeric(TPM)),
      treatment = recode_treatment(treatment),
      place = factor(place, levels = sites),
      treatment = factor(treatment, levels = c("Unburned","Burned")),
      process = factor(process, levels = process_order[[cycle]]),
      log_TPM = log10(TPM + 1)
    ) %>%
    filter(is.finite(TPM), !is.na(place), !is.na(treatment), !is.na(process))
}

read_cycle_summary <- function(file, cycle) {
  read.table(
    file, header = TRUE, sep = "\t",
    check.names = FALSE, stringsAsFactors = FALSE
  ) %>%
    mutate(across(where(is.character), ~ trimws(gsub("\r","",.)))) %>%
    standardize_cycle(cycle)
}

site_wilcox <- function(d, cycle) {
  d %>%
    group_by(place, process) %>%
    group_modify(~{
      ub <- .x$TPM[.x$treatment == "Unburned"]
      b  <- .x$TPM[.x$treatment == "Burned"]
      p <- if (length(ub) >= 2 && length(b) >= 2) {
        suppressWarnings(
          wilcox.test(b, ub, paired = FALSE, exact = FALSE, correct = FALSE)$p.value
        )
      } else {
        NA_real_
      }

      tibble(
        Cycle = cycle,
        n_unburned = length(ub),
        n_burned = length(b),
        mean_unburned = if(length(ub)) mean(ub) else NA_real_,
        mean_burned = if(length(b)) mean(b) else NA_real_,
        p_value = p
      )
    }) %>%
    ungroup() %>%
    mutate(
      Direction = case_when(
        mean_burned > mean_unburned ~ "Increase",
        mean_burned < mean_unburned ~ "Decrease",
        TRUE ~ "No change"
      ),
      sig_star = sig_star(p_value)
    )
}

cycle_theme <- function(cycle) {
  theme_bw(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = line_w),
      strip.background = element_rect(
        fill = cycle_cols[[cycle]], color = "black", linewidth = 1
      ),
      strip.text = element_text(face = "bold")
    )
}

plot_cycle <- function(d, cycle, tests, site_h, supplement_tag) {
  sm <- d %>%
    group_by(place, process, treatment) %>%
    summarise(
      n = n(),
      mean_TPM = mean(TPM, na.rm = TRUE),
      se_TPM = if_else(n > 1, sd(TPM, na.rm = TRUE) / sqrt(n), 0),
      .groups = "drop"
    )

  ann <- sm %>%
    group_by(place, process) %>%
    summarise(
      y = max(mean_TPM + se_TPM, na.rm = TRUE) * 1.12,
      .groups = "drop"
    ) %>%
    left_join(
      tests %>% select(place, process, sig_star, Direction),
      by = c("place", "process")
    )

  dodge <- position_dodge(.76)

  p2 <- ggplot(sm, aes(place, mean_TPM, fill = treatment)) +
    geom_col(
      position = dodge, width = .68,
      color = "black", linewidth = line_w
    ) +
    geom_errorbar(
      aes(ymin = pmax(mean_TPM - se_TPM, 0), ymax = mean_TPM + se_TPM),
      position = dodge, width = .20, linewidth = line_w
    ) +
    geom_text(
      data = ann,
      aes(place, y, label = sig_star, color = Direction),
      inherit.aes = FALSE,
      size = 6, fontface = "bold"
    ) +
    facet_wrap(~ process, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = treat_cols, name = "Fire treatment") +
    scale_color_manual(
      values = c(Increase="#D53E4F", Decrease="#3288BD", `No change`="black"),
      guide = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, .18))) +
    labs(x = "Sampling site", y = "Abundance (TPM)") +
    cycle_theme(cycle) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )

  ggsave(
    subfigure_file(supplement_tag, paste0(cycle, "_Site_Barplot")),
    p2,
    width = 18, height = site_h,
    device = cairo_pdf, bg = "white"
  )
}

df_c <- read_cycle_summary("TPM_summary_C_ALL.txt", "C")
df_n <- read_cycle_summary("TPM_summary_N_ALL.txt", "N")
df_p <- read_cycle_summary("TPM_summary_P_ALL.txt", "P")
df_s <- read_cycle_summary("TPM_summary_S_ALL.txt", "S")

site_tests <- bind_rows(
  site_wilcox(df_c, "C"),
  site_wilcox(df_n, "N"),
  site_wilcox(df_p, "P"),
  site_wilcox(df_s, "S")
)

write.csv(
  site_tests,
  file.path(out_dir, "CNPS_11SITE_Wilcoxon_rawP.csv"),
  row.names = FALSE
)

plot_cycle(df_c, "C", filter(site_tests, Cycle == "C"), 16, "S5")
plot_cycle(df_n, "N", filter(site_tests, Cycle == "N"), 16, "S6")
plot_cycle(df_p, "P", filter(site_tests, Cycle == "P"), 20, "S7")
plot_cycle(df_s, "S", filter(site_tests, Cycle == "S"), 16, "S8")


# Generate final figure manifest
figure_manifest <- tibble::tribble(
  ~Figure, ~File, ~Description,
  "1a", basename(subfigure_file("1a", "Sampling_Sites_Map")), "National sampling sites",
  "1b", basename(subfigure_file("1b", "Environmental_Responses")), "Environmental LMM responses",
  "2a", basename(subfigure_file("2a", "CNPS_Multifunctionality")), "CNPS multifunctionality",
  "2b", basename(subfigure_file("2b", "Random_Forest_Importance")), "Random-forest importance",
  "2c", basename(subfigure_file("2c", "Environmental_Regressions_Core_Soil_Properties")), "Environmental regressions: core soil properties",
  "3a", basename(subfigure_file("3a", "CNPS_Ordination")), "CNPS ordination",
  "3b", basename(subfigure_file("3b", "Cycle_Functional_Indices")), "Cycle functional indices",
  "3c", basename(subfigure_file("3c", "Process_Level_LMM_Effects")), "Process-level LMM effects",
  "4",  basename(figure_file(4, "XC_CNPS_Temporal_Functions")), "XC temporal CNPS functional responses",
  "5a", basename(subfigure_file("5a", "Global_Network_Architecture")), "Global network architecture",
  "5b", basename(subfigure_file("5b", "Overall_Edge_Turnover")), "Overall edge turnover",
  "5c", basename(subfigure_file("5c", "ZiPi_Roles")), "Zi-Pi roles and role counts",
  "5d", basename(subfigure_file("5d", "Element_Pair_Rewiring")), "Element-pair rewiring",
  "5e", basename(subfigure_file("5e", "Mixing_DeltaZWE")), "Element-pair weighted mixing enrichment",
  "5f", basename(subfigure_file("5f", "Module_Overlap")), "Cross-treatment module overlap",
  "5g", basename(subfigure_file("5g", "Module_Preservation")), "Module preservation",
  "S1", basename(subfigure_file("S1", "Multifunctionality_Correlations_pH_TP_DOC_DON_AP")), "Multifunctionality relationships with pH, TP, DOC, DON and AP",
  "S2", basename(subfigure_file("S2", "Multifunctionality_Correlations_Bulk_Stoichiometry")), "Multifunctionality relationships with SOC:TN, TN:TP, SOC:TS, TP:TS and TN:TS",
  "S3", basename(subfigure_file("S3", "Multifunctionality_Correlations_Available_Stoichiometry")), "Multifunctionality relationships with DOC:AN, DOC:AP, DOC:AS, AN:AS, AP:AS and AN:AP",
  "S4", basename(subfigure_file("S4", "Multifunctionality_Correlations_Microbial_Properties")), "Multifunctionality relationships with CO2, MBC and qCO2",
  "S5", basename(subfigure_file("S5", "C_Site_Barplot")), "Site-level C-process responses",
  "S6", basename(subfigure_file("S6", "N_Site_Barplot")), "Site-level N-process responses",
  "S7", basename(subfigure_file("S7", "P_Site_Barplot")), "Site-level P-process responses",
  "S8", basename(subfigure_file("S8", "S_Site_Barplot")), "Site-level S-process responses",
  "S9", basename(subfigure_file("S9", "Module_Functional_Composition")), "Module functional composition",
  "S11", basename(subfigure_file("S11", "XC_CNPS_PCoA")), "XC CNPS PCoA",
  "S12", basename(subfigure_file("S12", "Multifunctionality_Correlation_Heatmaps")), "Multifunctionality correlation heatmaps",
  "S13", basename(subfigure_file("S13", "CNPS_Environment_Correlation_Heatmap")), "CNPS-environment correlation heatmap",
  "S14", basename(subfigure_file("S14", "Module_Alluvial")), "Module splitting and fusion"
)

readr::write_csv(
  figure_manifest,
  file.path(results_dir, "Figure_Manifest.csv")
)

expected_figure_paths <- file.path(figures_dir, figure_manifest$File)
missing_figure_outputs <- expected_figure_paths[!file.exists(expected_figure_paths)]

if (length(missing_figure_outputs)) {
  stop(
    "The pipeline finished with missing figure outputs:
",
    paste0("  - ", basename(missing_figure_outputs), collapse = "
")
  )
}

message(
  "Integrated analysis complete. ",
  length(expected_figure_paths),
  " PDF figure files were exported successfully."
)