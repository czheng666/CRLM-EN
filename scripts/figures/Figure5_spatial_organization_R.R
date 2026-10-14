## ========================================================================
## Figure 5 — Spatial organization of frozen P1 in CRLM
##
## R-only, plotting-only implementation.
## The script reads the authoritative spatial addendum tables bundled with the
## current Figure 2 source-data package. The addendum name follows the final
## main-figure numbering; the spatial content itself is unchanged.
## It does not refit models or recompute Moran/LISA statistics.
##
## Public plotting patterns consulted:
## - 10XGenomics/HumanColonCancer_VisiumHD:
##   direct spatial coordinates, restrained scales and patchwork-style assembly
## - alberto-valdeolivas/ST_CRC_CMS:
##   shared spatial colour scales and compact multi-panel layouts
## - cnio-bu/breast-bcspatial:
##   reproducible figure scripts and consistent spatial plot geometry
## ========================================================================

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

required <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "ggrepel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required package(s): ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggrepel)
})

options(stringsAsFactors = FALSE)

FREEZE_ROOT <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
ADD_DIR <- file.path(
  FREEZE_ROOT,
  "Figure2_source_data_addendum_20260912",
  "spatial"
)
OUT_DIR <- file.path(CRLM_OUTPUT_ROOT, "Figure5_spatial_organization")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source_file <- function(name) {
  candidates <- c(
    file.path(ADD_DIR, name),
    file.path(
      FREEZE_ROOT,
      "complete_run_archive_20260912",
      "run_output",
      "results",
      "tables",
      name
    )
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("Could not locate authoritative frozen spatial source: ", name)
  }
  hit[1]
}

RESIDUAL_PATH <- source_file("spatial_P1_residuals.csv.gz")
MORAN_PATH <- source_file("31_spatial_global_Moran_999.csv")
LOCAL_PATH <- source_file("32_spatial_local_high_high_999.csv.gz")
BURDEN_PATH <- source_file("33_spatial_hotspot_burden_999.csv")
EXCLUDED_PATH <- source_file("82_GSE225857_marker_excluded_spatial_sensitivity.csv")

residuals_raw <- read.csv(gzfile(RESIDUAL_PATH), check.names = FALSE)
residuals <- residuals_raw
moran <- read.csv(MORAN_PATH, check.names = FALSE)
local_hh <- read.csv(gzfile(LOCAL_PATH), check.names = FALSE)
burden <- read.csv(BURDEN_PATH, check.names = FALSE)
excluded <- read.csv(EXCLUDED_PATH, check.names = FALSE)

required_residual <- c("dataset", "section", "pxl_row", "pxl_col", "P1_residual")
required_moran <- c("dataset", "section", "graph", "Moran_I", "empirical_p", "BH_q", "permutations")
required_local <- c("dataset", "section", "graph", "HH_component", "HH_component_size", "robust_HH_ge3")
required_burden <- c("dataset", "section", "graph", "robust_HH_ge3", "n_spots", "robust_HH_ge3_fraction")
required_excluded <- c("dataset", "section", "score_version", "graph", "Moran_I", "BH_q")

stopifnot(
  all(required_residual %in% names(residuals)),
  all(required_moran %in% names(moran)),
  all(required_local %in% names(local_hh)),
  all(required_burden %in% names(burden)),
  all(required_excluded %in% names(excluded))
)

## ------------------------------------------------------------------------
## Shared style
## ------------------------------------------------------------------------

COL_DARK <- "#22282E"
COL_GREY <- "#71808D"
COL_GRID <- "#E6EBEF"
COL_LIGHT <- "#D8DEE5"
COL_DISC <- "#4C78A8"
COL_REP <- "#4E9A82"
COL_LOW <- "#4C78A8"
COL_MID <- "#F7F6F4"
COL_HIGH <- "#B85F55"
COL_HH <- "#B5483B"
COL_TISSUE <- "#E7ECEF"
COL_EXCL <- "#D17A67"

dataset_cols <- c(
  "GSE217414" = COL_DISC,
  "GSE225857" = COL_REP
)

theme_crlm <- function(base_size = 7) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      text = element_text(colour = COL_DARK),
      plot.title = element_text(
        size = 8.7, face = "bold", colour = COL_DARK,
        margin = margin(b = 0.8, unit = "mm")
      ),
      plot.subtitle = element_text(
        size = 5.8, colour = COL_GREY,
        margin = margin(b = 1.3, unit = "mm")
      ),
      axis.title = element_text(size = 6.3, colour = COL_DARK),
      axis.text = element_text(size = 5.8, colour = COL_DARK),
      legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 5.5),
      plot.margin = margin(2, 3, 2, 3, unit = "mm")
    )
}

section_label <- function(x) as.character(x)

## Primary frozen graph for the main spatial summaries.
moran_knn6 <- moran %>%
  filter(graph == "knn6") %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    observed_Moran_I = as.numeric(Moran_I),
    empirical_p = as.numeric(empirical_p),
    BH_q = as.numeric(BH_q),
    permutations = as.numeric(permutations)
  )

if (nrow(moran_knn6) != 10L) {
  stop("Expected 10 frozen kNN6 section-level Moran summaries; found ", nrow(moran_knn6))
}

## ========================================================================
## A — Representative residual P1 maps
## ========================================================================

representative <- moran_knn6 %>%
  group_by(dataset) %>%
  mutate(
    cohort_median = median(observed_Moran_I, na.rm = TRUE),
    distance_to_median = abs(observed_Moran_I - cohort_median)
  ) %>%
  arrange(distance_to_median, section, .by_group = TRUE) %>%
  slice_head(n = 1L) %>%
  ungroup()

## The frozen residual table stores spot-centre coordinates as pxl_col/pxl_row.
residuals <- residuals %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    barcode = as.character(barcode),
    x = as.numeric(pxl_col),
    y = as.numeric(pxl_row),
    pixel_size_um = as.numeric(pixel_size_um),
    P1_residual = as.numeric(P1_residual)
  )

## Common display limits are based on the frozen residual distribution.  This
## improves visual contrast while preserving one scale across both maps.
residual_limits <- as.numeric(
  quantile(residuals$P1_residual, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
)
if (!all(is.finite(residual_limits)) || residual_limits[1] >= 0 || residual_limits[2] <= 0) {
  stop("Frozen P1_residual quantile limits do not span zero.")
}
residual_breaks <- pretty(residual_limits, n = 5)

residual_fill_scale <- scale_fill_gradient2(
  low = COL_LOW,
  mid = "#FAFAF8",
  high = COL_HIGH,
  midpoint = 0,
  limits = residual_limits,
  oob = scales::squish,
  breaks = residual_breaks,
  name = "Frozen P1 residual",
  guide = guide_colourbar(
    direction = "horizontal",
    title.position = "bottom",
    title.hjust = 0.5,
    barwidth = unit(36, "mm"),
    barheight = unit(1.8, "mm"),
    ticks = FALSE
  )
)

local_knn6 <- local_hh %>%
  filter(graph == "knn6") %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    barcode = as.character(barcode),
    x = as.numeric(pxl_col),
    y = as.numeric(pxl_row),
    HH_component = as.character(HH_component),
    robust_HH_ge3 = robust_HH_ge3 %in% TRUE
  )

BASE_SPOT_SIZE <- 0.72

make_residual_map <- function(ds, sec, moran_value) {
  current <- residuals %>%
    filter(
      dataset == ds,
      section == sec,
      is.finite(P1_residual),
      is.finite(x),
      is.finite(y)
      )
  ggplot(current, aes(x = x, y = y, fill = P1_residual)) +
    geom_point(
      shape = 21,
      size = BASE_SPOT_SIZE,
      stroke = 0,
      colour = "transparent",
      alpha = 1
    ) +
    residual_fill_scale +
    coord_fixed(expand = FALSE) +
    scale_y_reverse(expand = expansion(mult = c(0.005, 0.005))) +
    labs(
      title = paste0(ds, " / ", section_label(sec)),
      subtitle = paste0("Moran's I = ", sprintf("%.3f", moran_value))
    ) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 7.0, face = "bold", hjust = 0.5, colour = COL_DARK),
      plot.subtitle = element_text(size = 5.3, hjust = 0.5, colour = COL_GREY),
      legend.position = "bottom",
      legend.title = element_text(size = 5.0, colour = COL_DARK),
      legend.text = element_text(size = 4.8, colour = COL_DARK),
      legend.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "mm"),
      legend.box.spacing = unit(0, "mm"),
      plot.margin = margin(1, 1, 1, 1, unit = "mm")
    )
}

rep217 <- representative %>% filter(dataset == "GSE217414")
rep225 <- representative %>% filter(dataset == "GSE225857")

A217 <- make_residual_map("GSE217414", rep217$section, rep217$observed_Moran_I)
A225 <- make_residual_map("GSE225857", rep225$section, rep225$observed_Moran_I)

pA_header <- ggplot() +
  annotate(
    "text", x = 0, y = 1.02,
    label = "A  Spatial organization of residual P1 activity",
    hjust = 0, vjust = 1,
    family = "Arial", size = 8.7 / ggplot2:::.pt,
    fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 0, y = 0.46,
    label = "Representative sections; residual values shown on a common display scale",
    hjust = 0, vjust = 1,
    family = "Arial", size = 5.7 / ggplot2:::.pt,
    colour = COL_GREY
  ) +
  xlim(0, 1) + ylim(0, 1.05) +
  theme_void(base_family = "Arial")

A_maps <- (A217 + A225) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

panel_A <- pA_header / A_maps +
  plot_layout(heights = c(0.10, 1))

## ========================================================================
## B — Global spatial autocorrelation
## ========================================================================

moran_plot <- moran_knn6 %>%
  mutate(
    section_short = section_label(section),
    dataset = factor(dataset, levels = c("GSE217414", "GSE225857"))
  ) %>%
  arrange(dataset, desc(observed_Moran_I)) %>%
  mutate(section_short = factor(section_short, levels = rev(section_short)))

cohort_medians <- moran_plot %>%
  group_by(dataset) %>%
  summarise(med = median(observed_Moran_I, na.rm = TRUE), .groups = "drop")

pB <- ggplot(
  moran_plot,
  aes(x = observed_Moran_I, y = section_short, colour = dataset)
) +
  geom_segment(aes(x = 0, xend = observed_Moran_I, yend = section_short), linewidth = 0.65, alpha = 0.75) +
  geom_point(size = 2.35, stroke = 0.25) +
  geom_vline(
    data = cohort_medians,
    aes(xintercept = med, colour = dataset),
    linetype = "dashed", linewidth = 0.55, alpha = 0.55
  ) +
  scale_colour_manual(values = dataset_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 0.56), breaks = seq(0, 0.5, 0.1)) +
  labs(
    title = "B  Reproducible global spatial autocorrelation",
    subtitle = "kNN6; 10/10 sections positive; 999 permutations; all q = 0.001",
    x = "Residual Moran's I", y = NULL
  ) +
  theme_crlm() +
  theme(
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.35)
  ) +
  annotate(
    "text", x = 0.54, y = 9.55,
    label = "GSE217414", hjust = 1, size = 2.2,
    colour = scales::alpha(COL_DISC, 0.72), family = "Arial"
  ) +
  annotate(
    "text", x = 0.54, y = 4.55,
    label = "GSE225857", hjust = 1, size = 2.2,
    colour = scales::alpha(COL_REP, 0.72), family = "Arial"
  )

## ========================================================================
## C — Local High–High hotspot architecture
## ========================================================================

burden_knn6 <- burden %>%
  filter(graph == "knn6") %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    n_spots = as.numeric(n_spots),
    robust_HH_spots = as.numeric(robust_HH_ge3),
    robust_HH_fraction = as.numeric(robust_HH_ge3_fraction)
  )

component_summary <- local_hh %>%
  filter(graph == "knn6", robust_HH_ge3 %in% TRUE, is.finite(HH_component_size)) %>%
  group_by(dataset, section) %>%
  summarise(
    n_domains = dplyr::n_distinct(HH_component[!is.na(HH_component)]),
    largest_domain_spots = max(HH_component_size, na.rm = TRUE),
    local_robust_HH_spots = dplyr::n(),
    .groups = "drop"
  )

architecture <- burden_knn6 %>%
  left_join(component_summary, by = c("dataset", "section")) %>%
  mutate(
    n_domains = if_else(is.na(n_domains), 0, n_domains),
    largest_domain_spots = if_else(is.na(largest_domain_spots), 0, largest_domain_spots),
    extent_pct = 100 * robust_HH_fraction,
    dominance_pct = if_else(
      robust_HH_spots > 0,
      100 * largest_domain_spots / robust_HH_spots,
      0
    ),
    section_short = section_label(section)
  )

c_rep <- architecture %>%
  filter(dataset == "GSE225857") %>%
  arrange(desc(robust_HH_fraction), section) %>%
  slice_head(n = 1L)

c_current <- residuals %>%
  filter(
    dataset == c_rep$dataset[1],
    section == c_rep$section[1],
    is.finite(P1_residual),
    is.finite(x),
    is.finite(y)
  )
c_hh <- local_knn6 %>%
  filter(dataset == c_rep$dataset[1], section == c_rep$section[1], robust_HH_ge3) %>%
  distinct(barcode, .keep_all = TRUE)
c_map_data <- c_current %>%
  left_join(c_hh %>% select(barcode, HH_component) %>% mutate(is_HH = TRUE), by = "barcode") %>%
  mutate(is_HH = if_else(is.na(is_HH), FALSE, is_HH))

c_map <- ggplot(c_map_data, aes(x = x, y = y)) +
  geom_point(
    data = c_map_data %>% filter(!is_HH),
    shape = 21, size = 0.72, fill = COL_TISSUE, colour = "transparent", alpha = 1, stroke = 0
  ) +
  geom_point(
    data = c_map_data %>% filter(is_HH),
    shape = 21, size = 0.76, fill = COL_HH, colour = "transparent", alpha = 0.78, stroke = 0
  ) +
  coord_fixed(expand = FALSE) +
  scale_y_reverse(expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title = paste0(c_rep$dataset[1], " / ", c_rep$section[1]),
    subtitle = "Highest robust High-High burden"
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 6.8, face = "bold", hjust = 0.5, colour = COL_DARK),
    plot.subtitle = element_text(size = 5.2, hjust = 0.5, colour = COL_GREY),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

section_order_C <- c("19G081", "19G0619", "19G0635", "19G02977", "C2", "C1", "C4", "L1", "C3", "L2")

## Same frozen architecture metrics, shown as one burden-versus-concentration
## scatter so that the local domain structure is not presented as a QC table.
architecture_C <- architecture %>%
  filter(section_short %in% section_order_C) %>%
  mutate(
    section_short = factor(section_short, levels = section_order_C)
  ) %>%
  arrange(section_short)

extent_median <- median(architecture_C$extent_pct, na.rm = TRUE)
dominance_median <- median(architecture_C$dominance_pct, na.rm = TRUE)

pC_summary <- ggplot(
  architecture_C,
  aes(x = extent_pct, y = dominance_pct, colour = dataset)
) +
  geom_vline(
    xintercept = extent_median,
    linetype = "dashed",
    linewidth = 0.3,
    colour = COL_LIGHT
  ) +
  geom_hline(
    yintercept = dominance_median,
    linetype = "dashed",
    linewidth = 0.3,
    colour = COL_LIGHT
  ) +
  geom_point(size = 2.15, alpha = 0.95) +
  geom_text_repel(
    aes(label = section_short),
    hjust = 0,
    direction = "both",
    box.padding = grid::unit(0.14, "lines"),
    point.padding = grid::unit(0.12, "lines"),
    force = 0.35,
    max.overlaps = Inf,
    min.segment.length = Inf,
    seed = 20260913,
    size = 1.75,
    colour = COL_DARK,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = dataset_cols, name = NULL) +
  scale_x_continuous(
    limits = c(0, 12),
    breaks = seq(0, 10, 2),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Robust HH burden (%)",
    y = "Largest HH domain share (%)"
  ) +
  theme_crlm(base_size = 6.2) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 4.8, colour = COL_DARK),
    legend.key.size = unit(3.0, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    legend.box.spacing = unit(0, "mm"),
    axis.title = element_text(size = 5.8, colour = COL_DARK),
    axis.text = element_text(size = 5.2, colour = COL_DARK),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_line(colour = COL_GRID, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 2, 1, 1, unit = "mm")
  )

pC_header <- ggplot() +
  annotate(
    "text", x = 0, y = 1.05,
    label = "C  Local High-High domain architecture",
    hjust = 0, vjust = 1,
    family = "Arial", size = 8.7 / ggplot2:::.pt,
    fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 0, y = 0.52,
    label = "Representative domain map and section metrics",
    hjust = 0, vjust = 1,
    family = "Arial", size = 5.8 / ggplot2:::.pt,
    colour = COL_GREY
  ) +
  xlim(0, 1) + ylim(0, 1.1) +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(0, 1, 0, 1, unit = "mm"))

panel_C <- (pC_header / (c_map | pC_summary)) +
  plot_layout(heights = c(0.14, 1), widths = c(0.90, 1.35))

## ========================================================================
## D — Marker-exclusion robustness
## ========================================================================

marker_excl <- excluded %>%
  filter(graph == "knn6", score_version == "P1_marker_excluded") %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    excluded_Moran_I = as.numeric(Moran_I),
    excluded_BH_q = as.numeric(BH_q)
  )

marker <- moran_knn6 %>%
  filter(dataset == "GSE225857") %>%
  left_join(marker_excl, by = c("dataset", "section")) %>%
  mutate(
    retention = excluded_Moran_I / observed_Moran_I,
    section_short = factor(section_label(section), levels = section_label(section[order(observed_Moran_I)]))
  ) %>%
  arrange(desc(observed_Moran_I))

if (nrow(marker) != 6L || any(!is.finite(marker$retention))) {
  stop("GSE225857 marker-exclusion source does not contain 6 complete kNN6 sections.")
}

median_retention <- median(marker$retention, na.rm = TRUE)
marker_long <- marker %>%
  pivot_longer(
    cols = c(observed_Moran_I, excluded_Moran_I),
    names_to = "score", values_to = "moran"
  ) %>%
  mutate(score = recode(score, observed_Moran_I = "Full P1", excluded_Moran_I = "SPP1/APOC1/APOE removed"))

pD <- ggplot() +
  geom_segment(
    data = marker,
    aes(x = observed_Moran_I, xend = excluded_Moran_I, y = section_short, yend = section_short),
    colour = "#B0B7BF", linewidth = 0.7
  ) +
  geom_point(
    data = marker_long,
    aes(x = moran, y = section_short, shape = score, colour = score),
    size = 2.35
  ) +
  geom_text(
    data = marker,
    aes(x = 0.55, y = section_short, label = sprintf("%.0f%%", 100 * retention)),
    hjust = 0, size = 2.25, colour = COL_DARK
  ) +
  scale_shape_manual(
    values = c("Full P1" = 16, "SPP1/APOC1/APOE removed" = 15), name = NULL
  ) +
  scale_colour_manual(
    values = c("Full P1" = COL_DISC, "SPP1/APOC1/APOE removed" = COL_EXCL), name = NULL
  ) +
  scale_x_continuous(limits = c(0.05, 0.64), breaks = seq(0.1, 0.6, 0.1)) +
  scale_y_discrete(expand = expansion(add = c(0.25, 0.72))) +
  labs(
    title = "D  Robustness to dominant-marker exclusion",
    subtitle = paste0(
      "SPP1/APOC1/APOE removed; 6/6 sections; median Moran retention = ",
      round(100 * median_retention), "%"
    ),
    x = "Residual Moran's I", y = NULL
  ) +
  annotate(
    "text", x = 0.55, y = length(levels(marker$section_short)) + 0.42,
    label = "Retention", hjust = 0, size = 2.35, fontface = "bold", colour = COL_DARK
  ) +
  theme_crlm() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 5.4),
    legend.key.size = unit(3.5, "mm"),
    legend.spacing.x = unit(1.5, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = COL_GRID, linewidth = 0.35)
  )

## ========================================================================
## Final assembly and exports
## ========================================================================

top_row <-
  (patchwork::wrap_elements(full = panel_A) | pB) +
  plot_layout(widths = c(1.35, 0.85))

bottom_row <-
  (patchwork::wrap_elements(full = panel_C) | pD) +
  plot_layout(widths = c(1.15, 0.85))

Figure5 <- top_row / bottom_row +
  plot_layout(heights = c(1.18, 0.82))

png_path <- file.path(OUT_DIR, "Figure5_spatial_organization.png")
pdf_path <- file.path(OUT_DIR, "Figure5_spatial_organization.pdf")
svg_path <- file.path(OUT_DIR, "Figure5_spatial_organization.svg")
tiff_path <- file.path(OUT_DIR, "Figure5_spatial_organization.tiff")

ggsave(png_path, Figure5, width = 180, height = 170, units = "mm", dpi = 600, bg = "white")
ggsave(pdf_path, Figure5, width = 180, height = 170, units = "mm", device = grDevices::cairo_pdf, bg = "white")
ggsave(
  tiff_path, Figure5, width = 180, height = 170, units = "mm", dpi = 600,
  device = ragg::agg_tiff, compression = "lzw", bg = "white"
)
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(svg_path, width = 180 / 25.4, height = 170 / 25.4)
  print(Figure5)
  dev.off()
}

write.csv(representative, file.path(OUT_DIR, "Figure5A_representative_sections_USED.csv"), row.names = FALSE)
write.csv(c_rep, file.path(OUT_DIR, "Figure5C_representative_HH_section_USED.csv"), row.names = FALSE)
write.csv(moran_plot, file.path(OUT_DIR, "Figure5B_Moran_kNN6_USED.csv"), row.names = FALSE)
write.csv(architecture, file.path(OUT_DIR, "Figure5C_hotspot_architecture_kNN6_USED.csv"), row.names = FALSE)
write.csv(architecture_C, file.path(OUT_DIR, "Figure5C_summary_all_sections_USED.csv"), row.names = FALSE)
write.csv(marker, file.path(OUT_DIR, "Figure5D_marker_exclusion_kNN6_USED.csv"), row.names = FALSE)

legend_lines <- c(
  "Figure 5. Spatial organization and robustness of the frozen P1 myeloid program.", "",
  "A, Representative residual P1 maps. Whole-section maps are shown for one section per",
  "cohort selected closest to the cohort median frozen kNN6 Moran's I. Spot fill encodes",
  "the frozen P1_residual on a common display scale.",
  "Maps use frozen P1_residual values and stored spot-centre coordinates; no spatial statistic",
  "was recomputed for these displays.", "",
  "B, Reproducible global spatial autocorrelation across all 10 sections using the frozen",
  "kNN6 graph. Points show observed Moran's I; dashed vertical lines indicate cohort-specific",
  "median Moran's I values; cohort labels are shown within the panel; the subtitle reports the",
  "frozen permutation/q-value result.", "",
  "C, Local High-High domain architecture. The map shows the GSE225857 section with the highest",
  "frozen robust_HH_ge3 burden; red spots are frozen robust High-High spots and pale spots are",
  "other stored spots. The adjacent scatter shows frozen robust HH burden versus largest-domain",
  "share for all 10 sections, with colour indicating cohort. Dashed lines mark the medians of",
  "the corresponding metrics across all 10 sections; largest-domain share is the largest frozen",
  "robust High-High component divided by the section's robust HH spots; no local statistic was",
  "recomputed.", "",
  "D, Robustness to dominant-marker exclusion in GSE225857. Full P1 Moran's I is paired with the frozen",
  "P1_marker_excluded result after removing SPP1, APOC1 and APOE from the score. Retention",
  "is the excluded/full Moran's I ratio × 100%; values above 100% indicate slightly greater",
  "spatial autocorrelation after exclusion; no confidence interval is implied by this display.", "",
  "All panels are plotting-only displays of frozen spatial source tables. Primary spatial",
  "graph: kNN6. Residual P1 values are shown on a common display scale clipped at the pooled",
  "2.5th and 97.5th percentiles; no score",
  "standardization, spatial statistic refitting or permutation recomputation was performed."
)
writeLines(legend_lines, file.path(OUT_DIR, "Figure5_spatial_organization_legend.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))

message("Figure 5 exported")
message("PNG: ", png_path)
message("PDF: ", pdf_path)
message("SVG: ", svg_path)
message("Representative sections: ", paste(representative$dataset, representative$section, collapse = "; "))
message("Moran sections: ", nrow(moran_knn6), "; positive: ", sum(moran_knn6$observed_Moran_I > 0))
message("Marker-exclusion sections retained: ", nrow(marker))
