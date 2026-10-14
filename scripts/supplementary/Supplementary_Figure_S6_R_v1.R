## ============================================================================
## Supplementary Figure S6
## Spatial robustness across neighborhood definitions
##
## PLOTTING ONLY:
##   - no adjacency reconstruction
##   - no Moran's I refitting
##   - no local HH recomputation
##   - no marker-exclusion panel (the frozen table is incomplete for that use)
##
## A. Frozen robust local HH maps for the two Figure 5 representative sections
## B. Frozen global Moran's I across 10 sections and 4 graph definitions
## C. Frozen robust HH extent across the same 10 x 4 combinations
## ============================================================================

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(ragg)
  library(svglite)
})

options(stringsAsFactors = FALSE)

## ------------------------------------------------------------------------
## 0. Paths
## ------------------------------------------------------------------------

FREEZE_ROOT <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
SPATIAL_DIR <- file.path(
  FREEZE_ROOT,
  "Figure2_source_data_addendum_20260912",
  "spatial"
)

path_moran <- file.path(SPATIAL_DIR, "31_spatial_global_Moran_999.csv")
path_local_hh <- file.path(SPATIAL_DIR, "32_spatial_local_high_high_999.csv.gz")
path_hh_burden <- file.path(SPATIAL_DIR, "33_spatial_hotspot_burden_999.csv")
path_coords <- file.path(SPATIAL_DIR, "spatial_P1_residuals.csv.gz")

figure_dir <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S6_R_v1", "figures")
source_dir <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S6_R_v1", "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(path_moran),
  file.exists(path_local_hh),
  file.exists(path_hh_burden),
  file.exists(path_coords)
)

## ------------------------------------------------------------------------
## 1. Frozen display definitions
## ------------------------------------------------------------------------

rep_sections <- tibble::tribble(
  ~dataset,     ~section, 
  "GSE217414", "19G0635",
  "GSE225857", "L1"
)

graph_order <- c("knn6", "radius_100um", "radius_200um", "radius_300um")
um_symbol <- intToUtf8(0x00B5)
en_dash <- intToUtf8(0x2013)
times_symbol <- intToUtf8(0x00D7)
middle_dot <- intToUtf8(0x00B7)
graph_labels <- c(
  knn6 = "kNN6",
  radius_100um = paste0("100 ", um_symbol, "m"),
  radius_200um = paste0("200 ", um_symbol, "m"),
  radius_300um = paste0("300 ", um_symbol, "m")
)

section_order <- tibble::tribble(
  ~dataset,     ~section,
  "GSE217414", "19G081",
  "GSE217414", "19G0619",
  "GSE217414", "19G0635",
  "GSE217414", "19G02977",
  "GSE225857", "C1",
  "GSE225857", "C2",
  "GSE225857", "C3",
  "GSE225857", "C4",
  "GSE225857", "L1",
  "GSE225857", "L2"
)

## ------------------------------------------------------------------------
## 2. Shared style and helpers
## ------------------------------------------------------------------------

COL_DARK <- "#22282E"
COL_GREY <- "#71808D"
COL_GRID <- "#E6EBEF"
COL_HH <- "#B6534E"
COL_OTHER <- "#CCD5DB"
COL_MORAN_LOW <- "#EAF2F6"
COL_MORAN_HIGH <- "#3E739F"
COL_BURDEN_LOW <- "#F7EFED"
COL_BURDEN_HIGH <- "#B6534E"

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes")
}

theme_s6 <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_DARK),
      axis.title = element_text(size = 6.1),
      axis.text = element_text(size = 5.3, colour = COL_DARK),
      axis.line = element_line(linewidth = 0.3, colour = COL_DARK),
      axis.ticks = element_line(linewidth = 0.3, colour = COL_DARK),
      plot.title = element_text(
        size = 7.5, face = "bold", hjust = 0,
        margin = margin(0, 0, 1.5, 0)
      ),
      plot.subtitle = element_text(
        size = 5.5, colour = COL_GREY, hjust = 0,
        margin = margin(0, 0, 2.0, 0)
      ),
      strip.background = element_blank(),
      strip.text = element_text(size = 5.6, face = "bold", colour = COL_DARK),
      panel.spacing = grid::unit(1.4, "mm"),
      legend.title = element_text(size = 5.5),
      legend.text = element_text(size = 5.0),
      plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm")
    )
}

write_source_csv <- function(df, filename) {
  utils::write.csv(
    df,
    file.path(source_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

read_frozen_csv <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    utils::read.csv(gzfile(path), check.names = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE)
  }
}

## ------------------------------------------------------------------------
## 3. Read and standardise the authoritative frozen tables
## ------------------------------------------------------------------------

coords_raw <- read_frozen_csv(path_coords)
local_raw <- read_frozen_csv(path_local_hh)
moran_raw <- read_frozen_csv(path_moran)
burden_raw <- read_frozen_csv(path_hh_burden)

required_coords <- c("dataset", "section", "barcode", "pxl_row", "pxl_col")
required_local <- c("dataset", "section", "barcode", "graph", "robust_HH_ge3")
required_moran <- c("dataset", "section", "graph", "Moran_I", "empirical_p", "BH_q", "permutations")
required_burden <- c("dataset", "section", "graph", "robust_HH_ge3", "n_spots", "robust_HH_ge3_fraction")

stopifnot(
  all(required_coords %in% names(coords_raw)),
  all(required_local %in% names(local_raw)),
  all(required_moran %in% names(moran_raw)),
  all(required_burden %in% names(burden_raw))
)

coords <- coords_raw %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    barcode = as.character(barcode),
    pxl_col = as.numeric(pxl_col),
    pxl_row = as.numeric(pxl_row)
  ) %>%
  distinct()

local_hh <- local_raw %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    barcode = as.character(barcode),
    graph = as.character(graph),
    robust_HH_ge3 = as_bool(robust_HH_ge3)
  )

moran <- moran_raw %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    graph = as.character(graph),
    Moran_I = as.numeric(Moran_I),
    empirical_p = as.numeric(empirical_p),
    BH_q = as.numeric(BH_q),
    permutations = as.numeric(permutations)
  )

burden <- burden_raw %>%
  transmute(
    dataset = as.character(dataset),
    section = as.character(section),
    graph = as.character(graph),
    robust_HH_ge3 = as_bool(robust_HH_ge3),
    n_spots = as.numeric(n_spots),
    robust_HH_ge3_fraction = as.numeric(robust_HH_ge3_fraction)
  )

## ------------------------------------------------------------------------
## 4. Frozen-table audit
## ------------------------------------------------------------------------

expected_key <- section_order %>%
  tidyr::crossing(graph = graph_order)

key_moran <- moran %>% distinct(dataset, section, graph)
key_burden <- burden %>% distinct(dataset, section, graph)
key_local <- local_hh %>% distinct(dataset, section, graph)

stopifnot(
  nrow(moran) == 40L,
  nrow(burden) == 40L,
  nrow(key_moran) == 40L,
  nrow(key_burden) == 40L,
  nrow(key_local) == 40L,
  setequal(key_moran, expected_key),
  setequal(key_burden, expected_key),
  setequal(key_local, expected_key),
  all(moran$Moran_I > 0, na.rm = TRUE),
  all(abs(moran$empirical_p - 0.001) < 1e-12),
  all(abs(moran$BH_q - 0.001) < 1e-12),
  all(moran$permutations == 999)
)

stopifnot(!anyDuplicated(local_hh[c("dataset", "section", "barcode", "graph")]))
stopifnot(!anyDuplicated(coords[c("dataset", "section", "barcode")]))

local_keys <- unique(paste(local_hh$dataset, local_hh$section, local_hh$barcode, sep = "|"))
coord_keys <- unique(paste(coords$dataset, coords$section, coords$barcode, sep = "|"))
stopifnot(all(local_keys %in% coord_keys))

## ------------------------------------------------------------------------
## 5. Panel A - representative local HH architecture
## ------------------------------------------------------------------------

A_plot <- local_hh %>%
  semi_join(rep_sections, by = c("dataset", "section")) %>%
  filter(graph %in% graph_order) %>%
  left_join(coords, by = c("dataset", "section", "barcode")) %>%
  mutate(
    HH_status = factor(
      if_else(robust_HH_ge3, "Robust HH", "Other spots"),
      levels = c("Other spots", "Robust HH")
    )
  )

stopifnot(
  nrow(A_plot) > 0L,
  nrow(A_plot %>% distinct(dataset, section, graph)) == 8L,
  all(is.finite(A_plot$pxl_col)),
  all(is.finite(A_plot$pxl_row))
)

write_source_csv(A_plot, "S6A_representative_local_HH_plot_data.csv")

make_A_map <- function(ds, sec, graph_name) {
  d <- A_plot %>%
    filter(dataset == ds, section == sec, graph == graph_name)
  stopifnot(nrow(d) > 0L)

  ggplot(d, aes(pxl_col, pxl_row)) +
    geom_point(
      data = d %>% filter(!robust_HH_ge3),
      colour = COL_OTHER,
      size = 0.21,
      alpha = 0.75,
      stroke = 0
    ) +
    geom_point(
      data = d %>% filter(robust_HH_ge3),
      colour = COL_HH,
      size = 0.32,
      alpha = 0.95,
      stroke = 0
    ) +
    coord_fixed() +
    scale_y_reverse() +
    labs(title = unname(graph_labels[graph_name])) +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(
        size = 5.8, face = "bold", hjust = 0.5, colour = COL_DARK,
        margin = margin(b = 0.8, unit = "mm")
      ),
      legend.position = "none",
      plot.margin = margin(0.5, 0.7, 0.5, 0.7, unit = "mm")
    )
}

make_A_card <- function(ds, sec, card_label) {
  card_header <- ggplot() +
    annotate(
      "text", x = 0.5, y = 0.5, label = card_label,
      family = "Arial", fontface = "bold", hjust = 0.5,
      size = 5.5 / ggplot2:::.pt, colour = COL_DARK
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0.5, 0, 0.5, 0, unit = "mm"))

  maps <- patchwork::wrap_plots(
    lapply(
      graph_order,
      function(graph_name) make_A_map(ds, sec, graph_name)
    ),
    ncol = 2
  )

  card_header / maps + patchwork::plot_layout(heights = c(0.12, 1))
}

A_card_217414 <- make_A_card(
  "GSE217414", "19G0635", paste0("GSE217414  ", middle_dot, "  19G0635")
)
A_card_225857 <- make_A_card(
  "GSE225857", "L1", paste0("GSE225857  ", middle_dot, "  L1")
)

pA_header <- ggplot() +
  annotate(
    "text", x = 0, y = 1.00,
    label = paste0("A  Representative local high", en_dash, "high architecture across graph definitions"),
    hjust = 0, vjust = 1, family = "Arial", fontface = "bold",
    size = 7.5 / ggplot2:::.pt, colour = COL_DARK
  ) +
  annotate(
    "text", x = 0, y = 0.28,
    label = "Representative sections fixed from Fig. 5; background spots are shown in grey and robust HH spots in red",
    hjust = 0, vjust = 1, family = "Arial",
    size = 5.5 / ggplot2:::.pt, colour = COL_GREY
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))

pA_legend <- ggplot(
  data.frame(HH_status = factor(c("Other spots", "Robust HH"), levels = c("Other spots", "Robust HH"))),
  aes(0, 0, colour = HH_status)
) +
  geom_point(size = 1.4, alpha = 0) +
  scale_colour_manual(
    values = c("Other spots" = COL_OTHER, "Robust HH" = COL_HH),
    name = NULL
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.0, colour = COL_DARK),
    legend.key.width = grid::unit(4, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "mm"),
    plot.margin = margin(0, 0, 0, 0, unit = "mm")
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = 1, size = 1.4)
    )
  )

A_cards <- (
  patchwork::wrap_elements(full = A_card_217414) |
    patchwork::wrap_elements(full = A_card_225857)
) /
  (
    patchwork::plot_spacer() |
      patchwork::wrap_elements(full = pA_legend)
  ) +
  patchwork::plot_layout(heights = c(1, 0.08))

pA <- pA_header /
  patchwork::wrap_elements(full = A_cards) +
  patchwork::plot_layout(heights = c(0.12, 1))

## ------------------------------------------------------------------------
## 6. Panels B and C - aligned 10 x 4 heatmaps
## ------------------------------------------------------------------------

section_levels <- rev(section_order$section)

B_plot <- moran %>%
  mutate(
    dataset = factor(dataset, levels = c("GSE217414", "GSE225857")),
    section_plot = factor(section, levels = section_levels),
    graph_plot = factor(graph, levels = graph_order, labels = unname(graph_labels[graph_order]))
  ) %>%
  arrange(dataset, section_plot, graph_plot)

C_plot <- burden %>%
  mutate(
    dataset = factor(dataset, levels = c("GSE217414", "GSE225857")),
    section_plot = factor(section, levels = section_levels),
    graph_plot = factor(graph, levels = graph_order, labels = unname(graph_labels[graph_order]))
  ) %>%
  arrange(dataset, section_plot, graph_plot)

write_source_csv(B_plot, "S6B_global_Moran_heatmap_data.csv")
write_source_csv(C_plot, "S6C_local_HH_extent_heatmap_data.csv")

pB <- ggplot(B_plot, aes(graph_plot, section_plot, fill = Moran_I)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  facet_grid(
    dataset ~ ., scales = "free_y", space = "free_y", drop = TRUE,
    switch = "y"
  ) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = COL_MORAN_LOW, high = COL_MORAN_HIGH, name = "Moran's I") +
  labs(
    title = "B  Global spatial organization across neighborhood definitions",
    subtitle = paste0("40/40 section ", times_symbol, " graph combinations showed positive Moran's I; empirical P = 0.001 throughout"),
    x = NULL,
    y = NULL
  ) +
  theme_s6() +
  theme(
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    legend.key.height = grid::unit(3.6, "mm"),
    legend.key.width = grid::unit(3.0, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      barheight = grid::unit(28, "mm"),
      barwidth = grid::unit(3.0, "mm")
    )
  )

pC <- ggplot(C_plot, aes(graph_plot, section_plot, fill = robust_HH_ge3_fraction)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  facet_grid(
    dataset ~ ., scales = "free_y", space = "free_y", drop = TRUE,
    switch = "y"
  ) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(
    low = COL_BURDEN_LOW,
    high = COL_BURDEN_HIGH,
    name = "Robust HH\nfraction",
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = paste0("C  Scale dependence of local high", en_dash, "high extent"),
    subtitle = paste0("Frozen robust HH fraction; same 10 sections ", times_symbol, " 4 graph definitions"),
    x = NULL,
    y = NULL
  ) +
  theme_s6() +
  theme(
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    legend.key.height = grid::unit(3.6, "mm"),
    legend.key.width = grid::unit(3.0, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      barheight = grid::unit(28, "mm"),
      barwidth = grid::unit(3.0, "mm")
    )
  )

figure_s6 <- patchwork::wrap_elements(full = pA) /
  patchwork::wrap_elements(full = (pB | pC)) +
  patchwork::plot_layout(heights = c(0.95, 1.05))

## ------------------------------------------------------------------------
## 7. Export and provenance
## ------------------------------------------------------------------------

figure_base <- file.path(figure_dir, "Supplementary_Figure_S6_R_v1")
figure_width_mm <- 183
figure_height_mm <- 175
figure_dpi <- 600

ragg::agg_png(
  paste0(figure_base, ".png"),
  width = figure_width_mm, height = figure_height_mm,
  units = "mm", res = figure_dpi, background = "white"
)
print(figure_s6)
dev.off()

svglite::svglite(
  paste0(figure_base, ".svg"),
  width = figure_width_mm / 25.4, height = figure_height_mm / 25.4,
  bg = "white"
)
print(figure_s6)
dev.off()

grDevices::cairo_pdf(
  paste0(figure_base, ".pdf"),
  width = figure_width_mm / 25.4, height = figure_height_mm / 25.4,
  family = "Arial"
)
print(figure_s6)
dev.off()

ragg::agg_tiff(
  paste0(figure_base, ".tiff"),
  width = figure_width_mm, height = figure_height_mm,
  units = "mm", res = figure_dpi,
  compression = "lzw", background = "white"
)
print(figure_s6)
dev.off()

legend_lines <- c(
  "Supplementary Figure S6 | Spatial robustness across neighborhood definitions.",
  paste0("A, Frozen robust high", en_dash, "high spot classifications for the two representative sections used in Fig. 5,"),
  paste0("shown across kNN6 and 100, 200 and 300 ", um_symbol, "m radius graphs. Coordinates were joined by dataset, section"),
  "and barcode; no adjacency or local statistic was recomputed.",
  "B, Frozen global Moran's I for all 10 sections across the four graph definitions. All 40 combinations",
  "were positive and used 999 permutations with empirical P = 0.001.",
  paste0("C, Frozen robust HH fraction for the same 10 section ", times_symbol, " graph combinations. No hotspot burden or"),
  "confidence interval was recomputed for this display.",
  "All panels are plotting-only displays of frozen source tables."
)
writeLines(legend_lines, file.path(figure_dir, "Supplementary_Figure_S6_legend.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(figure_dir, "sessionInfo.txt"))

message("Supplementary Figure S6 exported")
message("PNG: ", paste0(figure_base, ".png"))
message("PDF: ", paste0(figure_base, ".pdf"))
message("SVG: ", paste0(figure_base, ".svg"))
message("TIFF: ", paste0(figure_base, ".tiff"))
