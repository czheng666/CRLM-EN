## ========================================================================
## Supplementary Figure S1
## Discovery preprocessing, atlas and final NMF input
##
## R-only plotting, export and visual-QA workflow.
## This script performs no refitting and introduces no new estimand. It reads
## the current formal rerun outputs and exports figure-source tables so every
## plotted quantity remains traceable.
## ========================================================================

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

required <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "svglite", "ragg", "Seurat", "ggrepel")
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
})

options(stringsAsFactors = FALSE)
set.seed(20260914)

## ------------------------------------------------------------------------
## Frozen inputs and output bundle
## ------------------------------------------------------------------------

FREEZE_ROOT <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
BASE_DIR <- file.path(
  FREEZE_ROOT,
  "complete_run_archive_20260912",
  "run_output",
  "results"
)
TABLE_DIR <- file.path(BASE_DIR, "tables")
AUDIT_DIR <- file.path(BASE_DIR, "audit")

## The formal rerun Seurat object is retained in the run cache. It provides
## the global UMAP and its cell-level frozen broad-lineage annotation; no
## embedding or annotation is recomputed here.
GLOBAL_RDS <- file.path(CRLM_PRIMARY_RERUN_ROOT, "cache/annotation/06_GSE178318_global_broad_lineage.rds")

## The five-state field is read directly from the frozen Figure 1 source
## table. In 36b_GSE178318_myeloid_UMAP_P1.csv.gz, the similarly named
## `fixed_state` column is a retention/QC status (retained versus
## flagged_or_excluded), not the five-state biological annotation.
FIGURE1_UMAP_SOURCE <- file.path(CRLM_OUTPUT_ROOT, "Figure1_P1_actual_v18", "source_data", "Figure1D_GSE178318_UMAP_P1.csv.gz")

OUT_DIR <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S1_R_v1")
FIG_DIR <- file.path(OUT_DIR, "figures")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(name, dir = TABLE_DIR) {
  path <- file.path(dir, name)
  if (!file.exists(path)) stop("Missing authoritative source: ", path)
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    read.csv(gzfile(path), check.names = FALSE)
  } else {
    read.csv(path, check.names = FALSE)
  }
}

write_csv_gz <- function(x, name) {
  path <- file.path(SOURCE_DIR, name)
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, na = "")
  path
}

write_csv_plain <- function(x, name) {
  path <- file.path(SOURCE_DIR, name)
  write.csv(x, path, row.names = FALSE, na = "")
  path
}

## ------------------------------------------------------------------------
## Figure contract
## ------------------------------------------------------------------------

## Core conclusion: the frozen discovery input is transparent, QC-audited,
## not dominated by a single library, and reduces to 12,433 myeloid cells for
## de novo P1 discovery.
## Archetype: asymmetric mixed-modality figure (quantitative audit + atlas).
## Export: 183 mm wide, editable SVG/PDF, 600-dpi TIFF/PNG, white background.

COL_DARK <- "#22282E"
COL_GREY <- "#71808D"
COL_GRID <- "#E6EBEF"
COL_BLUE <- "#4C78A8"
COL_TEAL <- "#4E9A82"
COL_RED <- "#B85F55"
COL_ORANGE <- "#D17A67"
COL_PURPLE <- "#756BB1"
COL_LIGHT <- "#D8DEE5"

lineage_cols <- c(
  T_NK = "#4E79A7",
  B = "#F28E2B",
  Plasma = "#E3B23C",
  Epithelial = "#E15759",
  Fibroblast = "#59A14F",
  Endothelial = "#76B7B2",
  Mast = "#B07AA1",
  Myeloid_DC = "#9C755F"
)

state_levels <- c("Mono", "Inflammatory", "C1Q/MHC-II", "P1-state", "Other")
state_cols <- c(
  Mono = "#4F94BE",
  Inflammatory = "#D58A2D",
  `C1Q/MHC-II` = "#438D75",
  `P1-state` = "#A63A45",
  Other = "#85858B"
)

## Slightly deeper display palette for Panel D only; Panel E retains state_cols.
state_cols_D <- c(
  Mono = "#3785B2",
  Inflammatory = "#C4771C",
  `C1Q/MHC-II` = "#357A64",
  `P1-state` = "#92303B",
  Other = "#6F7076"
)

tissue_cols <- c(CRC = "#4C78A8", LM = "#B65C5C", PBMC = "#9A9A9A")

theme_crlm <- function(base_size = 7) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      text = element_text(colour = COL_DARK),
      plot.title = element_text(size = 8.3, face = "bold", colour = COL_DARK,
                                margin = margin(b = 1.0, unit = "mm")),
      plot.subtitle = element_text(size = 5.8, colour = COL_GREY,
                                   margin = margin(b = 1.2, unit = "mm")),
      axis.title = element_text(size = 6.2, colour = COL_DARK),
      axis.text = element_text(size = 5.5, colour = COL_DARK),
      legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 5.5),
      panel.grid.major.y = element_line(colour = COL_GRID, linewidth = 0.25),
      panel.grid.minor = element_blank(),
      plot.margin = margin(2, 3, 2, 3, unit = "mm")
    )
}

theme_void_crlm <- function() {
  theme_void(base_family = "Arial", base_size = 7) +
    theme(
      plot.title = element_text(size = 8.0, face = "bold", hjust = 0,
                                colour = COL_DARK),
      plot.subtitle = element_text(size = 5.5, hjust = 0, colour = COL_GREY),
      legend.text = element_text(size = 5.5, colour = COL_DARK),
      legend.key.height = grid::unit(2.5, "mm"),
      plot.margin = margin(2, 3, 2, 3, unit = "mm")
    )
}

save_pub_r <- function(plot, stem, width_mm = 180, height_mm = 205, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(file.path(FIG_DIR, paste0(stem, ".svg")), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(file.path(FIG_DIR, paste0(stem, ".pdf")), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(FIG_DIR, paste0(stem, ".tiff")), width = w, height = h,
                 units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(FIG_DIR, paste0(stem, ".png")), width = w, height = h,
                units = "in", res = dpi)
  print(plot)
  dev.off()
}

## ------------------------------------------------------------------------
## Load current frozen tables
## ------------------------------------------------------------------------

qc <- read_csv("02_qc_retention_by_sample.csv")
scd <- read_csv("05_scdblfinder_by_sample.csv")
nmf_umap <- read_csv("36b_GSE178318_myeloid_UMAP_P1.csv.gz")
if (!file.exists(FIGURE1_UMAP_SOURCE)) stop("Missing frozen Figure 1 state source: ", FIGURE1_UMAP_SOURCE)
figure1_umap <- read.csv(gzfile(FIGURE1_UMAP_SOURCE), check.names = FALSE)
lineage_summary <- read_csv("17_lineage_subcluster_summary.csv")
singlet_audit <- read_csv("02_singlet_count_audit.csv", AUDIT_DIR)
atlas_audit <- read_csv("04_global_analysis_cell_count_audit.csv", AUDIT_DIR)

stopifnot(
  nrow(qc) == 15L,
  nrow(scd) == 15L,
  nrow(nmf_umap) == 12433L,
  nrow(figure1_umap) == 12433L,
  all(c("barcode", "state") %in% names(figure1_umap)),
  singlet_audit$n_singlets[1] == 108353L,
  singlet_audit$n_doublets[1] == 11047L,
  atlas_audit$n_analysis_cells[1] == 105267L,
  lineage_summary$n_cells[lineage_summary$lineage == "Myeloid_DC"] == 12433L
)

## ------------------------------------------------------------------------
## Extract the current global atlas coordinates from the formal rerun object
## ------------------------------------------------------------------------

if (!file.exists(GLOBAL_RDS)) stop("Missing formal rerun global atlas object: ", GLOBAL_RDS)
global_obj <- readRDS(GLOBAL_RDS)
global_meta <- global_obj@meta.data
global_emb <- as.data.frame(Seurat::Embeddings(global_obj, reduction = "umap.patient"))
if (!all(c("global_broad_lineage", "global_lineage_status", "patient_id", "tissue") %in% colnames(global_meta))) {
  stop("Formal rerun object does not contain the required global atlas metadata columns.")
}
global_meta <- global_meta[rownames(global_emb), , drop = FALSE]
global_umap <- data.frame(
  cell_id = rownames(global_emb),
  UMAP_1 = as.numeric(global_emb[[1]]),
  UMAP_2 = as.numeric(global_emb[[2]]),
  broad_lineage = as.character(global_meta$global_broad_lineage),
  annotation_status = as.character(global_meta$global_lineage_status),
  patient = as.character(global_meta$patient_id),
  tissue = as.character(global_meta$tissue),
  stringsAsFactors = FALSE
)
global_umap <- global_umap %>%
  filter(annotation_status == "retained", is.finite(UMAP_1), is.finite(UMAP_2))
rm(global_obj, global_meta, global_emb)
gc()

if (nrow(global_umap) != 105267L) {
  stop("Expected 105,267 retained final-atlas cells for the global UMAP; found ", nrow(global_umap))
}

## ------------------------------------------------------------------------
## Build R-generated S1 source tables
## ------------------------------------------------------------------------

flow <- data.frame(
  step = c("Raw GEO barcodes", "Primary QC", "Singlets", "Final atlas", "Myeloid NMF input"),
  n = c(
    sum(as.numeric(qc$n_raw)),
    sum(as.numeric(qc$n_primary_qc)),
    as.numeric(singlet_audit$n_singlets[1]),
    as.numeric(atlas_audit$n_analysis_cells[1]),
    12433L
  ),
  stringsAsFactors = FALSE
)

library_qc <- qc %>%
  select(sample_id, patient_id, tissue, n_raw, n_primary_qc, retention) %>%
  left_join(
    scd %>% select(sample_id, n_qc, n_singlet, n_doublet, singlet_rate, expected_doublet_rate),
    by = "sample_id"
  ) %>%
  mutate(
    qc_retention = n_primary_qc / n_raw,
    doublet_fraction = n_doublet / n_qc
  )

state_lookup <- figure1_umap %>%
  transmute(
    cell_id = as.character(barcode),
    frozen_broad_state = as.character(state)
  )

stopifnot(
  nrow(state_lookup) == 12433L,
  !anyDuplicated(state_lookup$cell_id),
  setequal(state_lookup$cell_id, as.character(nmf_umap$barcode)),
  all(state_lookup$frozen_broad_state %in% state_levels)
)

nmf_source <- nmf_umap %>%
  transmute(
    cell_id = barcode,
    patient,
    sample,
    tissue,
    fine_state,
    retention_status = fixed_state,
    P1_usage,
    P1_top15_score,
    UMAP_1,
    UMAP_2
  ) %>%
  left_join(state_lookup, by = "cell_id") %>%
  select(
    cell_id, patient, sample, tissue, fine_state, frozen_broad_state,
    retention_status, P1_usage, P1_top15_score, UMAP_1, UMAP_2
  )

stopifnot(
  nrow(nmf_source) == 12433L,
  !anyNA(nmf_source$frozen_broad_state)
)

state_audit <- data.frame(
  source_table = FIGURE1_UMAP_SOURCE,
  n_cells = nrow(nmf_source),
  state_field = "state",
  all_barcodes_matched = TRUE,
  all_states_in_frozen_levels = all(nmf_source$frozen_broad_state %in% state_levels),
  regex_remapping_used = FALSE,
  stringsAsFactors = FALSE
)

composition <- nmf_source %>%
  count(patient, tissue, frozen_broad_state, name = "n_cells") %>%
  group_by(patient, tissue) %>%
  mutate(fraction = n_cells / sum(n_cells)) %>%
  ungroup()

write_csv_plain(flow, "S1A_preprocessing_flow.csv")
write_csv_plain(library_qc, "S1B_library_qc_audit.csv")
write_csv_gz(global_umap, "S1C_global_atlas_coordinates.csv.gz")
write_csv_gz(nmf_source, "S1D_myeloid_atlas_coordinates.csv.gz")
write_csv_plain(composition, "S1E_patient_tissue_state_composition.csv")
write_csv_plain(state_audit, "S1D_state_definition_audit.csv")

## ------------------------------------------------------------------------
## Panel A — compact discovery flow
## ------------------------------------------------------------------------

flow_plot <- flow %>%
  mutate(x = seq_len(n()))

pA <- ggplot(flow_plot, aes(x = x, y = 1)) +
  geom_segment(
    data = flow_plot[-nrow(flow_plot), , drop = FALSE],
    aes(x = x + 0.22, xend = x + 0.78, y = 1, yend = 1),
    linewidth = 0.45, colour = "#AAB4BC",
    arrow = grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")
  ) +
  geom_label(
    aes(label = paste0(step, "\n", comma(n))),
    size = 2.65, family = "Arial", linewidth = 0.25,
    label.padding = grid::unit(1.5, "mm"), fill = "white", colour = COL_DARK
  ) +
  annotate("text", x = 0.55, y = 1.43,
           label = "A  Discovery cohort and NMF-input definition",
           hjust = 0, family = "Arial", fontface = "bold", size = 3.0, colour = COL_DARK) +
  coord_cartesian(xlim = c(0.5, 5.5), ylim = c(0.70, 1.50), clip = "off") +
  theme_void()

## ------------------------------------------------------------------------
## Panel B — 15-library QC retention and doublet audit
## ------------------------------------------------------------------------

library_plot <- library_qc %>%
  mutate(library = factor(sample_id, levels = rev(sample_id)))

pB1 <- ggplot(library_plot, aes(qc_retention, library, colour = tissue)) +
  geom_segment(aes(x = 0, xend = qc_retention, yend = library),
               colour = "#E2E7EA", linewidth = 0.45) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = tissue_cols, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(title = "B  Library-level QC audit", x = "QC retention", y = NULL) +
  theme_crlm() +
  theme(legend.position = "none", axis.text.y = element_text(size = 4.7))

pB2 <- ggplot(library_plot, aes(doublet_fraction, library, colour = tissue)) +
  geom_segment(aes(x = 0, xend = doublet_fraction, yend = library),
               colour = "#E2E7EA", linewidth = 0.45) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = tissue_cols, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, max(library_plot$doublet_fraction, na.rm = TRUE) * 1.12)) +
  labs(title = "Doublet calls", x = "Doublet fraction", y = NULL, colour = NULL) +
  theme_crlm() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")

panel_B <- pB1 | pB2

## ------------------------------------------------------------------------
## Panel C — global atlas UMAP
## ------------------------------------------------------------------------

set.seed(20260914)
global_umap_plot <- global_umap[sample.int(nrow(global_umap)), , drop = FALSE]
nmf_umap_plot <- nmf_source[sample.int(nrow(nmf_source)), , drop = FALSE]

## Keep only lineages actually displayed in retained atlas.
global_umap_plot$broad_lineage <- factor(global_umap_plot$broad_lineage, levels = names(lineage_cols))
present_lineages <- intersect(
  names(lineage_cols),
  unique(as.character(global_umap_plot$broad_lineage))
)
global_umap_plot$broad_lineage <- factor(
  global_umap_plot$broad_lineage,
  levels = present_lineages
)

## Robust label position: median coordinates are less affected by odd-shaped
## cluster tails than mean coordinates.
label_df <- global_umap_plot |>
  dplyr::filter(!is.na(broad_lineage)) |>
  dplyr::group_by(broad_lineage) |>
  dplyr::summarise(
    UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
    UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    label = dplyr::recode(
      as.character(broad_lineage),
      T_NK = "T/NK",
      B = "B",
      Plasma = "Plasma",
      Epithelial = "Epithelial",
      Fibroblast = "Fibroblast",
      Endothelial = "Endothelial",
      Mast = "Mast",
      Myeloid_DC = "Myeloid/DC"
    )
  )

## All retained lineage labels are placed outside the point clouds. Thin
## leader lines connect each label to the corresponding median island.
label_layout <- label_df |>
  dplyr::mutate(
    label_x = dplyr::case_when(
      broad_lineage == "B" ~ UMAP_1 + 1.3,
      broad_lineage == "Mast" ~ UMAP_1 - 2.1,
      broad_lineage == "Myeloid_DC" ~ UMAP_1 - 2.4,
      broad_lineage == "T_NK" ~ UMAP_1 + 3.2,
      broad_lineage == "Plasma" ~ UMAP_1 - 1.5,
      broad_lineage == "Endothelial" ~ UMAP_1 - 3.0,
      broad_lineage == "Fibroblast" ~ UMAP_1 - 0.2,
      broad_lineage == "Epithelial" ~ UMAP_1 + 3.0,
      TRUE ~ UMAP_1
    ),
    label_y = dplyr::case_when(
      broad_lineage == "B" ~ UMAP_2 + 0.6,
      broad_lineage == "Mast" ~ UMAP_2 + 1.0,
      broad_lineage == "Myeloid_DC" ~ UMAP_2 + 3.3,
      broad_lineage == "T_NK" ~ UMAP_2 + 4.5,
      broad_lineage == "Plasma" ~ UMAP_2 + 1.3,
      broad_lineage == "Endothelial" ~ UMAP_2 + 0.8,
      broad_lineage == "Fibroblast" ~ UMAP_2 - 1.2,
      broad_lineage == "Epithelial" ~ UMAP_2 + 2.0,
      TRUE ~ UMAP_2
    ),
    label_hjust = dplyr::case_when(
      broad_lineage %in% c("Fibroblast") ~ 0.5,
      broad_lineage %in% c("T_NK", "Epithelial") ~ 1.0,
      TRUE ~ 0.0
    )
  )

pC <- ggplot(global_umap_plot, aes(UMAP_1, UMAP_2, colour = broad_lineage)) +
  geom_point(size = 0.070, alpha = 0.95, stroke = 0) +
  geom_segment(
    data = label_layout,
    aes(x = UMAP_1, y = UMAP_2, xend = label_x, yend = label_y),
    inherit.aes = FALSE,
    linewidth = 0.22,
    colour = "#8A9299",
    lineend = "round"
  ) +
  geom_text(
    data = label_layout,
    aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "plain",
    size = 2.35,
    colour = "#22282E"
  ) +
  scale_colour_manual(values = lineage_cols[present_lineages], drop = TRUE) +
  ## Preserve the full global atlas while adding a small visual margin so
  ## boundary cells and labels are not clipped at the panel edge.
  coord_equal(expand = TRUE) +
  labs(title = "C  Global single-cell atlas", colour = NULL) +
  theme_void_crlm() +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = grid::unit(0.8, "mm"),
    legend.key.height = grid::unit(2.7, "mm"),
    legend.key.width = grid::unit(2.7, "mm"),
    legend.text = element_text(size = 5.5, colour = COL_DARK),
    plot.margin = margin(2, 1.0, 2, 3, unit = "mm")
  )

## ------------------------------------------------------------------------
## Panel D — myeloid-only atlas UMAP
## ------------------------------------------------------------------------

nmf_umap_plot$frozen_broad_state <- factor(nmf_umap_plot$frozen_broad_state, levels = state_levels)
d_xlim <- quantile(nmf_umap_plot$UMAP_1, probs = c(0.001, 0.999), na.rm = TRUE, names = FALSE)
d_ylim <- quantile(nmf_umap_plot$UMAP_2, probs = c(0.001, 0.999), na.rm = TRUE, names = FALSE)
pD <- ggplot(nmf_umap_plot, aes(UMAP_1, UMAP_2, colour = frozen_broad_state)) +
  geom_point(size = 0.105, alpha = 0.95, stroke = 0) +
  scale_colour_manual(values = state_cols_D, drop = FALSE) +
  coord_fixed(xlim = d_xlim, ylim = d_ylim, expand = FALSE, clip = "on") +
  labs(title = "D  Myeloid atlas entering P1 discovery", colour = NULL) +
  theme_void_crlm() +
  theme(legend.position = "right") +
  guides(
    colour = guide_legend(
      override.aes = list(size = 2.1, alpha = 1)
    )
  )

## ------------------------------------------------------------------------
## Panel E — patient/tissue contribution to the 12,433-cell NMF input
## ------------------------------------------------------------------------

pE <- ggplot(composition, aes(patient, fraction, fill = frozen_broad_state)) +
  geom_col(width = 0.74) +
  scale_fill_manual(values = state_cols, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
  labs(title = "E  Patient and tissue contributions to the 12,433-cell NMF input",
       x = NULL, y = "Myeloid-state composition", fill = NULL) +
  theme_crlm() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6.5, face = "bold", colour = COL_DARK),
    axis.text.x = element_text(size = 5.7, angle = 45, hjust = 1),
    legend.position = "bottom"
  )

## ------------------------------------------------------------------------
## Assemble and export
## ------------------------------------------------------------------------

fig_S1 <-
  pA /
  panel_B /
  ((pC | pD) + plot_layout(widths = c(0.85, 1.15))) /
  pE +
  plot_layout(heights = c(0.58, 1.00, 1.25, 1.00))

save_pub_r(fig_S1, "Supplementary_Figure_S1", width_mm = 180, height_mm = 205)

figure_manifest <- data.frame(
  file = list.files(FIG_DIR, full.names = FALSE),
  bytes = file.info(file.path(FIG_DIR, list.files(FIG_DIR, full.names = FALSE)))$size,
  stringsAsFactors = FALSE
)
write_csv_plain(figure_manifest, "S1_figure_manifest.csv")

source_manifest <- data.frame(
  source_file = c(
    "02_qc_retention_by_sample.csv",
    "05_scdblfinder_by_sample.csv",
    "02_singlet_count_audit.csv",
    "04_global_analysis_cell_count_audit.csv",
    "17_lineage_subcluster_summary.csv",
    "36b_GSE178318_myeloid_UMAP_P1.csv.gz",
    "Figure1D_GSE178318_UMAP_P1.csv.gz",
    "06_GSE178318_global_broad_lineage.rds"
  ),
  role = c(
    "15-library raw/QC counts",
    "15-library singlet/doublet calls",
    "cohort-level singlet audit",
    "final atlas count audit",
    "lineage-level NMF input audit",
    "myeloid UMAP and frozen P1 source table",
    "authoritative frozen five-state source for S1D/S1E",
    "global UMAP and broad-lineage metadata"
  ),
  stringsAsFactors = FALSE
)
write_csv_plain(source_manifest, "S1_source_manifest.csv")

writeLines(c(
  "Supplementary Figure S1 — R-only rebuild",
  paste0("R version: ", R.version.string),
  paste0("Generated: ", format(Sys.time(), tz = "Asia/Shanghai")),
  "Core conclusion: the discovery input is transparent, QC-audited and not library-dominated.",
  "Counts: 140,281 raw barcodes; 119,400 QC cells; 108,353 singlets; 105,267 final-atlas cells; 12,433 NMF myeloid cells.",
  "C/D UMAP point labels are intentionally omitted; legends are kept separate from the embedding field.",
  "Panel D uses the frozen Figure 1 state field directly; the 36b fixed_state field is a retention/QC status and is exported as retention_status.",
  "Panel D axes are restricted to the 0.1st–99.9th percentile coordinate range for display only; all 12,433 cells remain in the source data and analyses.",
  "Exports: SVG/PDF plus 600-dpi TIFF/PNG.",
  paste0("Output directory: ", OUT_DIR)
), file.path(OUT_DIR, "README.txt"))

message("S1 R rebuild complete: ", OUT_DIR)
