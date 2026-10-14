## Figure 2 v7: external P1 portability with a real Seurat-based Figure 2A.
## Plotting-only. All plotted values are read from frozen/addendum source files;
## no model is refitted and no inferential result is changed.

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

required <- c("Seurat", "ggplot2", "patchwork", "dplyr", "scales", "ComplexHeatmap", "circlize", "cowplot", "tibble", "tidyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) stop("Missing required R package(s): ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(cowplot)
  library(tibble)
  library(tidyr)
})
options(stringsAsFactors = FALSE)

freeze <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
add <- file.path(freeze, "Figure2_source_data_addendum_20260912")
out_dir <- file.path(CRLM_OUTPUT_ROOT, "Figure2_external_scRNA_portability_v7")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(add, "external_expression/GSE245552_external_P1_top15_logCP10K_expression.csv.gz")
weights_path <- file.path(add, "external_expression/GSE245552_frozen_P1_top15_weights_used.csv")
rep231_candidates <- file.path(add, "external_scRNA", c(
  "62_GSE231559_directional_deltas.csv",
  "62_GSE231559_replication_deltas.csv"
))
rep231_path <- rep231_candidates[file.exists(rep231_candidates)][1]
cohort245_path <- file.path(add, "external_scRNA/59_GSE245552_P1_Pk_state_pair_deltas.csv")
seurat_path <- file.path(add, "external_Seurat/GSE245552_frozen.rds")
legend_path <- file.path(out_dir, "Figure2_external_scRNA_portability_legend.txt")
script_path <- file.path(CRLM_CODE_RELEASE_ROOT, "scripts", "figures", "Figure2_external_scRNA_portability_v7.R")

required_paths <- c(expr_path, weights_path, rep231_path, cohort245_path, seurat_path)
if (!all(file.exists(required_paths))) {
  stop("Required Figure 2 source file(s) not found: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}
read_source_table <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rt", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

expr24552 <- read_source_table(expr_path)
weights <- read_source_table(weights_path)
rep231559 <- read_source_table(rep231_path)
pairs24552 <- read_source_table(cohort245_path)
obj24552 <- readRDS(seurat_path)
if (!inherits(obj24552, "Seurat")) stop("Serialized object is not a Seurat object: ", seurat_path)
if (!all(c("broad_state", "frozen_P1_score") %in% colnames(obj24552[[]]))) {
  stop("Serialized GSE245552 object is missing broad_state or frozen_P1_score metadata.")
}
if (!"umap.external" %in% names(obj24552@reductions)) {
  stop("Serialized GSE245552 object is missing the umap.external reduction.")
}

required_expr <- c("cell_id", "sample", "patient", "tissue", "state", "P1_full")
required_rep <- c("score_version", "patient", "state", "delta_T_minus_N")
required_pairs <- c("program", "patient", "state", "delta_LM_minus_PT")
if (!all(required_expr %in% names(expr24552))) stop("Expression export is missing required metadata columns.")
if (!all(required_rep %in% names(rep231559))) stop("GSE231559 directional-concordance table is missing required columns.")
if (!all(required_pairs %in% names(pairs24552))) stop("GSE245552 pair table is missing required columns.")
if (!all(c("gene", "weight_norm", "rank") %in% names(weights))) stop("Frozen P1 weight table is missing required columns.")

weights <- weights[
  weights$program == "P1" & is.finite(weights$rank),
  , drop = FALSE
]
weights <- weights[order(weights$rank), , drop = FALSE]
top_genes <- as.character(weights$gene)
if (length(top_genes) != 15L || !all(top_genes %in% names(expr24552))) {
  stop("The expression export and authoritative frozen P1 top-15 table do not match.")
}
message("Frozen P1 gene order: ", paste(top_genes, collapse = " > "))

tissue_levels <- c("Primary tumor", "Liver metastasis", "Colon-adjacent tissue", "Liver-adjacent tissue")
tissue_cols <- c(
  "Primary tumor" = "#3D6F95", "Liver metastasis" = "#B05A5F",
  "Colon-adjacent tissue" = "#D39A4A", "Liver-adjacent tissue" = "#4D8D82"
)
state_levels <- c("Mono", "Inflammatory", "C1Q/MHC-II", "P1-state", "Other")
state_levels_common <- c("Mono", "Inflammatory", "C1Q/MHC-II", "Other")
state_cols <- c(
  "Mono" = "#5A9FC6", "Inflammatory" = "#D99A42", "C1Q/MHC-II" = "#4E9A82",
  "P1-state" = "#A23E48", "Other" = "#96949A"
)
COL_P1 <- "#A23E48"
COL_P1_DARK <- "#8E343E"
COL_DARK <- "#22282E"
COL_TEXT2 <- "#71808D"
MINUS <- intToUtf8(0x2212)
DELTA <- intToUtf8(0x0394)
PAT_STATE <- "patient-by-state"

map_state <- function(x) recode(
  x,
  P2_inflammatory = "Inflammatory",
  P3_C1Q_MHCII = "C1Q/MHC-II",
  P1_state = "P1-state",
  .default = x
)

expr24552 <- expr24552 %>%
  mutate(
    tissue_label = recode(
      tissue,
      primary_tumor = "Primary tumor",
      liver_metastasis = "Liver metastasis",
      colon_adjacent_tissue = "Colon-adjacent tissue",
      liver_adjacent_tissue = "Liver-adjacent tissue"
    ),
    state_label = map_state(state)
  ) %>%
  filter(state_label %in% state_levels, is.finite(P1_full))
rep231559 <- rep231559 %>% filter(score_version == "P1_full", is.finite(delta_T_minus_N)) %>%
  mutate(state_label = map_state(state))
pairs24552 <- pairs24552 %>% filter(program == "P1", is.finite(delta_LM_minus_PT)) %>%
  mutate(state_label = map_state(state))

expr24552 <- expr24552[order(
  factor(expr24552$state_label, levels = state_levels),
  factor(expr24552$tissue_label, levels = tissue_levels),
  -expr24552$P1_full
), , drop = FALSE]
expr_mat <- t(as.matrix(expr24552[, top_genes, drop = FALSE]))
storage.mode(expr_mat) <- "numeric"
expr_z <- t(apply(expr_mat, 1L, function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(v)))
  as.numeric((v - mean(v, na.rm = TRUE)) / s)
}))
expr_z <- pmax(pmin(expr_z, 2.5), -2.5)
rownames(expr_z) <- top_genes
colnames(expr_z) <- expr24552$cell_id

p1_lim <- quantile(expr24552$P1_full, probs = c(0.02, 0.98), na.rm = TRUE)
p1_col_fun <- circlize::colorRamp2(p1_lim, c("#F4F1EF", COL_P1))
expr_col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#557C9B", "#F7F6F4", "#B85F55"))

## A: raw-vector frozen P1 projection from the serialized object.
## The two plots use exactly the same umap.external embedding.
## No rasterization, ggrastr, or Seurat automatic plotting is used for A.
emb_A <- as.data.frame(Seurat::Embeddings(obj24552, reduction = "umap.external"))
colnames(emb_A)[1:2] <- c("UMAP_1", "UMAP_2")
emb_A$cell_key <- rownames(emb_A)
meta_A <- obj24552[[]] %>%
  dplyr::select(broad_state, frozen_P1_score) %>%
  tibble::rownames_to_column("cell_key")
plot_A <- emb_A %>%
  dplyr::left_join(meta_A, by = "cell_key") %>%
  dplyr::mutate(broad_state = factor(broad_state, levels = state_levels)) %>%
  dplyr::filter(!is.na(broad_state), is.finite(frozen_P1_score))

set.seed(20260913)
df_state_A <- plot_A %>% dplyr::slice_sample(prop = 1)

pA1 <- ggplot(df_state_A, aes(UMAP_1, UMAP_2, colour = broad_state)) +
  geom_point(size = 0.19, alpha = 1, stroke = 0) +
  scale_colour_manual(values = state_cols, drop = FALSE) +
  coord_equal(expand = FALSE) +
  labs(title = "Myeloid state", colour = NULL) +
  guides(
    colour = guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(size = 1.7, alpha = 1)
    )
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 7.8, face = "bold", colour = COL_DARK,
      hjust = 0, margin = margin(b = 0.8, unit = "mm")),
    legend.position = "bottom",
    legend.text = element_text(size = 5.6, colour = COL_DARK),
    legend.key.width = grid::unit(2.6, "mm"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "mm")
  )

df_p1_A <- plot_A %>% dplyr::arrange(frozen_P1_score)
p1_limits_A <- quantile(df_p1_A$frozen_P1_score, c(0.02, 0.98), na.rm = TRUE)
pA2 <- ggplot(df_p1_A, aes(UMAP_1, UMAP_2, colour = frozen_P1_score)) +
  geom_point(size = 0.16, alpha = 1, stroke = 0) +
  scale_colour_gradientn(
    colours = c(
      "#CBD6DC", "#E5E3E1", "#E8C6BE", "#DD9688", "#CB685B", "#A54149"
    ),
    limits = p1_limits_A, oob = scales::squish, name = NULL,
    guide = guide_colourbar(
      direction = "horizontal", barwidth = grid::unit(24, "mm"),
      barheight = grid::unit(1.8, "mm"), ticks = TRUE
    )
  ) +
  coord_equal(expand = FALSE) +
  labs(title = "Frozen P1 score") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 7.8, face = "bold", colour = COL_DARK,
      hjust = 0, margin = margin(b = 0.8, unit = "mm")),
    legend.position = "bottom",
    legend.text = element_text(size = 5.6, colour = COL_DARK),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, unit = "mm")
  )

pA_body <- pA1 | pA2
pA_body <- pA_body + patchwork::plot_layout(widths = c(1, 1))
pA <- pA_body + patchwork::plot_annotation(
  title = "A  Frozen P1 projection in GSE245552",
  theme = theme(
    plot.title = element_text(
      family = "Arial", face = "bold", size = 9.5, colour = COL_DARK,
      hjust = 0, margin = margin(b = 0.8, unit = "mm")
    )
  )
)

## B: dense cell-level expression heatmap.
expr_state <- factor(expr24552$state_label, levels = state_levels)
expr_tissue <- factor(expr24552$tissue_label, levels = tissue_levels)
ha_expr <- ComplexHeatmap::HeatmapAnnotation(
  State = ComplexHeatmap::anno_simple(as.character(expr_state), col = state_cols, border = FALSE),
  `Tissue origin` = ComplexHeatmap::anno_simple(as.character(expr_tissue), col = tissue_cols, border = FALSE),
  `Frozen P1 score` = ComplexHeatmap::anno_simple(as.numeric(expr24552$P1_full), col = p1_col_fun, border = FALSE),
  annotation_name_gp = grid::gpar(fontsize = 5.3, fontface = "plain", col = "#66727B"),
  simple_anno_size = grid::unit(1.45, "mm"),
  gap = grid::unit(0.35, "mm")
)
ht_expr <- ComplexHeatmap::Heatmap(
  expr_z, name = "Expression", col = expr_col_fun, top_annotation = ha_expr,
  column_split = expr_state, column_gap = grid::unit(0.9, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE, cluster_column_slices = FALSE,
  show_column_names = FALSE, show_row_names = TRUE,
  row_names_side = "left", row_names_gp = grid::gpar(fontsize = 6.25, fontface = "italic", col = COL_DARK),
  column_title_gp = grid::gpar(fontsize = 7.4, fontface = "bold", col = COL_DARK),
  column_title_rot = 0, row_title = NULL, border = FALSE, use_raster = TRUE,
  raster_quality = 3, raster_device = "png",
  heatmap_legend_param = list(
    title = "Expression\n(row z-score)", at = c(-2.5, -1.25, 0, 1.25, 2.5),
    labels = c("-2.5", "-1.25", "0", "1.25", "2.5"), labels_gp = grid::gpar(fontsize = 5.7),
    title_gp = grid::gpar(fontsize = 6.3, fontface = "bold"), direction = "horizontal",
    legend_width = grid::unit(28, "mm")
  )
)
grob_expr <- grid::grid.grabExpr(
  ComplexHeatmap::draw(
    ht_expr, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
    annotation_legend_list = list(
      ComplexHeatmap::Legend(
        labels = names(tissue_cols), title = "Tissue origin", nrow = 1,
        legend_gp = grid::gpar(fill = unname(tissue_cols), col = NA),
        labels_gp = grid::gpar(fontsize = 5.5, col = COL_DARK),
        title_gp = grid::gpar(fontsize = 6.1, fontface = "bold", col = COL_DARK),
        grid_width = grid::unit(2.4, "mm"), grid_height = grid::unit(2.4, "mm"),
        gap = grid::unit(1.8, "mm")
      ),
      ComplexHeatmap::Legend(
        col_fun = p1_col_fun, title = "Frozen P1 score", direction = "horizontal",
        labels_gp = grid::gpar(fontsize = 5.5, col = COL_DARK),
        title_gp = grid::gpar(fontsize = 6.1, fontface = "bold", col = COL_DARK),
        legend_width = grid::unit(21, "mm")
      )
    ),
    merge_legends = TRUE
  )
)
title_expr <- cowplot::ggdraw() +
  cowplot::draw_label("B  Frozen P1 gene program across external myeloid cells", x = 0, y = 0.82,
    hjust = 0, vjust = 1, fontfamily = "Arial", fontface = "bold", size = 8.8) +
  cowplot::draw_label("GSE245552; cells ordered by state, tissue origin and frozen P1 score",
    x = 0, y = 0.13, hjust = 0, vjust = 0, fontfamily = "Arial", colour = COL_TEXT2, size = 6.5)
panel_expr <- cowplot::plot_grid(title_expr, grob_expr, ncol = 1, rel_heights = c(0.13, 0.87))

## Standalone B export for focused visual QA; this does not alter the assembled figure.
b_only_png <- file.path(out_dir, "Figure2B_external_P1_program_standalone.png")
b_only_pdf <- file.path(out_dir, "Figure2B_external_P1_program_standalone.pdf")
ggsave(
  b_only_png, panel_expr, width = 180, height = 85, units = "mm", dpi = 600, bg = "white"
)
ggsave(
  b_only_pdf, panel_expr, width = 180, height = 85, units = "mm", device = cairo_pdf, bg = "white"
)

## C: patient-level dot plot for independent GSE231559 directional concordance.
## Only evaluable patient-by-state contrasts are plotted; there are no NA tiles.
state_order_C <- c("Mono", "Inflammatory", "C1Q/MHC-II", "Other")
state_cols_C <- state_cols[state_order_C]
state_cols_C["Other"] <- "#85878D"
Cdat <- rep231559 %>%
  dplyr::filter(
    score_version == "P1_full",
    is.finite(delta_T_minus_N),
    state_label %in% state_order_C
  ) %>%
  dplyr::mutate(state = factor(state_label, levels = state_order_C))
Csum <- Cdat %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(
    n_pos = sum(delta_T_minus_N > 0),
    n_eval = dplyr::n(),
    median_delta = median(delta_T_minus_N),
    .groups = "drop"
  ) %>%
  dplyr::right_join(
    data.frame(state = factor(state_order_C, levels = state_order_C)),
    by = "state"
  ) %>%
  dplyr::mutate(
    n_pos = ifelse(is.na(n_pos), 0L, n_pos),
    n_eval = ifelse(is.na(n_eval), 0L, n_eval),
    label = ifelse(n_eval > 0, paste0(n_pos, "/", n_eval), ""),
    x = seq_along(state_order_C)
  )
yr_C <- range(Cdat$delta_T_minus_N, na.rm = TRUE)
yr_span_C <- diff(yr_C)
if (!is.finite(yr_span_C) || yr_span_C == 0) yr_span_C <- 1
pad_bottom_C <- max(0.18, 0.23 * yr_span_C)
pad_top_C <- max(0.10, 0.10 * yr_span_C)
# Keep the positive/evaluable counts clearly above the highest point.
label_y_C <- yr_C[2] + 0.85 * pad_top_C
ylow_C <- yr_C[1] - pad_bottom_C
yhigh_C <- yr_C[2] + 1.6 * pad_top_C
y_grid_C <- pretty(c(ylow_C, yhigh_C), n = 5)
y_grid_C <- y_grid_C[y_grid_C >= ylow_C & y_grid_C <= yhigh_C]
pC <- ggplot(Cdat, aes(state, delta_T_minus_N)) +
  geom_hline(
    yintercept = y_grid_C,
    linewidth = 0.25, colour = "#E5EBEF"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "#8E98A1") +
  geom_point(
    aes(fill = state), shape = 21, colour = "white", stroke = 0.25,
    size = 2.2,
    position = position_jitter(width = 0.085, height = 0, seed = 20260913),
    alpha = 0.95
  ) +
  geom_segment(
    data = Csum %>% dplyr::filter(n_eval > 0),
    aes(x = x - 0.17, xend = x + 0.17, y = median_delta, yend = median_delta),
    inherit.aes = FALSE,
    linewidth = 0.65, colour = "#59636B", lineend = "round"
  ) +
  geom_text(
    data = Csum,
    aes(x = state, y = label_y_C, label = label),
    inherit.aes = FALSE, size = 2.5, fontface = "bold", colour = COL_TEXT2
  ) +
  scale_fill_manual(values = state_cols_C, guide = "none") +
  scale_y_continuous(limits = c(ylow_C, yhigh_C), expand = expansion(mult = c(0, 0.01))) +
  labs(
    title = "C  Directional concordance in GSE231559",
    subtitle = paste0(
      paste0("Paired T", MINUS, "N differences across evaluable patient-state pairs")
    ),
    x = NULL, y = paste0(DELTA, " frozen P1 score (T", MINUS, "N)")
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 8.5, face = "bold", colour = COL_DARK, margin = margin(b = 0.8, unit = "mm")),
    plot.subtitle = element_text(size = 6.2, colour = COL_TEXT2, margin = margin(b = 1.4, unit = "mm")),
    axis.title.y = element_text(size = 6.6, face = "plain", colour = COL_DARK),
    axis.text.x = element_text(size = 6.5, face = "bold", colour = COL_DARK),
    axis.text.y = element_text(size = 6.0, colour = COL_DARK),
    axis.line.x = element_blank(), axis.line.y = element_line(linewidth = 0.45, colour = "#8E98A1"),
    axis.ticks.x = element_blank(), axis.ticks.y = element_line(linewidth = 0.35, colour = "#8E98A1"),
    panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    plot.margin = margin(1, 1, 1, 1, unit = "mm")
  )

## D: cohort-level dumbbell comparison of positive patient-by-state fractions.
summarise_D <- function(df, delta_col, cohort_name) {
  dat <- df %>%
    dplyr::filter(state_label %in% state_levels_common, is.finite(.data[[delta_col]])) %>%
    dplyr::mutate(delta = as.numeric(.data[[delta_col]]))
  by_state <- dat %>%
    dplyr::group_by(state_label) %>%
    dplyr::summarise(evaluable = dplyr::n(), positive = sum(delta > 0), .groups = "drop") %>%
    dplyr::mutate(fraction = positive / evaluable, cohort = cohort_name)
  overall <- dat %>%
    dplyr::summarise(state_label = "Overall (shared states)", evaluable = dplyr::n(), positive = sum(delta > 0)) %>%
    dplyr::mutate(fraction = positive / evaluable, cohort = cohort_name)
  out <- dplyr::bind_rows(by_state, overall)
  out$count <- paste0(out$positive, "/", out$evaluable)
  out
}
summary_24552_D <- summarise_D(pairs24552, "delta_LM_minus_PT", "GSE245552")
summary_231559_D <- summarise_D(rep231559, "delta_T_minus_N", "GSE231559")
state_order_D <- c("Mono", "Inflammatory", "C1Q/MHC-II", "Other", "Overall (shared states)")
Ddat <- dplyr::bind_rows(summary_24552_D, summary_231559_D) %>%
  dplyr::mutate(
    state = factor(state_label, levels = rev(state_order_D)),
    cohort = factor(cohort, levels = c("GSE245552", "GSE231559")),
    y_base = setNames(seq_along(rev(state_order_D)), rev(state_order_D))[as.character(state)],
    y = y_base + ifelse(cohort == "GSE245552", 0.28, -0.28),
    label_x = 1.045,
    label_hjust = 0
  )
cohort_cols_D <- c("GSE245552" = "#B64B3A", "GSE231559" = "#7A8FA8")
overall_row_D <- data.frame(ymin = 1 - 0.38, ymax = 1 + 0.38, xmin = 0, xmax = 1.28)
pD <- ggplot() +
  geom_rect(data = overall_row_D, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), inherit.aes = FALSE, fill = "#F6F8FA", colour = NA) +
  geom_hline(yintercept = seq_along(rev(state_order_D)), linewidth = 0.25, colour = "#E6EBEF") +
  annotate("point", x = 0.70, y = 5.92, colour = cohort_cols_D["GSE245552"], size = 1.75) +
  annotate("text", x = 0.74, y = 5.92, label = "GSE245552", hjust = 0, vjust = 0.5,
    size = 2.00, colour = COL_DARK, family = "Arial") +
  annotate("point", x = 0.70, y = 5.58, colour = cohort_cols_D["GSE231559"], size = 1.75) +
  annotate("text", x = 0.74, y = 5.58, label = "GSE231559", hjust = 0, vjust = 0.5,
    size = 2.00, colour = COL_DARK, family = "Arial") +
  
  geom_point(data = Ddat, aes(x = fraction, y = y, fill = cohort), shape = 21, size = 2.35, stroke = 0.25, colour = "white") +
  geom_text(data = Ddat, aes(x = label_x, y = y, label = count, colour = cohort, hjust = label_hjust), size = 2.1, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = cohort_cols_D, guide = "none") +
  scale_colour_manual(values = cohort_cols_D, guide = "none") +
  scale_x_continuous(limits = c(0, 1.28), breaks = c(0, 0.25, 0.5, 0.75, 1.0), labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(breaks = seq_along(rev(state_order_D)), labels = rev(state_order_D), limits = c(0.5, length(state_order_D) + 1.0), expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "D  Cross-cohort directional concordance",
    subtitle = paste0("Fraction of evaluable patient-by-state contrasts with ", DELTA, "P1 > 0"),
    x = "Positive contrasts (%)", y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_family = "Arial") +
  theme(
    plot.title = element_text(size = 7.8, face = "bold", colour = COL_DARK, margin = margin(b = 0.8, unit = "mm")),
    plot.subtitle = element_text(size = 5.8, colour = COL_TEXT2, margin = margin(b = 1.0, unit = "mm")),
    axis.title.x = element_text(size = 6.0, face = "plain", colour = COL_DARK),
    axis.text.x = element_text(size = 5.5, colour = COL_DARK),
    axis.text.y = element_text(size = 5.9, face = "bold", colour = COL_DARK),
    axis.line.y = element_blank(), axis.line.x = element_line(linewidth = 0.45, colour = "#9AA4AE"),
    axis.ticks.y = element_blank(), axis.ticks.x = element_line(linewidth = 0.4, colour = "#9AA4AE"),
    legend.position = "none", panel.grid = element_blank(), plot.margin = margin(1.5, 2, 2.5, 2, unit = "mm")
  )

## Standalone C/D export for focused visual QA; the assembled figure is unchanged.
cd_only <- cowplot::plot_grid(pC, pD, ncol = 2, rel_widths = c(1.22, 1.00), align = "h")
cd_only_png <- file.path(out_dir, "Figure2C_D_directional_concordance_standalone.png")
cd_only_pdf <- file.path(out_dir, "Figure2C_D_directional_concordance_standalone.pdf")
ggplot2::ggsave(cd_only_png, cd_only, width = 180, height = 82, units = "mm", dpi = 600, bg = "white")
ggplot2::ggsave(cd_only_pdf, cd_only, width = 180, height = 82, units = "mm", device = grDevices::cairo_pdf, bg = "white")

pA_grob <- grid::grid.grabExpr(print(pA))
bottom_row <- cd_only
fig4 <- cowplot::plot_grid(pA_grob, panel_expr, bottom_row, ncol = 1, rel_heights = c(0.18, 0.45, 0.37))

png_path <- file.path(out_dir, "Figure2_external_scRNA_portability.png")
pdf_path <- file.path(out_dir, "Figure2_external_scRNA_portability.pdf")
svg_path <- file.path(out_dir, "Figure2_external_scRNA_portability.svg")
tiff_path <- file.path(out_dir, "Figure2_external_scRNA_portability.tiff")
fig_width_mm <- 180
fig_height_mm <- 200
png(filename = png_path, width = round(fig_width_mm / 25.4 * 600), height = round(fig_height_mm / 25.4 * 600), res = 600, type = "cairo")
print(fig4)
dev.off()
cairo_pdf(filename = pdf_path, width = fig_width_mm / 25.4, height = fig_height_mm / 25.4, family = "Arial")
print(fig4)
dev.off()
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(svg_path, width = fig_width_mm / 25.4, height = fig_height_mm / 25.4)
  print(fig4)
  dev.off()
}
tiff(filename = tiff_path, width = round(fig_width_mm / 25.4 * 600), height = round(fig_height_mm / 25.4 * 600),
  res = 600, compression = "lzw", type = "cairo")
print(fig4)
dev.off()

legend_lines <- c(
  "Figure 2. External portability of frozen P1 across single-cell cohorts.", "",
  "A, Frozen P1 projection in the independent GSE245552 myeloid cell set. The left",
  "plot shows broad myeloid states and the right plot shows the continuous frozen P1",
  "score on the identical serialized umap.external embedding. Both maps use raw",
  "coordinates and vector ggplot2 point layers; scores were clipped to the 2nd-98th",
  "percentiles for visualization only and source values are unchanged.", "",
  "B, Cell-level expression view of frozen P1 in GSE245552. The 24,231 matched cells",
  "are ordered by conventional state, tissue and descending frozen P1 score. Rows are",
  "the top-15 genes read directly from the authoritative frozen P1 normalized loading",
  "weight table. Expression is log1p counts per 10,000 followed by row-wise z-scoring",
  "for display. This is a visualization-only derived export.", "",
  "C, Directional concordance in GSE231559. Each point represents one evaluable",
  paste0("patient-by-state T", MINUS, "N contrast from the frozen P1_full score. A patient-by-state pair"),
  "is evaluable only when at least 20 cells are present on each side (N and T).",
  "Point colour denotes the shared myeloid state and black bars denote state medians.",
  "Counts above points are positive/evaluable pairs. P1-state yielded no evaluable",
  "pair under this rule and is therefore omitted from C and from the D shared-state",
  "summary.", "",
  "D, Cross-cohort directional concordance across shared states. Points show the",
  "fraction of evaluable patient-by-state contrasts; labels show positive/evaluable",
  "counts. Overall (shared states) combines the four states shown; it excludes",
  "P1-state because GSE231559 had no evaluable pair for that state. GSE245552 uses",
  paste0("LM", MINUS, "PT and GSE231559 uses T", MINUS, "N contrasts; absolute score magnitudes"),
  "were not compared.", "",
  "All panels are plotting-only displays of frozen/addendum source tables; no model was",
  "refitted for this figure. The new expression and UMAP files are explicitly authorized",
  "visualization-only derived files in the Figure 2 addendum.", "",
  paste0("P1 top-15 genes (authoritative descending normalized loading): ", paste(top_genes, collapse = ", ")), "",
  "Source script: R figure-generation script (v7)"
)
## Preserve UTF-8 minus/delta symbols in the exported legend on Windows.
writeLines(legend_lines, legend_path, useBytes = TRUE)

message("Figure 2 v7 exported")
message("External expression cells: ", nrow(expr24552), "; genes: ", length(top_genes))
message("GSE245552 pairs: ", sum(pairs24552$delta_LM_minus_PT > 0), "/", nrow(pairs24552), " positive")
message("GSE231559 pairs: ", sum(rep231559$delta_T_minus_N > 0), "/", nrow(rep231559), " positive")
message("PNG: ", png_path)
message("PDF: ", pdf_path)
message("SVG: ", svg_path)
message("TIFF: ", tiff_path)
