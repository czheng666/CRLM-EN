# Figure 1 | CRLM P1 discovery and cross-state characterisation
# R-only, submission-grade draft using the current authoritative analysis outputs.
# This script reads frozen/current result tables and writes a new figure bundle.

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
  library(grid)
  library(svglite)
  library(ragg)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

results_root <- crlm_require_dir(CRLM_RESULTS_ROOT, "CRLM_RESULTS_ROOT")
freeze_dir <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
project_root <- CRLM_PROJECT_ROOT
out_dir <- file.path(CRLM_OUTPUT_ROOT, "Figure1_P1_actual_v18")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

table_dir <- file.path(results_root, "results", "tables")
need_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path)
  path
}
need_cols <- function(x, cols, label) {
  missing_cols <- setdiff(cols, names(x))
  if (length(missing_cols)) {
    stop(label, " is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  invisible(x)
}
read_result_table <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    tibble::as_tibble(utils::read.csv(gzfile(path), check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""))
  } else {
    tibble::as_tibble(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, comment.char = ""))
  }
}

umap_path <- need_file(file.path(table_dir, "36b_GSE178318_myeloid_UMAP_P1.csv.gz"))
usage_path <- need_file(file.path(table_dir, "36_discovery_program_usage.csv.gz"))
rank_path <- need_file(file.path(table_dir, "28_NMF_rank_summary.csv"))
selected_rank_path <- need_file(file.path(table_dir, "31_NMF_selected_rank.csv"))
genes_path <- need_file(file.path(table_dir, "34_frozen_P1_top15_weights.csv"))
lineage_path <- need_file(file.path(table_dir, "17_lineage_subcluster_summary.csv"))
benchmark_path <- need_file(file.path(table_dir, "58_GSE245552_P1_Pk_specificity_summary.csv"))

message("Resolved results root: ", results_root)
message("Resolved table directory: ", table_dir)

umap_df <- read_result_table(umap_path)
usage_df <- read_result_table(usage_path)
rank_df <- read_result_table(rank_path)
selected_df <- read_result_table(selected_rank_path)
genes_df <- read_result_table(genes_path)
lineage_df <- read_result_table(lineage_path)
benchmark_df <- read_result_table(benchmark_path)

need_cols(
  umap_df,
  c("barcode", "patient", "sample", "fine_state", "P1_top15_score", "UMAP_1", "UMAP_2"),
  "UMAP/P1 table"
)
need_cols(
  usage_df,
  c("barcode", "sample_id", "patient_id", "fine_label", "P1", "P2", "P3", "P4", "P5", "P6"),
  "program usage table"
)
need_cols(
  rank_df,
  c("rank", "n_runs.x", "median_pairwise_basis_spearman", "cophenetic_coefficient", "median_pairwise_top15_overlap"),
  "NMF rank audit table"
)
need_cols(selected_df, c("selected_rank"), "selected-rank table")
need_cols(genes_df, c("program", "rank", "gene"), "frozen P1 gene table")
need_cols(benchmark_df, c("program", "median_effect", "ci_low", "ci_high", "n_positive_patients", "positive_state_pair_fraction"), "P1-P6 benchmark table")

selected_rank <- as.integer(selected_df$selected_rank[[1]])
if (!identical(selected_rank, 6L)) {
  stop("The current selected rank is not 6; refusing to draw a rank-6 frozen figure.")
}

state_levels <- c("Mono", "Inflammatory", "C1Q/MHC-II", "P1-state", "Other")
state_cols <- c(
  "Mono" = "#4A8DB7",
  "Inflammatory" = "#C9862C",
  "C1Q/MHC-II" = "#3D8B74",
  "P1-state" = "#B64A3A",
  "Other" = "#8C8791"
)
ink <- "#20262E"
mid_ink <- "#6F7D8A"
rule <- "#CCD3D9"
p1_accent <- "#B64A3A"
heat_low <- "#F3F5F6"
heat_mid <- "#E8B7A9"
heat_high <- "#9F302F"

map_state <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  out <- dplyr::case_when(
    grepl("S100A8|S100A9|inflamm|IFN|FCN1|IL1B|TNF", x, ignore.case = TRUE) ~ "Inflammatory",
    grepl("Mono", x, ignore.case = TRUE) ~ "Mono",
    grepl("C1Q|FOLR2|MHC|HLA|CD74|FCER1A", x, ignore.case = TRUE) ~ "C1Q/MHC-II",
    grepl("SPP1|GPNMB|MARCO|APOC1|CCL18|P1", x, ignore.case = TRUE) ~ "P1-state",
    TRUE ~ "Other"
  )
  factor(out, levels = state_levels)
}

umap_df <- umap_df %>%
  mutate(
    state = map_state(fine_state),
    fine_state = as.character(fine_state)
  )
usage_df <- usage_df %>%
  mutate(state = map_state(fine_label))

current_n <- nrow(usage_df)
umap_n <- nrow(umap_df)
if (current_n != 12433L || umap_n != 12433L) {
  stop("Current authoritative Figure 1 discovery set must contain 12,433 rows; found usage=",
       current_n, ", UMAP=", umap_n)
}

# Freeze-package reconciliation is checked when present, but never edited.
reconciliation_path <- file.path(freeze_dir, "tables", "discovery_cell_set_reconciliation_summary.csv")
reconciliation_current_n <- NA_integer_
if (file.exists(reconciliation_path)) {
  reconciliation_df <- read_result_table(reconciliation_path)
  if (all(c("current_n") %in% names(reconciliation_df))) {
    reconciliation_current_n <- as.integer(reconciliation_df$current_n[[1]])
  }
}
if (!is.na(reconciliation_current_n) && reconciliation_current_n != 12433L) {
  stop("Freeze-package reconciliation does not report 12,433 current cells.")
}

if (dplyr::n_distinct(usage_df$patient_id) != 6L || dplyr::n_distinct(usage_df$sample_id) != 15L) {
  stop("Current discovery population does not reconcile to 6 patients / 15 samples.")
}
myeloid_dc_annotation_n <- if (all(c("n_clusters", "lineage") %in% names(lineage_df))) {
  as.integer(lineage_df$n_clusters[lineage_df$lineage == "Myeloid_DC"][[1]])
} else {
  NA_integer_
}
if (is.na(myeloid_dc_annotation_n) || myeloid_dc_annotation_n != 14L) {
  stop("Current discovery population does not reconcile to 14 annotated fine states.")
}
rank_audit_rows <- rank_df %>% filter(rank %in% 4:8)
if (nrow(rank_audit_rows) != 5L || !all(rank_audit_rows$n_runs.x == 30L)) {
  stop("NMF rank 4-8 audit does not reconcile to 30 runs per rank.")
}
rank6_audit_row <- rank_audit_rows %>% filter(as.integer(rank) == selected_rank)
if (nrow(rank6_audit_row) != 1L || round(as.numeric(rank6_audit_row$median_pairwise_basis_spearman), 3) != 0.972 ||
    round(as.numeric(rank6_audit_row$cophenetic_coefficient), 3) != 0.985 ||
    round(as.numeric(rank6_audit_row$median_pairwise_top15_overlap) * 15) != 14L) {
  stop("Rank-6 metrics do not reconcile to 0.972 / 0.985 / 14/15.")
}

rank_plot_df <- rank_df %>%
  filter(rank %in% 4:8) %>%
  transmute(
    rank = as.numeric(rank),
    `Basis stability` = median_pairwise_basis_spearman,
    `Cophenetic coefficient` = cophenetic_coefficient,
    `Top-15 gene overlap` = median_pairwise_top15_overlap
  ) %>%
  pivot_longer(-rank, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Basis stability", "Cophenetic coefficient", "Top-15 gene overlap")))

gene_value_candidates <- intersect(c("weight_normalized", "weight_norm", "weight_norm_within_top15", "weight"), names(genes_df))
if (length(gene_value_candidates) == 0) {
  stop("Frozen P1 gene table does not contain a normalized weight column.")
}
gene_value_col <- gene_value_candidates[[1]]
p1_genes <- genes_df %>%
  filter(toupper(program) == "P1") %>%
  mutate(gene_value = as.numeric(.data[[gene_value_col]])) %>%
  arrange(desc(gene_value)) %>%
  slice_head(n = 15)
if (nrow(p1_genes) != 15L) stop("Frozen P1 gene table does not contain exactly 15 P1 genes.")
expected_p1_genes <- c(
  "FTL", "APOC1", "SPP1", "APOE", "CTSD", "CSTB", "LGALS1", "ACTB", "CTSB", "VIM", "LYZ", "CD74", "CD63", "PSAP", "HLA-A"
)
gene_set_warning <- if (!setequal(p1_genes$gene, expected_p1_genes)) {
  paste0("WARNING: actual frozen P1 gene set differs from the supplied expected list: ", paste(setdiff(expected_p1_genes, p1_genes$gene), collapse = ", "))
} else {
  "Frozen P1 top-15 gene set matches the supplied current list."
}

p1_genes <- p1_genes %>%
  mutate(gene = factor(gene, levels = rev(gene)))
p1_gene_order <- as.character(p1_genes$gene)
if (length(p1_gene_order) != 15L || !all(diff(p1_genes$gene_value) <= 1e-12)) {
  stop("Frozen P1 top-15 genes are not sorted by descending normalized loading weight.")
}

rank6_labels <- rank_plot_df %>%
  group_by(metric) %>%
  mutate(label_y = value + 0.06 * (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))) %>%
  ungroup() %>%
  filter(rank == selected_rank) %>%
  mutate(
    value_label = dplyr::case_when(
      as.character(metric) == "Top-15 gene overlap" ~ paste0(round(value * 15), "/15"),
      TRUE ~ formatC(value, format = "f", digits = 3)
    )
  )

fine_state_summary <- umap_df %>%
  filter(!is.na(P1_top15_score)) %>%
  group_by(state, fine_state) %>%
  summarise(
    mean_P1_score = mean(P1_top15_score, na.rm = TRUE),
    median_P1_score = median(P1_top15_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  ) %>%
  mutate(state = factor(state, levels = state_levels)) %>%
  arrange(state, desc(mean_P1_score))
fine_state_summary <- fine_state_summary %>%
  mutate(
    point_size = scales::rescale(sqrt(n_cells), to = c(2.4, 3.8)),
    point_alpha = scales::rescale(log1p(n_cells), to = c(0.60, 0.90))
  )
state_medians <- fine_state_summary %>%
  group_by(state) %>%
  summarise(state_median = median(mean_P1_score, na.rm = TRUE), .groups = "drop") %>%
  mutate(x = as.numeric(state), xmin = x - 0.23, xmax = x + 0.23)

benchmark_denominator <- 50L
benchmark_order <- benchmark_df %>%
  arrange(desc(median_effect)) %>%
  pull(program) %>%
  as.character()
benchmark_plot_df <- benchmark_df %>%
  transmute(
    program = factor(as.character(program), levels = rev(benchmark_order)),
    median_effect = as.numeric(median_effect),
    ci_low = as.numeric(ci_low),
    ci_high = as.numeric(ci_high),
    positive_patients = paste0(as.integer(n_positive_patients), "/13"),
    positive_pairs = paste0(round(as.numeric(positive_state_pair_fraction) * benchmark_denominator), "/", benchmark_denominator),
    is_p1 = as.character(program) == "P1"
  )
benchmark_label_x_patients <- 1.00
benchmark_label_x_pairs <- 1.45
benchmark_header_y <- length(benchmark_order) + 0.50

p1_benchmark_row <- benchmark_df %>% filter(as.character(program) == "P1")
if (nrow(p1_benchmark_row) != 1L || round(as.numeric(p1_benchmark_row$median_effect), 3) != 0.608 ||
    as.integer(p1_benchmark_row$n_positive_patients) != 12L ||
    round(as.numeric(p1_benchmark_row$positive_state_pair_fraction) * benchmark_denominator) != 45L) {
  stop("Frozen P1 benchmark does not reconcile to 0.608 / 12/13 / 45/50.")
}

figure1_source_audit <- tibble::tibble(
  check = c(
    "authoritative discovery cells", "patients / samples", "annotated fine states",
    "NMF ranks / runs", "selected rank", "rank-6 metrics", "P1 top-15 order", "P1 benchmark"
  ),
  observed = c(
    as.character(current_n),
    paste0(dplyr::n_distinct(usage_df$patient_id), " / ", dplyr::n_distinct(usage_df$sample_id)),
    as.character(myeloid_dc_annotation_n),
    paste0("4-8 / ", paste(unique(rank_audit_rows$n_runs.x), collapse = ", ")),
    as.character(selected_rank),
    paste0(formatC(as.numeric(rank6_audit_row$median_pairwise_basis_spearman), format = "f", digits = 3),
           " / ", formatC(as.numeric(rank6_audit_row$cophenetic_coefficient), format = "f", digits = 3),
           " / ", round(as.numeric(rank6_audit_row$median_pairwise_top15_overlap) * 15), "/15"),
    paste(p1_gene_order, collapse = ", "),
    paste0(formatC(as.numeric(p1_benchmark_row$median_effect), format = "f", digits = 3),
           " / ", as.integer(p1_benchmark_row$n_positive_patients), "/13 / ",
           round(as.numeric(p1_benchmark_row$positive_state_pair_fraction) * benchmark_denominator), "/",
           benchmark_denominator)
  ),
  expected = c(
    "12,433", "6 / 15", "14", "4-8 / 30", "6", "0.972 / 0.985 / 14/15",
    "descending normalized frozen weights", "0.608 / 12/13 / 45/50"
  ),
  status = "PASS"
)

population_counts <- tibble::tibble(
  item = c("authoritative_discovery_cells", "patients", "samples", "annotation_clusters", "selected_rank"),
  value = c(
    nrow(usage_df),
    dplyr::n_distinct(usage_df$patient_id),
    dplyr::n_distinct(usage_df$sample_id),
    if (all(c("n_clusters", "lineage") %in% names(lineage_df))) lineage_df$n_clusters[lineage_df$lineage == "Myeloid_DC"][[1]] else NA_real_,
    selected_rank
  )
)

write_csv(umap_df %>% select(barcode, patient, sample, fine_state, state, P1_usage, P1_top15_score, UMAP_1, UMAP_2),
          file.path(source_dir, "Figure1D_GSE178318_UMAP_P1.csv.gz"))
write_csv(rank_plot_df, file.path(source_dir, "Figure1B_rank_stability_audit.csv"))
write_csv(p1_genes %>% mutate(gene = as.character(gene)) %>% select(-gene_value),
          file.path(source_dir, "Figure1C_frozen_P1_top15_genes.csv"))
write_csv(fine_state_summary, file.path(source_dir, "Figure1D_P1_by_fine_state.csv"))
write_csv(benchmark_plot_df %>% mutate(program = as.character(program)),
          file.path(source_dir, "Figure1E_P1_P6_prioritization_benchmark.csv"))
write_csv(population_counts, file.path(source_dir, "Figure1A_B_population_counts.csv"))
write_csv(figure1_source_audit, file.path(out_dir, "Figure1_source_data_audit.csv"))

base_theme <- theme_minimal(base_family = "Arial", base_size = 8.5) +
  theme(
    plot.title = element_text(size = 8.8, face = "bold", colour = ink, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 6.5, colour = mid_ink, margin = margin(b = 5)),
    axis.title = element_text(size = 8, colour = ink),
    axis.text = element_text(size = 7.2, colour = ink),
    axis.ticks = element_line(colour = rule, linewidth = 0.25),
    panel.grid.major = element_line(colour = "#E7EBEE", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 7.5, colour = ink),
    legend.text = element_text(size = 7, colour = ink),
    plot.margin = margin(4, 5, 4, 5)
  )

flow_panel <- function(title, subtitle, boxes, arrows = TRUE) {
  boxes <- as.data.frame(boxes, stringsAsFactors = FALSE)
  p <- ggplot() +
    geom_rect(
      data = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#F7F9FA", colour = rule, linewidth = 0.45
    ) +
    geom_text(
      data = boxes,
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
      family = "Arial", size = 2.55, lineheight = 0.92, colour = ink
    ) +
    labs(title = title, subtitle = subtitle) +
    coord_cartesian(xlim = c(0, max(boxes$xmax) + 0.08), ylim = c(0, 1), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 8.8, face = "bold", colour = ink, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 6.5, colour = mid_ink, margin = margin(b = 5)),
      plot.margin = margin(4, 4, 2, 4)
    )
  if (arrows && nrow(boxes) > 1) {
    arrow_df <- tibble::tibble(
      x = boxes$xmax[-nrow(boxes)] + 0.035,
      xend = boxes$xmin[-1] - 0.035,
      y = rep(0.5, nrow(boxes) - 1),
      yend = rep(0.5, nrow(boxes) - 1)
    )
    p <- p + geom_segment(
      data = arrow_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.45, colour = mid_ink,
      arrow = arrow(length = unit(0.12, "cm"), type = "closed")
    )
  }
  p
}

evidence_chain_panel <- function(title, nodes) {
  nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  arrows <- tibble::tibble(
    x = nodes$x[-nrow(nodes)] + 0.08,
    xend = nodes$x[-1] - 0.08,
    y = 0.44,
    yend = 0.44
  )
  ggplot(nodes, aes(x = x, y = 0.44)) +
    geom_segment(
      data = arrows,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE, colour = rule, linewidth = 0.65,
      arrow = arrow(length = unit(0.13, "cm"), type = "closed")
    ) +
    geom_point(
      aes(fill = is_p1, size = is_p1), shape = 21, stroke = 0.8, colour = rule
    ) +
    geom_text(
      aes(y = 0.80, label = label), family = "Arial", size = 2.75,
      fontface = "bold", colour = ink, lineheight = 0.92
    ) +
    geom_text(
      aes(y = 0.14, label = detail), family = "Arial", size = 2.0,
      colour = mid_ink, lineheight = 0.92
    ) +
    scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = p1_accent), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 3.6, `TRUE` = 4.0), guide = "none") +
    labs(title = title) +
    coord_cartesian(xlim = c(0.55, max(nodes$x) + 0.45), ylim = c(0.02, 0.94), clip = "off") +
    theme_void(base_family = "Arial") +
    theme(
      plot.title = element_text(size = 8.8, face = "bold", colour = ink, margin = margin(b = 4)),
      plot.margin = margin(3, 4, 2, 4)
    )
}

pA <- evidence_chain_panel(
  "A  Discovery and P1 definition workflow",
  tibble::tibble(
    x = 1:5,
    label = c("GSE178318\ndiscovery cohort", "QC and doublet\nfiltering", "12,433 myeloid/DC\ncells", "de novo NMF,\nranks 4\u20138", "Frozen P1\nprogram"),
    detail = c("6 patients / 15 samples", "fixed QC; sample-wise\ndoublet calling", "14 annotated fine states", "30 runs per rank", "rank-6 P1 definition"),
    is_p1 = c(FALSE, FALSE, FALSE, FALSE, TRUE)
  )
)

pB_rank <- ggplot(rank_plot_df, aes(x = rank, y = value, group = 1)) +
  annotate("rect", xmin = 5.72, xmax = 6.28, ymin = -Inf, ymax = Inf, fill = "#F8E9E5", colour = NA) +
  geom_line(colour = ink, linewidth = 0.55) +
  geom_point(aes(fill = rank == 6), shape = 21, size = 2.25, stroke = 0.45, colour = ink) +
  geom_text(
    data = rank6_labels,
    aes(x = 6.32, y = label_y, label = value_label),
    inherit.aes = FALSE, hjust = 0, vjust = 0,
    colour = p1_accent, family = "Arial", fontface = "bold", size = 2.0
  ) +
  scale_fill_manual(values = c(`TRUE` = p1_accent, `FALSE` = "white"), guide = "none") +
  scale_x_continuous(limits = c(4, 8.75), breaks = 4:8, minor_breaks = NULL) +
  facet_wrap(~metric, nrow = 1, scales = "free_y", labeller = as_labeller(c(
    `Basis stability` = "Basis stability", `Cophenetic coefficient` = "Cophenetic coefficient", `Top-15 gene overlap` = "Top-15 gene overlap"
  ))) +
  labs(title = "B  NMF rank selection and stability", subtitle = "Rank 6 was selected for downstream projection", x = "NMF rank", y = NULL) +
  base_theme +
  theme(
    strip.text = element_text(size = 5.5, face = "bold", colour = ink),
    strip.background = element_blank(),
    strip.clip = "off",
    panel.spacing.x = unit(0.4, "lines"),
    axis.text.x = element_text(size = 6.2),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(colour = rule, linewidth = 0.35),
    axis.line.y = element_blank()
  )

pC <- ggplot(p1_genes, aes(x = gene_value, y = gene)) +
  geom_col(fill = p1_accent, width = 0.46) +
  scale_x_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02), labels = label_percent(accuracy = 1), expand = c(0, 0)) +
  labs(title = "C  Frozen P1 gene weights", subtitle = "Top-15 genes ranked by normalized loading weight", x = "Normalized loading weight", y = NULL) +
  base_theme +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 5.9, face = "italic", margin = margin(r = 4)),
    axis.text.x = element_text(size = 6.7)
  )

pD <- ggplot(fine_state_summary, aes(x = as.numeric(state), y = mean_P1_score, colour = state)) +
  geom_segment(
    data = state_medians,
    aes(x = xmin, xend = xmax, y = state_median, yend = state_median),
    inherit.aes = FALSE, colour = ink, linewidth = 0.38
  ) +
  geom_point(aes(size = n_cells), alpha = 0.82, position = position_jitter(width = 0.07, height = 0, seed = 11)) +
  scale_colour_manual(values = state_cols, drop = FALSE, guide = "none") +
  scale_size_continuous(
    name = "Cells", trans = "sqrt", breaks = c(500, 1500, 2500),
    labels = scales::comma_format(), range = c(2.4, 4.6)
  ) +
  scale_x_continuous(breaks = seq_along(state_levels), labels = state_levels, limits = c(0.55, 5.45)) +
  labs(
    title = "D  P1 activity across myeloid states",
    subtitle = "Fine-state estimates with broad-state medians",
    x = NULL, y = "Mean frozen P1 weighted score"
  ) +
  base_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 6.8),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = "left",
    legend.title = element_text(size = 6.2, colour = mid_ink),
    legend.text = element_text(size = 5.8, colour = mid_ink),
    legend.key.width = unit(0.34, "cm"),
    legend.key.height = unit(0.25, "cm"),
    legend.spacing.x = unit(0.08, "cm"),
    legend.margin = margin(t = -5, r = 0, b = 0, l = 0),
    legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0)
  ) +
  guides(size = guide_legend(override.aes = list(alpha = 0.82, colour = mid_ink), nrow = 1))

pE <- ggplot(benchmark_plot_df, aes(y = program, colour = is_p1)) +
  geom_vline(xintercept = 0, colour = mid_ink, linewidth = 0.45) +
  geom_segment(aes(x = ci_low, xend = ci_high, y = program, yend = program), linewidth = 1.0) +
  geom_point(aes(x = median_effect), size = 2.5) +
  geom_text(aes(y = program, label = positive_patients), x = benchmark_label_x_patients,
            inherit.aes = FALSE, colour = ink, size = 2.7, family = "Arial") +
  geom_text(aes(y = program, label = positive_pairs), x = benchmark_label_x_pairs,
            inherit.aes = FALSE, colour = ink, size = 2.7, family = "Arial") +
  annotate("text", x = 0.24, y = benchmark_header_y, label = "Effect (95% CI)", hjust = 0.5, vjust = 0.5,
           colour = mid_ink, size = 2.1, family = "Arial", fontface = "bold") +
  annotate("text", x = benchmark_label_x_patients, y = benchmark_header_y, label = "Positive patients", hjust = 0.5, vjust = 0.5,
           colour = mid_ink, size = 1.65, family = "Arial", fontface = "bold") +
  annotate("text", x = benchmark_label_x_pairs, y = benchmark_header_y, label = "Positive patient\u00d7state\npairs", hjust = 0.5, vjust = 0.5,
           lineheight = 0.85, colour = mid_ink, size = 1.55, family = "Arial", fontface = "bold") +
  scale_colour_manual(values = c(`FALSE` = "#A3ACB4", `TRUE` = p1_accent), guide = "none") +
  scale_x_continuous(
    limits = c(-0.42, 1.90), breaks = c(-0.4, 0, 0.4, 0.8),
    labels = label_number(accuracy = 0.1), expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev(benchmark_order), expand = expansion(add = c(0.35, 0.85))) +
  labs(
    title = "E  P1\u2013P6 prioritization benchmark",
    subtitle = "P1 ranks first by effect magnitude and directional consistency",
    x = "Median LM \u2212 PT effect", y = NULL
  ) +
  base_theme +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(4, 4, 4, 5)
  )

pBC <- (pB_rank | pC) + plot_layout(widths = c(0.98, 1.02))
pDE <- (pD | pE) + plot_layout(widths = c(0.98, 1.02))
fig1 <- (pA / pBC / pDE) +
  plot_layout(heights = c(0.58, 1.12, 1.35)) +
  plot_annotation(theme = theme(plot.margin = margin(4, 5, 4, 5)))

width_mm <- 183
height_mm <- 148
png_path <- file.path(out_dir, "Figure1_CRLM_P1_final.png")
svg_path <- file.path(out_dir, "Figure1_CRLM_P1_final.svg")
pdf_path <- file.path(out_dir, "Figure1_CRLM_P1_final.pdf")
tiff_path <- file.path(out_dir, "Figure1_CRLM_P1_final.tiff")

ggsave(png_path, fig1, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_png, bg = "white")
ggsave(svg_path, fig1, width = width_mm, height = height_mm, units = "mm", device = svglite::svglite, bg = "white")
ggsave(pdf_path, fig1, width = width_mm, height = height_mm, units = "mm", device = grDevices::cairo_pdf, bg = "white")
ggsave(tiff_path, fig1, width = width_mm, height = height_mm, units = "mm", dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")

qa_notes <- c(
  "Figure 1 generated from current R analysis outputs; formal result tables were not modified.",
  "Final layout: A discovery and P1 definition workflow; B NMF rank selection and stability; C frozen P1 gene weights; D P1 activity across myeloid states; E P1\u2013P6 prioritization benchmark.",
  paste0("Current discovery rows: ", current_n, "; UMAP rows: ", umap_n, "; selected rank: ", selected_rank, "."),
  paste0("Freeze-package reconciliation current_n: ", ifelse(is.na(reconciliation_current_n), "not available", reconciliation_current_n), "."),
  "Panel B rank evidence: median basis stability, cophenetic coefficient and median pairwise top-15 gene overlap across ranks 4-8.",
  paste0("P1 top-15 genes plotted from: ", genes_path, "; sorted descending by ", gene_value_col, "."),
  gene_set_warning,
  "Panel D: each point represents one annotated fine state; black bars indicate the median within each conventional analysis state; point size denotes the number of cells.",
  paste0("Panel E: frozen GSE245552 P1\u2013P6 specificity benchmark; error bars are ci_low/ci_high; positive-pair denominator = ", benchmark_denominator, "."),
  "Figure 1 source-data audit: 12,433 cells; 6 patients / 15 samples; 14 fine states; ranks 4-8 with 30 runs per rank; rank 6; P1 top-15 order from descending normalized frozen weights; P1 benchmark 0.608 / 12/13 / 45/50.",
  "Outputs: 600-dpi PNG, SVG, PDF and 600-dpi TIFF."
)
writeLines(qa_notes, file.path(out_dir, "Figure1_QA_notes.txt"), useBytes = TRUE)

metric_definitions <- c(
  "Figure 1 | Panel B metric definitions for the manuscript legend",
  "Basis stability: median pairwise Spearman correlation of the NMF basis/loadings across 30 independent Brunet runs at each rank.",
  "Cophenetic coefficient: Pearson correlation between the consensus distance (1 - consensus matrix) and its average-linkage cophenetic distance, computed on the fixed balanced 999-cell audit subset.",
  "Top-15 gene overlap: median pairwise top-15 gene-overlap fraction, defined as |A intersection B| / 15; this is not a Jaccard index.",
  "Panel D legend: Each point represents one annotated fine state, point size indicates the number of contributing cells, and black horizontal bars denote medians across the corresponding conventional analysis states.",
  "Panel E ordering: Programs were ordered by the median patient-level LM\u2212PT effect; positive-patient and positive patient\u00d7state-pair fractions are shown as independent consistency metrics and were not combined into a composite score. This benchmark was exploratory and used only for prioritization, not independent confirmatory validation.",
  "Panel E benchmark metrics: patient-level median LM\u2212PT effect; fraction of patients with positive LM\u2212PT effect; and fraction of patient\u00d7state pairs with positive LM\u2212PT effect.",
  "Panel E intervals: ci_low and ci_high from the frozen GSE245552 P1\u2013P6 specificity benchmark."
)
writeLines(metric_definitions, file.path(out_dir, "Figure1_metric_definitions_for_legend.txt"), useBytes = TRUE)

message("Figure 1 written to: ", out_dir)
message("Current discovery cells: ", current_n)
message("Selected NMF rank: ", selected_rank)
