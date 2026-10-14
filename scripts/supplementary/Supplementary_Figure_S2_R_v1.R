## Supplementary Figure S2 — NMF rank, reproducibility and robustness audits
## Plotting-only audit of frozen rank-selection and repeated-fit outputs.
## No NMF fitting, re-annotation or optimization is performed here.

options(stringsAsFactors = FALSE)
set.seed(20260914)

required_packages <- c("ggplot2", "patchwork", "ragg", "svglite", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------------------
## Paths and output contract
## ------------------------------------------------------------------------

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

frozen_root <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
frozen_results <- file.path(frozen_root, "complete_run_archive_20260912", "run_output", "results")
source_tables <- file.path(frozen_results, "tables")
cache_root <- file.path(crlm_require_dir(CRLM_PRIMARY_RERUN_ROOT, "CRLM_PRIMARY_RERUN_ROOT"), "cache")
rank6_fit_dir <- file.path(cache_root, "nmf_fits", "rank_6")
rank6_consensus_file <- file.path(cache_root, "nmf_alignment", "rank6_consensus.rds")

output_root <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S2_R_v1")
figure_dir <- file.path(output_root, "figures")
source_dir <- file.path(output_root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(filename) {
  path <- file.path(source_tables, filename)
  if (!file.exists(path)) stop("Missing frozen table: ", filename)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_source <- function(x, filename) {
  write.csv(x, file.path(source_dir, filename), row.names = FALSE, na = "")
}

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(object_name, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

## ------------------------------------------------------------------------
## Frozen rank-selection inputs
## ------------------------------------------------------------------------

rank_summary <- read_table("28_NMF_rank_summary.csv")
selected_rank <- read_table("31_NMF_selected_rank.csv")
cophenetic <- read_table("31b_NMF_cophenetic_coefficients.csv")

assert_columns(
  rank_summary,
  c("rank", "median_pairwise_basis_spearman", "median_pairwise_usage_spearman", "consensus_rss"),
  "rank_summary"
)
assert_columns(selected_rank, c("selected_rank", "rss_elbow_rank"), "selected_rank")

rank_levels <- 4:8
rank_df <- rank_summary[rank_summary$rank %in% rank_levels, , drop = FALSE]
rank_df <- rank_df[order(rank_df$rank), , drop = FALSE]
if (!identical(as.integer(rank_df$rank), rank_levels)) {
  stop("Rank-selection table does not contain exactly ranks 4–8.")
}

k_selected <- as.integer(selected_rank$selected_rank[1])
k_elbow <- as.integer(selected_rank$rss_elbow_rank[1])
if (k_selected != 6L || k_elbow != 6L) {
  stop("Frozen rank-selection result is not K=6; current values are selected=", k_selected,
       ", RSS elbow=", k_elbow)
}

rank_diag <- data.frame(
  rank = as.integer(rank_df$rank),
  basis_spearman = as.numeric(rank_df$median_pairwise_basis_spearman),
  usage_spearman = as.numeric(rank_df$median_pairwise_usage_spearman),
  consensus_rss_million = as.numeric(rank_df$consensus_rss) / 1e6,
  stringsAsFactors = FALSE
)
write_source(rank_diag, "S2A_rank_selection_diagnostics.csv")
write_source(selected_rank, "S2A_selected_rank_audit.csv")
write_source(cophenetic[cophenetic$rank %in% rank_levels, , drop = FALSE], "S2A_cophenetic_coefficients.csv")

rank_stability_long <- rbind(
  data.frame(rank = rank_diag$rank, metric = "Basis Spearman", value = rank_diag$basis_spearman),
  data.frame(rank = rank_diag$rank, metric = "Usage Spearman", value = rank_diag$usage_spearman)
)
rank_stability_long$metric <- factor(rank_stability_long$metric, levels = c("Basis Spearman", "Usage Spearman"))

COL_DARK <- "#22282E"
COL_GRID <- "#D9E0E5"
COL_BLUE <- "#3E739F"
COL_GREEN <- "#4E9A80"
COL_ACCENT <- "#B6534E"

theme_s2 <- function(base_size = 6.4) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_DARK),
      axis.title = element_text(size = 6.0, colour = COL_DARK),
      axis.text = element_text(size = 5.3, colour = COL_DARK),
      axis.line = element_line(linewidth = 0.3, colour = COL_DARK),
      axis.ticks = element_line(linewidth = 0.3, colour = COL_DARK),
      panel.grid.major.y = element_line(linewidth = 0.22, colour = COL_GRID),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 5.7, colour = COL_DARK),
      legend.text = element_text(size = 5.2, colour = COL_DARK),
      legend.key.height = grid::unit(2.7, "mm"),
      legend.key.width = grid::unit(3.6, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(1.0, "mm"),
      plot.title = element_text(size = 7.4, face = "bold", hjust = 0, colour = COL_DARK,
                                margin = margin(1.2, 0, 2.0, 0)),
      plot.subtitle = element_text(size = 5.5, colour = "#59636C", hjust = 0),
      plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 5.6, face = "bold", colour = COL_DARK)
    )
}

pA_stability <- ggplot(rank_stability_long, aes(rank, value, colour = metric, group = metric)) +
  geom_vline(xintercept = k_selected, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.45) +
  scale_colour_manual(values = c("Basis Spearman" = COL_BLUE, "Usage Spearman" = COL_GREEN)) +
  scale_x_continuous(breaks = rank_levels, limits = c(3.8, 8.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.86, 1.005), breaks = c(0.88, 0.92, 0.96, 1.00), expand = c(0.01, 0)) +
  labs(title = "A  Rank-selection diagnostics", x = "NMF rank", y = "Median pairwise Spearman", colour = NULL) +
  theme_s2() +
  theme(legend.position = "bottom", legend.direction = "horizontal")

pA_rss <- ggplot(rank_diag, aes(rank, consensus_rss_million)) +
  geom_vline(xintercept = k_selected, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  geom_line(linewidth = 0.55, colour = COL_ACCENT) +
  geom_point(size = 1.45, colour = COL_ACCENT) +
  scale_x_continuous(breaks = rank_levels, limits = c(3.8, 8.2), expand = c(0, 0)) +
  scale_y_continuous(labels = function(x) sprintf("%.2f", x), expand = expansion(mult = c(0.03, 0.08))) +
  labs(x = "NMF rank", y = "Consensus RSS (×10^6)") +
  theme_s2() +
  theme(legend.position = "none")

pA <- pA_stability | pA_rss

## ------------------------------------------------------------------------
## Hungarian-aligned factor reproducibility across 30 rank-6 fits
## ------------------------------------------------------------------------

alignment <- read_table("25_NMF_factor_alignment.csv")
assert_columns(alignment, c("rank", "run", "factor", "matched_factor_in_run", "loading_spearman"), "alignment")
alignment_k6 <- alignment[alignment$rank == k_selected, , drop = FALSE]
alignment_k6 <- alignment_k6[order(alignment_k6$factor, alignment_k6$run), , drop = FALSE]
if (nrow(alignment_k6) != 180L || length(unique(alignment_k6$run)) != 30L) {
  stop("Rank-6 alignment table is not a complete 6 × 30 audit.")
}

consensus <- readRDS(rank6_consensus_file)
if (is.null(consensus$W) || is.null(consensus$H)) stop("Frozen rank-6 consensus object lacks W/H matrices.")

all_program_top15 <- read_table("35_frozen_all_program_top15_weights.csv")
assert_columns(all_program_top15, c("program", "rank", "gene"), "all_program_top15")
program_levels <- paste0("P", 1:6)
if (!setequal(unique(all_program_top15$program), program_levels)) {
  stop("Frozen program table does not contain exactly P1–P6.")
}

consensus_W <- consensus$W
consensus_genes <- rownames(consensus_W)
if (is.null(consensus_genes)) consensus_genes <- as.character(consensus$genes)
if (length(consensus_genes) != nrow(consensus_W)) stop("Consensus W gene names are unavailable.")

factor_program <- setNames(rep(NA_character_, 6), as.character(1:6))
factor_program_overlap <- numeric(6)
program_genes <- split(all_program_top15$gene, all_program_top15$program)
for (k in 1:6) {
  top_genes <- consensus_genes[order(consensus_W[, k], decreasing = TRUE, na.last = NA)[seq_len(15)]]
  overlaps <- vapply(program_genes[program_levels], function(g) sum(top_genes %in% g), numeric(1))
  best <- which.max(overlaps)
  if (overlaps[best] != 15) {
    stop("Consensus factor ", k, " does not map uniquely to a frozen program at top-15.")
  }
  factor_program[k] <- program_levels[best]
  factor_program_overlap[k] <- overlaps[best]
}
if (anyDuplicated(factor_program)) stop("Consensus-factor to frozen-program mapping is not one-to-one.")

factor_map_table <- data.frame(
  factor = 1:6,
  program = unname(factor_program),
  top15_overlap = factor_program_overlap,
  stringsAsFactors = FALSE
)
write_source(factor_map_table, "S2B_consensus_factor_to_program_map.csv")

alignment_plot <- alignment_k6
alignment_plot$program <- factor(factor_program[as.character(alignment_plot$factor)], levels = rev(program_levels))
alignment_plot$run_label <- sprintf("run%02d", alignment_plot$run)
write_source(alignment_plot, "S2B_hungarian_aligned_factor_reproducibility.csv")

pB <- ggplot(alignment_plot, aes(x = run, y = program, fill = loading_spearman)) +
  geom_tile(colour = "white", linewidth = 0.16) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30),
                     labels = c("1", "5", "10", "15", "20", "25", "30"),
                     expand = c(0, 0)) +
  scale_fill_viridis_c(option = "C", limits = c(0.68, 1.00),
                       breaks = c(0.7, 0.8, 0.9, 1.0),
                       name = "loading\nSpearman", oob = scales::squish,
                       guide = guide_colorbar(barheight = grid::unit(16, "mm"),
                                              barwidth = grid::unit(2.5, "mm"))) +
  labs(title = "B  Hungarian-aligned factor reproducibility",
       x = "Repeated rank-6 fit", y = "Frozen program") +
  theme_s2() +
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(), legend.position = "right")

## ------------------------------------------------------------------------
## P1-specific reproducibility across all 30 runs
## All three metrics use the frozen consensus P1 (factor 2) after the
## frozen Hungarian matching. The legacy marker correspondence table is
## retained only as a comparison audit and is not used for the figure.
## ------------------------------------------------------------------------

mapping_p1 <- alignment_k6[alignment_k6$factor == 2, c("run", "matched_factor_in_run")]
mapping_p1 <- mapping_p1[order(mapping_p1$run), , drop = FALSE]
if (nrow(mapping_p1) != 30L || !identical(as.integer(mapping_p1$run), 1:30)) {
  stop("Frozen Hungarian alignment does not contain a complete P1 mapping for runs 1–30.")
}

legacy_p1_audit <- read_table("43_selected_rank_P1_30run_audit.csv")
assert_columns(legacy_p1_audit, c("run", "matched_factor"), "legacy_p1_audit")
legacy_comparison <- merge(
  mapping_p1,
  legacy_p1_audit[, c("run", "matched_factor")],
  by = "run", suffixes = c("_hungarian", "_legacy_marker"), sort = TRUE
)
legacy_comparison$mapping_agrees <- legacy_comparison$matched_factor_in_run == legacy_comparison$matched_factor
write_source(legacy_comparison, "S2C_hungarian_vs_legacy_P1_mapping_audit.csv")

get_matrix_names <- function(x, component, fallback_name) {
  mat <- x[[component]]
  nm <- if (component == "W") rownames(mat) else colnames(mat)
  if (is.null(nm)) nm <- as.character(x[[fallback_name]])
  if (is.null(nm) || length(nm) != if (component == "W") nrow(mat) else ncol(mat)) {
    stop("Missing names for ", component, " in frozen fit object.")
  }
  as.character(nm)
}

consensus_cells <- get_matrix_names(consensus, "H", "cells")
consensus_genes <- get_matrix_names(consensus, "W", "genes")
consensus_h_cells <- colnames(consensus$H)
if (is.null(consensus_h_cells)) consensus_h_cells <- consensus_cells
consensus_p1_usage <- as.numeric(consensus$H[2, match(consensus_cells, consensus_h_cells)])
consensus_p1_loading <- as.numeric(consensus$W[, 2])
p1_genes <- as.character(read_table("34_frozen_P1_top15_weights.csv")$gene)
consensus_p1_top15 <- consensus_genes[order(consensus_p1_loading, decreasing = TRUE, na.last = NA)[seq_len(15)]]
if (!setequal(p1_genes, consensus_p1_top15)) {
  stop("Frozen P1 top-15 table does not match consensus factor 2 top-15 genes.")
}

p1_usage <- data.frame(
  run = integer(0),
  matched_factor = integer(0),
  usage_rho_vs_final_consensus_P1 = numeric(0),
  loading_rho_vs_final_consensus_P1 = numeric(0),
  top15_overlap_vs_final = numeric(0),
  stringsAsFactors = FALSE
)
p1_run_ranks <- vector("list", 30)
p1_gene_rank <- read_table("34_frozen_P1_top15_weights.csv")[, c("rank", "gene"), drop = FALSE]

for (rr in 1:30) {
  fit_file <- file.path(rank6_fit_dir, sprintf("rank6_run%02d.rds", rr))
  if (!file.exists(fit_file)) stop("Missing frozen rank-6 fit: run ", rr)
  fit <- readRDS(fit_file)
  fit_cells <- get_matrix_names(fit, "H", "cells")
  fit_genes <- get_matrix_names(fit, "W", "genes")
  if (!setequal(fit_cells, consensus_cells)) stop("Cell names differ in rank-6 run ", rr)
  if (!setequal(fit_genes, consensus_genes)) stop("Gene names differ in rank-6 run ", rr)

  mf <- mapping_p1$matched_factor_in_run[mapping_p1$run == rr]
  fit_h <- as.numeric(fit$H[mf, match(consensus_cells, fit_cells)])
  usage_rho <- suppressWarnings(cor(fit_h, consensus_p1_usage, method = "spearman", use = "pairwise.complete.obs"))
  fit_w <- as.numeric(fit$W[match(consensus_genes, fit_genes), mf])
  loading_rho <- suppressWarnings(cor(fit_w, consensus_p1_loading, method = "spearman", use = "pairwise.complete.obs"))
  top15_overlap <- length(intersect(
    consensus_genes[order(fit_w, decreasing = TRUE, na.last = NA)[seq_len(15)]],
    consensus_p1_top15
  )) / 15
  p1_usage <- rbind(p1_usage, data.frame(
    run = rr,
    matched_factor = mf,
    usage_rho_vs_final_consensus_P1 = usage_rho,
    loading_rho_vs_final_consensus_P1 = loading_rho,
    top15_overlap_vs_final = top15_overlap,
    stringsAsFactors = FALSE
  ))

  top_idx <- order(fit_w, decreasing = TRUE, na.last = NA)[seq_len(15)]
  top_run_genes <- consensus_genes[top_idx]
  rank_lookup <- setNames(seq_along(top_run_genes), top_run_genes)
  p1_run_ranks[[rr]] <- data.frame(
    run = rr,
    run_label = sprintf("run%02d", rr),
    gene = p1_genes,
    frozen_p1_rank = p1_gene_rank$rank[match(p1_genes, p1_gene_rank$gene)],
    run_top15_rank = unname(rank_lookup[p1_genes]),
    present_in_run_top15 = !is.na(unname(rank_lookup[p1_genes])),
    stringsAsFactors = FALSE
  )
}

p1_c <- p1_usage
p1_c <- p1_c[order(p1_c$run), , drop = FALSE]
p1_c$run_label <- sprintf("run%02d", p1_c$run)
write_source(p1_c, "S2C_P1_30run_reproducibility.csv")

p1_c_long <- rbind(
  data.frame(run = p1_c$run, metric = "Loading rho", value = p1_c$loading_rho_vs_final_consensus_P1),
  data.frame(run = p1_c$run, metric = "Usage rho", value = p1_c$usage_rho_vs_final_consensus_P1),
  data.frame(run = p1_c$run, metric = "Top-15 overlap", value = p1_c$top15_overlap_vs_final)
)
p1_metric_levels <- c("Loading rho", "Usage rho", "Top-15 overlap")
p1_c_long$metric <- factor(p1_c_long$metric, levels = p1_metric_levels)
p1_c_medians <- aggregate(value ~ metric, p1_c_long, median, na.rm = TRUE)

pC <- ggplot(p1_c_long, aes(x = run, y = value)) +
  geom_hline(data = p1_c_medians, aes(yintercept = value), linewidth = 0.45,
             linetype = 2, colour = COL_ACCENT, inherit.aes = FALSE) +
  geom_point(size = 1.25, colour = COL_BLUE, alpha = 0.90) +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_continuous(breaks = c(1, 10, 20, 30), limits = c(0.2, 30.8), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.55, 1.01), breaks = c(0.6, 0.8, 1.0), expand = c(0.01, 0)) +
  labs(title = "C  P1 reproducibility across 30 rank-6 fits",
       subtitle = "Hungarian-aligned to frozen consensus P1; dashed line: median",
       x = "Run", y = "Agreement with frozen P1") +
  theme_s2() +
  theme(panel.grid = element_blank(), strip.placement = "outside", strip.text = element_text(size = 5.7, face = "bold"))

## ------------------------------------------------------------------------
## Increased optimization-depth sensitivity
## The frozen audit contains three selected repeats (runs 1, 15 and 30).
## ------------------------------------------------------------------------

depth_summary <- read_table("NMF_iteration_depth_sensitivity_rank6_2000/rank6_500_1000_2000_summary.csv")
depth_run <- read_table("NMF_iteration_depth_sensitivity_rank6_2000/rank6_500_1000_2000_run_level.csv")
assert_columns(depth_summary, c("rank", "n_fits", "median_loading_spearman_1000_vs_2000", "median_usage_spearman_1000_vs_2000", "median_top15_overlap_1000_vs_2000"), "depth_summary")
assert_columns(depth_run, c("rank", "run", "p1_factor_in_2000", "loading_spearman_1000_vs_2000_median", "usage_spearman_1000_vs_2000_median", "top15_overlap_1000_vs_2000_median"), "depth_run")
depth_run <- depth_run[depth_run$rank == k_selected, , drop = FALSE]
if (nrow(depth_run) != 3L) stop("Expected three frozen rank-6 depth-sensitivity repeats.")

checkpoint_expected_p1 <- mapping_p1$matched_factor_in_run[match(depth_run$run, mapping_p1$run)]
checkpoint_identity <- data.frame(
  run = as.integer(depth_run$run),
  p1_factor_in_2000 = as.integer(depth_run$p1_factor_in_2000),
  hungarian_p1_factor = as.integer(checkpoint_expected_p1),
  identity_agrees = as.integer(depth_run$p1_factor_in_2000) == as.integer(checkpoint_expected_p1),
  stringsAsFactors = FALSE
)
if (anyNA(checkpoint_expected_p1) || !all(checkpoint_identity$identity_agrees)) {
  stop("2000-iteration checkpoint P1 identity does not agree with frozen Hungarian mapping.")
}
write_source(checkpoint_identity, "S2D_checkpoint_P1_identity_audit.csv")

depth_plot <- data.frame(
  run = depth_run$run,
  loading_rho = depth_run$loading_spearman_1000_vs_2000_median,
  usage_rho = depth_run$usage_spearman_1000_vs_2000_median,
  top15_overlap = depth_run$top15_overlap_1000_vs_2000_median,
  stringsAsFactors = FALSE
)
write_source(depth_summary, "S2D_optimization_depth_summary.csv")
write_source(depth_run, "S2D_optimization_depth_run_level.csv")

depth_long <- rbind(
  data.frame(run = depth_plot$run, metric = "Loading rho", value = depth_plot$loading_rho),
  data.frame(run = depth_plot$run, metric = "Usage rho", value = depth_plot$usage_rho),
  data.frame(run = depth_plot$run, metric = "Top-15 overlap", value = depth_plot$top15_overlap)
)
depth_long$metric <- factor(depth_long$metric, levels = p1_metric_levels)

pD <- ggplot(depth_long, aes(x = 1, y = value)) +
  geom_jitter(width = 0.055, height = 0, size = 1.75, colour = COL_BLUE, alpha = 0.90) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6.0,
               colour = COL_ACCENT) +
  facet_wrap(~ metric, nrow = 1) +
  scale_x_continuous(breaks = NULL, limits = c(0.78, 1.22), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.975, 1.002), breaks = c(0.98, 0.99, 1.00), expand = c(0.01, 0)) +
  labs(title = "D  Robustness to increased optimization depth",
       subtitle = "Primary 500; sensitivity 1,000/2,000 iterations; n = 3",
       x = NULL, y = "Agreement") +
  theme_s2() +
  theme(panel.grid = element_blank(), strip.placement = "outside")

## ------------------------------------------------------------------------
## Stability of leading P1 genes across 30 fits
## Heatmap fill is the run-specific top-15 rank for the frozen P1 gene set;
## blank cells are genes not present in that run's top 15.
## ------------------------------------------------------------------------

p1_gene_runs <- do.call(rbind, p1_run_ranks)
p1_gene_runs$gene <- factor(p1_gene_runs$gene, levels = rev(p1_genes))
p1_frequency <- aggregate(present_in_run_top15 ~ gene, p1_gene_runs, sum)
names(p1_frequency)[names(p1_frequency) == "present_in_run_top15"] <- "n_runs_top15"
p1_frequency$gene <- factor(as.character(p1_frequency$gene), levels = rev(p1_genes))
p1_frequency <- p1_frequency[order(p1_frequency$gene), , drop = FALSE]

write_source(p1_gene_runs, "S2E_P1_top15_gene_rank_by_run.csv")
write_source(p1_frequency, "S2E_P1_top15_gene_membership_frequency.csv")

pE_heat <- ggplot(p1_gene_runs, aes(x = run, y = gene, fill = run_top15_rank)) +
  geom_tile(colour = "white", linewidth = 0.13) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30),
                     labels = c("1", "5", "10", "15", "20", "25", "30"),
                     expand = c(0, 0)) +
  scale_fill_gradient(low = "#3C78A8", high = "#D8E4EC", na.value = "#F2F4F6",
                      limits = c(1, 15), breaks = c(1, 5, 10, 15),
                      name = "run top-15\nrank",
                      guide = guide_colorbar(barheight = grid::unit(14, "mm"),
                                             barwidth = grid::unit(2.3, "mm"))) +
  labs(x = "Repeated rank-6 fit", y = "Frozen P1 top genes") +
  theme_s2() +
  theme(panel.grid = element_blank(), axis.ticks.y = element_blank(), legend.position = "right")

pE_freq <- ggplot(p1_frequency, aes(x = n_runs_top15, y = gene)) +
  geom_col(width = 0.72, fill = COL_ACCENT) +
  geom_text(aes(x = pmax(n_runs_top15 - 0.65, 1), label = paste0(n_runs_top15, "/30")),
            hjust = 1, size = 1.75, colour = "white", fontface = "bold") +
  scale_x_continuous(limits = c(0, 35), breaks = c(0, 15, 30), expand = c(0, 0)) +
  labs(x = "Membership", y = NULL) +
  theme_s2() +
  theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), legend.position = "none")

pE_heat <- pE_heat + labs(title = "E  Stability of leading P1 genes")
pE <- (pE_heat | patchwork::plot_spacer() | pE_freq) +
  plot_layout(widths = c(1.42, 0.12, 0.74))

## ------------------------------------------------------------------------
## Assemble and export
## ------------------------------------------------------------------------

figure_s2 <- ((pA | pB) / pC / (pD | pE)) +
  plot_layout(heights = c(1.08, 0.92, 1.10), widths = c(1, 1))

figure_base <- file.path(figure_dir, "Supplementary_Figure_S2_R_v1")
figure_width_mm <- 183
figure_height_mm <- 190
figure_dpi <- 600

ragg::agg_png(paste0(figure_base, ".png"), width = figure_width_mm, height = figure_height_mm,
              units = "mm", res = figure_dpi, background = "white")
print(figure_s2)
dev.off()

svglite::svglite(paste0(figure_base, ".svg"), width = figure_width_mm / 25.4,
                height = figure_height_mm / 25.4, bg = "white")
print(figure_s2)
dev.off()

grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = figure_width_mm / 25.4,
                     height = figure_height_mm / 25.4, family = "Arial")
print(figure_s2)
dev.off()

ragg::agg_tiff(paste0(figure_base, ".tiff"), width = figure_width_mm, height = figure_height_mm,
               units = "mm", res = figure_dpi, compression = "lzw", background = "white")
print(figure_s2)
dev.off()

figure_manifest <- data.frame(
  panel = c("A", "B", "C", "D", "E"),
  evidence = c(
    "Ranks 4–8: pairwise stability and consensus RSS; frozen K=6 selection; cophenetic maximum also at K=6",
    "Six frozen consensus factors aligned to 30 rank-6 runs by Hungarian matching",
    "P1 loading rho, usage rho and top-15 overlap across all 30 rank-6 runs",
    "Rank-6 1,000-versus-2,000 iteration sensitivity across three selected repeats",
    "Run-specific top-15 rank of the frozen P1 gene set across 30 rank-6 fits"
  ),
  n = c("5 ranks", "6 × 30", "30 runs/metric", "3 repeats/metric", "15 genes × 30 runs"),
  source_audit = c(
    "Frozen rank summary and selected-rank audit",
    "Frozen alignment table; factor-to-program map validated at top-15 overlap 15/15",
    "Frozen Hungarian alignment plus loading/usage from frozen W/H matrices",
    "Frozen iteration-depth summary and run-level audit",
    "Frozen P1 top-15 table plus frozen rank-6 W matrices"
  ),
  stringsAsFactors = FALSE
)
write_source(figure_manifest, "S2_figure_manifest.csv")

writeLines(c(
  "Supplementary Figure S2 source-data note",
  "",
  "All plotted values were read from frozen rank-selection, alignment, audit and iteration-depth outputs.",
  "No NMF factorization, annotation, re-ranking of K, or optimization was performed by this plotting script.",
  "Panels B, C and E use the frozen consensus factor 2 (P1) and the frozen Hungarian alignment across all 30 rank-6 runs.",
  "The legacy marker-correspondence audit is exported for comparison only; it is not used to define the P1 run-level factor.",
  "Panel C displays all 30 rank-6 runs; loading and usage rho were calculated from the frozen run-level W/H matrices after Hungarian matching.",
  "K=6 also had the highest cophenetic coefficient among ranks 4–8; this statistic is retained in the source data but not added as a third plot in panel A.",
  "Panel D is a small depth-sensitivity audit with three selected repeats (runs 1, 15 and 30), and is not presented as a 30-run estimate.",
  "Panel E reports run-specific top-15 membership of the frozen P1 gene set; blank cells indicate absence from that run's top 15."
), file.path(output_root, "S2_source_data_note.txt"))

message("S2 R build complete: ", normalizePath(output_root, winslash = "/", mustWork = FALSE))
