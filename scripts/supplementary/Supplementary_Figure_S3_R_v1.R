## Supplementary Figure S3 — P1 relationship to literature-derived SPP1-TAM markers
## Plotting-only display of frozen comparator, external-state and residual audits.
## No new model fitting, hypothesis testing or score construction is performed.

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

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

main_tables <- file.path(crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR"), "tables")
extension_tables <- file.path(crlm_require_dir(CRLM_SUPPLEMENTARY_RESULTS_ROOT, "CRLM_SUPPLEMENTARY_RESULTS_ROOT"), "tables")

output_root <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S3_R_v1")
figure_dir <- file.path(output_root, "figures")
source_dir <- file.path(output_root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(directory, filename) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) stop("Missing frozen table: ", filename)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

write_source <- function(x, filename) {
  write.csv(x, file.path(source_dir, filename), row.names = FALSE, na = "")
}

write_source_gz <- function(x, filename) {
  con <- gzfile(file.path(source_dir, filename), open = "wt")
  write.csv(x, con, row.names = FALSE, na = "")
  close(con)
}

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(object_name, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

assert_near <- function(observed, expected, tolerance = 1e-6, label = "") {
  if (length(observed) != 1L || is.na(observed) ||
      abs(as.numeric(observed) - as.numeric(expected)) > tolerance) {
    stop(
      "Frozen value mismatch: ", label,
      ". observed=", observed,
      ", expected=", expected
    )
  }
}

COL_DARK <- "#22282E"
COL_BLUE <- "#3E739F"
COL_BLUE_GREY <- "#6F8796"
COL_GREY <- "#9EA5AA"
COL_GREEN <- "#4E9A80"
COL_ACCENT <- "#B6534E"
COL_GRID <- "#D9E0E5"

theme_s3 <- function(base_size = 6.4) {
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

## ------------------------------------------------------------------------
## Panel A — comparator landscape in discovery cells
## ------------------------------------------------------------------------

correlations <- read_table(extension_tables, "published_signature_discovery_cell_correlations.csv")
assert_columns(correlations, c("signature", "n_cells", "spearman_rho"), "correlations")

signature_order <- c(
  "SPP1_TAM_published_2025",
  "C1QC_TAM_published_2025",
  "SPP1_single_gene",
  "APOC1_single_gene",
  "P1_no_SPP1",
  "P1_no_SPP1_APOC1"
)
signature_labels <- c(
  SPP1_TAM_published_2025 = "Literature-derived SPP1-TAM marker score",
  C1QC_TAM_published_2025 = "C1QC-TAM",
  SPP1_single_gene = "SPP1",
  APOC1_single_gene = "APOC1",
  P1_no_SPP1 = "P1 excluding SPP1",
  P1_no_SPP1_APOC1 = "P1 excluding SPP1/APOC1"
)
correlations <- correlations[match(signature_order, correlations$signature), , drop = FALSE]
if (anyNA(correlations$signature)) stop("Panel A is missing one or more frozen comparator signatures.")

a_df <- data.frame(
  signature = signature_order,
  label = unname(signature_labels[signature_order]),
  rho = as.numeric(correlations$spearman_rho),
  n_cells = as.integer(correlations$n_cells),
  group = c(rep("Literature-derived comparator", 2), rep("Single-gene comparators", 2), rep("Marker-exclusion sensitivity", 2)),
  stringsAsFactors = FALSE
)

## Hard provenance check against the frozen discovery comparator output.
## This guards against accidentally substituting a historical comparator table.
a_expected_rho <- c(
  SPP1_TAM_published_2025 = 0.211,
  C1QC_TAM_published_2025 = 0.544,
  SPP1_single_gene = 0.623,
  APOC1_single_gene = 0.752,
  P1_no_SPP1 = 0.818,
  P1_no_SPP1_APOC1 = 0.796
)

for (ii in seq_along(signature_order)) {
  assert_near(
    a_df$rho[ii],
    a_expected_rho[[signature_order[ii]]],
    tolerance = 0.005,
    label = paste0("Panel A ", signature_order[ii], " Spearman rho")
  )
}

if (any(a_df$n_cells != 12433L)) {
  stop("Panel A comparator correlations must be based on 12,433 discovery cells.")
}

a_df$label <- factor(a_df$label, levels = rev(a_df$label))
a_df$group <- factor(a_df$group, levels = c("Literature-derived comparator", "Single-gene comparators", "Marker-exclusion sensitivity"))
write_source(a_df, "S3A_comparator_relationships.csv")

a_cols <- c(
  "Literature-derived comparator" = COL_BLUE_GREY,
  "Single-gene comparators" = COL_GREY,
  "Marker-exclusion sensitivity" = COL_GREEN
)

pA <- ggplot(a_df, aes(y = label, colour = group)) +
  geom_hline(yintercept = 2.5, linewidth = 0.3, linetype = 2, colour = "#AEB7BD") +
  geom_segment(aes(x = 0, xend = rho, yend = label), linewidth = 0.8) +
  geom_point(aes(x = rho), size = 1.8) +
  geom_text(aes(x = rho + 0.035, label = sprintf("rho == %.3f", rho)),
            hjust = 0, size = 2.05, colour = COL_DARK,
            show.legend = FALSE, parse = TRUE) +
  scale_colour_manual(values = a_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1.04), breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     expand = c(0, 0)) +
  labs(title = "A  Relationship of P1 to macrophage comparators",
       subtitle = "Discovery-cell descriptive Spearman associations; colours indicate comparator class",
       x = expression("Spearman " * rho * " with full P1"), y = NULL) +
  theme_s3() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_cartesian(clip = "off")

## ------------------------------------------------------------------------
## Panel B — state-stratified external cell-level relationship
## ------------------------------------------------------------------------

b_df <- read_table(main_tables, "GSE245552_cell_level_P1_vs_published_SPP1_TAM.csv.gz")
assert_columns(b_df, c("state", "P1_full", "SPP1_TAM_published_2025"), "external cell-level scores")

state_order <- c("Mono", "Other", "P1_state", "P2_inflammatory", "P3_C1Q_MHCII")
state_labels <- c(
  Mono = "Mono",
  Other = "Other",
  P1_state = "P1-state",
  P2_inflammatory = "Inflammatory",
  P3_C1Q_MHCII = "C1Q/MHC-II"
)
state_labels_b <- c(
  Mono = "Mono",
  Other = "Other",
  P1_state = "P1-state",
  P2_inflammatory = "Inflammatory",
  P3_C1Q_MHCII = "C1Q/\nMHC-II"
)
b_df <- b_df[b_df$state %in% state_order & is.finite(b_df$P1_full) & is.finite(b_df$SPP1_TAM_published_2025), , drop = FALSE]
b_df$state_label <- factor(unname(state_labels_b[b_df$state]), levels = unname(state_labels_b[state_order]))
if (nrow(b_df) != 24231L || !setequal(unique(b_df$state), state_order)) {
  stop("Panel B does not contain the expected 24,231 cells and five frozen states.")
}
write_source_gz(b_df, "S3B_external_cell_level_P1_vs_SPP1_TAM.csv.gz")

pB <- ggplot(b_df, aes(x = SPP1_TAM_published_2025, y = P1_full)) +
  geom_point(size = 0.14, alpha = 0.16, colour = COL_BLUE, stroke = 0) +
  facet_wrap(~ state_label, nrow = 1, scales = "fixed") +
  labs(title = expression("B  P1" - "SPP1-TAM relationship across myeloid states"),
       subtitle = "External GSE245552 cell-level scores stratified by frozen myeloid state",
       x = "Literature-derived SPP1-TAM marker score", y = "Full P1 score") +
  theme_s3() +
  theme(panel.grid = element_blank(), strip.placement = "outside",
        strip.text = element_text(size = 5.3, face = "bold"))

## ------------------------------------------------------------------------
## Corrected state-only residual audit used by Panels C and D
## ------------------------------------------------------------------------

residual_table <- read_table(main_tables, "GSE245552_patient_tissue_state_P1_residualized_SPP1_TAM_corrected_state_only.csv")
assert_columns(residual_table, c("patient", "tissue", "state", "P1_residual_state_only"), "residual_table")
residual_table <- residual_table[residual_table$tissue %in% c("primary_tumor", "liver_metastasis") &
                                   residual_table$state %in% state_order &
                                   is.finite(residual_table$P1_residual_state_only), , drop = FALSE]

pt <- residual_table[residual_table$tissue == "primary_tumor", c("patient", "state", "P1_residual_state_only"), drop = FALSE]
lm <- residual_table[residual_table$tissue == "liver_metastasis", c("patient", "state", "P1_residual_state_only"), drop = FALSE]
pair_df <- merge(pt, lm, by = c("patient", "state"), suffixes = c("_PT", "_LM"), sort = FALSE)
pair_df$delta_LM_minus_PT <- pair_df$P1_residual_state_only_LM - pair_df$P1_residual_state_only_PT
pair_df$state_label <- factor(unname(state_labels_b[pair_df$state]), levels = unname(state_labels_b[state_order]))
pair_df <- pair_df[order(pair_df$patient, pair_df$state_label), , drop = FALSE]
if (nrow(pair_df) != 57L || length(unique(pair_df$patient)) != 13L) {
  stop("Corrected residual table does not contain the expected 57 patient×state pairs and 13 patients.")
}
write_source(pair_df, "S3D_patient_state_residual_LM_minus_PT.csv")

## Authoritative frozen patient-level rule: mean across the available
## state-specific LM-PT deltas for each patient. Do not replace with a median
## or a newly chosen weighting rule.
patient_df <- aggregate(delta_LM_minus_PT ~ patient, pair_df, mean)
names(patient_df)[names(patient_df) == "delta_LM_minus_PT"] <- "patient_delta_LM_minus_PT"
patient_df <- patient_df[order(patient_df$patient_delta_LM_minus_PT), , drop = FALSE]
patient_df$patient <- factor(patient_df$patient, levels = patient_df$patient)

residual_summary <- read_table(main_tables, "GSE245552_patient_level_SPP1_residual_comparison_corrected_state_only.csv")
assert_columns(residual_summary, c("residual_definition", "n_pairs", "n_patients", "median_delta", "ci_low", "ci_high", "n_positive", "signflip_p_two_sided"), "residual_summary")
summary_row <- residual_summary[grepl("P1_residual_state_only", residual_summary$residual_definition, fixed = TRUE), , drop = FALSE]
if (nrow(summary_row) != 1L) stop("Corrected state-only patient summary is not uniquely identified.")
summary_row$residual_definition <- sub(
  "Literature-derived SPP1-TAM marker score",
  "literature-derived SPP1-TAM marker score",
  summary_row$residual_definition,
  fixed = TRUE
)
if (abs(median(patient_df$patient_delta_LM_minus_PT) - summary_row$median_delta) > 1e-8 ||
    sum(patient_df$patient_delta_LM_minus_PT > 0) != summary_row$n_positive) {
  stop("Patient-level residual deltas do not reproduce the frozen corrected summary.")
}

if (as.integer(summary_row$n_pairs) != 57L ||
    as.integer(summary_row$n_patients) != 13L) {
  stop("Frozen corrected summary does not report the expected 57 pairs and 13 patients.")
}

assert_near(summary_row$median_delta, 0.42987174, tolerance = 1e-6,
            label = "Panel C frozen median delta")
assert_near(summary_row$ci_low, 0.2490485, tolerance = 1e-6,
            label = "Panel C frozen bootstrap CI lower")
assert_near(summary_row$ci_high, 0.5328399, tolerance = 1e-6,
            label = "Panel C frozen bootstrap CI upper")
assert_near(summary_row$signflip_p_two_sided, 0.01659834, tolerance = 1e-6,
            label = "Panel C frozen sign-flip P")
write_source(patient_df, "S3C_patient_residual_LM_minus_PT.csv")
write_source(summary_row, "S3C_patient_residual_summary_corrected_state_only.csv")

patient_levels <- c(as.character(patient_df$patient), "Overall median")
patient_df$patient_plot <- factor(as.character(patient_df$patient), levels = patient_levels)
summary_patient <- factor("Overall median", levels = patient_levels)
summary_df <- data.frame(
  patient_plot = summary_patient,
  median_delta = as.numeric(summary_row$median_delta),
  ci_low = as.numeric(summary_row$ci_low),
  ci_high = as.numeric(summary_row$ci_high),
  n_positive = as.integer(summary_row$n_positive),
  signflip_p = as.numeric(summary_row$signflip_p_two_sided),
  stringsAsFactors = FALSE
)
C_summary_label <- bquote(
  "Median " * Delta * " = " * .(sprintf("%.3f", as.numeric(summary_row$median_delta))) *
    " [" * .(sprintf("%.3f", as.numeric(summary_row$ci_low))) *
    ", " * .(sprintf("%.3f", as.numeric(summary_row$ci_high))) *
    "]; " * .(as.integer(summary_row$n_positive)) * "/" *
    .(as.integer(summary_row$n_patients)) *
    " positive; sign-flip P = " *
    .(sprintf("%.4f", as.numeric(summary_row$signflip_p_two_sided)))
)

pC <- ggplot(patient_df, aes(x = patient_delta_LM_minus_PT, y = patient_plot)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  geom_hline(yintercept = length(patient_levels) - 0.5,
             linewidth = 0.28, linetype = 3, colour = "#AEB7BD") +
  geom_point(size = 1.65, colour = COL_BLUE) +
  geom_segment(data = summary_df,
               aes(x = ci_low, xend = ci_high, y = patient_plot, yend = patient_plot),
               inherit.aes = FALSE, linewidth = 0.85, colour = COL_ACCENT) +
  geom_point(data = summary_df,
             aes(x = median_delta, y = patient_plot), inherit.aes = FALSE,
             shape = 18, size = 2.6, colour = COL_ACCENT) +
  scale_x_continuous(limits = c(min(c(patient_df$patient_delta_LM_minus_PT, 0), na.rm = TRUE) - 0.10,
                                max(patient_df$patient_delta_LM_minus_PT, na.rm = TRUE) + 0.12),
                     breaks = scales::breaks_pretty(n = 6), expand = c(0, 0)) +
  labs(title = "C  CRLM-associated residual P1 elevation across patients",
       subtitle = bquote(
         atop(
           "Adjusted for literature-derived SPP1-TAM marker score and conventional myeloid state",
           .(C_summary_label)
         )
       ),
       x = expression("Residual LM " - " PT difference in P1"), y = NULL) +
  theme_s3() +
  theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_cartesian(clip = "off")

## ------------------------------------------------------------------------
## Panel D — state-stratified residual LM−PT deltas
## ------------------------------------------------------------------------

state_medians <- aggregate(delta_LM_minus_PT ~ state_label, pair_df, median)
state_medians$x <- match(state_medians$state_label, levels(pair_df$state_label))
state_medians$state_label <- factor(state_medians$state_label, levels = levels(pair_df$state_label))

pD <- ggplot(pair_df, aes(x = state_label, y = delta_LM_minus_PT)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  geom_point(position = position_jitter(width = 0.13, height = 0, seed = 20260914),
             size = 1.0, alpha = 0.62, colour = COL_BLUE) +
  geom_segment(data = state_medians,
               aes(x = x - 0.22, xend = x + 0.22, y = delta_LM_minus_PT, yend = delta_LM_minus_PT),
               inherit.aes = FALSE, linewidth = 0.9, colour = COL_ACCENT) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(title = "D  Residual P1 elevation across conventional myeloid states",
       subtitle = expression("57 evaluable patient " %*% " state pairs; bars indicate within-state medians"),
       x = NULL, y = expression("Residual LM " - " PT difference in P1")) +
  theme_s3() +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(linewidth = 0.22, colour = COL_GRID),
        axis.text.x = element_text(size = 5.0, angle = 0, hjust = 0.5, vjust = 0.5,
                                   lineheight = 0.9))

## ------------------------------------------------------------------------
## Assemble and export
## ------------------------------------------------------------------------

figure_s3 <- ((pA | pB) / pC / pD) +
  plot_layout(heights = c(1.00, 1.12, 1.12))

figure_base <- file.path(figure_dir, "Supplementary_Figure_S3_R_v1")
figure_width_mm <- 183
figure_height_mm <- 175
figure_dpi <- 600

ragg::agg_png(paste0(figure_base, ".png"), width = figure_width_mm, height = figure_height_mm,
              units = "mm", res = figure_dpi, background = "white")
print(figure_s3)
dev.off()

svglite::svglite(paste0(figure_base, ".svg"), width = figure_width_mm / 25.4,
                height = figure_height_mm / 25.4, bg = "white")
print(figure_s3)
dev.off()

grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = figure_width_mm / 25.4,
                     height = figure_height_mm / 25.4, family = "Arial")
print(figure_s3)
dev.off()

ragg::agg_tiff(paste0(figure_base, ".tiff"), width = figure_width_mm, height = figure_height_mm,
               units = "mm", res = figure_dpi, compression = "lzw", background = "white")
print(figure_s3)
dev.off()

figure_manifest <- data.frame(
  panel = c("A", "B", "C", "D"),
  evidence = c(
    "Discovery-cell descriptive relationship of P1 to literature-derived comparators, single genes and marker-excluded P1 scores",
    "External GSE245552 cell-level P1 versus literature-derived SPP1-TAM marker score by frozen five-state annotation",
    "Patient-level corrected state-only residual LM−PT delta, with frozen bootstrap/sign-flip summary",
    "Corrected state-only residual LM−PT delta across 57 evaluable patient×state pairs"
  ),
  n = c("12,433 discovery cells/signature", "24,231 external cells", "13 patients; 57 pairs", "57 patient×state pairs"),
  source_audit = c(
    "Literature-derived comparator discovery cell correlations; cells are descriptive units",
    "Frozen external cell-level score table; no inferential statistics shown",
    "Frozen corrected state-only residual table and patient-level summary",
    "Frozen corrected state-only residual table"
  ),
  stringsAsFactors = FALSE
)
write_source(figure_manifest, "S3_figure_manifest.csv")

writeLines(c(
  "Supplementary Figure S3 source-data note",
  "",
  "Panel A is a discovery-cell descriptive comparator audit (12,433 cells); these cells are not independent biological replicates.",
  "Panels B–D use the external GSE245552 frozen score/state table and corrected state-only residualization outputs.",
  "Panel C uses the corrected model P1_full ~ literature-derived SPP1-TAM marker score + state, with tissue retained for the LM−PT contrast.",
  "For Panel C, the frozen patient-level effect is the mean of the available state-specific LM−PT deltas within each patient.",
  "Panel C reports the frozen patient-level bootstrap interval and sign-flip P value; no new inferential statistic is calculated by this plotting script.",
  "Panel D shows all 57 evaluable patient×state residual deltas and overlays only the within-state median.",
  "Panel B contains no fitted line and no cell-level inferential test.",
  "The prior tissue-plus-state residual output is legacy audit material and is not used for the main S3 panels."
), file.path(output_root, "S3_source_data_note.txt"))

message("S3 R build complete: ", normalizePath(output_root, winslash = "/", mustWork = FALSE))
