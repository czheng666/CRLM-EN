## Supplementary Figure S4
## Composition-activity decomposition and sensitivity analyses
##
## PLOTTING-ONLY SCRIPT
## No model fitting, re-estimation, bootstrap, CI calculation or new P values.

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

options(stringsAsFactors = FALSE)
set.seed(20260914)

required_packages <- c("ggplot2", "patchwork", "ragg", "svglite", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ------------------------------------------------------------------------
## Paths
## ------------------------------------------------------------------------

freeze_dir <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
main_tables <- file.path(freeze_dir, "tables")
results_tables <- paste0(
  freeze_dir, "/",
  "complete_run_archive_20260912/run_output/results/tables"
)

output_root <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S4_R_v1")
figure_dir <- file.path(output_root, "figures")
source_dir <- file.path(output_root, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(directory, filename) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) stop("Missing frozen table: ", path)
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

assert_near <- function(observed, expected, tolerance = 1e-6, label = "") {
  if (length(observed) != 1L || is.na(observed) ||
      abs(as.numeric(observed) - as.numeric(expected)) > tolerance) {
    stop("Frozen value mismatch: ", label,
         ". observed=", observed, ", expected=", expected)
  }
}

## ------------------------------------------------------------------------
## Style
## ------------------------------------------------------------------------

COL_DARK <- "#22282E"
COL_GRID <- "#D9E0E5"
COL_COMP <- "#7A94A3"
COL_ACTIVITY <- "#4E9A80"
COL_ACCENT <- "#B6534E"
COL_NONP1 <- "#3E739F"
COL_LINK <- "#C8D0D6"

theme_s4 <- function(base_size = 6.4) {
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
      legend.box.spacing = grid::unit(0.8, "mm"),
      plot.title = element_text(size = 7.4, face = "bold", hjust = 0,
                                colour = COL_DARK,
                                margin = margin(1.2, 0, 2.0, 0)),
      plot.subtitle = element_text(size = 5.5, colour = "#59636C", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = 5.6, face = "bold", colour = COL_DARK),
      plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm")
    )
}

component_cols <- c(
  "Composition" = COL_COMP,
  "Within-state activity" = COL_ACTIVITY
)

## ------------------------------------------------------------------------
## Panel A — overall decomposition
## ------------------------------------------------------------------------

A_raw <- read_table(main_tables, "composition_activity_contribution_summary.csv")
assert_columns(A_raw, c("estimand", "component", "estimate", "ci_low", "ci_high", "n_patients"), "A_raw")

A <- A_raw[
  A_raw$estimand == "ratio_of_sums_signed_contribution" &
    A_raw$component %in% c("composition", "within_state_activity"),
  c("component", "estimate", "ci_low", "ci_high", "n_patients"),
  drop = FALSE
]
A$component[A$component == "composition"] <- "Composition"
A$component[A$component == "within_state_activity"] <- "Within-state activity"
A$component <- factor(A$component, levels = names(component_cols))
A$estimate <- as.numeric(A$estimate)
A$ci_low <- as.numeric(A$ci_low)
A$ci_high <- as.numeric(A$ci_high)
A$percent <- 100 * A$estimate
A$ci_low_percent <- 100 * A$ci_low
A$ci_high_percent <- 100 * A$ci_high

if (nrow(A) != 2L || length(unique(A$n_patients)) != 1L || unique(A$n_patients) != 12L) {
  stop("Panel A must contain two ratio-of-sums components across 12 patients.")
}
assert_near(A$estimate[A$component == "Composition"], 0.1915785, 5e-4,
            "Panel A composition fraction")
assert_near(A$estimate[A$component == "Within-state activity"], 0.8084215, 5e-4,
            "Panel A within-state activity fraction")
assert_near(sum(A$estimate), 1, 1e-6, "Panel A decomposition closure")
write_source(A, "S4A_overall_composition_activity_decomposition.csv")

A_plot <- data.frame(
  x = "Overall",
  component = A$component,
  percent = A$percent,
  stringsAsFactors = FALSE
)

pA <- ggplot(A_plot, aes(y = x, x = percent, fill = component)) +
  geom_col(width = 0.45, colour = "white", linewidth = 0.35,
           position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = sprintf("%.1f%%", percent)),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 2.35, colour = "white", fontface = "bold") +
  scale_fill_manual(values = component_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(title = paste0("A  Overall composition", intToUtf8(0x2013),
                      "activity decomposition"),
       subtitle = "Ratio-of-sums signed contribution across 12 evaluable patients",
       x = "Contribution to overall P1 elevation", y = NULL) +
  theme_s4() +
  theme(legend.position = "bottom", axis.text.y = element_blank(),
        plot.title = element_text(family = "DejaVu Sans", size = 7.4,
                                  face = "bold", hjust = 0, colour = COL_DARK,
                                  margin = margin(1.2, 0, 2.0, 0)),
        axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())

## ------------------------------------------------------------------------
## Panel B — patient-level paired decomposition
## ------------------------------------------------------------------------

B_raw <- read_table(main_tables, "composition_activity_contribution_by_patient.csv")
assert_columns(B_raw, c("score_version", "patient", "overall_delta",
                        "state_composition", "within_state_activity", "closure_error"), "B_raw")
B <- B_raw[B_raw$score_version == "P1_full", c(
  "patient", "overall_delta", "state_composition", "within_state_activity", "closure_error"
), drop = FALSE]
names(B) <- c("patient", "overall_delta", "composition", "activity", "closure_error")
B[] <- lapply(B, function(x) if (is.character(x)) x else as.numeric(x))
B$patient <- as.character(B$patient)

if (nrow(B) != 12L) stop("Panel B must contain 12 frozen evaluable patients.")
if (sum(B$overall_delta > 0) != 12L || sum(B$activity > B$composition) != 11L ||
    sum(B$activity > 0) != 12L) {
  stop("Panel B qualitative frozen audit does not reproduce 12/12, 11/12 and 12/12.")
}
max_closure <- max(abs(B$closure_error), na.rm = TRUE)
if (max_closure > 1e-12) stop("Panel B closure error exceeds tolerance: ", max_closure)
B <- B[order(B$activity - B$composition), , drop = FALSE]
B$patient <- factor(B$patient, levels = B$patient)
write_source(B, "S4B_patient_level_composition_activity.csv")

pB <- ggplot(B) +
  geom_segment(aes(x = composition, xend = activity, y = patient, yend = patient),
               linewidth = 0.65, colour = COL_LINK) +
  geom_point(aes(x = composition, y = patient, colour = "Composition"), size = 1.75) +
  geom_point(aes(x = activity, y = patient, colour = "Within-state activity"), size = 1.75) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  scale_colour_manual(values = component_cols, name = NULL) +
  scale_x_continuous(limits = c(min(B$composition, na.rm = TRUE) - 0.10,
                                max(B$activity, na.rm = TRUE) + 0.14),
                     breaks = scales::breaks_pretty(n = 5), expand = c(0, 0)) +
  labs(title = "B  Patient-level composition and within-state contributions",
       subtitle = "12 evaluable patients; activity > composition in 11/12; activity positive in 12/12",
       x = "Signed contribution to P1 elevation", y = NULL) +
  theme_s4() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

## ------------------------------------------------------------------------
## Panel C — cell-count threshold sensitivity
## ------------------------------------------------------------------------

C_raw <- read_table(results_tables, "52_GSE245552_composition_activity_summary_threshold_LOPO.csv")
assert_columns(C_raw, c("score_version", "threshold", "component", "n_patients",
                        "median", "ci_low", "ci_high"), "C_raw")
C <- C_raw[C_raw$score_version == "P1_full" & C_raw$threshold %in% c(1, 5, 10) &
             C_raw$component %in% c("state_composition", "within_state_activity"),
           c("threshold", "component", "n_patients", "median", "ci_low", "ci_high"),
           drop = FALSE]
C$component[C$component == "state_composition"] <- "Composition"
C$component[C$component == "within_state_activity"] <- "Within-state activity"
C$component <- factor(C$component, levels = names(component_cols))
C$threshold <- as.integer(C$threshold)
C$threshold_label <- factor(paste0("threshold_", C$threshold),
                            levels = c("threshold_10", "threshold_5", "threshold_1"))
C$median <- as.numeric(C$median)
C$ci_low <- as.numeric(C$ci_low)
C$ci_high <- as.numeric(C$ci_high)

if (nrow(C) != 6L || any(C$ci_low > C$median | C$ci_high < C$median)) {
  stop("Panel C must contain six valid frozen median/CI rows.")
}
for (thr in c(1, 5, 10)) {
  if (sum(C$threshold == thr) != 2L) stop("Panel C threshold is incomplete: ", thr)
}
assert_near(C$median[C$threshold == 1 & C$component == "Composition"], 0.1000156, 5e-4,
            "Panel C threshold 1 composition")
assert_near(C$median[C$threshold == 1 & C$component == "Within-state activity"], 0.7361997, 5e-4,
            "Panel C threshold 1 activity")
assert_near(C$median[C$threshold == 10 & C$component == "Within-state activity"], 0.8133281, 5e-4,
            "Panel C threshold 10 activity")
write_source(C, "S4C_cell_count_threshold_sensitivity.csv")

pC <- ggplot(C, aes(x = threshold_label, y = median, colour = component)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.12, linewidth = 0.65,
                position = position_dodge(width = 0.36)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.36)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  scale_colour_manual(values = component_cols, name = NULL) +
  labs(title = "C  Sensitivity to cell-count threshold",
       subtitle = "Frozen median contributions with 95% bootstrap confidence intervals",
       x = NULL, y = "Signed contribution to P1 elevation") +
  coord_flip() +
  scale_x_discrete(
    labels = c(
      threshold_1 = expression(phantom(1) >= "1 cell"),
      threshold_5 = expression(phantom(1) >= "5 cells"),
      threshold_10 = expression(phantom(1) >= "10 cells")
    )
  ) +
  theme_s4() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

## ------------------------------------------------------------------------
## Panel D — P1-state inclusion sensitivity
## ------------------------------------------------------------------------

D_patient <- read_table(results_tables, "53_GSE245552_five_state_by_patient.csv")
assert_columns(D_patient, c("score_version", "patient", "overall_delta",
                            "non_P1_activity", "P1_state_contribution", "closure_error"), "D_patient")
D_patient <- D_patient[D_patient$score_version == "P1_full", , drop = FALSE]
if (nrow(D_patient) != 12L || max(abs(as.numeric(D_patient$closure_error)), na.rm = TRUE) > 1e-12) {
  stop("Panel D patient-level audit does not reproduce 12 rows with closed decomposition.")
}
write_source(D_patient, "S4D_five_state_by_patient_audit.csv")

D_raw <- read_table(results_tables, "54_GSE245552_five_state_presence_sensitivity_summary.csv")
assert_columns(D_raw, c("score_version", "subset", "component", "n_patients",
                        "median", "ci_low", "ci_high"), "D_raw")
D <- D_raw[D_raw$score_version == "P1_full" &
             D_raw$subset %in% c("all_complete_four_state", "both_endpoints_positive") &
             D_raw$component %in% c("non_P1_activity", "P1_state_contribution"),
           c("subset", "component", "n_patients", "median", "ci_low", "ci_high"),
           drop = FALSE]
D$subset_label <- ifelse(D$subset == "all_complete_four_state",
                         "Complete four-state set\n(n = 12)",
                         "Both endpoints positive\n(n = 6)")
D$component[D$component == "non_P1_activity"] <- "Non-P1 activity"
D$component[D$component == "P1_state_contribution"] <- "P1-state contribution"
D$subset_label <- factor(D$subset_label,
                         levels = c("Both endpoints positive\n(n = 6)",
                                    "Complete four-state set\n(n = 12)"))
D$component <- factor(D$component, levels = c("Non-P1 activity", "P1-state contribution"))
D$median <- as.numeric(D$median)
D$ci_low <- as.numeric(D$ci_low)
D$ci_high <- as.numeric(D$ci_high)

if (nrow(D) != 4L || any(D$ci_low > D$median | D$ci_high < D$median)) {
  stop("Panel D must contain four valid frozen median/CI rows.")
}
if (any(D$n_patients[D$subset == "all_complete_four_state"] != 12L) ||
    any(D$n_patients[D$subset == "both_endpoints_positive"] != 6L)) {
  stop("Panel D patient counts do not match the frozen analysis sets.")
}
assert_near(D$median[D$subset == "all_complete_four_state" & D$component == "Non-P1 activity"],
            0.73099565, 5e-4, "Panel D all-complete non-P1 activity")
assert_near(D$median[D$subset == "all_complete_four_state" & D$component == "P1-state contribution"],
            0.01327755, 5e-4, "Panel D all-complete P1-state contribution")
write_source(D, "S4D_P1_state_inclusion_sensitivity.csv")

D_cols <- c("Non-P1 activity" = COL_NONP1, "P1-state contribution" = COL_ACCENT)
pD <- ggplot(D, aes(x = subset_label, y = median, colour = component)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.12, linewidth = 0.65,
                position = position_dodge(width = 0.36)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.36)) +
  geom_hline(yintercept = 0, linewidth = 0.35, linetype = 2, colour = "#7D8A94") +
  scale_colour_manual(values = D_cols, name = NULL) +
  labs(title = "D  Sensitivity to P1-state inclusion",
       subtitle = "Frozen median contributions with 95% bootstrap confidence intervals",
       x = NULL, y = "Signed contribution to P1 elevation") +
  coord_flip() +
  theme_s4() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank())

## ------------------------------------------------------------------------
## Assemble and export
## ------------------------------------------------------------------------

figure_s4 <- ((pA | pC) / (pB | pD)) +
  plot_layout(heights = c(0.90, 1.15), widths = c(1.05, 0.95))

figure_base <- file.path(figure_dir, "Supplementary_Figure_S4_R_v1")
figure_width_mm <- 183
figure_height_mm <- 145
figure_dpi <- 600

ragg::agg_png(paste0(figure_base, ".png"), width = figure_width_mm,
              height = figure_height_mm, units = "mm", res = figure_dpi,
              background = "white")
print(figure_s4)
dev.off()

svglite::svglite(paste0(figure_base, ".svg"), width = figure_width_mm / 25.4,
                height = figure_height_mm / 25.4, bg = "white")
print(figure_s4)
dev.off()

grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = figure_width_mm / 25.4,
                     height = figure_height_mm / 25.4, family = "Arial")
print(figure_s4)
dev.off()

ragg::agg_tiff(paste0(figure_base, ".tiff"), width = figure_width_mm,
               height = figure_height_mm, units = "mm", res = figure_dpi,
               compression = "lzw", background = "white")
print(figure_s4)
dev.off()

figure_manifest <- data.frame(
  panel = c("A", "B", "C", "D"),
  evidence = c(
    "Overall ratio-of-sums composition versus within-state activity decomposition",
    "Frozen patient-level paired composition and within-state contributions",
    "Frozen contribution medians and bootstrap CIs across cell-count thresholds",
    "Frozen non-P1 activity versus P1-state contribution across inclusion sets"
  ),
  n = c("2 components; 12 patients", "12 patients", "3 thresholds × 2 components",
        "2 analysis sets × 2 components"),
  source_audit = c(
    "composition_activity_contribution_summary.csv; ratio_of_sums_signed_contribution only",
    "composition_activity_contribution_by_patient.csv",
    "52_GSE245552_composition_activity_summary_threshold_LOPO.csv",
    "53_GSE245552_five_state_by_patient.csv and 54_GSE245552_five_state_presence_sensitivity_summary.csv"
  ),
  stringsAsFactors = FALSE
)
write_source(figure_manifest, "S4_figure_manifest.csv")

writeLines(c(
  "Supplementary Figure S4 source-data note",
  "",
  "This is a plotting-only script.",
  "No model fitting, decomposition re-estimation, bootstrap resampling, CI calculation, hypothesis testing or new P-value calculation was performed.",
  "Panel A uses the authoritative ratio_of_sums_signed_contribution estimand. The median_patient_ratio_signed_contribution estimand was not plotted.",
  "Panel A frozen values are Composition = 19.16% and Within-state activity = 80.84% across 12 evaluable patients.",
  "Panel B displays all 12 frozen patient-level paired contributions. Overall delta was positive in 12/12 patients, activity exceeded composition in 11/12 patients, and activity was positive in 12/12 patients.",
  paste0("Maximum absolute Panel B closure error = ", format(max_closure, scientific = TRUE, digits = 3), "."),
  "Panel C displays the frozen contribution medians and 95% bootstrap CIs at thresholds of 1, 5 and 10 cells per patient-state unit; no new attributable fraction was calculated.",
  "Panel D is a P1-state inclusion sensitivity analysis comparing non-P1 activity with P1-state contribution in the frozen all-complete-four-state (n = 12) and both-endpoints-positive (n = 6) sets.",
  "State-specific P1 shifts were not interpreted or plotted as decomposition contributions."
), file.path(output_root, "S4_source_data_note.txt"))

capture.output(sessionInfo(), file = file.path(output_root, "S4_sessionInfo.txt"))
message("Supplementary Figure S4 plotting-only build complete: ",
        normalizePath(output_root, winslash = "/", mustWork = FALSE))
