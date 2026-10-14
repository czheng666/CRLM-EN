## ============================================================================
## Supplementary Figure S5 v5
## Directional portability and pair-level distribution
##
## Unified ggplot2 + patchwork redesign:
##   A. State-level evaluability and directional portability
##   C. Compact external-cohort directional concordance summary
##   B. Pair-level GSE231559 T-N distribution
##
## Plotting-only / deterministic summaries.
## No model fitting, score construction, hypothesis testing, new P values,
## effect-size pooling, or re-estimation of frozen quantities.
##
## GSE231559 contrast = T - N.
## GSE245552 contrast = LM - PT.
## P1-state has source observations but zero evaluable N/T pairs; it is NE,
## not a negative directional-concordance result.
## ============================================================================

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

options(stringsAsFactors = FALSE)
set.seed(20260914)

required_packages <- c(
  "ggplot2",
  "patchwork",
  "ggbeeswarm",
  "scales",
  "ragg",
  "svglite"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

## ============================================================================
## 0. Paths
## ============================================================================

frozen_root <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")

results_tables <- file.path(
  frozen_root,
  "complete_run_archive_20260912",
  "run_output",
  "results",
  "tables"
)

external_scRNA <- file.path(
  frozen_root,
  "Figure2_source_data_addendum_20260912",
  "external_scRNA"
)

output_root <- file.path(CRLM_OUTPUT_ROOT, "Supplementary_Figure_S5_R_v5")
figure_dir <- file.path(output_root, "figures")
source_dir <- file.path(output_root, "source_data")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

file_projection <- file.path(
  external_scRNA,
  "GSE231559_frozen_projection.csv"
)

file_delta_candidates <- file.path(
  results_tables,
  c("62_GSE231559_directional_deltas.csv", "62_GSE231559_replication_deltas.csv")
)
file_delta <- file_delta_candidates[file.exists(file_delta_candidates)][1]

file_summary_candidates <- file.path(
  results_tables,
  c("63_GSE231559_directional_portability_summary.csv", "63_GSE231559_replication_summary.csv")
)
file_summary <- file_summary_candidates[file.exists(file_summary_candidates)][1]

file_gse245552 <- file.path(
  results_tables,
  "59_GSE245552_P1_Pk_state_pair_deltas.csv"
)

## ============================================================================
## 1. Helpers, symbols, and style
## ============================================================================

read_csv_safe <- function(path) {
  if (!file.exists(path)) {
    stop("Missing frozen input: ", path)
  }

  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

write_source <- function(x, filename) {
  write.csv(
    x,
    file.path(source_dir, filename),
    row.names = FALSE,
    na = ""
  )
}

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))

  if (length(missing) > 0L) {
    stop(
      object_name,
      " missing columns: ",
      paste(missing, collapse = ", ")
    )
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
COL_MUTED <- "#59636C"
COL_GRID <- "#D9E0E5"
COL_NEUTRAL <- "#C9D1D7"
COL_POS <- "#4E9A80"
COL_NEG <- "#B6534E"
COL_NONPOS <- "#DDE3E7"
COL_ZERO <- "#7A8790"
COL_EXT1 <- "#3E739F"
COL_EXT2 <- "#4E9A80"

## Shared positive-direction colour across A, B, and C.

## Runtime-generated symbols avoid encoding/font boxes in exported figures.
minus_sign <- intToUtf8(0x2212)
times_sign <- intToUtf8(0x00D7)
geq_sign <- intToUtf8(0x2265)
delta_sign <- intToUtf8(0x0394)

theme_s5 <- function(base_size = 6.5) {
  theme_classic(
    base_size = base_size,
    base_family = "Arial"
  ) +
    theme(
      text = element_text(colour = COL_DARK),
      axis.title = element_text(size = 6.1, colour = COL_DARK),
      axis.text = element_text(size = 5.5, colour = COL_DARK),
      axis.line = element_line(linewidth = 0.3, colour = COL_DARK),
      axis.ticks = element_line(linewidth = 0.3, colour = COL_DARK),
      plot.title = element_text(
        size = 7.5,
        face = "bold",
        hjust = 0,
        colour = COL_DARK,
        margin = margin(0, 0, 1.5, 0)
      ),
      plot.subtitle = element_text(
        size = 5.5,
        colour = COL_MUTED,
        hjust = 0,
        margin = margin(0, 0, 2.0, 0)
      ),
      legend.title = element_text(size = 5.5, colour = COL_DARK),
      legend.text = element_text(size = 5.2, colour = COL_DARK),
      legend.key.height = grid::unit(2.7, "mm"),
      legend.key.width = grid::unit(3.5, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(0.8, "mm"),
      plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm")
    )
}

state_order <- c(
  "Mono",
  "Other",
  "P1_state",
  "P2_inflammatory",
  "P3_C1Q_MHCII"
)

state_labels <- c(
  Mono = "Mono",
  Other = "Other",
  P1_state = "P1-state",
  P2_inflammatory = "Inflammatory",
  P3_C1Q_MHCII = "C1Q/MHC-II"
)

## ============================================================================
## 2. Frozen evaluability source and state-level directional concordance summary
## ============================================================================

projection <- read_csv_safe(file_projection)

assert_columns(
  projection,
  c(
    "score_version",
    "sample",
    "patient",
    "tissue",
    "reference_class",
    "state",
    "n_cells"
  ),
  "GSE231559 projection"
)

projection <- projection[
  projection$score_version == "P1_full" &
    projection$tissue %in% c("N", "T") &
    projection$state %in% state_order,
  ,
  drop = FALSE
]

projection$n_cells <- as.numeric(projection$n_cells)

if (nrow(projection) != 130L ||
    length(unique(projection$sample)) != 26L ||
    length(unique(projection$patient)) != 18L ||
    !setequal(unique(projection$state), state_order)) {
  stop("Frozen GSE231559 projection dimensions or state definition failed.")
}

count_df <- aggregate(
  n_cells ~ patient + state + tissue,
  data = projection[, c("patient", "state", "tissue", "n_cells")],
  FUN = sum,
  na.rm = TRUE
)

patients <- sort(unique(count_df$patient))

eval_df <- merge(
  expand.grid(
    patient = patients,
    state = state_order,
    stringsAsFactors = FALSE
  ),
  count_df[
    count_df$tissue == "N",
    c("patient", "state", "n_cells")
  ],
  by = c("patient", "state"),
  all.x = TRUE,
  sort = FALSE
)
names(eval_df)[names(eval_df) == "n_cells"] <- "n_N"

eval_df <- merge(
  eval_df,
  count_df[
    count_df$tissue == "T",
    c("patient", "state", "n_cells")
  ],
  by = c("patient", "state"),
  all.x = TRUE,
  sort = FALSE
)
names(eval_df)[names(eval_df) == "n_cells"] <- "n_T"

eval_df$n_N[is.na(eval_df$n_N)] <- 0
eval_df$n_T[is.na(eval_df$n_T)] <- 0
eval_df$n_N <- as.integer(round(eval_df$n_N))
eval_df$n_T <- as.integer(round(eval_df$n_T))
eval_df$evaluable <- eval_df$n_N >= 20L & eval_df$n_T >= 20L
eval_df$state_label <- unname(state_labels[eval_df$state])

expected_evaluable <- c(
  Mono = 4L,
  Other = 6L,
  P1_state = 0L,
  P2_inflammatory = 1L,
  P3_C1Q_MHCII = 7L
)

observed_evaluable <- tapply(eval_df$evaluable, eval_df$state, sum)

if (sum(eval_df$evaluable) != 18L ||
    any(observed_evaluable[state_order] != expected_evaluable[state_order])) {
  stop("Frozen evaluability audit does not reproduce the expected 18 pairs.")
}

write_source(
  eval_df[order(match(eval_df$patient, patients), match(eval_df$state, state_order)), ],
  "S5_evaluability_audit_source.csv"
)

delta <- read_csv_safe(file_delta)

assert_columns(
  delta,
  c("score_version", "patient", "state", "delta_T_minus_N"),
  "GSE231559 directional concordance delta"
)

delta <- delta[
  delta$score_version == "P1_full" &
    delta$state %in% state_order &
    is.finite(as.numeric(delta$delta_T_minus_N)),
  ,
  drop = FALSE
]
delta$delta_T_minus_N <- as.numeric(delta$delta_T_minus_N)

eval_key <- paste(
  eval_df$patient[eval_df$evaluable],
  eval_df$state[eval_df$evaluable],
  sep = "||"
)
delta_key <- paste(delta$patient, delta$state, sep = "||")

if (nrow(delta) != 18L ||
    length(unique(delta_key)) != 18L ||
    !setequal(eval_key, delta_key)) {
  stop("Frozen T-N delta table does not match the 18 evaluable pairs.")
}

if (sum(delta$delta_T_minus_N > 0) != 13L ||
    any(delta$state == "P1_state")) {
  stop("Frozen T-N directional concordance or P1-state exclusion failed.")
}

summary_231 <- read_csv_safe(file_summary)
assert_columns(
  summary_231,
  c("score_version", "n_pairs", "n_positive", "median_delta"),
  "GSE231559 directional concordance summary"
)

summary_231_p1 <- summary_231[
  summary_231$score_version == "P1_full",
  ,
  drop = FALSE
]

if (nrow(summary_231_p1) != 1L ||
    as.integer(summary_231_p1$n_pairs) != 18L ||
    as.integer(summary_231_p1$n_positive) != 13L) {
  stop("Frozen GSE231559 summary does not reproduce 18/13.")
}

assert_near(
  median(delta$delta_T_minus_N),
  as.numeric(summary_231_p1$median_delta),
  tolerance = 5e-10,
  label = "GSE231559 median T-N delta"
)

state_summary <- data.frame(
  state = state_order,
  positive = vapply(
    state_order,
    function(ss) {
      n_ss <- sum(eval_df$evaluable[eval_df$state == ss])
      if (n_ss == 0L) NA_integer_ else sum(delta$delta_T_minus_N[delta$state == ss] > 0)
    },
    integer(1)
  ),
  evaluable = as.integer(expected_evaluable[state_order]),
  stringsAsFactors = FALSE
)

state_summary$nonpositive <- ifelse(
  state_summary$evaluable > 0,
  state_summary$evaluable - state_summary$positive,
  NA_integer_
)
state_summary$state_label <- unname(state_labels[state_summary$state])

if (!identical(
  as.integer(state_summary$positive[match(state_order, state_summary$state)]),
  c(4L, 2L, NA_integer_, 0L, 7L)
)) {
  stop("State-level positive counts do not reproduce the frozen audit.")
}

## ============================================================================
## 3. PANEL A — state-level evaluability and directional portability
## ============================================================================

A_state_order <- c(
  "P3_C1Q_MHCII",
  "Mono",
  "Other",
  "P2_inflammatory",
  "P1_state"
)

A_state <- state_summary
A_state$state_plot <- factor(
  A_state$state,
  levels = rev(A_state_order),
  labels = unname(state_labels[rev(A_state_order)])
)

A_state_rect <- rbind(
  data.frame(
    state_plot = A_state$state_plot,
    direction = "Positive",
    xmin = 0,
    xmax = A_state$positive,
    stringsAsFactors = FALSE
  ),
  data.frame(
    state_plot = A_state$state_plot,
    direction = "Non-positive",
    xmin = A_state$positive,
    xmax = A_state$evaluable,
    stringsAsFactors = FALSE
  )
)
A_state_rect <- A_state_rect[
  is.finite(A_state_rect$xmin) & is.finite(A_state_rect$xmax),
  ,
  drop = FALSE
]
A_state_rect$direction <- factor(
  A_state_rect$direction,
  levels = c("Positive", "Non-positive")
)

state_labels_df <- data.frame(
  state_plot = A_state$state_plot,
  x = ifelse(A_state$evaluable > 0, A_state$evaluable + 0.18, 0.30),
  label = ifelse(
    A_state$evaluable == 0,
    "NE (0 evaluable)",
    paste0(A_state$positive, "/", A_state$evaluable)
  ),
  evaluable = A_state$evaluable,
  stringsAsFactors = FALSE
)

pA_state <- ggplot() +
  geom_rect(
    data = A_state_rect,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = as.numeric(state_plot) - 0.25,
      ymax = as.numeric(state_plot) + 0.25,
      fill = direction
    ),
    colour = "white",
    linewidth = 0.35
  ) +
  geom_text(
    data = state_labels_df[state_labels_df$evaluable > 0, , drop = FALSE],
    aes(x = x, y = state_plot, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.15,
    fontface = "bold",
    colour = COL_DARK
  ) +
  geom_text(
    data = state_labels_df[state_labels_df$evaluable == 0, , drop = FALSE],
    aes(x = x, y = state_plot, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.05,
    fontface = "italic",
    colour = COL_MUTED
  ) +
  scale_fill_manual(
    values = c(
      "Positive" = COL_POS,
      "Non-positive" = COL_NEG
    ),
    name = NULL,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 8.3),
    breaks = c(0, 2, 4, 6, 8),
    expand = c(0, 0)
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "A  State-level evaluability and directional portability",
    subtitle = "Evaluable pairs split by directional concordance",
    x = "Number of evaluable pairs",
    y = NULL
  ) +
  theme_s5() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      linewidth = 0.14,
      colour = "#E4E8EB"
    ),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 5.0),
    plot.margin = margin(2.5, 2.5, 2.2, 2.5, "mm")
  )

## ============================================================================
## 4. PANEL A2 — compact external-cohort summary
## ============================================================================

ext245 <- read_csv_safe(file_gse245552)
assert_columns(
  ext245,
  c("program", "delta_LM_minus_PT"),
  "GSE245552 directional pairs"
)

ext245 <- ext245[
  ext245$program == "P1" &
    is.finite(as.numeric(ext245$delta_LM_minus_PT)),
  ,
  drop = FALSE
]
ext245$delta_LM_minus_PT <- as.numeric(ext245$delta_LM_minus_PT)

if (nrow(ext245) != 50L ||
    sum(ext245$delta_LM_minus_PT > 0) != 45L) {
  stop("Expected 45/50 positive GSE245552 LM-PT pairs.")
}

cohort_df <- data.frame(
  cohort = c("GSE245552", "GSE231559"),
  contrast_ascii = c("LM-PT", "T-N"),
  positive = c(45L, 13L),
  total = c(50L, 18L),
  stringsAsFactors = FALSE
)
cohort_df$nonpositive <- cohort_df$total - cohort_df$positive
cohort_df$fraction <- cohort_df$positive / cohort_df$total
cohort_df$cohort_label <- c(
  paste0("GSE245552  LM ", minus_sign, " PT"),
  paste0("GSE231559  T ", minus_sign, " N")
)
cohort_df$cohort_plot <- factor(
  cohort_df$cohort_label,
  levels = rev(cohort_df$cohort_label)
)
cohort_df$y_num <- c(1.65, 1.00)
cohort_df$display <- paste0(
  cohort_df$positive,
  "/",
  cohort_df$total,
  " (",
  sprintf("%.0f%%", 100 * cohort_df$fraction),
  ")"
)

cohort_rect <- rbind(
  data.frame(
    cohort_plot = cohort_df$cohort_plot,
    y_num = cohort_df$y_num,
    direction = "Positive",
    xmin = 0,
    xmax = cohort_df$fraction,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort_plot = cohort_df$cohort_plot,
    y_num = cohort_df$y_num,
    direction = "Non-positive",
    xmin = cohort_df$fraction,
    xmax = 1,
    stringsAsFactors = FALSE
  )
)
cohort_rect$direction <- factor(
  cohort_rect$direction,
  levels = c("Positive", "Non-positive")
)

cohort_label_df <- data.frame(
  cohort_plot = cohort_df$cohort_plot,
  y_num = cohort_df$y_num,
  x = cohort_df$fraction / 2,
  label = cohort_df$display,
  stringsAsFactors = FALSE
)

pA_ext <- ggplot() +
  annotate(
    "segment",
    x = 0,
    xend = 1,
    y = 1.94,
    yend = 1.94,
    linewidth = 0.30,
    colour = "#E4E8EB"
  ) +
  geom_rect(
    data = cohort_rect,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = y_num - 0.17,
      ymax = y_num + 0.17,
      fill = direction
    ),
    colour = "white",
    linewidth = 0.30
  ) +
  geom_text(
    data = cohort_label_df,
    aes(x = x, y = y_num, label = label),
    inherit.aes = FALSE,
    colour = "white",
    size = 2.0,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Positive" = COL_POS,
      "Non-positive" = COL_NONPOS
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = cohort_df$y_num,
    labels = cohort_df$cohort_label,
    limits = c(0.76, 2.08),
    expand = c(0, 0)
  ) +
  labs(
    title = "C  External-cohort directional concordance",
    x = "Positive fraction",
    y = NULL
  ) +
  theme_s5() +
  theme(
    plot.title = element_text(size = 6.5, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(
      linewidth = 0.14,
      colour = "#E4E8EB"
    ),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 4.75, lineheight = 0.88),
    axis.text.x = element_text(size = 4.9),
    plot.margin = margin(2.0, 2.5, 2.5, 2.5, "mm")
  )

pA <- (
  pA_state /
    pA_ext
) +
  plot_layout(heights = c(1.55, 0.65))

write_source(
  state_summary[, c(
    "state",
    "state_label",
    "positive",
    "nonpositive",
    "evaluable"
  )],
  "S5A_state_level_directional_concordance_summary.csv"
)

write_source(
  cohort_df[, c(
    "cohort",
    "contrast_ascii",
    "positive",
    "nonpositive",
    "total",
    "fraction",
    "display"
  )],
  "S5A_external_cohort_summary.csv"
)

## ============================================================================
## 5. PANEL B — horizontal quasirandom distribution
## ============================================================================

## The authoritative pair-level object is reconstructed directly from the frozen
## GSE231559 directional-deltas table, preserving the user-specified
## patient/state/delta mapping without fitting or transformation.
B <- data.frame(
  patient = as.character(delta$patient),
  state = as.character(delta$state),
  delta = as.numeric(delta$delta_T_minus_N),
  stringsAsFactors = FALSE
)

B$state_label <- unname(state_labels[B$state])
B_state_order <- c(
  "P3_C1Q_MHCII",
  "Mono",
  "Other",
  "P2_inflammatory"
)

state_axis_labels_B <- c(
  P3_C1Q_MHCII = paste0("C1Q/MHC-II\n(7/7 positive)"),
  Mono = "Mono (4/4 positive)",
  Other = "Other (2/6 positive)",
  P2_inflammatory = paste0("Inflammatory\n(0/1 positive)")
)

B$state_plot <- factor(
  B$state,
  levels = rev(B_state_order),
  labels = unname(state_axis_labels_B[rev(B_state_order)])
)
B$direction <- ifelse(B$delta > 0, "Positive", "Negative")

write_source(
  B[order(match(B$state, B_state_order), B$delta),
    c("patient", "state", "state_label", "delta", "direction")],
  "S5B_GSE231559_pair_level_quasirandom_source.csv"
)

pB <- ggplot(
  B,
  aes(
    x = delta,
    y = state_plot
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.40,
    colour = COL_ZERO
  ) +
  ggbeeswarm::geom_quasirandom(
    aes(colour = direction),
    groupOnX = FALSE,
    orientation = "x",
    width = 0.18,
    size = 1.95,
    alpha = 0.95
  ) +
  scale_colour_manual(
    values = c(
      "Negative" = COL_NEG,
      "Positive" = COL_POS
    ),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-0.85, 1.20),
    breaks = c(-0.5, 0, 0.5, 1.0),
    expand = c(0, 0)
  ) +
  labs(
    title = "B  Pair-level P1 directional concordance in GSE231559",
    subtitle = paste0(
      "18 evaluable patient ",
      times_sign,
      " state pairs; all pair-level effects shown"
    ),
    x = paste0("P1 ", delta_sign, " (T ", minus_sign, " N)"),
    y = NULL
  ) +
  theme_s5() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(
      linewidth = 0.20,
      colour = COL_GRID
    ),
    panel.grid.major.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 4.75, lineheight = 0.78),
    plot.margin = margin(2.5, 2.5, 2.5, 2.5, "mm")
  )

## ============================================================================
## 6. Final layout and export
## ============================================================================

figure_s5 <- (
  pA |
    pB
) +
  plot_layout(widths = c(0.82, 1.38))

figure_base <- file.path(
  figure_dir,
  "Supplementary_Figure_S5_R_v5"
)

figure_width_mm <- 183
figure_height_mm <- 105
figure_dpi <- 600

ragg::agg_png(
  paste0(figure_base, ".png"),
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  res = figure_dpi,
  background = "white"
)
print(figure_s5)
dev.off()

svglite::svglite(
  paste0(figure_base, ".svg"),
  width = figure_width_mm / 25.4,
  height = figure_height_mm / 25.4,
  bg = "white"
)
print(figure_s5)
dev.off()

grDevices::cairo_pdf(
  paste0(figure_base, ".pdf"),
  width = figure_width_mm / 25.4,
  height = figure_height_mm / 25.4,
  family = "Arial"
)
print(figure_s5)
dev.off()

ragg::agg_tiff(
  paste0(figure_base, ".tiff"),
  width = figure_width_mm,
  height = figure_height_mm,
  units = "mm",
  res = figure_dpi,
  compression = "lzw",
  background = "white"
)
print(figure_s5)
dev.off()

## ============================================================================
## 7. Manifest and provenance
## ============================================================================

figure_manifest <- data.frame(
  panel = c("A", "C", "B"),
  evidence = c(
    "State-level GSE231559 P1 T-N evaluability and directional concordance",
    "External directional concordance summary for GSE245552 LM-PT and GSE231559 T-N",
    "Pair-level GSE231559 P1 T-N distribution with all raw effects shown"
  ),
  n = c(
    "Mono 4/4; Other 2/6; Inflammatory 0/1; C1Q/MHC-II 7/7; P1-state NE",
    "GSE245552 45/50; GSE231559 13/18",
    "18 evaluable patient x state pairs; 13 positive"
  ),
  stringsAsFactors = FALSE
)

write_source(figure_manifest, "S5_figure_manifest.csv")

writeLines(
  c(
    "Supplementary Figure S5 v5 source-data note",
    "",
    "This is a plotting-only / deterministic-summary script.",
    "No model fitting, score construction, hypothesis testing, new P-value calculation, effect-size pooling, or re-estimation of frozen quantities was performed.",
    "",
    "Panel A summarizes state-level GSE231559 T-N directional concordance using the frozen evaluability criterion. Bars show evaluable patient-state pairs split into positive and non-positive differences. P1-state has source observations but zero evaluable N/T pairs and is shown as NE (0 evaluable), not as a negative directional-concordance result.",
    "Panel C summarizes directional concordance only: GSE245552 uses LM-PT (45/50 positive), whereas GSE231559 uses T-N (13/18 positive). These contrasts are distinct and are not pooled.",
    "",
    "Panel B displays all 18 evaluable GSE231559 patient-state T-N effects using a horizontal quasirandom distribution. No median overlay or state-specific inferential summary is shown. Patient IDs are retained in the source-data table but are intentionally omitted from the figure.",
    "",
    "The complete patient x state N/T cell-count audit is retained in S5_evaluability_audit_source.csv even though the full matrix is not shown in this v5 visual design."
  ),
  file.path(output_root, "S5_source_data_note.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(output_root, "S5_sessionInfo.txt")
)

message(
  "S5 v5 plotting-only build complete: ",
  normalizePath(output_root, winslash = "/", mustWork = FALSE)
)
