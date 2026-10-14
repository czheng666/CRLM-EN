## ============================================================
## CRLM Figure 4
## Cross-state P1 elevation and SPP1-TAM comparator analysis
##
## Plotting from the 20260912 frozen result package only.
## No model fitting, re-estimation, or new statistical inference.
## ============================================================

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(scales)
})

## ------------------------------------------------------------
## 0. Paths and output
## ------------------------------------------------------------

freeze_dir <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
table_dir  <- file.path(freeze_dir, "tables")

out_dir <- file.path(CRLM_OUTPUT_ROOT, "Figure4_CRLM_cross_state_SPP1_v1")
source_dir <- file.path(out_dir, "source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
message("Reading frozen tables from: ", table_dir)

read_frozen <- function(file) {
  path <- file.path(table_dir, file)
  message("  input: ", file, " -> ", path)
  if (!file.exists(path)) {
    stop("Missing frozen input: ", path)
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

require_cols <- function(dat, cols, label) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "))
  }
}

## ------------------------------------------------------------
## 1. Frozen source tables
## ------------------------------------------------------------

state_map <- c(
  "Mono"            = "Mono",
  "P2_inflammatory" = "Inflammatory",
  "P3_C1Q_MHCII"    = "C1Q/MHC-II",
  "P1_state"        = "P1-state",
  "Other"           = "Other"
)

state_order <- c(
  "Mono",
  "Inflammatory",
  "C1Q/MHC-II",
  "P1-state",
  "Other"
)

## A: authoritative state-specific mixed-model contrasts + BH FDR
A_raw <- read_frozen("patient_state_mixed_model_state_contrasts_BH.csv")
require_cols(
  A_raw,
  c("state", "estimate_LM_minus_PT", "ci_low", "ci_high", "BH_FDR"),
  "A frozen mixed-model table"
)

A_dat <- A_raw %>%
  transmute(
    state = unname(state_map[as.character(state)]),
    estimate = estimate_LM_minus_PT,
    ci_low = ci_low,
    ci_high = ci_high,
    BH_FDR = BH_FDR
  )

if (anyNA(A_dat$state)) stop("A contains an unmapped state label.")

## B: frozen LM-PT deltas. The current frozen table contains more than
## one paired row for some patient-state combinations. For a single
## patient x state heatmap cell, retain the frozen delta values and take
## their arithmetic mean within patient x state; this is a display-level
## collapse, not a new inferential analysis.
B_raw <- read_frozen("GSE245552_patient_state_P1_vs_published_SPP1_TAM_deltas.csv")
require_cols(B_raw, c("patient", "state", "score_version", "delta"), "B frozen delta table")

B_dat <- B_raw %>%
  filter(score_version == "P1_full") %>%
  transmute(
    patient = as.character(patient),
    state = unname(state_map[as.character(state)]),
    delta_P1 = as.numeric(delta)
  ) %>%
  filter(!is.na(state), is.finite(delta_P1)) %>%
  group_by(patient, state) %>%
  summarise(delta_P1 = mean(delta_P1), .groups = "drop")

if (nrow(B_dat) == 0) stop("B has no P1_full rows.")
if (anyNA(B_dat$state)) stop("B contains an unmapped state label.")

patient_order <- unique(
  B_raw$patient[B_raw$score_version == "P1_full"]
)

## C: the first four rows are the frozen patient-level pooled summaries,
## with state missing. They are the prespecified comparator summaries.
C_raw <- read_frozen("GSE245552_patient_level_P1_vs_published_SPP1_TAM_summary_primary.csv")
require_cols(
  C_raw,
  c("score_version", "median_delta", "ci_low", "ci_high", "n_positive", "n_patients", "state"),
  "C frozen comparator summary"
)

sig_order <- c(
  "P1_full",
  "P1_no_SPP1",
  "P1_no_SPP1_APOC1",
  "published_SPP1_TAM"
)

C_source_order <- c(
  "P1_full",
  "P1_no_SPP1",
  "P1_no_SPP1_APOC1",
  "SPP1_TAM_published_2025"
)

C_dat <- C_raw %>%
  filter(
    (is.na(state) | state == ""),
    score_version %in% C_source_order
  ) %>%
  transmute(
    signature = recode(
      score_version,
      P1_full = "P1_full",
      P1_no_SPP1 = "P1_no_SPP1",
      P1_no_SPP1_APOC1 = "P1_no_SPP1_APOC1",
      SPP1_TAM_published_2025 = "published_SPP1_TAM"
    ),
    median_effect = median_delta,
    ci_low = ci_low,
    ci_high = ci_high,
    positive_patients = n_positive,
    total_patients = n_patients
  )

if (!setequal(C_dat$signature, c("P1_full", "P1_no_SPP1", "P1_no_SPP1_APOC1", "published_SPP1_TAM"))) {
  stop("C does not contain exactly the four frozen comparator summaries.")
}
C_dat <- C_dat[match(sig_order, C_dat$signature), , drop = FALSE]

## D: frozen corrected residual scores. The authoritative summary defines
## patient-level residual LM-PT effects by pairing LM and PT corrected
## residual scores within state and then averaging available state deltas
## within patient. No residualization is performed here.
D_score_raw <- read_frozen("GSE245552_patient_tissue_state_P1_residualized_SPP1_TAM_corrected_state_only.csv")
require_cols(
  D_score_raw,
  c("patient", "tissue", "state", "P1_residual_state_only"),
  "D frozen corrected residual table"
)

D_summary_raw <- read_frozen("GSE245552_patient_level_SPP1_residual_comparison_corrected_state_only.csv")
require_cols(
  D_summary_raw,
  c("residual_definition", "median_delta", "ci_low", "ci_high", "n_positive", "n_patients", "signflip_p_two_sided"),
  "D frozen residual summary"
)

D_summary <- D_summary_raw %>%
  filter(residual_definition == "P1_residual_state_only: P1_full ~ published SPP1-TAM + state (corrected)")
if (nrow(D_summary) != 1) stop("D corrected state-only summary is not unique.")

paired_patients <- patient_order

D_pair <- D_score_raw %>%
  filter(patient %in% paired_patients, tissue %in% c("primary_tumor", "liver_metastasis")) %>%
  transmute(
    patient = as.character(patient),
    state = unname(state_map[as.character(state)]),
    tissue_group = if_else(tissue == "primary_tumor", "PT", "LM"),
    corrected_residual = as.numeric(P1_residual_state_only)
  ) %>%
  filter(!is.na(state), is.finite(corrected_residual)) %>%
  group_by(patient, state, tissue_group) %>%
  summarise(corrected_residual = mean(corrected_residual), .groups = "drop") %>%
  pivot_wider(
    names_from = tissue_group,
    values_from = corrected_residual,
    names_prefix = "residual_"
  ) %>%
  mutate(delta_state = residual_LM - residual_PT) %>%
  filter(is.finite(delta_state))

D_dat <- D_pair %>%
  group_by(patient) %>%
  summarise(delta_residual = mean(delta_state), .groups = "drop")

if (nrow(D_dat) != as.integer(D_summary$n_patients)) {
  stop("D patient-level display data does not match frozen n_patients.")
}

## Write the exact plotting inputs used below, for source-data traceability.
utils::write.csv(A_dat, file.path(source_dir, "Figure4A_state_mixed_model_contrasts.csv"), row.names = FALSE)
utils::write.csv(B_dat, file.path(source_dir, "Figure4B_patient_state_delta.csv"), row.names = FALSE)
utils::write.csv(C_dat, file.path(source_dir, "Figure4C_P1_SPP1_comparator.csv"), row.names = FALSE)
utils::write.csv(D_dat, file.path(source_dir, "Figure4D_corrected_residual_patient_delta.csv"), row.names = FALSE)
utils::write.csv(
  D_summary %>% transmute(
    median_delta,
    ci_low,
    ci_high,
    n_positive,
    n_patients,
    p_signflip = signflip_p_two_sided
  ),
  file.path(source_dir, "Figure4D_corrected_residual_summary.csv"), row.names = FALSE
)

## C score-scale audit. This is a provenance statement, not a new analysis.
## All four scores use the same log1p(CP10K) expression input, but P1 and
## leave-out scores are weighted projections with renormalized weights,
## whereas the literature-derived SPP1-TAM marker score is an unweighted module mean. No common
## post-hoc standardization was applied.
C_scale_audit <- data.frame(
  scope = "C four score definitions",
  expression_normalization = "log1p(CP10K)",
  P1_score_construction = "weighted projection using frozen P1 weights; weights renormalized after gene exclusion",
  published_comparator_construction = "unweighted mean of the literature-derived 8-gene SPP1-TAM marker panel",
  common_posthoc_standardization = FALSE,
  raw_cross_score_magnitude_comparable = FALSE,
  interpretation = "Interpret direction, positive-patient consistency, and leave-out persistence; do not compare raw effect magnitudes across score definitions.",
  stringsAsFactors = FALSE
)
utils::write.csv(
  C_scale_audit,
  file.path(source_dir, "Figure4C_score_scale_audit.csv"),
  row.names = FALSE
)

## ------------------------------------------------------------
## 2. Visual language
## ------------------------------------------------------------

COL_P1    <- "#BE4D3B"
COL_RED2  <- "#D17A6B"
COL_RED3  <- "#E2AAA0"
COL_GREY  <- "#A4AEB6"
COL_GREY2 <- "#D6DCE1"
COL_DARK  <- "#22282E"
COL_TEXT2 <- "#71808D"
COL_BLUE  <- "#688DA6"

state_cols <- c(
  "Mono"         = "#5AA0C8",
  "Inflammatory" = "#D99A42",
  "C1Q/MHC-II"   = "#509A83",
  "P1-state"     = "#C45B4A",
  "Other"        = "#96949A"
)

sig_labels <- c(
  "P1_full" = "Frozen P1",
  "P1_no_SPP1" = "P1 excluding SPP1",
  "P1_no_SPP1_APOC1" = "P1 excluding SPP1 and APOC1",
  "published_SPP1_TAM" = "Literature-derived\nSPP1-TAM marker score"
)

sig_cols <- c(
  "P1_full" = COL_P1,
  "P1_no_SPP1" = COL_RED2,
  "P1_no_SPP1_APOC1" = COL_RED3,
  "published_SPP1_TAM" = COL_GREY
)

theme_fig <- theme_classic(base_size = 9) +
  theme(
    text = element_text(family = "Arial", colour = COL_DARK),
    plot.title = element_text(face = "bold", size = 10.5, hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 8, colour = COL_TEXT2, hjust = 0, margin = margin(b = 8)),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, colour = COL_DARK),
    axis.line = element_line(colour = COL_GREY2, linewidth = 0.45),
    axis.ticks = element_line(colour = COL_GREY2, linewidth = 0.4),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(5, 6, 5, 6)
  )

panel_tag_theme <- theme(
  plot.tag = element_text(size = 11, face = "bold", colour = COL_DARK),
  plot.tag.position = c(0, 1)
)

## ------------------------------------------------------------
## 3. Panel A: formal five-state mixed-model forest
## ------------------------------------------------------------

A <- A_dat %>%
  mutate(
    state = factor(state, levels = rev(state_order)),
    y = as.numeric(state),
    FDR_label = if_else(BH_FDR < 0.001, "<0.001", sprintf("%.3f", BH_FDR))
  )

A_range <- range(c(A$ci_low, A$ci_high), na.rm = TRUE)
A_pad <- diff(A_range) * 0.10
A_text_x <- max(A$ci_high, na.rm = TRUE) + A_pad
A_xlim <- c(min(A$ci_low, 0, na.rm = TRUE) - diff(A_range) * 0.04, A_text_x + diff(A_range) * 0.32)

pA <- ggplot(A, aes(colour = state)) +
  geom_vline(xintercept = 0, colour = "#687682", linewidth = 0.55) +
  geom_segment(aes(x = ci_low, xend = ci_high, y = y, yend = y), linewidth = 0.9) +
  geom_point(aes(x = estimate, y = y), size = 3.2) +
  geom_text(
    aes(x = A_text_x, y = y, label = FDR_label),
    hjust = 0, colour = COL_TEXT2, size = 2.8, inherit.aes = FALSE
  ) +
  annotate(
    "text", x = A_text_x, y = 5.58, label = "BH-FDR",
    hjust = 0, size = 2.8, fontface = "bold", colour = COL_TEXT2
  ) +
  scale_colour_manual(values = state_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(rev(state_order)),
    labels = rev(state_order),
    expand = expansion(mult = c(0.12, 0.22))
  ) +
  scale_x_continuous(limits = A_xlim, labels = label_number(accuracy = 0.1)) +
  labs(
    tag = "A",
    title = "Cross-state LM-associated elevation of P1 activity",
    subtitle = "State-specific LM\u2212PT contrasts; 16 patients, 136 patient\u00d7tissue\u00d7state units",
    x = "Estimated LM\u2212PT difference in P1 score",
    y = NULL
  ) +
  theme_fig + panel_tag_theme +
  theme(plot.margin = margin(5, 32, 8, 6))

## ------------------------------------------------------------
## 4. Panel B: patient x state consistency heatmap
## ------------------------------------------------------------

B_grid <- expand_grid(
  patient = patient_order,
  state = state_order
) %>%
  left_join(B_dat, by = c("patient", "state")) %>%
  mutate(
    patient = factor(patient, levels = rev(patient_order)),
    state = factor(state, levels = state_order)
  )

B_limit <- max(abs(B_grid$delta_P1), na.rm = TRUE)

pB <- ggplot(B_grid, aes(x = state, y = patient, fill = delta_P1)) +
  geom_tile(width = 0.94, height = 0.94, colour = "white", linewidth = 0.45) +
  scale_fill_gradient2(
    low = COL_BLUE, mid = "white", high = COL_P1, midpoint = 0,
    limits = c(-B_limit, B_limit), oob = squish, na.value = COL_GREY2,
    name = "LM\u2212PT P1 difference"
  ) +
  labs(
    tag = "B",
    title = "Patient-level consistency\nacross states",
    subtitle = "Paired LM\u2212PT differences in P1 score",
    x = NULL, y = NULL
  ) +
  theme_fig + panel_tag_theme +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text.y = element_text(size = 7.5),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(8, "mm"),
    plot.margin = margin(5, 6, 5, 6)
  )

## ------------------------------------------------------------
## 5. Panel C: P1 versus SPP1-associated comparator forest
## ------------------------------------------------------------

C <- C_dat %>%
  mutate(
    signature = factor(signature, levels = rev(sig_order)),
    y = as.numeric(signature),
    positive_label = paste0(positive_patients, "/", total_patients)
  )

C_range <- range(c(C$ci_low, C$ci_high), na.rm = TRUE)
C_pad <- diff(C_range) * 0.10
C_text_x <- max(C$ci_high, na.rm = TRUE) + C_pad
C_xlim <- c(min(C$ci_low, 0, na.rm = TRUE) - diff(C_range) * 0.05, C_text_x + diff(C_range) * 0.36)

pC <- ggplot(C, aes(colour = signature)) +
  geom_vline(xintercept = 0, colour = "#687682", linewidth = 0.55) +
  geom_segment(aes(x = ci_low, xend = ci_high, y = y, yend = y), linewidth = 0.9) +
  geom_point(aes(x = median_effect, y = y), size = 3.1) +
  geom_text(
    aes(x = C_text_x, y = y, label = positive_label),
    hjust = 0, size = 2.8, colour = COL_DARK, inherit.aes = FALSE
  ) +
  annotate(
    "text", x = C_text_x + diff(C_range) * 0.08, y = 4.52, label = "Positive patients",
    hjust = 1, size = 2.6, fontface = "bold", colour = COL_TEXT2
  ) +
  scale_colour_manual(values = sig_cols, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(rev(sig_order)),
    labels = rev(sig_labels[sig_order]),
    expand = expansion(mult = c(0.16, 0.16))
  ) +
  scale_x_continuous(limits = C_xlim, labels = label_number(accuracy = 0.1)) +
  labs(
    tag = "C",
    title = "P1 and literature-derived\nSPP1-TAM markers",
    subtitle = "GSE245552 patient-level comparison",
    x = "Median LM\u2212PT effect",
    y = NULL
  ) +
  theme_fig + panel_tag_theme +
  theme(
    plot.title = element_text(size = 8.7, face = "bold", hjust = 0, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.5, colour = COL_TEXT2, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(5, 32, 5, 6)
  )

## ------------------------------------------------------------
## 6. Panel D: ordered patient residual effect plot
## ------------------------------------------------------------

D <- D_dat %>%
  arrange(delta_residual) %>%
  mutate(
    patient = factor(patient, levels = patient),
    y = seq_len(n()) + 2,
    direction = if_else(delta_residual > 0, "Positive", "Negative")
  )

summary_med <- D_summary$median_delta[[1]]
summary_lo  <- D_summary$ci_low[[1]]
summary_hi  <- D_summary$ci_high[[1]]

summary_text <- paste0(
  D_summary$n_positive[[1]], "/", D_summary$n_patients[[1]],
  " positive   |   two-sided sign-flip P = ",
  sprintf("%.4f", D_summary$signflip_p_two_sided[[1]])
)

D_range <- range(c(D$delta_residual, summary_lo, summary_hi, 0), na.rm = TRUE)
D_pad <- diff(D_range) * 0.12
D_xlim <- c(min(D_range[1] - D_pad, -0.55), D_range[2] + D_pad)

summary_y <- 1.25

D_summary_dat <- tibble(
  y = summary_y,
  x0 = summary_lo,
  x1 = summary_hi,
  x = summary_med
)

pD <- ggplot(D, aes(x = delta_residual, y = y)) +
  geom_vline(xintercept = 0, colour = "#687682", linewidth = 0.55) +
  geom_segment(
    aes(x = 0, xend = delta_residual, yend = y, colour = direction),
    linewidth = 0.28, alpha = 0.48
  ) +
  geom_point(aes(colour = direction), size = 2.8) +
  geom_hline(yintercept = 2.40, colour = COL_GREY2, linewidth = 0.6) +
  geom_errorbarh(
    data = D_summary_dat,
    aes(y = y, xmin = x0, xmax = x1),
    inherit.aes = FALSE,
    height = 0.12,
    linewidth = 1.0,
    colour = COL_P1
  ) +
  geom_point(
    data = D_summary_dat,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    shape = 23,
    size = 3.6,
    fill = COL_P1,
    colour = COL_P1
  ) +
  annotate(
    "text",
    x = -0.02,
    y = summary_y,
    label = paste0(
      "Median \u0394 = ", sprintf("%.3f", summary_med),
      " (95% CI ", sprintf("%.3f", summary_lo),
      "\u2013", sprintf("%.3f", summary_hi), ")"
    ),
    hjust = 1,
    size = 2.60,
    fontface = "bold",
    colour = COL_P1
  ) +
  scale_colour_manual(
    values = c("Positive" = COL_P1, "Negative" = COL_BLUE),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = c(summary_y, D$y),
    labels = c("Median", as.character(D$patient)),
    limits = c(0.15, nrow(D) + 2.55),
    expand = c(0, 0)
  ) +
  scale_x_continuous(limits = D_xlim, labels = label_number(accuracy = 0.1)) +
  labs(
    x = "Residual LM\u2212PT difference in P1 score",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_fig +
  theme(
    axis.text.y = element_text(size = 7.5),
    plot.margin = margin(5, 6, 4, 6)
  )

pD_header <- ggplot() +
  annotate(
    "text", x = 0.00, y = 1.00, label = "D",
    hjust = 0, vjust = 1, size = 4.0, fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 0.035, y = 1.00,
    label = "Residual P1 signal after SPP1-TAM marker score\nand state adjustment",
    hjust = 0, vjust = 1, size = 3.55, lineheight = 0.92,
    fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 0.995, y = 0.32,
    label = summary_text,
    hjust = 1, vjust = 1, size = 2.25, colour = COL_TEXT2
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(5, 6, 4, 6))

pD_final <- pD_header / pD + plot_layout(heights = c(0.25, 1))

## ------------------------------------------------------------
## 7. Assemble and export
## ------------------------------------------------------------

fig2 <- pA + pB + pC + pD_final + plot_layout(
  design = "AAAA\nBBCC\nDDDD",
  heights = c(0.95, 1.10, 1.15)
)

ggsave(
  file.path(out_dir, "Figure4_CRLM_cross_state_SPP1.png"), fig2,
  width = 180, height = 225, units = "mm", dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure4_CRLM_cross_state_SPP1.pdf"), fig2,
  width = 180, height = 225, units = "mm", device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure4_CRLM_cross_state_SPP1.svg"), fig2,
  width = 180, height = 225, units = "mm", bg = "white"
)

## High-resolution TIFF for journals that require raster input.
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    file.path(out_dir, "Figure4_CRLM_cross_state_SPP1.tiff"), fig2,
    width = 180, height = 225, units = "mm", dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", background = "white"
  )
}

## Final numerical audit for this plotted version.
audit <- tibble(
  item = c(
    "A states", "B patients", "B patient-state cells", "C signatures",
    "D patients", "D median_delta", "D ci_low", "D ci_high",
    "D n_positive", "D n_patients", "D p_signflip"
  ),
  value = c(
    paste(state_order, collapse = " | "),
    length(patient_order), nrow(B_dat), paste(C$signature, collapse = " | "),
    nrow(D_dat), summary_med, summary_lo, summary_hi,
    D_summary$n_positive[[1]], D_summary$n_patients[[1]], D_summary$signflip_p_two_sided[[1]]
  )
)
utils::write.csv(audit, file.path(out_dir, "Figure4_source_data_audit.csv"), row.names = FALSE)

writeLines(c(
  "Figure 4 source-data provenance",
  "",
  "All numerical inputs originate from CRLM_final_frozen_20260912_R/tables.",
  "No model fitting, re-estimation, or new inferential analysis is performed.",
  "",
  "A: patient_state_mixed_model_state_contrasts_BH.csv.",
  "B: GSE245552_patient_state_P1_vs_published_SPP1_TAM_deltas.csv, score_version=P1_full.",
  "   Duplicate frozen rows within patient x state are collapsed by arithmetic mean for one heatmap tile.",
  "C: GSE245552_patient_level_P1_vs_published_SPP1_TAM_summary_primary.csv, pooled rows with missing state.",
  "   C is an independently reprocessed comparator framework; no common post-hoc z-standardization was applied.",
  "   P1/leave-out scores are weighted projections; the literature-derived SPP1-TAM marker score is an unweighted mean of SPP1, PPARG, ADM, MARCO, VEGFA, VCAN, CXCL8 and ANGPTL4.",
  "   Raw effect magnitudes across score definitions are therefore not interpreted as directly comparable.",
  "D: corrected state-only residual scores are paired LM vs PT within state and averaged within patient for display;",
  "   the residualization model itself is not rerun. Summary values come from the frozen summary table.",
  "",
  "D frozen summary target: median 0.42987174; CI 0.2490485\u20130.5328399; 11/13 positive; P=0.01659834.",
  "C-D positioning: prespecified comparator analyses using an independently reprocessed GSE245552 framework;",
  "   they do not redefine the frozen P1 program or primary cross-state inference."
), file.path(out_dir, "Figure4_PROVENANCE.txt"))

writeLines(c(
  "Figure 4 legend notes",
  "A, State-specific LM\u2212PT contrasts from the prespecified mixed model; model form: mean P1 score ~ tissue \u00d7 state + (1 | patient). The plotted units are 136 patient\u00d7tissue\u00d7state aggregated observation units from 16 patients. Points and horizontal bars show estimates and two-sided 95% confidence intervals; BH-FDR values are shown at right.",
  "B, Each tile is a patient\u00d7state LM\u2212PT difference in P1 score. Grey tiles indicate non-evaluable combinations and are not zero-filled.",
  "Panels C\u2013D present prespecified comparator analyses using an independently reprocessed GSE245552 framework and do not redefine the frozen P1 program or the primary cross-state inference.",
  "C, Patient-level LM\u2212PT effect summaries for the frozen P1 score, leave-out scores, and the literature-derived 8-gene SPP1-TAM marker score (SPP1, PPARG, ADM, MARCO, VEGFA, VCAN, CXCL8 and ANGPTL4). Scores in panel C are shown on their native scales and were not jointly standardized; comparisons therefore concern direction and patient-level consistency rather than the absolute magnitude of effects across score definitions.",
  "D, Patient-level residual LM\u2212PT differences after adjustment for the literature-derived 8-gene SPP1-TAM marker score (SPP1, PPARG, ADM, MARCO, VEGFA, VCAN, CXCL8 and ANGPTL4) and conventional myeloid state. Patients are ordered by residual LM\u2212PT difference. The bottom Median row shows the median difference and its 95% bootstrap confidence interval (B = 10,000); the P value is from a two-sided sign-flip test."
), file.path(out_dir, "Figure4_legend_notes.txt"), useBytes = TRUE)

message("Figure 4 rendered to: ", out_dir)
