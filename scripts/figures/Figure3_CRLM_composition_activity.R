## ============================================================
## CRLM Figure 3
## Composition versus within-state activity
##
## Plotting only from the frozen result package.
## No model fitting, re-estimation, or new inferential analysis.
## ============================================================

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(scales)

## ------------------------------------------------------------
## 0. Frozen inputs and output
## ------------------------------------------------------------

freeze_dir <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
table_dir  <- file.path(freeze_dir, "tables")
out_dir    <- file.path(CRLM_OUTPUT_ROOT, "Figure3_CRLM_composition_activity_v9")
source_dir <- file.path(out_dir, "source_data")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

read_frozen <- function(filename) {
  path <- file.path(table_dir, filename)
  if (!file.exists(path)) stop("Missing frozen input: ", path)
  message("Reading frozen table: ", path)
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

B_raw <- read_frozen("composition_activity_contribution_by_patient.csv")
C_raw <- read_frozen("composition_activity_contribution_summary.csv")

B_dat <- B_raw %>%
  filter(score_version == "P1_full") %>%
  distinct(patient, .keep_all = TRUE)

C_dat <- C_raw %>%
  filter(estimand == "ratio_of_sums_signed_contribution")

C_totals <- C_raw %>%
  filter(estimand == "observed_totals")

E_median <- C_raw %>%
  filter(estimand == "median_patient_ratio_signed_contribution") %>%
  mutate(
    component = recode(
      component,
      composition = "Composition",
      within_state_activity = "Within-state activity"
    ),
    component = factor(component, levels = c("Composition", "Within-state activity")),
    percent = estimate * 100
  )

required_B <- c(
  "patient", "overall_delta", "state_composition",
  "within_state_activity", "composition_fraction_signed",
  "within_state_fraction_signed"
)
required_C <- c("component", "estimate", "ci_low", "ci_high")
required_totals <- c("component", "estimate")

if (!all(required_B %in% names(B_dat))) {
  stop("Missing Figure 3 patient columns: ",
       paste(setdiff(required_B, names(B_dat)), collapse = ", "))
}
if (!all(required_C %in% names(C_dat))) {
  stop("Missing Figure 3 summary columns: ",
       paste(setdiff(required_C, names(C_dat)), collapse = ", "))
}
if (!all(required_totals %in% names(C_totals))) {
  stop("Missing Figure 3 observed-total columns: ",
       paste(setdiff(required_totals, names(C_totals)), collapse = ", "))
}

file.copy(
  file.path(table_dir, "composition_activity_contribution_by_patient.csv"),
  file.path(source_dir, "composition_activity_contribution_by_patient.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(table_dir, "composition_activity_contribution_summary.csv"),
  file.path(source_dir, "composition_activity_contribution_summary.csv"),
  overwrite = TRUE
)

## ------------------------------------------------------------
## 1. Figure contract
## ------------------------------------------------------------
## Core conclusion:
## Within-state P1 activity accounts for most of the aggregate
## metastasis-associated remodeling, with the same pattern visible
## across the patient-level signed decompositions.
##
## Archetype:
## Schematic-led quantitative composite with a quantitative hero panel.
##
## Palette:
## Composition = blue-grey; within-state activity = brick red.
## ------------------------------------------------------------

COL_ACTIVITY <- "#BE4D3B"
COL_COMP     <- "#6E8FA6"
COL_DARK     <- "#22282E"
COL_TEXT2    <- "#71808D"
COL_GREY2    <- "#D6DCE1"

theme_fig <- theme_classic(base_size = 8.5) +
  theme(
    text = element_text(family = "Arial", colour = COL_DARK),
    plot.title = element_text(
      face = "bold", size = 9.6, hjust = 0,
      margin = margin(b = 2)
    ),
    plot.subtitle = element_text(
      size = 7.2, colour = COL_TEXT2, hjust = 0,
      margin = margin(b = 6)
    ),
    axis.title = element_text(size = 8.2),
    axis.text = element_text(size = 7.2, colour = COL_DARK),
    axis.line = element_line(colour = COL_GREY2, linewidth = 0.42),
    axis.ticks = element_line(colour = COL_GREY2, linewidth = 0.36),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.8),
    panel.grid = element_blank(),
    plot.margin = margin(5, 6, 5, 6)
  )

panel_tag_theme <- theme(
  plot.tag = element_text(size = 10.5, face = "bold", colour = COL_DARK),
  plot.tag.position = c(0, 1)
)

## ------------------------------------------------------------
## 2. Panel A: decomposition logic
## ------------------------------------------------------------

pA <- ggplot() +
  annotate(
    "text", x = 1, y = 0.83,
    label = "\u0394P1 total",
    size = 3.0, fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 2, y = 0.83,
    label = "=",
    size = 3.35, fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 3.1, y = 0.83,
    label = "composition component",
    size = 2.85, fontface = "bold", colour = COL_COMP
  ) +
  annotate(
    "text", x = 4.2, y = 0.83,
    label = "+",
    size = 3.35, fontface = "bold", colour = COL_DARK
  ) +
  annotate(
    "text", x = 5.35, y = 0.83,
    label = "within-state activity component",
    size = 2.75, fontface = "bold", colour = COL_ACTIVITY
  ) +
  annotate(
    "text", x = 3.6, y = 0.45,
    label = "12 evaluable patients; 5 conventional states",
    size = 2.15, colour = COL_TEXT2
  ) +
  coord_cartesian(xlim = c(0.35, 6.25), ylim = c(0.32, 1.08)) +
  labs(tag = "A") +
  theme_void(base_family = "Arial") +
  panel_tag_theme +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.margin = margin(2, 6, 1, 6)
  )

## ------------------------------------------------------------
## 3. Panel B: signed patient-level decomposition
## ------------------------------------------------------------

patient_order_B <- B_dat %>%
  arrange(desc(overall_delta)) %>%
  pull(patient)

B <- B_dat %>%
  mutate(
    dominant_component = if_else(
      within_state_activity > state_composition,
      "Activity", "Composition"
    ),
    patient = factor(patient, levels = rev(patient_order_B))
  )

B_long <- B %>%
  select(patient, state_composition, within_state_activity) %>%
  pivot_longer(
    cols = c(state_composition, within_state_activity),
    names_to = "component",
    values_to = "contribution"
  ) %>%
  mutate(
    component = recode(
      component,
      state_composition = "Composition",
      within_state_activity = "Within-state activity"
    ),
    component = factor(
      component,
      levels = c("Composition", "Within-state activity")
    )
  )

B_total_max <- max(B$overall_delta, na.rm = TRUE)
B_separator_x <- B_total_max + 0.18
B_total_x <- B_total_max + 0.34
B_header_x <- B_total_x + 0.12
B_x_max <- B_total_max + 0.80

pB <- ggplot(
  B_long,
  aes(x = contribution, y = patient, fill = component)
) +
  geom_hline(
    yintercept = seq_len(nrow(B)), colour = "#EEF1F3",
    linewidth = 0.22
  ) +
  geom_vline(xintercept = 0, colour = "#687682", linewidth = 0.55) +
  geom_col(
    width = 0.54,
    position = position_stack(reverse = TRUE)
  ) +
  geom_point(
    data = B,
    aes(x = overall_delta, y = patient),
    inherit.aes = FALSE,
    shape = 21, size = 2.1, stroke = 0.45,
    fill = "white", colour = COL_DARK
  ) +
  geom_vline(
    xintercept = B_separator_x,
    colour = COL_GREY2, linewidth = 0.45
  ) +
  geom_text(
    data = B,
    aes(
      x = B_total_x,
      y = patient,
      label = sprintf("%.2f", overall_delta)
    ),
    inherit.aes = FALSE, hjust = 1, size = 2.15, colour = COL_DARK
  ) +
  annotate(
    "text", x = B_header_x,
    y = nrow(B) + 0.72, label = "Observed \u0394P1",
    hjust = 1, size = 2.05, fontface = "bold", colour = COL_TEXT2
  ) +
  scale_fill_manual(
    values = c(
      "Composition" = COL_COMP,
      "Within-state activity" = COL_ACTIVITY
    )
  ) +
  scale_x_continuous(
    limits = c(min(B_long$contribution, na.rm = TRUE) - 0.08,
               B_x_max),
    breaks = pretty(c(min(B_long$contribution, na.rm = TRUE),
                      max(B$overall_delta, na.rm = TRUE)), n = 5),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    ylim = c(0.5, nrow(B) + 1.05),
    clip = "off"
  ) +
  labs(
    tag = "B",
    title = "Patient-level decomposition",
    subtitle = "Additive components; open circles = observed total",
    x = "Contribution to LM\u2212PT P1 shift",
    y = NULL,
    fill = NULL
  ) +
  theme_fig +
  panel_tag_theme +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.x = element_text(size = 7.0),
    plot.margin = margin(5, 6, 5, 6)
  )

## ------------------------------------------------------------
## 4. Panel C: aggregate contribution hero bar
## ------------------------------------------------------------

component_order <- c("Composition", "Within-state activity")

C <- C_dat %>%
  mutate(
    component = recode(
      component,
      composition = "Composition",
      within_state_activity = "Within-state activity",
      within_state = "Within-state activity"
    ),
    component = factor(component, levels = component_order),
    percent = estimate * 100,
    ci_low_percent = ci_low * 100,
    ci_high_percent = ci_high * 100,
    label = paste0(
      sprintf("%.2f%%", percent),
      "\n[", sprintf("%.2f", ci_low_percent),
      "\u2013", sprintf("%.2f", ci_high_percent), "%]"
    )
  )

C_bar <- C %>%
  arrange(component) %>%
  mutate(
    start = lag(cumsum(percent), default = 0),
    end = cumsum(percent),
    midpoint = (start + end) / 2
  )

C_ci <- C %>%
  arrange(component) %>%
  mutate(
    track_y = c(0.66, 0.26),
    estimate_label = sprintf("%.2f%%", percent),
    ci_label = paste0(
      "95% CI ", sprintf("%.2f", ci_low_percent),
      "\u2013", sprintf("%.2f", ci_high_percent), "%"
    )
  )

C_total <- C_totals %>%
  mutate(
    component = recode(
      component,
      overall_delta = "Observed total",
      state_composition = "Composition",
      within_state_activity = "Within-state activity"
    )
  )

aggregate_total <- C_total$estimate[C_total$component == "Observed total"]
aggregate_composition <- C_total$estimate[C_total$component == "Composition"]
aggregate_activity <- C_total$estimate[C_total$component == "Within-state activity"]

C_abs <- tibble(
  component = factor(
    c("Composition", "Within-state activity"),
    levels = c("Composition", "Within-state activity")
  ),
  start = c(0, aggregate_composition),
  end = c(aggregate_composition, aggregate_total),
  midpoint = (start + end) / 2,
  absolute = c(aggregate_composition, aggregate_activity),
  percent = c(
    C$percent[C$component == "Composition"],
    C$percent[C$component == "Within-state activity"]
  )
)

pC_abs <- ggplot(C_abs) +
  geom_rect(
    aes(xmin = start, xmax = end, ymin = 0.76, ymax = 1.28, fill = component),
    colour = "white", linewidth = 0.7
  ) +
  geom_text(
    aes(x = midpoint, y = 1, label = paste0(
      sprintf("%.2f", absolute), "\n(", sprintf("%.2f%%", percent), ")"
    )),
    colour = "white", fontface = "bold", size = 2.35,
    lineheight = 0.9
  ) +
  annotate(
    "text", x = aggregate_total, y = 1.43,
    label = paste0("Summed patient-level P1 shift = ", sprintf("%.2f", aggregate_total)),
    hjust = 1, size = 2.35, fontface = "bold", colour = COL_DARK
  ) +
  scale_fill_manual(
    values = c(
      "Composition" = COL_COMP,
      "Within-state activity" = COL_ACTIVITY
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, aggregate_total * 1.08),
    breaks = c(0, 5, 10),
    labels = function(x) sprintf("%.0f", x),
    expand = c(0, 0)
  ) +
  labs(
    tag = "B",
    title = "Aggregate ratio-of-sums\ndecomposition",
    subtitle = "Across 12 evaluable patients",
    x = "P1 shift",
    y = NULL
  ) +
  theme_fig +
  panel_tag_theme +
  theme(
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    plot.margin = margin(4, 6, 1, 6),
    plot.title = element_text(size = 8.7, margin = margin(l = 8, b = 2))
  )

pC_ci_label <- ggplot(C_ci) +
  geom_text(
    aes(x = 1, y = track_y, label = component),
    hjust = 1, size = 2.30, fontface = "bold", colour = COL_TEXT2
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0.08, 0.84), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(1, 1, 3, 6)
  )

pC_ci_plot <- ggplot(C_ci) +
  geom_segment(
    aes(
      x = ci_low_percent, xend = ci_high_percent,
      y = track_y, yend = track_y, colour = component
    ), linewidth = 0.90
  ) +
  geom_segment(
    aes(
      x = ci_low_percent, xend = ci_low_percent,
      y = track_y - 0.065, yend = track_y + 0.065,
      colour = component
    ), linewidth = 0.65
  ) +
  geom_segment(
    aes(
      x = ci_high_percent, xend = ci_high_percent,
      y = track_y - 0.065, yend = track_y + 0.065,
      colour = component
    ), linewidth = 0.65
  ) +
  geom_point(aes(x = percent, y = track_y, colour = component), size = 2.20) +
  geom_text(
    aes(x = 105, y = track_y + 0.055, label = estimate_label, colour = component),
    hjust = 0, size = 2.60, fontface = "bold"
  ) +
  geom_text(
    aes(x = 105, y = track_y - 0.060, label = ci_label),
    hjust = 0, size = 2.50, colour = COL_TEXT2
  ) +
  scale_colour_manual(
    values = c(
      "Composition" = COL_COMP,
      "Within-state activity" = COL_ACTIVITY
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 200),
    breaks = c(0, 50, 100),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  labs(title = "Bootstrap uncertainty", x = NULL, y = NULL) +
  theme_fig +
  theme(
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(colour = COL_GREY2, linewidth = 0.4),
    plot.margin = margin(1, 6, 3, 6),
    plot.title = element_text(size = 6.3, colour = COL_TEXT2, margin = margin(b = 1))
  )

pC_ci <- pC_ci_label + pC_ci_plot + plot_layout(widths = c(1.05, 3.25))

pC <- pC_abs / pC_ci + plot_layout(heights = c(1.22, 0.90))

## ------------------------------------------------------------
## 5. Panel D: patient-level activity-dominance margin
## ------------------------------------------------------------

within_dominant_n <- B %>%
  summarise(n = sum(within_state_fraction_signed > composition_fraction_signed)) %>%
  pull(n)

D_margin <- B %>%
  transmute(
    patient,
    margin = within_state_activity - state_composition,
    positive = margin > 0
  ) %>%
  arrange(desc(margin)) %>%
  mutate(patient = factor(patient, levels = rev(patient)))

D_margin_min <- min(D_margin$margin)
D_margin_max <- max(D_margin$margin)
D_margin_lim <- c(
  floor(D_margin_min * 10) / 10 - 0.08,
  ceiling(D_margin_max * 10) / 10 + 0.25
)

pD <- ggplot(D_margin, aes(x = margin, y = patient, colour = positive)) +
  geom_vline(xintercept = 0, colour = "#687682", linewidth = 0.55) +
  geom_segment(
    aes(x = 0, xend = margin, y = patient, yend = patient),
    linewidth = 0.72, alpha = 0.20
  ) +
  geom_point(size = 2.9) +
  scale_colour_manual(
    values = c(`FALSE` = COL_COMP, `TRUE` = COL_ACTIVITY),
    guide = "none"
  ) +
  scale_x_continuous(
    limits = D_margin_lim,
    breaks = pretty(D_margin_lim, n = 5),
    labels = function(x) sprintf("%.1f", x),
    expand = c(0, 0)
  ) +
  labs(
    tag = "C",
    title = "Within-state activity exceeds composition in most patients",
    subtitle = "Patient-level component difference",
    x = "Within-state \u2212 composition contribution",
    y = NULL
  ) +
  theme_fig +
  panel_tag_theme +
  theme(
    axis.text.y = element_text(size = 7.1),
    plot.title = element_text(size = 7.6, margin = margin(b = 2)),
    plot.margin = margin(5, 6, 5, 6)
  )

## ------------------------------------------------------------
## 6. Panel E: patient-level consistency summary
## ------------------------------------------------------------

E_consistency <- tibble(
  metric = c(
    "Within-state activity > composition",
    "Within-state activity positive",
    "Composition positive"
  ),
  count = c(
    sum(B$within_state_activity > B$state_composition),
    sum(B$within_state_activity > 0),
    sum(B$state_composition > 0)
  ),
  component = c(
    "Within-state activity",
    "Within-state activity",
    "Composition"
  )
) %>%
  mutate(
    total = nrow(B),
    percent = count / total * 100,
    label = paste0(count, "/", total, " (", sprintf("%.1f%%", percent), ")"),
    metric = factor(metric, levels = rev(c(
      "Within-state activity > composition",
      "Within-state activity positive",
      "Composition positive"
    )))
  )

E_block <- E_consistency %>%
  mutate(
    y = c(3, 2, 1),
    metric_short = c(
      "activity-dominant",
      "activity-positive",
      "composition-positive"
    )
  )

pE <- ggplot(E_block) +
  geom_text(
    aes(x = 0.86, y = y, label = paste0(count, "/", total), colour = component),
    hjust = 1, size = 3.55, fontface = "bold"
  ) +
  geom_text(
    aes(x = 1.06, y = y, label = metric_short),
    hjust = 0, size = 2.40, colour = COL_TEXT2
  ) +
  scale_colour_manual(
    values = c(
      "Composition" = COL_COMP,
      "Within-state activity" = COL_ACTIVITY
    ),
    guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 2.65), ylim = c(0.5, 3.5), clip = "off") +
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = NULL
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.margin = margin(5, 1, 5, 1)
  )

## ------------------------------------------------------------
## 7. Assemble and export the three-panel main figure
## ------------------------------------------------------------

pA_patient <- pB +
  labs(
    tag = "A",
    title = "Patient-level decomposition of the LM\u2212PT P1 shift",
    subtitle = "Additive components; open circles denote observed LM\u2212PT P1 shift"
  ) +
  theme(
    plot.title = element_text(size = 8.6, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 6.8, margin = margin(b = 5))
  )

pB_aggregate <- pC

fig3 <- pA_patient + pB_aggregate + pD + pE +
  plot_layout(
    design = "AABB\nAABB\nAACD\nAACD",
    widths = c(1.3, 1.0, 0.95, 1.1),
    heights = c(1.0, 1.0, 1.0, 1.0)
  )

ggsave(
  file.path(out_dir, "Figure3_CRLM_composition_activity.png"),
  fig3, width = 180, height = 175, units = "mm", dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure3_CRLM_composition_activity.pdf"),
  fig3, width = 180, height = 175, units = "mm",
  device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure3_CRLM_composition_activity.svg"),
  fig3, width = 180, height = 175, units = "mm", bg = "white"
)

if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    file.path(out_dir, "Figure3_CRLM_composition_activity.tiff"),
    fig3, width = 180, height = 175, units = "mm", dpi = 600,
    device = ragg::agg_tiff, compression = "lzw", background = "white"
  )
}

audit <- tibble(
  item = c(
    "patient_count", "overall_delta_sum", "composition_sum",
    "within_state_activity_sum", "aggregate_composition_percent",
    "aggregate_composition_ci_low_percent",
    "aggregate_composition_ci_high_percent",
    "aggregate_within_state_percent",
    "aggregate_within_state_ci_low_percent",
    "aggregate_within_state_ci_high_percent",
    "within_state_dominant_patients",
    "median_patient_composition_percent",
    "median_patient_within_state_percent",
    "activity_dominance_margin_positive_patients",
    "within_state_activity_positive_patients"
  ),
  value = c(
    nrow(B),
    sum(B$overall_delta),
    sum(B$state_composition),
    sum(B$within_state_activity),
    C$percent[C$component == "Composition"],
    C$ci_low_percent[C$component == "Composition"],
    C$ci_high_percent[C$component == "Composition"],
    C$percent[C$component == "Within-state activity"],
    C$ci_low_percent[C$component == "Within-state activity"],
    C$ci_high_percent[C$component == "Within-state activity"],
    within_dominant_n,
    E_median$percent[E_median$component == "Composition"],
    E_median$percent[E_median$component == "Within-state activity"],
    sum(D_margin$positive),
    sum(B$within_state_activity > 0)
  )
)
utils::write.csv(
  audit,
  file.path(out_dir, "Figure3_source_data_audit.csv"),
  row.names = FALSE
)

writeLines(c(
  "Figure 3 legend notes",
  "",
  "All numerical values are read from the frozen composition/activity contribution tables.",
  "No model fitting, re-estimation, or new inferential analysis is performed.",
  "A, Total LM\u2212PT change in P1 is partitioned into a state-composition component and a within-state activity component across 12 evaluable patients and 5 conventional states.",
  "B, Bars show signed patient-level additive components; open circles show the observed total LM\u2212PT P1 shift for each patient. Patients are ordered by the observed overall shift. Right-hand annotations report the observed total shift and the larger component for each patient.",
  "C, The stacked bar shows the summed patient-level P1 shift (11.20) decomposed into composition (2.15; 19.16%) and within-state activity (9.05; 80.84%) across 12 evaluable patients. The lower tracks show frozen bootstrap 95% confidence intervals for the relative contribution estimates.",
  "The aggregate ratio-of-sums estimate is the prespecified aggregate estimand; patient-level components are shown as descriptive evidence of the decomposition."
), file.path(out_dir, "Figure3_legend_notes.txt"), useBytes = TRUE)

writeLines(c(
  "Figure 3 provenance",
  "",
  "Input package: CRLM_final_frozen_20260912_R.",
  "Input tables:",
  "  composition_activity_contribution_by_patient.csv",
  "  composition_activity_contribution_summary.csv",
  "",
  "Frozen aggregate estimates:",
  "  Composition: 19.16% [6.66-30.63%].",
  "  Within-state activity: 80.84% [68.75-93.53%].",
  "",
  "The Figure 3 script is plotting-only and does not alter the frozen analysis."
), file.path(out_dir, "Figure3_PROVENANCE.txt"))

message("Figure 3 rendered to: ", out_dir)
