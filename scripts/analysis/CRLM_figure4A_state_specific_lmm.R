## Reproduce Figure 4A state-specific LM-PT contrasts from frozen aggregate data
##
## This script is standalone and does not modify frozen outputs. It fits the
## prespecified REML model on patient x tissue x state aggregates, extracts
## state-specific LM-PT contrasts with emmeans, uses Satterthwaite inference,
## and applies BH correction across the five prespecified states.

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

if (!requireNamespace("lme4", quietly = TRUE) ||
    !requireNamespace("lmerTest", quietly = TRUE) ||
    !requireNamespace("emmeans", quietly = TRUE)) {
  stop("This script requires lme4, lmerTest, and emmeans.")
}

freeze_dir <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
output_root <- normalizePath(
  Sys.getenv("CRLM_OUTPUT_ROOT", unset = file.path(getwd(), "outputs")),
  winslash = "/", mustWork = FALSE
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(
  freeze_dir,
  "complete_run_archive_20260912",
  "run_output",
  "results",
  "tables",
  "48_GSE245552_patient_tissue_state_aggregated.csv"
)
input_file <- crlm_require_file(input_file, "frozen patient x tissue x state aggregate table")

dat <- utils::read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("score_version", "patient", "tissue", "state", "mean_score")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

dat <- dat[dat$score_version == "P1_full" & is.finite(dat$mean_score), , drop = FALSE]
tissue_key <- tolower(trimws(as.character(dat$tissue)))
lm_levels <- c("liver_metastasis", "lm", "metastasis")
pt_levels <- c("primary_tumor", "pt", "primary")
if (any(!tissue_key %in% c(lm_levels, pt_levels))) {
  stop("Unsupported tissue label in frozen aggregate table.")
}
dat$tissue_group <- factor(
  ifelse(tissue_key %in% lm_levels, "LM", "PT"),
  levels = c("PT", "LM")
)
dat$patient <- factor(dat$patient)
dat$state <- factor(dat$state)
dat <- dat[complete.cases(dat[, c("mean_score", "tissue_group", "patient", "state")]), , drop = FALSE]

fit <- lmerTest::lmer(
  mean_score ~ tissue_group * state + (1 | patient),
  data = dat,
  REML = TRUE
)

emm <- emmeans::emmeans(fit, ~ tissue_group | state)
con <- summary(
  emmeans::contrast(
    emm,
    method = list(LM_minus_PT = c(-1, 1)),
    by = "state",
    adjust = "none"
  ),
  infer = c(TRUE, TRUE)
)

out <- data.frame(
  state = as.character(con$state),
  estimate_LM_minus_PT = con$estimate,
  se = con$SE,
  ci_low = con$lower.CL,
  ci_high = con$upper.CL,
  df_formal = con$df,
  p_formal_two_sided = con$p.value,
  ci_low_formal = con$lower.CL,
  ci_high_formal = con$upper.CL,
  p_approx_two_sided = 2 * stats::pnorm(-abs(con$estimate / con$SE)),
  inference_method = "emmeans_lmerTest_Satterthwaite_same_fitted_model",
  n_units = nrow(dat),
  n_patients = length(unique(dat$patient)),
  stringsAsFactors = FALSE
)
out$BH_FDR <- stats::p.adjust(out$p_formal_two_sided, method = "BH")
out$BH_scope <- "five prespecified state-specific contrasts"

output_file <- file.path(output_root, "patient_state_mixed_model_state_contrasts_BH.csv")
utils::write.csv(out, output_file, row.names = FALSE, quote = TRUE)
message("Wrote: ", normalizePath(output_file, winslash = "/", mustWork = TRUE))
