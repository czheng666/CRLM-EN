## Build Supplementary Table S13 from the isolated 5-gene sensitivity output.
## This script does not modify any frozen Figure 2 or Supplementary package files.

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

options(stringsAsFactors = FALSE)

fivegene_root <- crlm_require_dir(CRLM_5GENE_ROOT, "CRLM_5GENE_ROOT")
source_csv <- file.path(fivegene_root, "tables", "8_vs_5_comparator_summary.csv")
out_dir <- fivegene_root
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_csv)) {
  stop("Authoritative 8-vs-5 comparator summary not found: ", source_csv)
}

dat <- read.csv(source_csv, check.names = FALSE, stringsAsFactors = FALSE)
expected <- c(
  "8-gene literature-derived marker score",
  "5-gene macrophage-specific SPP1-TAM score"
)
stopifnot(identical(dat$comparator, expected))
stopifnot(all(dat$discovery_n_cells == 12433L))
stopifnot(all(dat$discovery_spearman_rho >= 0 & dat$discovery_spearman_rho <= 1))
stopifnot(all(dat$residual_n_patients == 13L))
stopifnot(all(dat$residual_n_positive <= dat$residual_n_patients))

display <- data.frame(
  `Definition` = c(
    "Literature-derived 8-gene SPP1-TAM marker score",
    "Alternative macrophage-specific 5-gene SPP1-TAM marker score"
  ),
  `Included genes` = gsub(";", ", ", dat$genes, fixed = TRUE),
  `Score calculation` = c(
    "Unweighted mean of log1p(CP10K) expression; frozen primary comparator, original implementation retained",
    "Gene-wise z-score across GSE245552 cells, followed by an unweighted mean z-score"
  ),
  `Discovery-cohort correlation with P1 (Spearman ρ; n = 12,433 cells)` = sprintf("%.6f", dat$discovery_spearman_rho),
  `Median residual LM−PT difference` = sprintf("%.6f", dat$residual_median_delta),
  `95% CI` = sprintf("%.6f-%.6f", dat$residual_ci_low, dat$residual_ci_high),
  `Positive patients` = sprintf("%d/%d", dat$residual_n_positive, dat$residual_n_patients),
  `Sign-flip P` = sprintf("%.6f", dat$residual_signflip_p_two_sided),
  check.names = FALSE
)

out_csv <- file.path(out_dir, "Supplementary_Table_S13_alternative_SPP1_TAM_marker_definitions.csv")
write.csv(display, out_csv, row.names = FALSE, fileEncoding = "UTF-8")

## Keep exact numeric values and provenance alongside the reader-facing table.
exact_csv <- file.path(out_dir, "Supplementary_Table_S13_exact_source_values.csv")
write.csv(dat, exact_csv, row.names = FALSE, fileEncoding = "UTF-8")

note <- c(
  "Supplementary Table S13 presents a sensitivity analysis using alternative SPP1-TAM marker definitions.",
  "",
  "The eight-gene row is the frozen primary comparator, comprising SPP1, PPARG, ADM, MARCO, VEGFA, VCAN, CXCL8, and ANGPTL4; scores were calculated as the unweighted mean of log1p(CP10K) expression without post hoc standardization.",
  "The five-gene row is the alternative macrophage-specific definition, comprising PLAUR, SLC11A1, BRI3, FBP1, and C15orf48; each gene was z-standardized across GSE245552 cells before taking the unweighted mean z-score.",
  "Discovery-cohort correlations are descriptive results from 12,433 GSE178318 discovery-cohort myeloid cells. For each comparator, P1_full ~ comparator + conventional myeloid state was fitted on all 164 patient × tissue × state units, followed by within-state pairing of primary tumour and liver metastasis samples and averaging across available states within each patient.",
  "Confidence intervals used 10,000 patient bootstrap resamples; directionality used 10,000 two-sided patient-level Monte Carlo sign-flip resamples with a plus-one empirical P-value correction. Values are rounded to six decimal places; exact source values are stored in Supplementary_Table_S13_exact_source_values.csv.",
  "The five-gene branch is an alternative-definition sensitivity analysis and does not replace the frozen eight-gene comparator."
)
writeLines(note, file.path(out_dir, "Supplementary_Table_S13_note.txt"), useBytes = TRUE)

cat("Wrote:\n", out_csv, "\n", exact_csv, "\n", sep = "")
cat(paste(capture.output(print(display, row.names = FALSE)), collapse = "\n"), "\n")
