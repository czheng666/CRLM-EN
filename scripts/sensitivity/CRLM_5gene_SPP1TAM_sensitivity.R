###############################################################################

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))
# Alternative SPP1-TAM definition sensitivity analysis
#
# Purpose:
#   Evaluate the van Buuren et al. 2025 macrophage-specific five-gene SPP1-TAM
#   comparator without replacing the frozen eight-gene literature-derived
#   comparator or modifying any frozen analysis artifact.
#
# Scope:
#   1) Discovery-cell P1 versus five-gene descriptive Spearman association.
#   2) GSE245552 patient x tissue x state five-gene score aggregation.
#   3) State-only residualization: P1_full ~ comparator + state.
#   4) Patient-level LM-PT residual summaries using the frozen mean-across-
#      available-state rule.
#
# The primary five-gene score follows the cited paper's score concept:
# gene-wise z-standardization followed by the unweighted mean of the five
# standardized genes. For this sensitivity branch, z-standardization is
# performed within each dataset across all retained cells. The direct
# log1p(CP10K) mean is retained only as an auxiliary diagnostic.
###############################################################################

options(stringsAsFactors = FALSE, width = 200)

required_packages <- c("Seurat", "Matrix", "data.table", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

## ------------------------------------------------------------------------
## Frozen inputs and isolated output contract
## ------------------------------------------------------------------------

PRIMARY <- crlm_require_dir(CRLM_PRIMARY_RERUN_ROOT, "CRLM_PRIMARY_RERUN_ROOT")
FREEZE <- crlm_require_dir(CRLM_FREEZE_DIR, "CRLM_FREEZE_DIR")
RAW_GSE245552 <- crlm_require_dir(CRLM_RAW_GSE245552_DIR, "CRLM_RAW_GSE245552_DIR")

OUT <- if (nzchar(CRLM_5GENE_ROOT)) {
  CRLM_5GENE_ROOT
} else {
  file.path(CRLM_OUTPUT_ROOT, "supplementary_extension_5gene_SPP1TAM_20260915")
}
TABLES <- file.path(OUT, "tables")
AUDIT <- file.path(OUT, "audit")
MANIFEST <- file.path(OUT, "manifest")
dir.create(TABLES, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT, recursive = TRUE, showWarnings = FALSE)
dir.create(MANIFEST, recursive = TRUE, showWarnings = FALSE)

DISCOVERY_RDS <- file.path(PRIMARY, "cache/annotation/fine_myeloid_dc_annotated.rds")
FROZEN_P1_RDS <- file.path(PRIMARY, "cache/projection/frozen_P1_v4.rds")
GSE245_CELL_SCORE <- file.path(
  FREEZE, "tables/GSE245552_cell_level_P1_vs_published_SPP1_TAM.csv.gz"
)
GSE245_SAMPLE_SCORE_FROZEN <- file.path(
  FREEZE, "tables/GSE245552_patient_tissue_state_P1_vs_published_SPP1_TAM.csv"
)
GSE245_RESIDUAL_SUMMARY_FROZEN <- file.path(
  FREEZE, "tables/GSE245552_patient_level_SPP1_residual_comparison_corrected_state_only.csv"
)
DISCOVERY_8GENE_CORR <- file.path(
  crlm_require_dir(CRLM_SUPPLEMENTARY_RESULTS_ROOT, "CRLM_SUPPLEMENTARY_RESULTS_ROOT"),
  "tables",
  "published_signature_discovery_cell_correlations.csv"
)

required_files <- c(
  DISCOVERY_RDS, FROZEN_P1_RDS, GSE245_CELL_SCORE,
  GSE245_SAMPLE_SCORE_FROZEN, GSE245_RESIDUAL_SUMMARY_FROZEN,
  DISCOVERY_8GENE_CORR
)
if (!all(file.exists(required_files))) {
  stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}
if (!dir.exists(RAW_GSE245552)) stop("Missing GSE245552 raw directory: ", RAW_GSE245552)

FIVE_GENES <- c("PLAUR", "SLC11A1", "BRI3", "FBP1", "C15orf48")
OLD8_GENES <- c("SPP1", "PPARG", "ADM", "MARCO", "VEGFA", "VCAN", "CXCL8", "ANGPTL4")

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

write_csv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wt")
  utils::write.csv(x, con, row.names = FALSE, na = "")
  close(con)
}

sha256 <- function(path) digest::digest(file = path, algo = "sha256")

assert_columns <- function(x, required, label) {
  miss <- setdiff(required, names(x))
  if (length(miss) > 0L) stop(label, " missing columns: ", paste(miss, collapse = ", "))
}

assert_near <- function(x, y, tolerance, label) {
  if (length(x) != 1L || !is.finite(x) || abs(as.numeric(x) - as.numeric(y)) > tolerance) {
    stop(label, " mismatch: observed=", x, "; expected=", y)
  }
}

gene_z_score <- function(expr_genes_by_cells) {
  expr_genes_by_cells <- as.matrix(expr_genes_by_cells)
  z <- do.call(rbind, lapply(seq_len(nrow(expr_genes_by_cells)), function(i) {
    v <- as.numeric(expr_genes_by_cells[i, ])
    mu <- mean(v, na.rm = TRUE)
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, length(v)) else (v - mu) / s
  }))
  rownames(z) <- rownames(expr_genes_by_cells)
  z
}

mean_z_score <- function(expr_genes_by_cells) {
  z <- gene_z_score(expr_genes_by_cells)
  as.numeric(colMeans(z, na.rm = TRUE))
}

collapse_patient_tissue_state <- function(x, score_columns) {
  key <- interaction(x$patient, x$tissue, x$state, drop = TRUE)
  rows <- lapply(split(x, key), function(dd) {
    out <- data.frame(
      patient = as.character(dd$patient[[1]]),
      tissue = as.character(dd$tissue[[1]]),
      state = as.character(dd$state[[1]]),
      n_cells = sum(as.numeric(dd$n_cells), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (nm in score_columns) {
      out[[nm]] <- sum(as.numeric(dd[[nm]]) * as.numeric(dd$n_cells), na.rm = TRUE) /
        sum(as.numeric(dd$n_cells), na.rm = TRUE)
    }
    out
  })
  do.call(rbind, rows)
}

collapse_sample_tissue_state <- function(x, score_columns) {
  key <- interaction(x$sample, x$patient, x$tissue, x$state, drop = TRUE)
  rows <- lapply(split(x, key), function(dd) {
    out <- data.frame(
      sample = as.character(dd$sample[[1]]),
      patient = as.character(dd$patient[[1]]),
      tissue = as.character(dd$tissue[[1]]),
      state = as.character(dd$state[[1]]),
      n_cells = sum(as.numeric(dd$n_cells), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    for (nm in score_columns) {
      out[[nm]] <- sum(as.numeric(dd[[nm]]) * as.numeric(dd$n_cells), na.rm = TRUE) /
        sum(as.numeric(dd$n_cells), na.rm = TRUE)
    }
    out
  })
  do.call(rbind, rows)
}

## ------------------------------------------------------------------------
## 1. Discovery-cell five-gene score
## ------------------------------------------------------------------------

message("Reading frozen discovery objects")
discovery_obj <- readRDS(DISCOVERY_RDS)
frozen_p1 <- readRDS(FROZEN_P1_RDS)
if (is.null(frozen_p1$H) || !"P1" %in% rownames(frozen_p1$H)) {
  stop("Frozen P1 object does not contain H['P1', ].")
}

discovery_cells <- intersect(colnames(discovery_obj), colnames(frozen_p1$H))
if (length(discovery_cells) != 12433L) {
  stop("Discovery cell count mismatch: observed ", length(discovery_cells), "; expected 12433.")
}

discovery_expr <- Seurat::GetAssayData(
  discovery_obj, assay = "RNA", layer = "data"
)[, discovery_cells, drop = FALSE]
if (!all(FIVE_GENES %in% rownames(discovery_expr))) {
  stop("Discovery expression is missing five-gene member(s): ",
       paste(setdiff(FIVE_GENES, rownames(discovery_expr)), collapse = ", "))
}

discovery_five_expr <- discovery_expr[FIVE_GENES, discovery_cells, drop = FALSE]
discovery_five_z_matrix <- gene_z_score(discovery_five_expr)
discovery_five_z <- as.numeric(colMeans(discovery_five_z_matrix))
discovery_p1 <- as.numeric(frozen_p1$H["P1", discovery_cells])
discovery_ok <- is.finite(discovery_p1) & is.finite(discovery_five_z)

discovery_scores <- data.frame(
  cell_id = discovery_cells,
  P1_frozen_usage = discovery_p1,
  SPP1_TAM_5gene_macrophage_specific_z = discovery_five_z,
  stringsAsFactors = FALSE
)
for (g in FIVE_GENES) {
  discovery_scores[[paste0(g, "_z")]] <- discovery_five_z_matrix[g, ]
}
write_csv_gz(discovery_scores, file.path(TABLES, "discovery_5gene_cell_scores.csv.gz"))

discovery_rho5 <- suppressWarnings(stats::cor(
  discovery_p1[discovery_ok], discovery_five_z[discovery_ok], method = "spearman"
))
discovery_5gene_summary <- data.frame(
  comparator = "SPP1_TAM_5gene_macrophage_specific",
  genes = paste(FIVE_GENES, collapse = ";"),
  n_cells = sum(discovery_ok),
  spearman_rho_with_frozen_P1 = discovery_rho5,
  score_definition = "gene-wise z-standardization across discovery cells, followed by unweighted mean z-score",
  stringsAsFactors = FALSE
)
write_csv(discovery_5gene_summary, file.path(TABLES, "discovery_5gene_summary.csv"))

## ------------------------------------------------------------------------
## 2. GSE245552 GEO-supplied 10x-format extraction and five-gene score
## ------------------------------------------------------------------------

message("Reading authoritative GSE245552 cell-level score table")
external_score <- data.table::fread(GSE245_CELL_SCORE, data.table = FALSE, check.names = FALSE)
assert_columns(
  external_score,
  c("sample", "patient", "tissue", "state", "barcode", "P1_full",
    "SPP1_TAM_published_2025"),
  "GSE245552 cell-level score table"
)
external_score$sample <- as.character(external_score$sample)
external_score$patient <- as.character(external_score$patient)
external_score$tissue <- as.character(external_score$tissue)
external_score$state <- as.character(external_score$state)
external_score$barcode <- as.character(external_score$barcode)
external_score$P1_full <- as.numeric(external_score$P1_full)
external_score$SPP1_TAM_published_2025 <- as.numeric(external_score$SPP1_TAM_published_2025)

matrix_files <- sort(list.files(
  RAW_GSE245552, pattern = "_matrix\\.mtx\\.gz$", full.names = TRUE
))
if (length(matrix_files) == 0L) stop("No GSE245552 10x matrix files found.")

sample_from_matrix <- function(path) {
  stem <- sub("_matrix\\.mtx\\.gz$", "", basename(path))
  sub("^GSM[0-9]+_", "", stem)
}

raw_samples <- vapply(matrix_files, sample_from_matrix, character(1))
score_samples <- sort(unique(external_score$sample))
if (!setequal(raw_samples, score_samples)) {
  stop(
    "Raw/sample mismatch. Missing raw samples: ",
    paste(setdiff(score_samples, raw_samples), collapse = ", "),
    "; raw samples absent from score table: ",
    paste(setdiff(raw_samples, score_samples), collapse = ", ")
  )
}

score_by_sample <- split(external_score, external_score$sample)
meta_parts <- vector("list", length(matrix_files))
five_expr_parts <- vector("list", length(matrix_files))
eight_expr_parts <- vector("list", length(matrix_files))
old8_recomputed_parts <- vector("list", length(matrix_files))

for (ii in seq_along(matrix_files)) {
  matrix_file <- matrix_files[[ii]]
  sample_id <- raw_samples[[ii]]
  base <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
  barcode_file <- file.path(RAW_GSE245552, paste0(base, "_barcodes.tsv.gz"))
  feature_file <- file.path(RAW_GSE245552, paste0(base, "_features.tsv.gz"))
  if (!file.exists(barcode_file) || !file.exists(feature_file)) {
    stop("Missing barcode/features file for sample ", sample_id)
  }

  barcodes <- read.delim(gzfile(barcode_file), header = FALSE, stringsAsFactors = FALSE)[[1]]
  features <- read.delim(gzfile(feature_file), header = FALSE, stringsAsFactors = FALSE)
  if (ncol(features) < 2L) stop("Feature file has fewer than two columns: ", feature_file)
  genes <- make.unique(as.character(features[[2]]))
  if (length(genes) != nrow(features)) stop("Feature/matrix row mismatch in ", sample_id)

  sample_score <- score_by_sample[[sample_id]]
  raw_idx <- match(sample_score$barcode, barcodes)
  if (anyNA(raw_idx)) {
    stop("Not all authoritative score-table barcodes match raw data for ", sample_id,
         "; missing=", sum(is.na(raw_idx)))
  }
  if (anyDuplicated(raw_idx)) stop("Duplicate raw barcode match in ", sample_id)

  con <- gzfile(matrix_file, open = "rt")
  counts <- Matrix::readMM(con)
  close(con)
  counts <- as(counts, "dgCMatrix")
  if (nrow(counts) != length(genes)) stop("Matrix/features dimension mismatch in ", sample_id)
  counts <- counts[, raw_idx, drop = FALSE]

  five_idx <- match(FIVE_GENES, genes)
  old8_idx <- match(OLD8_GENES, genes)
  if (anyNA(five_idx)) stop("Five-gene member missing from raw matrix ", sample_id, ": ", paste(FIVE_GENES[is.na(five_idx)], collapse = ", "))
  if (anyNA(old8_idx)) stop("Eight-gene member missing from raw matrix ", sample_id, ": ", paste(OLD8_GENES[is.na(old8_idx)], collapse = ", "))

  library_size <- Matrix::colSums(counts)
  if (any(!is.finite(library_size) | library_size <= 0)) {
    stop("Non-positive library size in authoritative cells for ", sample_id)
  }

  five_counts <- as.matrix(counts[five_idx, , drop = FALSE])
  old8_counts <- as.matrix(counts[old8_idx, , drop = FALSE])
  five_expr <- log1p(sweep(five_counts, 2, library_size, "/") * 1e4)
  old8_expr <- log1p(sweep(old8_counts, 2, library_size, "/") * 1e4)

  meta_parts[[ii]] <- sample_score
  five_expr_parts[[ii]] <- five_expr
  eight_expr_parts[[ii]] <- old8_expr
  old8_recomputed_parts[[ii]] <- colMeans(old8_expr)

  rm(counts, five_counts, old8_counts, five_expr, old8_expr, library_size,
     raw_idx, barcodes, features, genes, sample_score)
  gc(verbose = FALSE)
}

external_meta <- do.call(rbind, meta_parts)
external_five_expr <- do.call(cbind, five_expr_parts)
external_old8_expr <- do.call(cbind, eight_expr_parts)
rownames(external_five_expr) <- FIVE_GENES
rownames(external_old8_expr) <- OLD8_GENES
external_old8_recomputed <- unlist(old8_recomputed_parts, use.names = FALSE)

if (nrow(external_meta) != 24231L) {
  stop("GSE245552 authoritative cell count mismatch: observed ", nrow(external_meta), "; expected 24231.")
}
if (ncol(external_five_expr) != nrow(external_meta) || ncol(external_old8_expr) != nrow(external_meta)) {
  stop("External expression/metadata cell order mismatch.")
}

old8_authoritative <- as.numeric(external_meta$SPP1_TAM_published_2025)
old8_delta <- external_old8_recomputed - old8_authoritative
old8_ok <- is.finite(old8_delta)
if (!all(old8_ok) || max(abs(old8_delta[old8_ok])) > 1e-8) {
  stop(
    "Raw-matrix recomputation of the frozen eight-gene score failed: max_abs_difference=",
    max(abs(old8_delta[old8_ok]), na.rm = TRUE)
  )
}

external_five_z_matrix <- gene_z_score(external_five_expr)
external_five_z <- as.numeric(colMeans(external_five_z_matrix))
external_five_log1p_mean <- as.numeric(colMeans(external_five_expr))

external_scores <- external_meta[, c("sample", "patient", "tissue", "state", "barcode"), drop = FALSE]
external_scores$P1_full <- as.numeric(external_meta$P1_full)
external_scores$SPP1_TAM_8gene_frozen <- old8_authoritative
external_scores$SPP1_TAM_8gene_recomputed <- external_old8_recomputed
external_scores$SPP1_TAM_5gene_log1p_mean_auxiliary <- external_five_log1p_mean
external_scores$SPP1_TAM_5gene_macrophage_specific_z <- external_five_z
for (g in FIVE_GENES) {
  external_scores[[paste0(g, "_z")]] <- external_five_z_matrix[g, ]
}
write_csv_gz(external_scores, file.path(TABLES, "GSE245552_5gene_cell_scores.csv.gz"))

## ------------------------------------------------------------------------
## 3. Patient x tissue x state aggregation and frozen-table reconciliation
## ------------------------------------------------------------------------

## Reconstruct all frozen sample x tissue x state rows first, including
## adjacent tissues.  The frozen sample-level score table contains these rows,
## so they must be included in the reconciliation gate.  The LM-PT estimand
## is restricted to primary_tumor and liver_metastasis only after this gate,
## in the residualization branch below.
external_for_aggregation <- external_scores[, c(
  "sample", "patient", "tissue", "state", "P1_full",
  "SPP1_TAM_8gene_frozen", "SPP1_TAM_5gene_macrophage_specific_z"
), drop = FALSE]
external_for_aggregation$n_cells <- 1L
external_sample_agg <- collapse_sample_tissue_state(
  external_for_aggregation,
  c("P1_full", "SPP1_TAM_8gene_frozen", "SPP1_TAM_5gene_macrophage_specific_z")
)
## Verify that cell-level reconstruction reproduces the frozen sample-level
## P1 and eight-gene comparator table used by the corrected S3 residual branch.
frozen_sample <- utils::read.csv(GSE245_SAMPLE_SCORE_FROZEN, check.names = FALSE, stringsAsFactors = FALSE)
assert_columns(frozen_sample, c("patient", "tissue", "state", "score_version", "n_cells", "mean_score"), "frozen sample score table")
p1_sample <- frozen_sample[frozen_sample$score_version == "P1_full",
                           c("sample", "patient", "tissue", "state", "n_cells", "mean_score"), drop = FALSE]
names(p1_sample)[names(p1_sample) == "mean_score"] <- "P1_frozen_mean_score"
names(p1_sample)[names(p1_sample) == "n_cells"] <- "P1_frozen_n_cells"
old8_sample <- frozen_sample[frozen_sample$score_version == "SPP1_TAM_published_2025",
                             c("sample", "patient", "tissue", "state", "n_cells", "mean_score"), drop = FALSE]
names(old8_sample)[names(old8_sample) == "mean_score"] <- "old8_frozen_mean_score"
names(old8_sample)[names(old8_sample) == "n_cells"] <- "old8_frozen_n_cells"
recon_sample <- merge(external_sample_agg, p1_sample,
                       by = c("sample", "patient", "tissue", "state"), all = TRUE, sort = FALSE)
recon_sample <- merge(recon_sample, old8_sample,
                      by = c("sample", "patient", "tissue", "state"), all = TRUE, sort = FALSE)
if (anyNA(recon_sample$P1_full) || anyNA(recon_sample$P1_frozen_mean_score) ||
    anyNA(recon_sample$SPP1_TAM_8gene_frozen) || anyNA(recon_sample$old8_frozen_mean_score)) {
  stop("Sample-level reconciliation has unmatched patient/tissue/state rows.")
}
recon_sample$P1_abs_difference <- abs(recon_sample$P1_full - recon_sample$P1_frozen_mean_score)
recon_sample$old8_abs_difference <- abs(recon_sample$SPP1_TAM_8gene_frozen - recon_sample$old8_frozen_mean_score)
if (max(recon_sample$P1_abs_difference, na.rm = TRUE) > 1e-8 ||
    max(recon_sample$old8_abs_difference, na.rm = TRUE) > 1e-8) {
  stop("Cell-level reconstruction does not reproduce frozen sample-level P1/8-gene table; max_abs_difference=",
       max(c(recon_sample$P1_abs_difference, recon_sample$old8_abs_difference), na.rm = TRUE))
}
write_csv(recon_sample, file.path(AUDIT, "GSE245552_sample_level_reconciliation.csv"))

## Collapse the reconciled sample/state rows to the patient/tissue/state unit
## used by the corrected S3 residualization.
external_agg <- collapse_patient_tissue_state(
  external_sample_agg[, c("patient", "tissue", "state", "n_cells", "P1_full",
                          "SPP1_TAM_8gene_frozen", "SPP1_TAM_5gene_macrophage_specific_z"), drop = FALSE],
  c("P1_full", "SPP1_TAM_8gene_frozen", "SPP1_TAM_5gene_macrophage_specific_z")
)
write_csv(external_agg, file.path(TABLES, "GSE245552_5gene_patient_tissue_state_aggregated.csv"))

## Descriptive LM-PT pattern of the alternative five-gene score itself.
## This is a raw patient x state contrast and is not used as a new inferential
## test; the adjusted patient-level result is calculated below.
lmpt_agg <- external_agg[
  external_agg$tissue %in% c("primary_tumor", "liver_metastasis"), , drop = FALSE
]
five_pt <- lmpt_agg[lmpt_agg$tissue == "primary_tumor",
                    c("patient", "state", "SPP1_TAM_5gene_macrophage_specific_z"), drop = FALSE]
five_lm <- lmpt_agg[lmpt_agg$tissue == "liver_metastasis",
                    c("patient", "state", "SPP1_TAM_5gene_macrophage_specific_z"), drop = FALSE]
five_lmpt <- merge(
  five_pt, five_lm,
  by = c("patient", "state"),
  suffixes = c("_PT", "_LM"),
  sort = FALSE
)
five_lmpt$delta_LM_minus_PT <- five_lmpt$SPP1_TAM_5gene_macrophage_specific_z_LM -
  five_lmpt$SPP1_TAM_5gene_macrophage_specific_z_PT
five_lmpt <- five_lmpt[order(five_lmpt$state, five_lmpt$patient), , drop = FALSE]
write_csv(five_lmpt, file.path(TABLES, "GSE245552_5gene_raw_LM_PT_patient_state_deltas.csv"))

five_state_pattern <- do.call(rbind, lapply(split(five_lmpt, five_lmpt$state), function(dd) {
  data.frame(
    state = as.character(dd$state[[1]]),
    n_pairs = nrow(dd),
    median_delta_LM_minus_PT = stats::median(dd$delta_LM_minus_PT),
    mean_delta_LM_minus_PT = mean(dd$delta_LM_minus_PT),
    n_positive = sum(dd$delta_LM_minus_PT > 0),
    positive_fraction = sum(dd$delta_LM_minus_PT > 0) / nrow(dd),
    stringsAsFactors = FALSE
  )
}))
rownames(five_state_pattern) <- NULL
write_csv(five_state_pattern, file.path(TABLES, "GSE245552_5gene_raw_LM_PT_state_pattern.csv"))

## ------------------------------------------------------------------------
## 4. State-only residualization and patient-level summaries
## ------------------------------------------------------------------------

long <- external_agg[, c(
  "patient", "tissue", "state", "P1_full",
  "SPP1_TAM_8gene_frozen", "SPP1_TAM_5gene_macrophage_specific_z"
), drop = FALSE]
## Fit the corrected state-only model on all 164 frozen patient x tissue x
## state units.  Tissue is retained only when forming the final LM-PT pairs.
long <- long[is.finite(long$P1_full) &
             is.finite(long$SPP1_TAM_8gene_frozen) &
             is.finite(long$SPP1_TAM_5gene_macrophage_specific_z), , drop = FALSE]

patient_summary <- function(long_data, comparator, comparator_label, seed) {
  model_data <- long_data[, c("patient", "tissue", "state", "P1_full", comparator), drop = FALSE]
  model_data <- model_data[complete.cases(model_data), , drop = FALSE]
  formula_text <- stats::as.formula(paste0("P1_full ~ ", comparator, " + state"))
  fit <- stats::lm(formula_text, data = model_data)
  model_data$residual <- stats::residuals(fit)

  contrast_data <- model_data[
    model_data$tissue %in% c("primary_tumor", "liver_metastasis"), , drop = FALSE
  ]
  pt <- contrast_data[contrast_data$tissue == "primary_tumor", c("patient", "state", "residual"), drop = FALSE]
  lm <- contrast_data[contrast_data$tissue == "liver_metastasis", c("patient", "state", "residual"), drop = FALSE]
  pair <- merge(pt, lm, by = c("patient", "state"), suffixes = c("_PT", "_LM"), sort = FALSE)
  pair$delta_LM_minus_PT <- pair$residual_LM - pair$residual_PT
  pair <- pair[order(pair$patient, pair$state), , drop = FALSE]

  if (nrow(pair) == 0L) stop("No paired residuals for ", comparator_label)
  patient_delta <- stats::aggregate(delta_LM_minus_PT ~ patient, pair, mean)
  x <- as.numeric(patient_delta$delta_LM_minus_PT)
  set.seed(seed)
  b_boot <- 10000L
  med_boot <- vapply(seq_len(b_boot), function(i) median(sample(x, length(x), replace = TRUE)), numeric(1))
  null <- vapply(seq_len(b_boot), function(i) median(sample(c(-1, 1), length(x), replace = TRUE) * x), numeric(1))

  summary <- data.frame(
    comparator = comparator_label,
    comparator_column = comparator,
    model = paste0("P1_full ~ ", comparator, " + state"),
    n_pairs = nrow(pair),
    n_patients = length(x),
    median_delta = median(x),
    ci_low = unname(stats::quantile(med_boot, 0.025)),
    ci_high = unname(stats::quantile(med_boot, 0.975)),
    n_positive = sum(x > 0),
    signflip_p_two_sided = (1 + sum(abs(null) >= abs(median(x)))) / (b_boot + 1L),
    bootstrap_B = b_boot,
    seed = seed,
    pooling = "mean across available state-specific LM-PT residual deltas within patient",
    stringsAsFactors = FALSE
  )
  list(summary = summary, pair = pair, patient_delta = patient_delta, fit = fit)
}

res8 <- patient_summary(
  long, "SPP1_TAM_8gene_frozen", "8-gene literature-derived marker score", 20260969L
)
res5 <- patient_summary(
  long, "SPP1_TAM_5gene_macrophage_specific_z", "5-gene macrophage-specific SPP1-TAM score", 20260970L
)

write_csv(res8$pair, file.path(TABLES, "GSE245552_8gene_residual_patient_state_deltas.csv"))
write_csv(res5$pair, file.path(TABLES, "GSE245552_5gene_residual_patient_state_deltas.csv"))
write_csv(res8$patient_delta, file.path(TABLES, "GSE245552_8gene_residual_patient_deltas.csv"))
write_csv(res5$patient_delta, file.path(TABLES, "GSE245552_5gene_residual_patient_deltas.csv"))
write_csv(rbind(res8$summary, res5$summary), file.path(TABLES, "GSE245552_8_vs_5_residual_summary.csv"))

## Reproduce the frozen eight-gene corrected result exactly before accepting the
## alternative branch. This guards against accidental changes to the patient
## aggregation rule or residual model.
frozen_residual_summary <- utils::read.csv(GSE245_RESIDUAL_SUMMARY_FROZEN, check.names = FALSE, stringsAsFactors = FALSE)
frozen_row <- frozen_residual_summary[
  grepl("P1_residual_state_only", frozen_residual_summary$residual_definition, fixed = TRUE), , drop = FALSE
]
if (nrow(frozen_row) != 1L) stop("Frozen eight-gene residual summary is not uniquely identified.")
assert_near(res8$summary$median_delta, frozen_row$median_delta, 1e-8, "8-gene residual median")
assert_near(res8$summary$ci_low, frozen_row$ci_low, 1e-8, "8-gene residual CI low")
assert_near(res8$summary$ci_high, frozen_row$ci_high, 1e-8, "8-gene residual CI high")
if (res8$summary$n_pairs != frozen_row$n_pairs || res8$summary$n_patients != frozen_row$n_patients ||
    res8$summary$n_positive != frozen_row$n_positive) {
  stop("8-gene reconstructed residual summary does not reproduce the frozen pair/patient counts.")
}

discovery_8 <- utils::read.csv(DISCOVERY_8GENE_CORR, check.names = FALSE, stringsAsFactors = FALSE)
assert_columns(discovery_8, c("signature", "n_cells", "spearman_rho"), "frozen eight-gene discovery correlation")
discovery_8_row <- discovery_8[discovery_8$signature == "SPP1_TAM_published_2025", , drop = FALSE]
if (nrow(discovery_8_row) != 1L) stop("Frozen eight-gene discovery correlation is not uniquely identified.")

comparison <- rbind(
  data.frame(
    comparator = "8-gene literature-derived marker score",
    genes = paste(OLD8_GENES, collapse = ";"),
    discovery_n_cells = as.integer(discovery_8_row$n_cells),
    discovery_spearman_rho = as.numeric(discovery_8_row$spearman_rho),
    residual_median_delta = res8$summary$median_delta,
    residual_ci_low = res8$summary$ci_low,
    residual_ci_high = res8$summary$ci_high,
    residual_n_positive = res8$summary$n_positive,
    residual_n_patients = res8$summary$n_patients,
    residual_signflip_p_two_sided = res8$summary$signflip_p_two_sided,
    score_definition = "frozen unweighted mean of log1p(CP10K) expression; comparator retained unchanged",
    source_status = "frozen comparator",
    stringsAsFactors = FALSE
  ),
  data.frame(
    comparator = "5-gene macrophage-specific SPP1-TAM score",
    genes = paste(FIVE_GENES, collapse = ";"),
    discovery_n_cells = discovery_5gene_summary$n_cells,
    discovery_spearman_rho = discovery_5gene_summary$spearman_rho_with_frozen_P1,
    residual_median_delta = res5$summary$median_delta,
    residual_ci_low = res5$summary$ci_low,
    residual_ci_high = res5$summary$ci_high,
    residual_n_positive = res5$summary$n_positive,
    residual_n_patients = res5$summary$n_patients,
    residual_signflip_p_two_sided = res5$summary$signflip_p_two_sided,
    score_definition = "gene-wise z-standardization across GSE245552 cells, followed by unweighted mean z-score",
    source_status = "new alternative-definition sensitivity",
    stringsAsFactors = FALSE
  )
)
write_csv(comparison, file.path(TABLES, "8_vs_5_comparator_summary.csv"))

signature_definitions <- data.frame(
  comparator = c("8-gene literature-derived marker score", "5-gene macrophage-specific SPP1-TAM score"),
  genes = c(paste(OLD8_GENES, collapse = ";"), paste(FIVE_GENES, collapse = ";")),
  n_genes = c(length(OLD8_GENES), length(FIVE_GENES)),
  source = c(
    "van Buuren et al. 2025 documentation of markers attributed to SPP1-TAM in prior literature; not a universal gold-standard signature",
    "van Buuren et al. 2025 macrophage-specific SPP1-TAM signature"
  ),
  doi = c("10.3389/fimmu.2025.1699237", "10.3389/fimmu.2025.1699237"),
  score_implementation = c(
    "frozen unweighted mean of log1p(CP10K) expression",
    "gene-wise z-standardization across dataset cells, followed by unweighted mean z-score"
  ),
  role = c("frozen primary comparator", "alternative-definition sensitivity only"),
  stringsAsFactors = FALSE
)
write_csv(signature_definitions, file.path(TABLES, "SPP1_TAM_comparator_definitions_8_vs_5.csv"))

## ------------------------------------------------------------------------
## 5. Audit and reproducibility metadata
## ------------------------------------------------------------------------

audit <- data.frame(
  check = c(
    "discovery_cells", "discovery_five_genes_present", "GSE245552_cells",
    "old8_raw_score_reproduction", "P1_aggregate_reproduction",
    "old8_aggregate_reproduction", "old8_residual_summary_reproduction",
    "new5_branch_status"
  ),
  observed = c(
    length(discovery_cells),
    all(FIVE_GENES %in% rownames(discovery_expr)),
    nrow(external_meta),
    max(abs(old8_delta), na.rm = TRUE),
    max(recon_sample$P1_abs_difference, na.rm = TRUE),
    max(recon_sample$old8_abs_difference, na.rm = TRUE),
    TRUE,
    "COMPLETED"
  ),
  expected_or_rule = c(
    "12433", "TRUE", "24231", "<=1e-8", "<=1e-8", "<=1e-8", "TRUE", "alternative sensitivity only"
  ),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "COMPLETED"),
  stringsAsFactors = FALSE
)
write_csv(audit, file.path(AUDIT, "5gene_sensitivity_audit.csv"))

package_names <- c("Seurat", "SeuratObject", "Matrix", "data.table", "digest", "NMF", "harmony", "lmerTest", "emmeans")
package_versions <- data.frame(
  package = package_names,
  version = vapply(package_names, function(pkg) {
    if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
write_csv(package_versions, file.path(MANIFEST, "package_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(MANIFEST, "sessionInfo.txt"))

writeLines(c(
  "Independent alternative-definition sensitivity analysis",
  paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai")),
  "primary_frozen_comparator_retained=TRUE",
  "frozen_outputs_overwritten=FALSE",
  paste0("discovery_cells=", length(discovery_cells)),
  paste0("GSE245552_cells=", nrow(external_meta)),
  "residual_model=P1_full ~ comparator + state",
  "patient_level_rule=mean across available state-specific LM-PT residual deltas",
  "five_gene_score=gene-wise z-standardization across dataset cells followed by mean",
  "source=van Buuren et al. 2025; DOI:10.3389/fimmu.2025.1699237"
), file.path(MANIFEST, "run_manifest.txt"))

all_files <- list.files(OUT, recursive = TRUE, full.names = TRUE)
all_files <- all_files[!dir.exists(all_files)]
hash_manifest <- data.frame(
  file = substring(all_files, nchar(OUT) + 2L),
  bytes = file.info(all_files)$size,
  sha256 = vapply(all_files, sha256, character(1)),
  stringsAsFactors = FALSE
)
write_csv(hash_manifest, file.path(MANIFEST, "SHA256_manifest.csv"))

message("COMPLETED independent 5-gene SPP1-TAM sensitivity analysis: ", OUT)
