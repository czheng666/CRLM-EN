################################################################################
# NMF iteration-depth sensitivity audit (500 -> 1000)
#
# This is an isolated reviewer-defense analysis.  It never overwrites the
# 150 primary fits, rank selection, or the frozen P1 object.  It reuses the
# exact rank/run/seed choices from the primary NMF manifest: runs 1, 15 and 30
# for each rank (4:8), and compares each 1,000-iteration fit with its matching
# 500-iteration fit by LSAP-aligned loading/usage correlations and top-gene
# overlap.  The script can use a small PSOCK cluster; if that fails, it falls
# back to serial execution with the same tasks and seeds.
################################################################################

options(stringsAsFactors = FALSE)

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args, value = TRUE)
.script <- if (length(.file_arg)) sub("^--file=", "", .file_arg[[1]]) else "NMF_iteration_depth_sensitivity.R"
.main <- file.path(CRLM_CODE_RELEASE_ROOT, "scripts", "analysis", "CRLM_methodologically_upgraded_rerun_v1.R")
if (!file.exists(.main)) stop("Main pipeline not found: ", .main)

# The caller may override these values.  Defaults point to the current
# project output and are intentionally separate from primary fit files.
if (!nzchar(Sys.getenv("CRLM_OUTPUT_DIR", unset = ""))) {
  Sys.setenv(CRLM_OUTPUT_DIR = CRLM_PRIMARY_RERUN_ROOT)
}
Sys.setenv(CRLM_SKIP_INITIAL_FUTURE_PLAN = "1")
source(.main, local = FALSE)

# The main pipeline initializes a multisession future plan for ordinary
# analyses. Reset it here: this isolated audit must not silently spawn four
# additional workers when workers=1, and should never duplicate the large NMF
# matrix behind the user's back.
if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

check_packages(c("NMF", "clue", "data.table"))
if (!"package:NMF" %in% search()) {
  suppressPackageStartupMessages(library(NMF))
}

input_file <- p(DIR$cache_nmf_matrix, "GSE178318_myeloid_NMF_input.rds")
manifest_file <- p(DIR$tables, "NMF_run_manifest.csv")
rank_selection_file <- p(DIR$cache_nmf_alignment, "rank_selection.rds")
frozen_file <- p(DIR$cache_projection, "frozen_P1_v4.rds")
if (!file.exists(input_file)) stop("NMF input cache not found: ", input_file)
if (!file.exists(manifest_file)) stop("NMF run manifest not found: ", manifest_file)
if (!file.exists(rank_selection_file)) stop("Rank-selection cache not found: ", rank_selection_file)
if (!file.exists(frozen_file)) stop("Frozen P1 cache not found: ", frozen_file)

dat <- readRDS(input_file)
X <- dat$X
if (is.null(dim(X)) || nrow(X) != 1500L || ncol(X) != 12433L) {
  stop("NMF input dimensions do not match the frozen input: ", paste(dim(X), collapse = " x "))
}
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
selected_runs <- c(1L, 15L, 30L)
tasks <- manifest[
  manifest$rank %in% 4:8 & manifest$run %in% selected_runs,
  c("rank", "run", "seed"), drop = FALSE
]
tasks <- tasks[order(tasks$rank, tasks$run), , drop = FALSE]
if (nrow(tasks) != 15L) stop("The fixed-seed audit does not contain 15 tasks: ", nrow(tasks))

out_dir <- file.path(DIR$tables, "NMF_iteration_depth_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fit_cache_dir <- file.path(out_dir, "fit_cache")
dir.create(fit_cache_dir, recursive = TRUE, showWarnings = FALSE)
run_out <- file.path(out_dir, "NMF_iteration_depth_sensitivity_run_level.csv")
summary_out <- file.path(out_dir, "NMF_iteration_depth_sensitivity_summary.csv")
manifest_out <- file.path(out_dir, "NMF_iteration_depth_sensitivity_manifest.csv")
session_out <- file.path(out_dir, "NMF_iteration_depth_sensitivity_sessionInfo.txt")
checkpoint_out <- file.path(out_dir, "NMF_iteration_depth_sensitivity_checkpoint.rds")

# Exporting the full 1,500 x 12,433 matrix to several PSOCK workers can exceed
# desktop RAM and terminate workers before a result is returned. Serial is the
# safe default; a small worker count can still be opted into explicitly.
workers <- suppressWarnings(as.integer(Sys.getenv("CRLM_NMF_DEPTH_WORKERS", unset = "1")))
if (!is.finite(workers) || workers < 1L) workers <- 1L
workers <- min(workers, 4L)
max_iter_new <- 1000L

fit_cache_file_for <- function(rank, run, max_iter = max_iter_new) {
  file.path(
    fit_cache_dir,
    sprintf("rank%d_run%02d_maxIter%d.rds", as.integer(rank), as.integer(run), as.integer(max_iter))
  )
}

save_fit_cache <- function(z, cache_file) {
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(cache_file, ".tmp-", Sys.getpid())
  saveRDS(z, tmp, compress = "xz")
  if (!file.rename(tmp, cache_file)) {
    file.copy(tmp, cache_file, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(cache_file)
}

valid_fit_cache <- function(cache_file, rank, run) {
  if (!file.exists(cache_file)) return(FALSE)
  z <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  is.list(z) &&
    identical(as.integer(z$rank), as.integer(rank)) &&
    identical(as.integer(z$run), as.integer(run)) &&
    is.matrix(z$W) && is.matrix(z$H) &&
    nrow(z$W) == nrow(X) && ncol(z$W) == as.integer(rank) &&
    nrow(z$H) == as.integer(rank) && ncol(z$H) == ncol(X) &&
    all(is.finite(z$W)) && all(is.finite(z$H))
}

fit_depth_worker <- function(task, Xmat, max_iter, cache_file = NULL) {
  set.seed(as.integer(task$seed))
  tm <- system.time({
    fit <- NMF::nmf(
      Xmat,
      rank = as.integer(task$rank),
      method = "brunet",
      nrun = 1L,
      maxIter = as.integer(max_iter),
      .options = ""
    )
  })
  z <- list(
    rank = as.integer(task$rank),
    run = as.integer(task$run),
    seed = as.integer(task$seed),
    W = NMF::basis(fit),
    H = NMF::coef(fit),
    n_iterations = tryCatch(as.integer(NMF::niter(fit)), error = function(e) NA_integer_),
    elapsed_seconds = unname(tm[["elapsed"]])
  )
  # Write the fit from the worker itself as soon as that individual run
  # returns.  Distinct tasks use distinct paths, so PSOCK workers do not
  # contend for the same file.  The main process also verifies the cache
  # after the worker result arrives.
  if (!is.null(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(cache_file, ".tmp-", Sys.getpid())
    saveRDS(z, tmp, compress = "xz")
    if (!file.rename(tmp, cache_file)) {
      file.copy(tmp, cache_file, overwrite = TRUE)
      unlink(tmp)
    }
  }
  z
}

fit_rank <- function(rank_tasks) {
  cache_files <- vapply(
    seq_len(nrow(rank_tasks)),
    function(i) fit_cache_file_for(rank_tasks$rank[[i]], rank_tasks$run[[i]], max_iter_new),
    character(1)
  )
  if (workers <= 1L) {
    return(lapply(seq_len(nrow(rank_tasks)), function(i)
      fit_depth_worker(rank_tasks[i, , drop = FALSE], X, max_iter_new, cache_files[[i]])))
  }
  cl <- NULL
  ans <- tryCatch({
    cl <- parallel::makeCluster(workers, type = "PSOCK")
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(NMF)))
    parallel::clusterExport(
      cl,
      c("X", "fit_depth_worker", "max_iter_new", "cache_files"),
      envir = .GlobalEnv
    )
    parallel::parLapplyLB(
      cl,
      seq_len(nrow(rank_tasks)),
      function(i) fit_depth_worker(
        rank_tasks[i, , drop = FALSE], X, max_iter_new, cache_files[[i]]
      )
    )
  }, error = function(e) {
    message("[NMF-depth] PSOCK execution failed; falling back to serial mode: ", conditionMessage(e))
    NULL
  }, finally = {
    if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
  })
  if (is.null(ans)) {
    lapply(seq_len(nrow(rank_tasks)), function(i)
      fit_depth_worker(rank_tasks[i, , drop = FALSE], X, max_iter_new, cache_files[[i]]))
  } else ans
}

top_overlap <- function(a, b, n = 15L) {
  n <- min(n, length(a), length(b))
  if (!n) return(NA_real_)
  length(intersect(a[seq_len(n)], b[seq_len(n)])) / n
}

rank_rows <- if (file.exists(checkpoint_out)) {
  readRDS(checkpoint_out)
} else {
  list()
}
if (!is.list(rank_rows)) rank_rows <- list()
rr <- length(rank_rows)
done_key <- if (rr) {
  vapply(rank_rows, function(z) paste(z$rank[[1]], z$run[[1]], z$seed[[1]], sep = ":"), character(1))
} else character()
message("[NMF-depth] start: 15 fits, maxIter=1000, workers=", workers)
for (k in 4:8) {
  tk <- tasks[tasks$rank == k, , drop = FALSE]
  if (length(done_key)) {
    tk_key <- paste(tk$rank, tk$run, tk$seed, sep = ":")
    tk <- tk[!tk_key %in% done_key, , drop = FALSE]
  }
  if (!nrow(tk)) next
  # Resume from an individual fit cache whenever possible.  Only missing or
  # invalid fits are sent to NMF; a completed fit is never recomputed merely
  # because the audit process was interrupted later.
  cache_files <- vapply(
    seq_len(nrow(tk)),
    function(i) fit_cache_file_for(tk$rank[[i]], tk$run[[i]], max_iter_new),
    character(1)
  )
  cached_ok <- vapply(
    seq_len(nrow(tk)),
    function(i) valid_fit_cache(cache_files[[i]], tk$rank[[i]], tk$run[[i]]),
    logical(1)
  )
  fits_new <- vector("list", nrow(tk))
  if (any(cached_ok)) {
    for (j in which(cached_ok)) {
      fits_new[[j]] <- readRDS(cache_files[[j]])
      message(sprintf("[NMF-depth] reuse fit cache: rank %d run %02d", tk$rank[[j]], tk$run[[j]]))
    }
  }
  if (any(!cached_ok)) {
    pending <- tk[!cached_ok, , drop = FALSE]
    pending_fits <- fit_rank(pending)
    for (j in seq_along(pending_fits)) {
      z_pending <- pending_fits[[j]]
      pending_cache <- fit_cache_file_for(z_pending$rank, z_pending$run, max_iter_new)
      if (!valid_fit_cache(pending_cache, z_pending$rank, z_pending$run)) {
        save_fit_cache(z_pending, pending_cache)
      }
      idx <- which(tk$run == z_pending$run & tk$rank == z_pending$rank)
      fits_new[[idx]] <- z_pending
      message(sprintf("[NMF-depth] fit cache written: rank %d run %02d -> %s", z_pending$rank, z_pending$run, pending_cache))
    }
  }
  for (j in seq_along(fits_new)) {
    z <- fits_new[[j]]
    if (is.null(z)) stop("NMF-depth fit result is missing: rank ", k, " index ", j)
    r <- as.integer(z$run)
    fit_cache_file <- fit_cache_file_for(k, r, max_iter_new)
    if (!valid_fit_cache(fit_cache_file, k, r)) {
      save_fit_cache(z, fit_cache_file)
    }
    old_file <- file.path(DIR$cache_nmf_fits, paste0("rank_", k), sprintf("rank%d_run%02d.rds", k, r))
    if (!file.exists(old_file)) stop("Corresponding 500-iteration fit is missing: ", old_file)
    old <- readRDS(old_file)
    old_norm <- normalize_nmf_scale(old$W, old$H)
    new_norm <- normalize_nmf_scale(z$W, z$H)
    W0 <- old_norm$W; H0 <- old_norm$H
    W1 <- new_norm$W; H1 <- new_norm$H
    rownames(W1) <- rownames(W0)
    colnames(W1) <- paste0("factor_", seq_len(k))
    rownames(H1) <- colnames(W1)
    colnames(H1) <- colnames(H0)
    mt <- match_factors_lsap(W0, W1)
    ord <- mt$factor_other[match(seq_len(k), mt$factor_ref)]
    W1a <- W1[, ord, drop = FALSE]
    H1a <- H1[ord, , drop = FALSE]
    loading_rho <- vapply(seq_len(k), function(i) safe_cor(W0[, i], W1a[, i]), numeric(1))
    usage_rho <- vapply(seq_len(k), function(i) safe_cor(H0[i, ], H1a[i, ]), numeric(1))
    top0 <- lapply(seq_len(k), function(i) rownames(W0)[order(W0[, i], decreasing = TRUE)][seq_len(min(15L, nrow(W0)))])
    top1 <- lapply(seq_len(k), function(i) rownames(W1a)[order(W1a[, i], decreasing = TRUE)][seq_len(min(15L, nrow(W1a)))])
    top15 <- vapply(seq_len(k), function(i) top_overlap(top0[[i]], top1[[i]]), numeric(1))
    rss500 <- as.numeric(old$rss)[1L]
    rss1000 <- compute_rss(X, W1a, H1a)
    p1_loading_rho <- NA_real_; p1_usage_rho <- NA_real_; p1_top15 <- NA_real_; p1_ref_factor <- NA_integer_
    if (k == 6L) {
      frozen <- readRDS(frozen_file)
      common <- intersect(rownames(frozen$W), rownames(W0))
      p1_vec <- frozen$W[common, "P1"]
      p1_ref_factor <- which.max(vapply(seq_len(k), function(i)
        safe_cor(W0[common, i], p1_vec), numeric(1)))
      p1_loading_rho <- safe_cor(W1a[common, p1_ref_factor], p1_vec)
      if ("P1" %in% rownames(frozen$H) && all(colnames(frozen$H) %in% colnames(H1a))) {
        p1_usage_rho <- safe_cor(
          H1a[p1_ref_factor, colnames(frozen$H)],
          frozen$H["P1", colnames(frozen$H)]
        )
      }
      p1_top15 <- top_overlap(
        rownames(W1a)[order(W1a[, p1_ref_factor], decreasing = TRUE)],
        as.character(frozen$P1_top15$gene),
        n = 15L
      )
    }
    rr <- rr + 1L
    rank_rows[[rr]] <- data.frame(
      rank = k, run = r, seed = z$seed, old_max_iter = old$max_iter,
      new_max_iter = max_iter_new, old_n_iterations = old$n_iterations,
      new_n_iterations = z$n_iterations,
      new_convergence_status = ifelse(!is.finite(z$n_iterations), "unknown",
        ifelse(z$n_iterations < max_iter_new, "criterion_met_before_cap",
               "reached_cap_stop_reason_unavailable")),
      basis_loading_median_spearman = median(loading_rho, na.rm = TRUE),
      basis_loading_min_spearman = min(loading_rho, na.rm = TRUE),
      usage_median_spearman = median(usage_rho, na.rm = TRUE),
      usage_min_spearman = min(usage_rho, na.rm = TRUE),
      top15_median_overlap = median(top15, na.rm = TRUE),
      top15_min_overlap = min(top15, na.rm = TRUE),
      rss_500 = rss500, rss_1000 = rss1000,
      rss_relative_change = if (is.finite(rss500) && rss500 > 0) (rss1000 - rss500) / rss500 else NA_real_,
      p1_ref_factor_in_500_fit = p1_ref_factor,
      p1_loading_spearman_vs_frozen = p1_loading_rho,
      p1_usage_spearman_vs_frozen = p1_usage_rho,
      p1_top15_overlap_vs_frozen = p1_top15,
      fit_cache_file = fit_cache_file,
      elapsed_seconds_1000 = z$elapsed_seconds,
      stringsAsFactors = FALSE
    )
    # Persist after every completed fit so an interrupted audit does not erase
    # earlier results and the run can be inspected incrementally.
    saveRDS(rank_rows[seq_len(rr)], checkpoint_out)
    write.csv(do.call(rbind, rank_rows[seq_len(rr)]), run_out, row.names = FALSE)
    done_key <- c(done_key, paste(k, r, z$seed, sep = ":"))
    message(sprintf("[NMF-depth] rank %d run %02d/30 complete; loading median=%.4f; top15=%.3f; %.1fs",
                    k, r, median(loading_rho, na.rm = TRUE), median(top15, na.rm = TRUE), z$elapsed_seconds))
  }
}

run_df <- do.call(rbind, rank_rows)
write.csv(run_df, run_out, row.names = FALSE)

summary_df <- do.call(rbind, lapply(split(run_df, run_df$rank), function(q) {
  finite_values <- function(x) {
    x <- as.numeric(x)
    x[is.finite(x)]
  }
  p1_loading <- finite_values(q$p1_loading_spearman_vs_frozen)
  p1_top15 <- finite_values(q$p1_top15_overlap_vs_frozen)
  data.frame(
    rank = q$rank[[1]], n_fits = nrow(q),
    median_loading_spearman = median(q$basis_loading_median_spearman, na.rm = TRUE),
    min_loading_spearman = min(q$basis_loading_min_spearman, na.rm = TRUE),
    median_usage_spearman = median(q$usage_median_spearman, na.rm = TRUE),
    min_usage_spearman = min(q$usage_min_spearman, na.rm = TRUE),
    median_top15_overlap = median(q$top15_median_overlap, na.rm = TRUE),
    min_top15_overlap = min(q$top15_min_overlap, na.rm = TRUE),
    median_rss_relative_change = median(q$rss_relative_change, na.rm = TRUE),
    median_p1_loading_spearman_vs_frozen = if (length(p1_loading)) median(p1_loading) else NA_real_,
    min_p1_loading_spearman_vs_frozen = if (length(p1_loading)) min(p1_loading) else NA_real_,
    median_p1_top15_overlap_vs_frozen = if (length(p1_top15)) median(p1_top15) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_df, summary_out, row.names = FALSE)

manifest_out_df <- data.frame(
  analysis = "NMF iteration-depth sensitivity",
  primary_max_iter = 500L,
  sensitivity_max_iter = max_iter_new,
  ranks = "4,5,6,7,8",
  runs_per_rank = 3L,
  selected_primary_runs = paste(selected_runs, collapse = ","),
  n_fits = nrow(run_df),
  alignment = "LSAP/Hungarian against matching 500-iteration fit",
  input_rds = input_file,
  input_sha256 = sha256_file(input_file),
  primary_manifest_sha256 = sha256_file(manifest_file),
  frozen_P1_sha256 = sha256_file(frozen_file),
  workers_requested = workers,
  R_version = R.version.string,
  NMF_version = as.character(packageVersion("NMF")),
  generated_at = format(Sys.time(), tz = "Asia/Shanghai"),
  stringsAsFactors = FALSE
)
write.csv(manifest_out_df, manifest_out, row.names = FALSE)
writeLines(capture.output(sessionInfo()), session_out, useBytes = TRUE)
message("[NMF-depth] DONE; outputs: ", out_dir)
