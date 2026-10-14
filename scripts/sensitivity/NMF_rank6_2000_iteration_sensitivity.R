#!/usr/bin/env Rscript
# Isolated reviewer-defense audit: rank 6, fixed primary seeds 1/15/30,
# maxIter 500 -> 1000 -> 2000. Does not overwrite primary NMF/P1.
options(stringsAsFactors = FALSE, width = 200)
suppressPackageStartupMessages(library(NMF))

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()), "config.R"))

PRIMARY <- crlm_require_dir(CRLM_PRIMARY_RERUN_ROOT, "CRLM_PRIMARY_RERUN_ROOT")
OUT <- file.path(PRIMARY, "results/tables/NMF_iteration_depth_sensitivity_rank6_2000")
CACHE <- file.path(OUT, "fit_cache")
dir.create(CACHE, recursive = TRUE, showWarnings = FALSE)
LOG <- file.path(OUT, "run.log")
logmsg <- function(...) {
  z <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "Asia/Shanghai"), " | ", paste0(..., collapse = ""))
  cat(z, "\n", file = LOG, append = TRUE)
  cat(z, "\n")
}
safe_write <- function(x, f) { data.table::fwrite(x, f, na = "NA") }

input_file <- file.path(PRIMARY, "cache/nmf_matrix/GSE178318_myeloid_NMF_input.rds")
manifest_file <- file.path(PRIMARY, "results/tables/NMF_run_manifest.csv")
frozen_file <- file.path(PRIMARY, "cache/projection/frozen_P1_v4.rds")
if (!file.exists(input_file) || !file.exists(manifest_file) || !file.exists(frozen_file)) stop("Required primary NMF files missing")
dat <- readRDS(input_file); X <- dat$X
if (!identical(dim(X), c(1500L, 12433L))) stop("Unexpected NMF matrix dimension: ", paste(dim(X), collapse = "x"))
mf <- read.csv(manifest_file, stringsAsFactors = FALSE)
tk <- mf[mf$rank == 6 & mf$run %in% c(1L, 15L, 30L), c("rank", "run", "seed"), drop = FALSE]
tk <- tk[order(tk$run), , drop = FALSE]
if (nrow(tk) != 3L) stop("Expected three rank-6 tasks, got ", nrow(tk))
frozen <- readRDS(frozen_file)

normalize_nmf_scale <- function(W, H) {
  sc <- colSums(W); sc[!is.finite(sc) | sc <= 0] <- 1
  list(W = sweep(W, 2, sc, "/"), H = sweep(H, 1, sc, "*"))
}
spearman <- function(a, b) suppressWarnings(stats::cor(as.numeric(a), as.numeric(b), method = "spearman", use = "pairwise.complete.obs"))
top_genes <- function(W, genes, factor, n = 15L) genes[order(W[, factor], decreasing = TRUE)][seq_len(min(n, nrow(W)))]
get_rss <- function(z) {
  if (!is.null(z$rss) && length(z$rss) == 1L && is.finite(z$rss)) return(as.numeric(z$rss))
  if (!is.null(z$residual) && length(z$residual) == 1L && is.finite(z$residual)) return(as.numeric(z$residual))
  # Older 1000-iteration caches did not retain RSS; compute it from the
  # normalized reconstruction so the depth comparison is complete.
  if (!is.null(z$W) && !is.null(z$H)) {
    zz <- normalize_nmf_scale(z$W, z$H)
    return(sum((X - zz$W %*% zz$H)^2))
  }
  NA_real_
}
fit_file <- function(run, iter) file.path(CACHE, sprintf("rank6_run%02d_maxIter%d.rds", run, iter))
valid_fit <- function(f, run, iter) {
  if (!file.exists(f)) return(FALSE)
  z <- tryCatch(readRDS(f), error = function(e) NULL)
  is.list(z) && identical(as.integer(z$run), as.integer(run)) && identical(as.integer(z$max_iter), as.integer(iter)) &&
    is.matrix(z$W) && is.matrix(z$H) && identical(dim(z$W), c(1500L, 6L)) && identical(dim(z$H), c(6L, 12433L)) &&
    all(is.finite(z$W)) && all(is.finite(z$H))
}
run_one <- function(task, iter) {
  f <- fit_file(task$run, iter)
  if (valid_fit(f, task$run, iter)) { logmsg("reuse cache rank6 run", task$run, " maxIter", iter); return(readRDS(f)) }
  logmsg("START fit rank6 run", task$run, " seed=", task$seed, " maxIter=", iter)
  set.seed(as.integer(task$seed))
  tm <- system.time({ fit <- NMF::nmf(X, rank = 6L, method = "brunet", nrun = 1L, maxIter = as.integer(iter), .options = "") })
  z <- list(rank = 6L, run = as.integer(task$run), seed = as.integer(task$seed), max_iter = as.integer(iter),
            W = NMF::basis(fit), H = NMF::coef(fit), rss = tryCatch(as.numeric(NMF::residuals(fit, X)), error = function(e) NA_real_),
            n_iterations = tryCatch(as.integer(NMF::niter(fit)), error = function(e) NA_integer_), elapsed_seconds = unname(tm[["elapsed"]]),
            stop_reason = "NMF package stop reason not retained in primary object")
  tmp <- paste0(f, ".tmp-", Sys.getpid()); saveRDS(z, tmp, compress = "xz"); if (!file.rename(tmp, f)) { file.copy(tmp, f, overwrite = TRUE); unlink(tmp) }
  logmsg("DONE fit rank6 run", task$run, " maxIter=", iter, " niter=", z$n_iterations, " elapsed_min=", round(z$elapsed_seconds / 60, 2))
  z
}

gene_names <- if (!is.null(dat$genes)) as.character(dat$genes) else rownames(X)
if (is.null(gene_names) || length(gene_names) != nrow(X)) stop("NMF gene names unavailable")
old500 <- lapply(seq_len(nrow(tk)), function(i) readRDS(file.path(PRIMARY, "cache/nmf_fits", "rank_6", sprintf("rank6_run%02d.rds", tk$run[i]))))
old1000 <- lapply(seq_len(nrow(tk)), function(i) {
  f <- file.path(PRIMARY, "results/tables/NMF_iteration_depth_sensitivity/fit_cache", sprintf("rank6_run%02d_maxIter1000.rds", tk$run[i]))
  if (!file.exists(f)) stop("Missing existing 1000-iteration cache: ", f)
  readRDS(f)
})
new2000 <- lapply(seq_len(nrow(tk)), function(i) run_one(tk[i, , drop = FALSE], 2000L))

rows <- lapply(seq_len(nrow(tk)), function(i) {
  o5 <- normalize_nmf_scale(old500[[i]]$W, old500[[i]]$H)
  o1 <- normalize_nmf_scale(old1000[[i]]$W, old1000[[i]]$H)
  o2 <- normalize_nmf_scale(new2000[[i]]$W, new2000[[i]]$H)
  # Pairwise factor alignment by absolute loading correlation (Hungarian LSAP).
  cc <- matrix(0, 6, 6); for (a in 1:6) for (b in 1:6) cc[a, b] <- spearman(o1$W[, a], o2$W[, b])
  if (requireNamespace("clue", quietly = TRUE)) {
    asn <- as.integer(clue::solve_LSAP(1 - abs(cc)))
  } else asn <- apply(abs(cc), 1, which.max)
  w_corr <- vapply(1:6, function(a) spearman(o1$W[, a], o2$W[, asn[a]]), numeric(1))
  h_corr <- vapply(1:6, function(a) spearman(o1$H[a, ], o2$H[asn[a], ]), numeric(1))
  overlaps <- vapply(1:6, function(a) length(intersect(top_genes(o1$W, gene_names, a), top_genes(o2$W, gene_names, asn[a]))) / 15, numeric(1))
  # Frozen P1 is row 1 in the frozen object; match new factors directly to frozen P1.
  fp <- normalize_nmf_scale(frozen$W, frozen$H)
  p1_cor <- vapply(1:6, function(a) spearman(fp$W[, 1], o2$W[, a]), numeric(1))
  p1_factor <- which.max(p1_cor)
  p1_overlap <- length(intersect(top_genes(fp$W, rownames(fp$W), 1), top_genes(o2$W, gene_names, p1_factor))) / 15
  data.frame(rank = 6L, run = tk$run[i], seed = tk$seed[i], maxIter_500_niter = old500[[i]]$n_iterations,
             maxIter_1000_niter = old1000[[i]]$n_iterations, maxIter_2000_niter = new2000[[i]]$n_iterations,
             loading_spearman_1000_vs_2000_median = stats::median(w_corr), loading_spearman_1000_vs_2000_min = min(w_corr),
             usage_spearman_1000_vs_2000_median = stats::median(h_corr), usage_spearman_1000_vs_2000_min = min(h_corr),
             top15_overlap_1000_vs_2000_median = stats::median(overlaps), top15_overlap_1000_vs_2000_min = min(overlaps),
             rss_500 = get_rss(old500[[i]]), rss_1000 = get_rss(old1000[[i]]), rss_2000 = get_rss(new2000[[i]]),
             p1_factor_in_2000 = p1_factor, p1_loading_spearman_vs_frozen = max(p1_cor), p1_top15_overlap_vs_frozen = p1_overlap,
             stop_reason = "algorithm-level stop reason unavailable at cap", stringsAsFactors = FALSE)
})
run_df <- do.call(rbind, rows)
safe_write(run_df, file.path(OUT, "rank6_500_1000_2000_run_level.csv"))
summary <- data.frame(rank = 6L, n_fits = nrow(run_df),
                      median_loading_spearman_1000_vs_2000 = median(run_df$loading_spearman_1000_vs_2000_median),
                      min_loading_spearman_1000_vs_2000 = min(run_df$loading_spearman_1000_vs_2000_min),
                      median_usage_spearman_1000_vs_2000 = median(run_df$usage_spearman_1000_vs_2000_median),
                      min_usage_spearman_1000_vs_2000 = min(run_df$usage_spearman_1000_vs_2000_min),
                      median_top15_overlap_1000_vs_2000 = median(run_df$top15_overlap_1000_vs_2000_median),
                      min_top15_overlap_1000_vs_2000 = min(run_df$top15_overlap_1000_vs_2000_min),
                      median_p1_loading_spearman_vs_frozen = median(run_df$p1_loading_spearman_vs_frozen),
                      min_p1_loading_spearman_vs_frozen = min(run_df$p1_loading_spearman_vs_frozen),
                      median_p1_top15_overlap_vs_frozen = median(run_df$p1_top15_overlap_vs_frozen),
                      stringsAsFactors = FALSE)
safe_write(summary, file.path(OUT, "rank6_500_1000_2000_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
logmsg("FINISHED rank6 2000 sensitivity; summary written")
print(run_df)
