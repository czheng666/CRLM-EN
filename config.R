## Shared configuration for the public R code release.
## Set these variables before running a script; no scientific default path is
## embedded in the release scripts.

options(stringsAsFactors = FALSE)

CRLM_CODE_RELEASE_ROOT <- normalizePath(
  Sys.getenv("CRLM_CODE_RELEASE_ROOT", unset = getwd()),
  winslash = "/", mustWork = FALSE
)
CRLM_PROJECT_ROOT <- normalizePath(
  Sys.getenv("CRLM_PROJECT_ROOT", unset = CRLM_CODE_RELEASE_ROOT),
  winslash = "/", mustWork = FALSE
)
CRLM_RESULTS_ROOT <- Sys.getenv("CRLM_RESULTS_ROOT", unset = "")
CRLM_FREEZE_DIR <- Sys.getenv("CRLM_FREEZE_DIR", unset = "")
CRLM_SUPPLEMENTARY_RESULTS_ROOT <- Sys.getenv(
  "CRLM_SUPPLEMENTARY_RESULTS_ROOT", unset = ""
)
CRLM_PRIMARY_RERUN_ROOT <- Sys.getenv("CRLM_PRIMARY_RERUN_ROOT", unset = CRLM_RESULTS_ROOT)
CRLM_5GENE_ROOT <- Sys.getenv("CRLM_5GENE_ROOT", unset = "")
CRLM_RAW_GSE245552_DIR <- Sys.getenv("CRLM_RAW_GSE245552_DIR", unset = "")
CRLM_OUTPUT_ROOT <- normalizePath(
  Sys.getenv("CRLM_OUTPUT_ROOT", unset = file.path(CRLM_CODE_RELEASE_ROOT, "outputs")),
  winslash = "/", mustWork = FALSE
)

crlm_require_dir <- function(path, env_name) {
  if (!nzchar(path) || !dir.exists(path)) {
    stop("Set ", env_name, " to an existing directory before running this script.")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

crlm_require_file <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
