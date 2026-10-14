################################################################################
################################################################################
# CRLM 方法学升级规范化重跑总脚本 v1.0
# Formal single-file rerun pipeline — full manuscript chain
#
# 目标：
#   在尽量保持原文章科学主线与实际方法一致的前提下，从 GSE178318 原始
#   10x 矩阵重新执行 QC、doublet、Harmony/聚类、marker-module 注释与 NMF；
#   对旧版最主要的方法学不足进行修复：
#       rank 4–8 全矩阵、每个 rank 30 次独立 Brunet NMF；
#       LSAP/Hungarian 因子对齐；
#       loading / usage / top-gene 稳定性；
#       基于“稳定性 + 重构误差拐点”的机械 rank 重选；
#       10→20→30 次运行收敛审计。
#
# 原研究保持不变的主要设定：
#   GSE178318 GEO-supplied 10x-format count matrix
#   nFeature_RNA >= 500
#   percent.mt <= 15%
#   每文库 log10(nFeature), log10(nCount) median ± 3 MAD
#   scDblFinder 按 15 个文库分别运行，dbr.per1k=0.008
#   全局 atlas: 3000 HVG, PCA 50, Harmony(patient), 邻居/UMAP PCs 1:30
#   lineage reclustering: 2000 HVG, PCA 30, Harmony(sample), PCs 1:20
#   NMF: log1p(CP10K)（与 Seurat LogNormalize scale.factor=10000 对齐），
#        1500 HVG，Brunet，rank 4:8
#   外部：冻结 top-15 权重，不重新 NMF
#   统计：患者为生物学重复；state-specific paired Wilcoxon 作为探索性分支；
#         Figure 4A 的主状态推断由独立的 lmerTest/emmeans 脚本完成；
#         patient bootstrap 10000；composition/activity 患者层分解
#   空间：kNN(k=6)为主分析；物理距离100/200/300 μm为敏感性；
#         残差 Moran I / local High-High；999 permutations
#
# 使用方法：
#   1) 通过 CRLM_RERUN_DIR 和 CRLM_GSE*_DIR 设置输入路径；不依赖旧版绝对路径
#   2) source("CRLM_methodologically_upgraded_rerun_v1.R")
#   3) 选择运行：
#        run_stage(1)
#        run_stage(2)
#        ...
#      或：
#        run_many(1:8)
#
# 注意：
#   - 本脚本不会自动一次性运行所有 stage。
#   - 每个耗时 stage 均带缓存/断点续跑。
#   - 不读取旧 P1/NMF 来决定 rank；旧 P1/annotation 只在新版结果完成后作为“对应关系审计”。
#   - 新版 broad/fine annotation 由预先定义的 marker modules 重算；旧 cluster map 不决定纳入或标签。
################################################################################


############################
# Stage 0：配置、目录、工具函数
############################
options(stringsAsFactors = FALSE)

library(future)

# Standalone audits (for example the NMF iteration-depth sensitivity) source
# this file only to reuse helper functions.  They must not inherit a new
# multisession pool and duplicate the large expression matrix.
if (identical(Sys.getenv("CRLM_SKIP_INITIAL_FUTURE_PLAN", unset = ""), "1")) {
  plan("sequential")
} else {
  plan("multisession", workers = 4)
}

options(
  future.globals.maxSize = 60 * 1024^3
)

PIPELINE_VERSION <- "CRLM_methodologically_upgraded_rerun_v1"
# CACHE_SIGNATURE is computed after PARAM and the required packages are
# available.  It invalidates caches when frozen parameters or key package
# versions change, even if PIPELINE_VERSION was not manually edited.
CACHE_SIGNATURE <- NULL
# ============================
# Step-level cache helpers
# ============================

cache_version_file <- function(path) {
  paste0(path, ".pipeline_version.txt")
}

cache_signature_file <- function(path) {
  paste0(path, ".cache_signature.txt")
}

cache_expected_signature <- function() {
  if (!is.null(CACHE_SIGNATURE) && length(CACHE_SIGNATURE) == 1L) {
    CACHE_SIGNATURE
  } else {
    PIPELINE_VERSION
  }
}

cache_stamp_object <- function(object) {
  # Keep version and signature inside the object as well as in small sidecars.
  # Sidecars let cache checks avoid deserializing large Seurat objects.
  signature <- cache_expected_signature()
  if (inherits(object, "Seurat")) {
    misc <- object@misc
    misc$pipeline_version <- PIPELINE_VERSION
    misc$cache_signature <- signature
    object@misc <- misc
  } else if (is.list(object) && !is.data.frame(object)) {
    object$pipeline_version <- PIPELINE_VERSION
    object$cache_signature <- signature
  }
  object
}

cache_has_current_version <- function(path, expected = NULL) {
  if (is.null(expected)) expected <- cache_expected_signature()
  sf <- cache_signature_file(path)
  if (!file.exists(path) || !file.exists(sf)) return(FALSE)
  got <- tryCatch(
    trimws(readLines(sf, n = 1L, warn = FALSE)),
    error = function(e) character()
  )
  if (length(got) != 1L) return(FALSE)
  if (identical(got[[1]], expected)) return(TRUE)
  # Optional, explicit one-run migration escape hatch.  This is deliberately
  # disabled by default; it is only for reusing caches produced immediately
  # before a code-only repair whose analytical inputs and parameters are
  # unchanged.  Set CRLM_CACHE_COMPAT_SIGNATURES to a comma-separated list
  # for that recovery invocation, then unset it for ordinary runs.
  compat <- trimws(strsplit(Sys.getenv("CRLM_CACHE_COMPAT_SIGNATURES", ""),
                            ",", fixed = TRUE)[[1]])
  length(compat) > 0L && nzchar(compat[[1]]) && got[[1]] %in% compat
}

save_step_cache <- function(object, path){
  object <- cache_stamp_object(object)
  
  dir.create(
    dirname(path),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  tmp <- paste0(path, ".tmp")
  
  saveRDS(
    object,
    tmp,
    compress = "gzip"
  )
  
  if(file.exists(path)){
    file.remove(path)
  }
  
  if (!file.rename(tmp, path)) {
    stop("无法完成缓存替换：", path)
  }
  
  # Tiny sidecars make future signature checks cheap and avoid loading a
  # potentially multi-GB Seurat object merely to inspect its metadata.
  vf <- cache_version_file(path)
  vf_tmp <- paste0(vf, ".tmp")
  writeLines(PIPELINE_VERSION, vf_tmp, useBytes = TRUE)
  if (file.exists(vf)) file.remove(vf)
  if (!file.rename(vf_tmp, vf)) {
    stop("无法写入缓存版本标记：", vf)
  }
  sf <- cache_signature_file(path)
  sf_tmp <- paste0(sf, ".tmp")
  writeLines(cache_expected_signature(), sf_tmp, useBytes = TRUE)
  if (file.exists(sf)) file.remove(sf)
  if (!file.rename(sf_tmp, sf)) {
    stop("无法写入缓存签名标记：", sf)
  }
}

load_step_cache <- function(path){
  if (cache_has_current_version(path)) {
    message("Loading current cache (signature matched): ", basename(path))
    return(readRDS(path))
  }
  if (file.exists(path)) {
    message("Cache exists but signature is missing or mismatched: ", basename(path))
  }
  NULL
}

must_load_step_cache <- function(path, label = basename(path)) {
  object <- load_step_cache(path)
  if (is.null(object)) {
    stop("缓存缺失或版本不匹配（", label, "）：", path)
  }
  object
}

cache_exists <- function(path){
  file.exists(path)
}
# ---- 路径：新版输出与旧版完全隔离 ----
# Prefer the actual script path both for Rscript and for source().  The
# CRLM_RERUN_DIR environment variable remains an explicit override.
.cmd_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .cmd_args, value = TRUE)
.source_file <- if (length(.file_arg)) {
  sub("^--file=", "", .file_arg[[1]])
} else {
  .frame_files <- vapply(
    sys.frames(),
    function(fr) {
      if (exists("ofile", envir = fr, inherits = FALSE)) {
        as.character(get("ofile", envir = fr, inherits = FALSE))
      } else {
        NA_character_
      }
    },
    character(1)
  )
  .frame_files <- .frame_files[!is.na(.frame_files) & nzchar(.frame_files)]
  if (length(.frame_files)) .frame_files[[1]] else NA_character_
}
.script_dir <- if (!is.na(.source_file) && nzchar(.source_file)) {
  dirname(normalizePath(.source_file, mustWork = FALSE))
} else {
  getwd()
}
PROJECT_ROOT <- normalizePath(
  Sys.getenv("CRLM_RERUN_DIR", unset = .script_dir),
  mustWork = FALSE
)
.base_env <- Sys.getenv("CRLM_OUTPUT_DIR", unset = "")
# On Windows, normalizePath() can expand a supplied 8.3 short path back to a
# non-ASCII long path.  Rscript may then fail to translate that path when
# reading large RDS files.  Preserve an explicit user-supplied path verbatim;
# the default remains normalized for ordinary interactive use.
BASE <- if (nzchar(.base_env)) {
  gsub("\\\\", "/", .base_env)
} else {
  # Portable default: keep all generated files below the project root.  A
  # machine-specific location is used only when CRLM_OUTPUT_DIR is supplied.
  file.path(PROJECT_ROOT, "outputs", PIPELINE_VERSION)
}
RAW_ROOT <- normalizePath(
  Sys.getenv("CRLM_RAW_DIR", unset = file.path(PROJECT_ROOT, "data", "GEO_raw")),
  mustWork = FALSE
)

# ---- GSE178318 输入路径 ----

GSE178318_DIR <- Sys.getenv(
  "CRLM_GSE178318_DIR",
  unset = file.path(RAW_ROOT, "GSE178318")
)

# Accept both the conventional GEO-export names and the user's leading-
# underscore names (_matrix.mtx.gz, _barcodes.tsv.gz, _genes.tsv.gz).
first_existing_input <- function(root, candidates) {
  paths <- file.path(root, candidates)
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else paths[[1]]
}
GSE178318_MATRIX <- first_existing_input(
  GSE178318_DIR, c("GSE178318_matrix.mtx.gz", "_matrix.mtx.gz")
)
GSE178318_BARCODES <- first_existing_input(
  GSE178318_DIR, c("GSE178318_barcodes.tsv.gz", "_barcodes.tsv.gz")
)
GSE178318_FEATURES <- first_existing_input(
  GSE178318_DIR, c("GSE178318_genes.tsv.gz", "GSE178318_features.tsv.gz", "_genes.tsv.gz", "_features.tsv.gz")
)

# 外部数据下载完成后放到下面目录；不存在时相应 Stage 会明确停止。
GSE245552_DIR <- Sys.getenv("CRLM_GSE245552_DIR", unset = file.path(RAW_ROOT, "GSE245552"))
GSE231559_DIR <- Sys.getenv("CRLM_GSE231559_DIR", unset = file.path(RAW_ROOT, "GSE231559"))
GSE115469_DIR <- Sys.getenv("CRLM_GSE115469_DIR", unset = file.path(RAW_ROOT, "GSE115469"))
GSE136103_DIR <- Sys.getenv("CRLM_GSE136103_DIR", unset = file.path(RAW_ROOT, "GSE136103"))
GSE235079_DIR <- Sys.getenv("CRLM_GSE235079_DIR", unset = file.path(RAW_ROOT, "GSE235079"))
GSE168933_DIR <- Sys.getenv("CRLM_GSE168933_DIR", unset = file.path(RAW_ROOT, "GSE168933"))

# GSE168933 在当前代码包中没有找到可执行的原始重建脚本。
# 因此 v3 不会猜测其映射算法；若要从头重建，请提供下面的“冻结映射输入”
# （必须是根据该数据集官方/作者标签预先得到、且不使用P1结果挑选的表）。
# 期望列：reference_label, marker_gene（可多行/label）。
GSE168933_REFERENCE_MARKER_FILE <- first_existing_input(
  GSE168933_DIR,
  c("GSE168933_reference_marker_map.csv",
    "02_reference_marker_genes_GSE168933.csv")
)

# 空间数据用 manifest 管理。Stage 12 首次运行会生成模板。
SPATIAL_MANIFEST <- file.path(BASE, "inputs", "spatial_manifest.csv")

# 可选：如果你有旧版完整 NMF RDS，可填写路径。
# 仅在“新 rank 已经确定之后”用于 P1 correspondence audit，不参与 rank 选择。
LEGACY_NMF_RDS <- ""

# ---- 冻结参数 ----
PARAM <- list(
  seed = 178318L,
  
  # QC
  min_features = 500L,
  max_percent_mt = 15,
  mad_multiplier = 3,
  scdblfinder_dbr_per1k = 0.008,
  scdblfinder_dbr_method = "fixed",
  
  # atlas
  global_hvg = 3000L,
  global_pcs = 50L,
  global_neighbor_dims = 1:30,
  global_primary_resolution = 0.6,
  
  # myeloid reclustering
  lineage_hvg = 2000L,
  lineage_pcs = 30L,
  lineage_neighbor_dims = 1:20,
  lineage_primary_resolution = 0.6,
  annotation_module_margin = 0.05,
  # Marker-audit fallback labels are assigned only when the corresponding
  # prespecified quality module is among the top two cluster-level scores and
  # exceeds the across-cluster median.  A failed guard stops the rerun rather
  # than silently carrying a legacy cluster-number override forward.
  marker_audit_top_rank = 2L,
  p1_anchor_min_margin = 0.15,
  
  # NMF
  nmf_hvg = 1500L,
  nmf_ranks = 4:8,
  nmf_runs_per_rank = 30L,
  nmf_method = "brunet",
  # Explicit 500-iteration ceiling for every primary independent fit.  The default in
  # NMF is implicit and can make a full 1,500 x 12,433 matrix run for hours
  # without a reproducible audit trail.  NMFfit does not expose a portable
  # package-level convergence flag across methods, so `converged` below is
  # retained only as an iteration-ceiling proxy.  Reaching this ceiling is
  # audited and reported, but is not by itself treated as a failed fit or as
  # biological evidence for rank selection.
  nmf_max_iter = 500L,
  nmf_loading_spearman_threshold = 0.70,
  nmf_usage_spearman_threshold = 0.70,
  factor_recurrence_threshold = 0.80,
  top15_overlap_threshold = 0.50,
  eligible_rank_stable_factor_fraction = 0.80,
  
  # patient statistics
  n_boot = 10000L,
  n_signflip = 10000L,
  
  # spatial
  n_spatial_perm = 999L,
  matched_null_n = 5000L,
  # Expression/detection-matched null: retain the ten closest candidates
  # after an explicit raw-scale caliper.  A target for which no unused gene
  # satisfies both calipers is retained as an auditable extreme target rather
  # than being silently drawn from a broad 50-gene neighbourhood.
  matched_null_knn = 10L,
  matched_null_expr_caliper = 0.5,
  matched_null_detection_caliper = 0.05,
  projection_normalization = "log1p_cp10k",
  gse217414_radius_px = 150,
  gse225857_radius_px = 250,
  spatial_knn_k = 6L,
  spatial_radius_um_sensitivity = c(100, 200, 300)
)

STRICT_EXPECTED_COUNTS <- TRUE
EXPECTED_GSE178318_GENE_COUNT <- 33694L
EXPECTED_GSE178318_CELL_COUNT <- 140281L
# Feature annotations can change slightly between GEO exports; cell count is
# kept exact while the feature count allows a small, explicit tolerance.
EXPECTED_GSE178318_GENE_COUNT_TOLERANCE <- 0.05

# P1选择锚点必须在发现队列内预先定义；不使用外部、空间或患者结局。
# Core anchors identify the intended axis; auxiliary genes prevent selection
# from being driven by only two transcripts.  These are factor-selection
# anchors, not the final P1 score gene set.
P1_CORE_ANCHOR_GENES <- c("SPP1", "APOC1")
P1_AUXILIARY_ANCHOR_GENES <- c("APOE", "GPNMB", "MARCO", "LPL", "CST3", "LGALS3")
P1_ANCHOR_GENES <- c(P1_CORE_ANCHOR_GENES, P1_AUXILIARY_ANCHOR_GENES)
# Provisional sensitivity definition.  This is intentionally not treated as
# frozen until the new 30-run NMF has produced the final P1 top-15 and the
# marker-exclusion rule is reviewed against that actual gene set.  Full-P1
# branches may run before confirmation; only marker-excluded branches are held.
# GPNMB is retained in the P1-state module but is not silently counted as
# deleted when it is absent from the final P1 scoring genes.
P1_MARKER_EXCLUDED_GENES <- c("SPP1", "APOC1", "APOE")
# Confirmed after inspection of the rerun-derived P1 top-15: all three genes
# participate in the actual score and are removed for the prespecified
# marker-exclusion sensitivity (GPNMB is not part of the final P1 top-15).
P1_MARKER_EXCLUSION_CONFIRMED <- TRUE

# 原文章冻结的 P1 top-15，仅在最终 rank 选定后做 correspondence audit。
LEGACY_P1_TOP15 <- data.frame(
  gene = c("FTL","APOC1","SPP1","APOE","CSTB","LGALS1","CTSD","ACTB",
           "VIM","CD63","CTSB","LYZ","HLA-A","PSAP","CAPG"),
  weight = c(123.199,113.953,110.295,104.873,96.243,93.224,90.478,85.968,
             79.371,71.778,71.539,69.206,68.240,64.538,64.144),
  stringsAsFactors = FALSE
)
LEGACY_P1_TOP15$weight_norm <- LEGACY_P1_TOP15$weight / sum(LEGACY_P1_TOP15$weight)

# Legacy maps are retained only for concordance auditing.  They are never used
# to assign the new rerun's primary annotations.
LEGACY_GLOBAL_MAP <- data.frame(
  cluster = as.character(0:22),
  broad_lineage = c(
    "T_NK","T_NK","Stress_mixed","T_NK","T_NK","T_NK","B","Myeloid_DC",
    "Myeloid_DC","T_NK","Epithelial","Epithelial","Cycling_mixed","T_NK",
    "Plasma","Myeloid_DC","Low_quality_high_mt","Fibroblast","Mast",
    "Stress_mixed","Rare_neuronal_like","Endothelial",
    "Mixed_epithelial_immunoglobulin"
  ),
  annotation_status = c(
    "retained","retained","excluded_from_lineage_primary","retained","retained",
    "retained","retained","retained","retained","retained","retained","retained",
    "excluded_from_lineage_primary","retained","retained","retained",
    "excluded_from_lineage_primary","retained","retained",
    "excluded_from_lineage_primary","excluded_from_lineage_primary","retained",
    "excluded_from_lineage_primary"
  ),
  stringsAsFactors = FALSE
)

# Legacy fine-state map, used only for comparison after new marker/module
# annotation.
LEGACY_MYELOID_MAP <- data.frame(
  cluster = as.character(0:13),
  fine_label = c(
    "Mac_FOLR2_C1Q",
    "Mono_S100A8_A12",
    "Mac_SPP1_MARCO",
    "Myeloid_TWISTNB_THBS1",
    "DC_CD1C_FCER1A",
    "Mac_SPP1_TREM2_GPNMB",
    "Mac_CCL18_APOC1",
    "Myeloid_stress",
    "Myeloid_IFN_CLEC5A",
    "Mixed_T_NK_myeloid",
    "Mono_ADGRE1",
    "Myeloid_cycling",
    "DC_CCL19_CCL22",
    "DC_CLEC9A_XCR1"
  ),
  annotation_status = c(
    "retained","retained","retained","exploratory","retained","retained",
    "retained","flagged_stress","exploratory","excluded_mixed","exploratory",
    "excluded_cycling","retained","retained"
  ),
  stringsAsFactors = FALSE
)


# Legacy T/NK fine-state map, used only for comparison.
LEGACY_TNK_MAP <- data.frame(
  cluster = as.character(0:11),
  fine_label = c(
    "CD8_ZNF683_checkpoint_cytotoxic",
    "CD8_IL7R_memory",
    "T_naive_CCR7_LEF1",
    "NK_FGFBP2_PRF1",
    "Treg_FOXP3",
    "CD8_CXCL13_checkpoint",
    "CD8_GZMK_memory",
    "T_gamma_delta",
    "NK_FCER1G_KLRC1",
    "T_gamma_delta_NCR3",
    "T_stress",
    "T_activated_TNFSF4"
  ),
  annotation_status = c(
    "retained","retained","retained","retained","retained","retained",
    "retained","retained","retained","retained","flagged_stress",
    "exploratory_low_coverage"
  ),
  stringsAsFactors = FALSE
)

LINEAGES_TO_SUBCLUSTER <- c(
  "T_NK","Myeloid_DC","B","Plasma","Epithelial","Fibroblast","Endothelial","Mast"
)

# New-rerun annotation modules.  These are the primary annotation rule; the
# legacy cluster maps above are retained only for post-hoc concordance tables.
BROAD_LINEAGE_MODULES <- list(
  T_NK = c("CD3D","CD3E","TRAC","NKG7","GNLY","KLRD1"),
  B = c("CD79A","MS4A1","CD37","CD74","CD79B","CD19"),
  Plasma = c("JCHAIN","MZB1","SDC1","IGHG1","DERL3"),
  Myeloid_DC = c("LYZ","LST1","TYROBP","FCER1G","CTSS","AIF1"),
  Epithelial = c("EPCAM","KRT8","KRT18","KRT19","KRT20","CEACAM5"),
  Fibroblast = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA"),
  Endothelial = c("PECAM1","VWF","EMCN","KDR","ESAM","RAMP2"),
  Mast = c("TPSAB1","TPSB2","KIT","MS4A2","CPA3"),
  Erythroid = c("HBA1","HBA2","HBB","ALAS2"),
  Cycling = c("MKI67","TOP2A","STMN1","TUBB"),
  Stress = c("JUN","FOS","HSPA1A","HSPA1B","DNAJB1")
)

FINE_MYELOID_MODULES <- list(
  Mac_FOLR2_C1Q = c("FOLR2","C1QA","C1QB","C1QC","SLC40A1","VSIG4"),
  Mono_S100A8_A12 = c("S100A8","S100A9","FCN1","VCAN","CTSD"),
  Mac_SPP1_MARCO = c("SPP1","MARCO","LPL","TREM2","LGALS3"),
  Myeloid_TWISTNB_THBS1 = c("TWISTNB","THBS1","CTSD","LILRB1"),
  DC_CD1C_FCER1A = c("CD1C","FCER1A","CLEC10A","CST3"),
  Mac_SPP1_MARCO_LPL = c("SPP1","TREM2","GPNMB","APOC1","LPL"),
  Mac_CCL18_APOC1 = c("CCL18","APOC1","GPNMB","MMP9","MMP12"),
  Myeloid_stress = c("JUN","FOS","HSPA1A","HSPA1B","DNAJB1"),
  Myeloid_IFN_CLEC5A = c("IFITM3","ISG15","IFI6","CLEC5A","IFITM1"),
  # Explicit mixed-lineage / low-confidence modules are quality overlays;
  # they are retained for transparent auditing and excluded from the primary
  # biological state set below.
  Mixed_TNK_contamination = c("LYZ","LST1","CD3D","TRAC","NKG7","KLRD1","GZMK"),
  Mixed_low_confidence = c("ADGRE1","CD79B","LYPD2","TPPP3","CDKN1C"),
  Mixed_T_NK_myeloid = c("LYZ","LST1","CD3D","TRAC","NKG7"),
  Mono_ADGRE1 = c("ADGRE1","FCN1","CTSD","LYZ","S100A10"),
  Myeloid_cycling = c("MKI67","TOP2A","STMN1","TUBB"),
  DC_CCL19_CCL22 = c("CCL19","CCL22","CD1C","FCER1A","HLA-DRA"),
  DC_CLEC9A_XCR1 = c("CLEC9A","XCR1","BATF3","IRF8","HLA-DRA")
)

FINE_TNK_MODULES <- list(
  CD8_ZNF683_checkpoint_cytotoxic = c("CD8A","CD8B","ZNF683","HAVCR2","LAG3","GZMB"),
  CD8_IL7R_memory = c("CD3D","IL7R","LTB","MALAT1","CCR7"),
  T_naive_CCR7_LEF1 = c("CCR7","LEF1","TCF7","MAL","LTB"),
  NK_FGFBP2_PRF1 = c("NKG7","GNLY","FGFBP2","PRF1","GZMB"),
  Treg_FOXP3 = c("FOXP3","IL2RA","CTLA4","TNFRSF18","IKZF2"),
  CD8_CXCL13_checkpoint = c("CD8A","CXCL13","PDCD1","CTLA4","TOX"),
  CD8_GZMK_memory = c("CD8A","GZMK","IL7R","LTB","MALAT1"),
  T_gamma_delta = c("TRDC","TRGC1","TRGC2","CD3D","CD3E"),
  NK_FCER1G_KLRC1 = c("FCER1G","KLRD1","KLRC1","TYROBP","NKG7"),
  T_gamma_delta_NCR3 = c("TRDC","TRGC1","NCR3","KLRB1","CD3D"),
  T_stress = c("JUN","FOS","HSPA1A","DNAJB1","DUSP1"),
  T_activated_TNFSF4 = c("TNFSF4","ICOS","CD4","IL7R","CD3D")
)

# 外部 frozen module definitions
STATE_MODULES <- list(
  P1_state = c("SPP1","APOC1","APOE","GPNMB","LPL","MARCO"),
  P2_inflammatory = c("CXCL8","CXCL2","IL1B","CCL3","CCL20"),
  P3_C1Q_MHCII = c("C1QA","C1QB","C1QC","CD74","HLA-DRA"),
  Mono = c("S100A8","S100A9","FCN1","CTSD"),
  Other = c("LYZ","LST1","TYROBP","FCER1G","CTSS")
)
MYELOID_MODULE <- c("LYZ","LST1","TYROBP","FCER1G","CTSS","AIF1")
EPITHELIAL_MODULE <- c("EPCAM","KRT8","KRT18","KRT19","CEACAM5","CEACAM6")

SPATIAL_COVARIATE_MODULES <- list(
  myeloid = c("TYROBP","LST1","FCER1G","CTSS"),
  epithelial = c("EPCAM","KRT8","KRT18","KRT19","CEACAM5","CEACAM6"),
  liver = c("ALB","APOA1","CYP3A4","TTR","HP","SERPINA1")
)

# ---- 目录 ----
DIR <- list(
  inputs = file.path(BASE, "inputs"),
  cache_input = file.path(BASE, "cache", "input"),
  cache_annotation = file.path(BASE, "cache", "annotation"),
  cache_nmf_matrix = file.path(BASE, "cache", "nmf_matrix"),
  cache_nmf_fits = file.path(BASE, "cache", "nmf_fits"),
  cache_nmf_alignment = file.path(BASE, "cache", "nmf_alignment"),
  cache_projection = file.path(BASE, "cache", "projection"),
  cache_spatial = file.path(BASE, "cache", "spatial"),
  tables = file.path(BASE, "results", "tables"),
  figures = file.path(BASE, "results", "figures"),
  audit = file.path(BASE, "results", "audit"),
  manifest = file.path(BASE, "manifest"),
  logs = file.path(BASE, "logs")
)
invisible(lapply(DIR, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- 包检查：不自动安装，避免悄悄改变环境 ----
REQUIRED_PACKAGES <- c(
  "Seurat","Matrix","SingleCellExperiment","SummarizedExperiment",
  "scDblFinder","BiocParallel","harmony","NMF","data.table","ggplot2",
  "clue","digest","sandwich","jsonlite"
)
OPTIONAL_SPATIAL_PACKAGES <- c("dbscan","FNN")

check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "缺少R包：", paste(missing, collapse = ", "),
      "\n请先安装后再运行。脚本不会自动修改你的R环境。"
    )
  }
  invisible(TRUE)
}

check_packages(REQUIRED_PACKAGES)

# Cache signature: changes to the frozen parameter bundle, R, or key analysis
# package versions invalidate existing step caches automatically.
.cache_signature_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "scDblFinder", "harmony", "NMF"
)
.cache_signature_versions <- vapply(
  .cache_signature_packages,
  function(pkg) as.character(utils::packageVersion(pkg)),
  character(1)
)
# Fingerprint the actual files that can affect a run.  The mandatory discovery
# matrix is always hashed.  Small, directly supplied external files are also
# included when present; for very large external objects, size and mtime are
# retained without forcing a multi-gigabyte hash during every pipeline start.
.cache_signature_file_candidates <- c(
  GSE178318_MATRIX = GSE178318_MATRIX,
  GSE178318_BARCODES = GSE178318_BARCODES,
  GSE178318_FEATURES = GSE178318_FEATURES,
  GSE168933_REFERENCE_MARKER_FILE = GSE168933_REFERENCE_MARKER_FILE
)
.cache_signature_external_dirs <- c(
  GSE245552_DIR, GSE231559_DIR, GSE115469_DIR, GSE136103_DIR,
  GSE235079_DIR, GSE168933_DIR
)
.cache_signature_external_files <- unique(unlist(lapply(
  .cache_signature_external_dirs,
  function(root) {
    if (!dir.exists(root)) return(character())
    fs <- list.files(root, full.names = TRUE, recursive = FALSE,
                     include.dirs = FALSE)
    fs[grepl("\\.(csv|csv\\.gz|rds|rds\\.gz|h5|h5ad|mtx|mtx\\.gz|tsv|tsv\\.gz|json)$",
             fs, ignore.case = TRUE)]
  }
)))
if (length(.cache_signature_external_files)) {
  .cache_signature_file_candidates <- c(
    .cache_signature_file_candidates,
    stats::setNames(.cache_signature_external_files,
                    paste0("external_", seq_along(.cache_signature_external_files)))
  )
}
.cache_signature_file_candidates <- .cache_signature_file_candidates[
  !is.na(.cache_signature_file_candidates) &
    nzchar(.cache_signature_file_candidates)
]
.cache_signature_hash_limit <- 2 * 1024^3
.cache_signature_script_hash <- if (!is.na(.source_file) &&
                                      nzchar(.source_file) &&
                                      file.exists(.source_file)) {
  digest::digest(file = .source_file, algo = "sha256")
} else NA_character_
.cache_signature_manifest <- do.call(rbind, lapply(
  seq_along(.cache_signature_file_candidates),
  function(ii) {
    fp <- unname(.cache_signature_file_candidates[[ii]])
    fi <- if (file.exists(fp)) file.info(fp) else NULL
    sz <- if (!is.null(fi)) as.numeric(fi$size) else NA_real_
    mt <- if (!is.null(fi)) as.character(fi$mtime) else NA_character_
    sha <- if (!is.null(fi) && is.finite(sz) && sz <= .cache_signature_hash_limit) {
      digest::digest(file = fp, algo = "sha256")
    } else NA_character_
    data.frame(
      input = names(.cache_signature_file_candidates)[ii],
      path = fp,
      exists = file.exists(fp),
      size_bytes = sz,
      modified = mt,
      sha256 = sha,
      hash_status = if (!file.exists(fp)) "missing" else if (is.na(sha)) "size_mtime_only" else "sha256",
      stringsAsFactors = FALSE
    )
  }
))
CACHE_SIGNATURE <- digest::digest(
  list(
    pipeline = PIPELINE_VERSION,
    script_sha256 = .cache_signature_script_hash,
    PARAM = PARAM,
    R = R.version.string,
    packages = .cache_signature_versions,
    input_manifest = .cache_signature_manifest,
    input_expectations = list(
      strict = STRICT_EXPECTED_COUNTS,
      gene_count = EXPECTED_GSE178318_GENE_COUNT,
      cell_count = EXPECTED_GSE178318_CELL_COUNT,
      gene_count_tolerance = EXPECTED_GSE178318_GENE_COUNT_TOLERANCE
    ),
    marker_definitions = list(
      P1_ANCHOR_GENES = P1_ANCHOR_GENES,
      P1_MARKER_EXCLUDED_GENES = P1_MARKER_EXCLUDED_GENES,
      STATE_MODULES = STATE_MODULES,
      BROAD_LINEAGE_MODULES = BROAD_LINEAGE_MODULES,
      FINE_MYELOID_MODULES = FINE_MYELOID_MODULES,
      FINE_TNK_MODULES = FINE_TNK_MODULES,
      SPATIAL_COVARIATE_MODULES = SPATIAL_COVARIATE_MODULES
    )
  ),
  algo = "sha256"
)
message(
  "CRLM pipeline loaded. GSE178318 directory: ", GSE178318_DIR,
  "\nUse run_stage(0) for input preflight; Stage 0 is not an analysis stage."
)

writeLines(PIPELINE_VERSION, file.path(DIR$manifest, "pipeline_version.txt"))
writeLines(CACHE_SIGNATURE, file.path(DIR$manifest, "cache_signature.txt"))
writeLines(ifelse(is.na(.cache_signature_script_hash),
                  "script_hash_unavailable", .cache_signature_script_hash),
           file.path(DIR$manifest, "script_sha256.txt"))
write.csv(.cache_signature_manifest,
          file.path(DIR$manifest, "cache_input_manifest.csv"),
          row.names = FALSE)
jsonlite::write_json(PARAM, file.path(DIR$manifest, "parameter_freeze.json"),
                     auto_unbox = TRUE, pretty = TRUE)
git_hash <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) character()
)
writeLines(
  if (length(git_hash)) git_hash[[1]] else "git_unavailable_or_not_a_repository",
  file.path(DIR$manifest, "git_hash.txt")
)

# ---- 通用工具函数 ----
p <- function(...) file.path(...)

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

write_session_info <- function(tag) {
  writeLines(
    capture.output(sessionInfo()),
    p(DIR$audit, paste0(tag, "_sessionInfo.txt"))
  )
}

STAGE_ROUGH_ESTIMATES <- c(
  `01_import_qc` = "2–10 min",
  `02_samplewise_scdblfinder` = "10–45 min",
  `03_global_atlas_and_sensitivities` = "20–90 min",
  `04_all_lineage_subclustering_and_fine_maps` = "20–120 min",
  `05_core_patient_statistics_tcell_and_NMF_input` = "5–30 min",
  `06_multiseed_nmf` = "2–12 h",
  `07_fullmatrix_NMF_consensus_stability_rank_selection` = "20–120 min",
  `08_post_rank_P1_correspondence_and_freeze` = "5–30 min",
  `09_GSE245552_frozen_projection_all_programs` = "10–60 min",
  `10_GSE245552_full_patient_level_statistics` = "5–30 min",
  `11_GSE231559_directional_portability` = "10–60 min",
  `12_spatial_projection` = "10–90 min",
  `13_spatial_statistics` = "20–180 min",
  `14_discovery_supplementary_robustness` = "10–90 min",
  `14b_complexity_HC3_audit` = "5–30 min",
  `15_matched_random_signature_null_5000` = "30 min–6 h",
  `16_nonmalignant_liver_boundary` = "10–90 min",
  `17_spatial_supplementary_audits` = "20–180 min",
  `18_final_analysis_audit` = "5–30 min"
)

stage_rough_estimate <- function(stage_name) {
  x <- unname(STAGE_ROUGH_ESTIMATES[stage_name])
  if (!length(x) || is.na(x)) "未设定" else x
}

format_elapsed <- function(started_at, finished_at = Sys.time()) {
  round(as.numeric(difftime(finished_at, started_at, units = "secs")), 1)
}

append_stage_completion <- function(stage_name, status, started_at,
                                    message = NA_character_,
                                    estimated = stage_rough_estimate(stage_name),
                                    finished_at = Sys.time()) {
  f <- p(DIR$audit, "stage_completion.csv")
  row <- data.frame(
    stage = stage_name,
    status = status,
    started_at = format(started_at, tz = "Asia/Shanghai"),
    finished_at = format(finished_at, tz = "Asia/Shanghai"),
    elapsed_seconds = format_elapsed(started_at, finished_at),
    rough_estimate = estimated,
    R_version = R.version.string,
    pipeline_version = PIPELINE_VERSION,
    cache_signature = CACHE_SIGNATURE,
    message = as.character(message),
    stringsAsFactors = FALSE
  )
  if (file.exists(f)) {
    write.table(row, f, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE, qmethod = "double")
  } else {
    write.csv(row, f, row.names = FALSE)
  }
}

write_pipeline_config <- function() {
  cfg <- data.frame(
    item = c(
      "pipeline_version","cache_signature","base","GSE178318_matrix","GSE178318_barcodes",
      "GSE178318_features","min_features","max_percent_mt","mad_multiplier",
      "strict_expected_counts","expected_gene_count","expected_cell_count",
      "expected_gene_count_tolerance",
      "scdblfinder_dbr_per1k","scdblfinder_dbr_method","global_hvg","global_pcs","global_primary_resolution",
      "lineage_hvg","lineage_pcs","lineage_primary_resolution",
      "annotation_module_margin","annotation_margin_definition",
      "p1_anchor_min_margin",
      "nmf_hvg","nmf_ranks","nmf_runs_per_rank","nmf_method","nmf_max_iter",
      "nmf_loading_spearman_threshold","nmf_usage_spearman_threshold",
      "factor_recurrence_threshold",
      "top15_overlap_threshold",
      "P1_anchor_genes",
      "P1_core_anchor_genes",
      "P1_auxiliary_anchor_genes",
      "P1_marker_excluded_genes_provisional",
      "P1_marker_exclusion_confirmed",
      "eligible_rank_stable_factor_fraction",
      "n_boot","n_signflip","n_spatial_perm",
      "matched_null_n","matched_null_knn","matched_null_expr_caliper",
      "matched_null_detection_caliper","projection_normalization",
      "spatial_primary_graph","spatial_radius_um_sensitivity",
      "GSE245552_dir","GSE231559_dir","GSE115469_dir",
      "GSE136103_dir","GSE235079_dir","GSE168933_dir",
      "R_version","NMF_version","Seurat_version"
    ),
    value = c(
      PIPELINE_VERSION, CACHE_SIGNATURE, BASE, GSE178318_MATRIX, GSE178318_BARCODES,
      GSE178318_FEATURES, PARAM$min_features, PARAM$max_percent_mt,
      PARAM$mad_multiplier, STRICT_EXPECTED_COUNTS,
      EXPECTED_GSE178318_GENE_COUNT, EXPECTED_GSE178318_CELL_COUNT,
      EXPECTED_GSE178318_GENE_COUNT_TOLERANCE,
      PARAM$scdblfinder_dbr_per1k, PARAM$scdblfinder_dbr_method,
      PARAM$global_hvg, PARAM$global_pcs, PARAM$global_primary_resolution,
      PARAM$lineage_hvg, PARAM$lineage_pcs, PARAM$lineage_primary_resolution,
      PARAM$annotation_module_margin,
      "winner_minus_runner_up_on_cluster_mean_log1p_normalized_module_score; below margin -> unresolved/excluded",
      PARAM$p1_anchor_min_margin,
      PARAM$nmf_hvg, paste(PARAM$nmf_ranks, collapse = ","),
      PARAM$nmf_runs_per_rank, PARAM$nmf_method, PARAM$nmf_max_iter,
      PARAM$nmf_loading_spearman_threshold,
      PARAM$nmf_usage_spearman_threshold,
      PARAM$factor_recurrence_threshold,
      PARAM$top15_overlap_threshold,
      paste(P1_ANCHOR_GENES, collapse = ";"),
      paste(P1_CORE_ANCHOR_GENES, collapse = ";"),
      paste(P1_AUXILIARY_ANCHOR_GENES, collapse = ";"),
      paste(P1_MARKER_EXCLUDED_GENES, collapse = ";"),
      P1_MARKER_EXCLUSION_CONFIRMED,
      PARAM$eligible_rank_stable_factor_fraction,
      PARAM$n_boot, PARAM$n_signflip, PARAM$n_spatial_perm,
      PARAM$matched_null_n, PARAM$matched_null_knn,
      PARAM$matched_null_expr_caliper, PARAM$matched_null_detection_caliper,
      PARAM$projection_normalization,
      paste0("knn", PARAM$spatial_knn_k),
      paste(PARAM$spatial_radius_um_sensitivity, collapse = ";"),
      GSE245552_DIR, GSE231559_DIR, GSE115469_DIR,
      GSE136103_DIR, GSE235079_DIR, GSE168933_DIR,
      R.version.string,
      if (requireNamespace("NMF", quietly = TRUE)) as.character(packageVersion("NMF")) else NA_character_,
      if (requireNamespace("Seurat", quietly = TRUE)) as.character(packageVersion("Seurat")) else NA_character_
    ),
    stringsAsFactors = FALSE
  )
  write.csv(cfg, p(DIR$audit, "00_pipeline_config.csv"), row.names = FALSE)
}


# ============================================================
# Step-level console monitor
# Provides:
#   - current step name
#   - start/end time
#   - duration
#   - explicit error location
# ============================================================
run_step <- function(step_name, expr, stage_name = NA_character_,
                     step_i = NA_integer_, step_n = NA_integer_,
                     rough_estimate = NA_character_) {
  
  t0 <- Sys.time()
  
  cat("\n")
  cat("------------------------------------------------------------\n")
  prefix <- if (is.na(stage_name)) "" else paste0("[", stage_name, "] ")
  position <- if (is.finite(step_i) && is.finite(step_n))
    paste0(" step ", step_i, "/", step_n) else ""
  estimate <- if (is.na(rough_estimate)) "未设定" else rough_estimate
  cat(prefix, "START", position, ": ", step_name, "\n", sep = "")
  cat("CURRENT TIME: ", format(t0, tz="Asia/Shanghai"), "\n", sep = "")
  cat("ROUGH ESTIMATE: ", estimate,
      "（粗略估计，会随硬件/数据量变化）\n", sep = "")
  cat("------------------------------------------------------------\n")
  flush.console()
  
  out <- tryCatch({
    
    ans <- force(expr)
    
    t1 <- Sys.time()
    
    cat("\n")
    cat("DONE STEP :", step_name, "\n")
    cat("END TIME  :", format(t1, tz="Asia/Shanghai"), "\n")
    cat(
      "DURATION  :",
      round(as.numeric(difftime(t1, t0, units="mins")), 2),
      "minutes\n"
    )
    cat("------------------------------------------------------------\n")
    flush.console()
    
    ans
    
  }, error=function(e){
    
    cat("\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cat("ERROR STEP:", step_name, "\n")
    cat("TIME:", format(Sys.time(), tz="Asia/Shanghai"), "\n")
    cat("MESSAGE:\n")
    cat(conditionMessage(e), "\n")
    cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    flush.console()
    
    stop(e)
  })
  
  out
}

# Unified stage sub-step monitor.  It deliberately only wraps evaluation and
# console/log output; it does not alter seeds, objects, or analysis settings.
stage_step <- function(stage_name, step_i, step_n, description, expr,
                       rough_estimate = "未设定") {
  run_step(
    description, expr, stage_name = stage_name, step_i = step_i,
    step_n = step_n, rough_estimate = rough_estimate
  )
}

stage_log <- function(stage_name, expr) {
  log_file <- p(DIR$logs, paste0(stage_name, ".log"))
  zz <- file(log_file, open = "at", encoding = "UTF-8")
  sink(zz, split = TRUE)
  sink(zz, type = "message")
  started_at <- Sys.time()
  rough_estimate <- stage_rough_estimate(stage_name)
  stage_status <- "running"
  completion_written <- FALSE
  append_stage_completion(stage_name, "started", started_at,
                          estimated = rough_estimate)
  on.exit({
    # `return()` inside the stage expression also exits stage_log directly.
    # Finalize that path here so cache hits receive the same DONE record as a
    # normal fall-through, while failed stages are not marked completed twice.
    if (!completion_written && identical(stage_status, "running")) {
      finished_at <- Sys.time()
      append_stage_completion(stage_name, "completed", started_at,
                              estimated = rough_estimate,
                              finished_at = finished_at)
      cat("[", stage_name, "] DONE | completed time: ",
          format(finished_at, tz = "Asia/Shanghai"),
          " | actual elapsed: ", format_elapsed(started_at, finished_at),
          " seconds | early return\n", sep = "")
      completion_written <- TRUE
    }
    sink(type = "message")
    sink()
    close(zz)
  }, add = TRUE)
  cat("\n============================================================\n")
  cat("[", stage_name, "] START | current time: ",
      format(started_at, tz = "Asia/Shanghai"),
      " | rough estimate: ", rough_estimate,
      "（粗略估计，会随硬件/数据量变化）\n", sep = "")
  ans <- tryCatch(
    force(expr),
    error = function(e) {
      cat("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
      failed_at <- Sys.time()
      cat("FAILED STAGE:", stage_name, "\n")
      cat("FAIL TIME:", format(failed_at, tz="Asia/Shanghai"), "\n")
      cat("ELAPSED:", format_elapsed(started_at, failed_at), "seconds\n")
      cat("MESSAGE:\n")
      cat(conditionMessage(e), "\n")
      cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
      stage_status <<- "failed"
      append_stage_completion(stage_name, "failed", started_at,
                              conditionMessage(e), rough_estimate, failed_at)
      completion_written <<- TRUE
      stop(e)
    }
  )
  finished_at <- Sys.time()
  append_stage_completion(stage_name, "completed", started_at,
                          estimated = rough_estimate,
                          finished_at = finished_at)
  completion_written <- TRUE
  cat("[", stage_name, "] DONE | completed time: ",
      format(finished_at, tz = "Asia/Shanghai"),
      " | actual elapsed: ", format_elapsed(started_at, finished_at),
      " seconds\n", sep = "")
  ans
}

assert_file <- function(path, label = path) {
  if (!file.exists(path)) stop("缺少文件：", label, "\n", path)
  invisible(TRUE)
}

# Explicit name: sparse log1p(CP10K), matching Seurat LogNormalize's
# scale.factor=10000 transformation while avoiding a hidden assay side effect.
normalize_log1p_cp10k_sparse <- function(X) {
  X <- methods::as(X, "dgCMatrix")
  lib <- Matrix::colSums(X)
  sf <- 1e4 / pmax(lib, 1)
  Y <- X %*% Matrix::Diagonal(x = sf)
  Y@x <- log1p(Y@x)
  Y
}

collapse_duplicate_genes <- function(X, genes) {
  keep <- !is.na(genes) & nzchar(genes)
  X <- X[keep, , drop = FALSE]
  genes <- genes[keep]
  lev <- unique(genes)
  if (!anyDuplicated(genes)) {
    rownames(X) <- genes
    return(X)
  }
  grp <- match(genes, lev)
  A <- Matrix::sparseMatrix(
    i = grp,
    j = seq_along(grp),
    x = 1,
    dims = c(length(lev), length(grp))
  )
  out <- A %*% X
  rownames(out) <- lev
  out
}

read_10x_triplet <- function(matrix_file, barcode_file, feature_file) {
  assert_file(matrix_file)
  assert_file(barcode_file)
  assert_file(feature_file)
  X <- Matrix::readMM(gzfile(matrix_file))
  X <- methods::as(X, "dgCMatrix")
  bar <- read.delim(gzfile(barcode_file), header = FALSE, stringsAsFactors = FALSE)
  feat <- read.delim(gzfile(feature_file), header = FALSE, stringsAsFactors = FALSE)
  genes <- if (ncol(feat) >= 2) feat[[2]] else feat[[1]]
  stopifnot(nrow(X) == length(genes), ncol(X) == nrow(bar))
  X <- collapse_duplicate_genes(X, as.character(genes))
  colnames(X) <- as.character(bar[[1]])
  X
}

module_score <- function(norm, genes) {
  g <- intersect(genes, rownames(norm))
  if (!length(g)) return(rep(0, ncol(norm)))
  Matrix::colMeans(norm[g, , drop = FALSE])
}

cluster_module_scores <- function(obj, cluster_col, modules) {
  md <- obj[[]]
  cl <- as.character(md[[cluster_col]])
  clusters <- sort(unique(cl))
  expr <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
  avg <- do.call(cbind, lapply(clusters, function(cc) {
    Matrix::rowMeans(expr[, cl == cc, drop = FALSE])
  }))
  colnames(avg) <- clusters
  scores <- do.call(cbind, lapply(modules, function(gs) module_score(avg, gs)))
  colnames(scores) <- names(modules)
  rownames(scores) <- clusters
  list(clusters = clusters, average = avg, scores = scores)
}

annotation_score_audit <- function(scores, assigned_label, annotation_status,
                                   level, lineage = NA_character_,
                                   margin_threshold = PARAM$annotation_module_margin) {
  stopifnot(nrow(scores) == length(assigned_label),
            nrow(scores) == length(annotation_status))
  rows <- lapply(seq_len(nrow(scores)), function(ii) {
    v <- scores[ii, ]
    ok <- is.finite(v)
    if (any(ok)) {
      vv <- v[ok]
      labs <- colnames(scores)[ok]
      oo <- order(vv, decreasing = TRUE)
      winner <- labs[oo[1L]]
      winner_score <- unname(vv[oo[1L]])
      if (length(oo) > 1L) {
        runner_up <- labs[oo[2L]]
        runner_up_score <- unname(vv[oo[2L]])
        margin <- winner_score - runner_up_score
      } else {
        runner_up <- NA_character_
        runner_up_score <- NA_real_
        margin <- NA_real_
      }
    } else {
      winner <- NA_character_
      runner_up <- NA_character_
      winner_score <- NA_real_
      runner_up_score <- NA_real_
      margin <- NA_real_
    }
    data.frame(
      level = level,
      lineage = lineage,
      cluster = rownames(scores)[ii],
      winner = winner,
      runner_up = runner_up,
      assigned_label = unname(assigned_label[ii]),
      winner_score = winner_score,
      runner_up_score = runner_up_score,
      margin = margin,
      margin_threshold = margin_threshold,
      annotation_status = unname(annotation_status[ii]),
      unresolved = unname(assigned_label[ii]) %in% c("Unresolved", "Unknown_mixed"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

weighted_score <- function(norm, weights_df, genes_to_use = NULL) {
  w <- weights_df
  if (!is.null(genes_to_use)) w <- w[w$gene %in% genes_to_use, , drop = FALSE]
  w <- w[w$gene %in% rownames(norm), , drop = FALSE]
  if (!nrow(w)) return(rep(NA_real_, ncol(norm)))
  ww <- w$weight
  ww <- ww / sum(ww)
  as.numeric(Matrix::t(norm[w$gene, , drop = FALSE]) %*% ww)
}

# P1 score-set factory.  Full P1 and the single-gene SPP1 sensitivity are
# available immediately after Stage 8.  The multi-gene marker-excluded set is
# deliberately withheld until the new P1 top-15 has been inspected and the
# exclusion rule has been explicitly confirmed.
make_p1_score_sets <- function(p1_genes) {
  sets <- list(
    P1_full = p1_genes,
    P1_no_SPP1 = setdiff(p1_genes, "SPP1")
  )
  if (isTRUE(P1_MARKER_EXCLUSION_CONFIRMED)) {
    sets$P1_marker_excluded <- setdiff(p1_genes, P1_MARKER_EXCLUDED_GENES)
  }
  sets
}

bootstrap_median_ci <- function(x, B = PARAM$n_boot, seed = PARAM$seed + 900000L) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  set.seed(seed)
  b <- replicate(B, median(sample(x, length(x), replace = TRUE)))
  as.numeric(stats::quantile(b, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

exact_sign_p_greater <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  stats::binom.test(sum(x > 0), length(x), p = 0.5, alternative = "greater")$p.value
}

signflip_two_sided <- function(x, B = PARAM$n_signflip, seed = PARAM$seed + 910000L) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(NA_real_)
  obs <- mean(x)
  # 小样本可完整枚举
  if (n <= 18L) {
    grid <- expand.grid(rep(list(c(-1, 1)), n))
    vals <- apply(as.matrix(grid), 1, function(s) mean(x * s))
    return(mean(abs(vals) >= abs(obs)))
  }
  set.seed(seed)
  vals <- replicate(B, mean(x * sample(c(-1,1), n, replace = TRUE)))
  (1 + sum(abs(vals) >= abs(obs))) / (1 + B)
}

# NMF 的 W/H 存在尺度不唯一性。
# 统一把每个 W column 归一化到和为1，同时相应放大 H row，保持 W%*%H 不变。

score_genes_z <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(rep(0, ncol(obj)))
  mat <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE]
  z <- do.call(rbind, lapply(seq_len(nrow(mat)), function(i) {
    v <- as.numeric(mat[i, ])
    s <- stats::sd(v)
    if (!is.finite(s) || s == 0) rep(0, length(v)) else (v - mean(v)) / s
  }))
  sc <- colMeans(z)
  ss <- stats::sd(sc)
  if (!is.finite(ss) || ss == 0) return(rep(0, length(sc)))
  as.numeric((sc - mean(sc)) / ss)
}

wilcox_signed_rank_auto <- function(delta, alternative = "two.sided") {
  delta <- delta[is.finite(delta)]
  if (length(delta) < 2L || all(delta == 0)) {
    return(list(statistic = 0, p.value = 1, exact_used = FALSE))
  }
  nz <- delta[delta != 0]
  # R 的 exact signed-rank 要求没有ties/zero。能exact时显式exact；
  # 否则按预设的 exact/渐近检验规则转为渐近法。
  can_exact <- length(nz) > 0L &&
    length(nz) < 50L &&
    !anyDuplicated(abs(nz))
  ans <- suppressWarnings(stats::wilcox.test(
    delta, mu = 0, alternative = alternative,
    exact = can_exact, correct = FALSE
  ))
  list(
    statistic = unname(ans$statistic),
    p.value = ans$p.value,
    exact_used = can_exact
  )
}

plus_one_p <- function(k, n) (1 + k) / (1 + n)

safe_cor <- function(x, y, method = "spearman") {
  suppressWarnings(stats::cor(x, y, method = method, use = "pairwise.complete.obs"))
}

normalize_nmf_scale <- function(W, H) {
  s <- colSums(W)
  s[!is.finite(s) | s <= 0] <- 1
  Wn <- sweep(W, 2, s, "/")
  Hn <- sweep(H, 1, s, "*")
  list(W = Wn, H = Hn)
}

spearman_matrix <- function(A, B) {
  # A/B: gene x factor
  out <- matrix(NA_real_, ncol(A), ncol(B))
  for (i in seq_len(ncol(A))) {
    for (j in seq_len(ncol(B))) {
      out[i, j] <- suppressWarnings(cor(A[, i], B[, j], method = "spearman"))
    }
  }
  out[!is.finite(out)] <- -1
  out
}

match_factors_lsap <- function(W_ref, W_other) {
  common <- intersect(rownames(W_ref), rownames(W_other))
  if (length(common) < 100) stop("LSAP因子匹配的共同基因过少。")
  S <- spearman_matrix(W_ref[common, , drop = FALSE], W_other[common, , drop = FALSE])
  # solve_LSAP在部分clue版本要求非负矩阵；整体平移不改变最大匹配解。
  S_lsap <- S - min(S, na.rm = TRUE)
  perm <- as.integer(clue::solve_LSAP(S_lsap, maximum = TRUE))
  data.frame(
    factor_ref = seq_len(ncol(W_ref)),
    factor_other = perm,
    basis_spearman = S[cbind(seq_len(ncol(W_ref)), perm)],
    stringsAsFactors = FALSE
  )
}

top_gene_overlap <- function(W1, f1, W2, f2, n = 15L) {
  g1 <- names(sort(W1[, f1], decreasing = TRUE))[seq_len(min(n, nrow(W1)))]
  g2 <- names(sort(W2[, f2], decreasing = TRUE))[seq_len(min(n, nrow(W2)))]
  length(intersect(g1, g2)) / n
}

compute_rss <- function(X, W, H) {
  sum((X - W %*% H)^2)
}

compute_elbow_rank <- function(rank, rss) {
  ok <- is.finite(rank) & is.finite(rss) & rss > 0
  rank <- rank[ok]; y <- log(rss[ok])
  if (length(rank) < 3) return(rank[which.min(rss[ok])])
  x1 <- rank[1]; y1 <- y[1]
  x2 <- rank[length(rank)]; y2 <- y[length(y)]
  denom <- sqrt((y2-y1)^2 + (x2-x1)^2)
  d <- abs((y2-y1)*rank - (x2-x1)*y + x2*y1 - y2*x1) / denom
  rank[which.max(d)]
}

select_consensus_cells <- function(metadata, target = 999L, seed = 257319L) {
  n <- nrow(metadata)
  if (!n) return(integer())
  target <- min(as.integer(target), n)
  strata <- if ("fine_label" %in% names(metadata))
    as.character(metadata$fine_label) else rep("all", n)
  strata[is.na(strata) | !nzchar(strata)] <- "unresolved"
  groups <- split(seq_len(n), strata)
  set.seed(seed)
  groups <- lapply(groups, function(ii) sample(ii, length(ii), replace = FALSE))
  take <- integer()
  pos <- setNames(rep(1L, length(groups)), names(groups))
  while (length(take) < target) {
    progressed <- FALSE
    for (nm in names(groups)) {
      jj <- pos[[nm]]
      if (jj <= length(groups[[nm]]) && length(take) < target) {
        take <- c(take, groups[[nm]][jj])
        pos[[nm]] <- jj + 1L
        progressed <- TRUE
      }
    }
    if (!progressed) break
  }
  sort(take)
}

write_pipeline_config()
write.csv(
  data.frame(
    stage = c("global_atlas", "lineage_reclustering", "NMF_rank_run", "NMF_screening",
              "GSE245552_bootstrap", "GSE231559_bootstrap", "spatial_permutation",
              "matched_random_null"),
    seed_rule = c(
      "PARAM$seed + 500",
      "PARAM$seed + 800",
      "PARAM$seed + 600000 + rank*1000 + run",
      "PARAM$seed + 1301",
      "20260813 + version_index",
      "20260815 + pair_index",
      "cohort-specific explicit seed passed to permutation function",
      "20260815"
    ),
    stringsAsFactors = FALSE
  ),
  p(DIR$audit, "random_seed_manifest.csv"),
  row.names = FALSE
)
write_session_info("00_initial")
write.csv(LEGACY_GLOBAL_MAP, p(DIR$audit, "00_legacy_global_map.csv"), row.names = FALSE)

write.csv(LEGACY_MYELOID_MAP, p(DIR$audit, "00_legacy_myeloid_map.csv"), row.names = FALSE)
write.csv(LEGACY_TNK_MAP, p(DIR$audit, "00_legacy_tnk_map.csv"), row.names = FALSE)
write.csv(LEGACY_P1_TOP15, p(DIR$audit, "00_legacy_P1_top15_reference.csv"), row.names = FALSE)
write.csv(
  do.call(rbind, lapply(names(STATE_MODULES), function(st) {
    data.frame(state = st, marker_gene = STATE_MODULES[[st]], stringsAsFactors = FALSE)
  })),
  p(DIR$audit, "00_frozen_state_marker_modules.csv"), row.names = FALSE
)
write.csv(
  do.call(rbind, lapply(names(SPATIAL_COVARIATE_MODULES), function(mod) {
    data.frame(module = mod, gene = SPATIAL_COVARIATE_MODULES[[mod]], stringsAsFactors = FALSE)
  })),
  p(DIR$audit, "00_spatial_covariate_modules.csv"), row.names = FALSE
)
write.csv(
  do.call(rbind, lapply(names(BROAD_LINEAGE_MODULES), function(mod) {
    data.frame(module_family = "broad_lineage", module = mod,
               marker_gene = BROAD_LINEAGE_MODULES[[mod]], stringsAsFactors = FALSE)
  })),
  p(DIR$audit, "00_new_annotation_modules_broad.csv"), row.names = FALSE
)
write.csv(
  rbind(
    do.call(rbind, lapply(names(FINE_MYELOID_MODULES), function(mod) {
      data.frame(module_family = "fine_myeloid", module = mod,
                 marker_gene = FINE_MYELOID_MODULES[[mod]], stringsAsFactors = FALSE)
    })),
    do.call(rbind, lapply(names(FINE_TNK_MODULES), function(mod) {
      data.frame(module_family = "fine_TNK", module = mod,
                 marker_gene = FINE_TNK_MODULES[[mod]], stringsAsFactors = FALSE)
    }))
  ),
  p(DIR$audit, "00_new_annotation_modules_fine.csv"), row.names = FALSE
)
write.csv(
  data.frame(
    projection_normalization = PARAM$projection_normalization,
    formula = "log1p(10000 * raw_counts / per_cell_total_counts)",
    score_formula = "sum(normalized_expression_i * frozen_weight_i) after renormalizing weights over available genes",
    external_retraining = FALSE,
    stringsAsFactors = FALSE
  ),
  p(DIR$audit, "00_projection_normalization_manifest.csv"), row.names = FALSE
)


############################
# Stage 1：GSE178318 原始导入 + 冻结 QC
############################

stage1_import_qc <- function(force = FALSE) {
  stage_log("01_import_qc", {
    stage_step("01_import_qc", 1L, 5L,
      "检查输入文件、参数包与输入清单",
      {
        check_packages(c("Seurat","Matrix","digest"))
        if (!dir.exists(GSE178318_DIR)) {
          stop("GSE178318目录不存在：", GSE178318_DIR)
        }
        invisible(lapply(
                   c(GSE178318_MATRIX, GSE178318_BARCODES, GSE178318_FEATURES),
                   assert_file
                 ))
                 message(
                   "Stage 1 input files:\n  matrix: ", GSE178318_MATRIX,
                   "\n  barcodes: ", GSE178318_BARCODES,
                   "\n  features: ", GSE178318_FEATURES
                 )
                 input_manifest <- data.frame(
                   file = c("matrix","barcodes","features"),
                   path = c(GSE178318_MATRIX,GSE178318_BARCODES,GSE178318_FEATURES),
                   size_bytes = file.info(c(GSE178318_MATRIX,GSE178318_BARCODES,GSE178318_FEATURES))$size,
                   sha256 = vapply(c(GSE178318_MATRIX, GSE178318_BARCODES, GSE178318_FEATURES),
                                   sha256_file, character(1)),
                   stringsAsFactors = FALSE
                 )
                 write.csv(input_manifest, p(DIR$audit, "01_GSE178318_input_manifest.csv"), row.names = FALSE)
               }, "<5 min")
    
    out_all <- p(DIR$cache_input, "01_GSE178318_all_140281_with_qc.rds")
    out_qc  <- p(DIR$cache_input, "02_GSE178318_primary_qc.rds")
    if (!force && cache_has_current_version(out_all) &&
        cache_has_current_version(out_qc)) {
      message("[01_import_qc] CACHE HIT：直接返回 QC 对象。")
      return(load_step_cache(out_qc))
    }
    
    counts <- stage_step(
      "01_import_qc", 2L, 5L,
      "读取 GSE178318 10x 矩阵并核对细胞/基因数",
      {
        set.seed(PARAM$seed)
        ans <- Seurat::ReadMtx(
          mtx = GSE178318_MATRIX,
          cells = GSE178318_BARCODES,
          features = GSE178318_FEATURES,
          feature.column = 2,
          cell.column = 1,
          unique.features = TRUE,
          strip.suffix = FALSE
        )
        if (anyNA(rownames(ans)) || anyDuplicated(rownames(ans))) {
          stop("GSE178318 feature names contain missing or duplicated entries。")
        }
        if (anyNA(colnames(ans)) || anyDuplicated(colnames(ans))) {
          stop("GSE178318 cell barcodes contain missing or duplicated entries。")
        }
        if (STRICT_EXPECTED_COUNTS) {
          if (ncol(ans) != EXPECTED_GSE178318_CELL_COUNT) {
            stop(
              "GSE178318细胞数不符合预期：observed=", ncol(ans),
              "; expected=", EXPECTED_GSE178318_CELL_COUNT,
              "。请确认输入矩阵未被替换。"
            )
          }
          gene_lower <- floor(
            EXPECTED_GSE178318_GENE_COUNT *
              (1 - EXPECTED_GSE178318_GENE_COUNT_TOLERANCE)
          )
          gene_upper <- ceiling(
            EXPECTED_GSE178318_GENE_COUNT *
              (1 + EXPECTED_GSE178318_GENE_COUNT_TOLERANCE)
          )
          if (nrow(ans) < gene_lower || nrow(ans) > gene_upper) {
            stop(
              "GSE178318 feature数超出允许范围：observed=", nrow(ans),
              "; expected=", EXPECTED_GSE178318_GENE_COUNT,
              " +/- ", EXPECTED_GSE178318_GENE_COUNT_TOLERANCE * 100,
              "%。请确认输入矩阵未被替换。"
            )
          }
        }
        ans
      }, "1–8 min")
    
    obj <- stage_step(
      "01_import_qc", 3L, 5L,
      "创建 Seurat 对象",
      {
        ans <- Seurat::CreateSeuratObject(
          counts = counts,
          assay = "RNA",
          project = "GSE178318_formal_rerun",
          min.cells = 0,
          min.features = 0
        )
        ans
      }, "1–5 min")
    rm(counts); invisible(gc())
    qc_step4_started_at <- Sys.time()
    cat("[01_import_qc] STEP 4/5 IN PROGRESS | 解析 barcode、计算 percent.mt/HB/ribo、应用 MAD/QC 规则 | current time: ",
        format(qc_step4_started_at, tz = "Asia/Shanghai"),
        " | rough estimate: 1–5 min（粗略估计，会随硬件/数据量变化）\n", sep = "")
    
    bc <- colnames(obj)
    rx <- "^([ACGTN]+)_((COL[0-9]+)_(CRC|LM|PBMC))$"
    if (!all(grepl(rx, bc))) stop("GSE178318 barcode格式与冻结解析规则不一致。")
    if (anyDuplicated(bc)) stop("发现重复barcode。")
    
    obj$sample_id  <- sub(rx, "\\2", bc)
    obj$patient_id <- sub(rx, "\\3", bc)
    obj$tissue     <- sub(rx, "\\4", bc)
    stopifnot(length(unique(obj$sample_id)) == 15L,
              length(unique(obj$patient_id)) == 6L)
    
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
    hb_genes <- grep("^(HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ)$",
                     rownames(obj), value = TRUE)
    ribo_genes <- grep("^(RPL|RPS)", rownames(obj), value = TRUE)
    obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, features = hb_genes)
    obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(obj, features = ribo_genes)
    
    meta <- obj[[]]
    meta$barcode <- rownames(meta)
    meta$log10_nFeature <- log10(meta$nFeature_RNA)
    meta$log10_nCount <- log10(meta$nCount_RNA)
    
    sample_levels <- sort(unique(meta$sample_id))
    bounds <- do.call(rbind, lapply(sample_levels, function(s) {
      d <- meta[meta$sample_id == s, , drop = FALSE]
      med_f <- median(d$log10_nFeature)
      mad_f <- mad(d$log10_nFeature)
      med_c <- median(d$log10_nCount)
      mad_c <- mad(d$log10_nCount)
      if (!is.finite(mad_f) || !is.finite(mad_c) || mad_f <= 0 || mad_c <= 0) {
        stop("样本 ", s, " 的MAD不可计算。")
      }
      data.frame(
        sample_id = s, patient_id = d$patient_id[1], tissue = d$tissue[1],
        n_raw = nrow(d),
        feature_median = med_f, feature_mad = mad_f,
        feature_lower = med_f - PARAM$mad_multiplier * mad_f,
        feature_upper = med_f + PARAM$mad_multiplier * mad_f,
        count_median = med_c, count_mad = mad_c,
        count_lower = med_c - PARAM$mad_multiplier * mad_c,
        count_upper = med_c + PARAM$mad_multiplier * mad_c
      )
    }))
    
    i <- match(meta$sample_id, bounds$sample_id)
    meta$pass_fixed_features <- meta$nFeature_RNA >= PARAM$min_features
    meta$pass_fixed_mt <- meta$percent.mt <= PARAM$max_percent_mt
    meta$pass_mad_features <- meta$log10_nFeature >= bounds$feature_lower[i] &
      meta$log10_nFeature <= bounds$feature_upper[i]
    meta$pass_mad_counts <- meta$log10_nCount >= bounds$count_lower[i] &
      meta$log10_nCount <= bounds$count_upper[i]
    meta$pass_primary_qc <- with(
      meta,
      pass_fixed_features & pass_fixed_mt & pass_mad_features & pass_mad_counts
    )
    
    for (nm in c("pass_fixed_features","pass_fixed_mt","pass_mad_features",
                 "pass_mad_counts","pass_primary_qc")) {
      obj[[nm]] <- meta[[nm]]
    }
    
    write.csv(bounds, p(DIR$tables, "01_sample_mad_boundaries.csv"), row.names = FALSE)
    write.csv(bounds, p(DIR$audit, "QC_boundary_manifest.csv"), row.names = FALSE)
    
    qc_summary <- do.call(rbind, lapply(sample_levels, function(s) {
      d <- meta[meta$sample_id == s, , drop = FALSE]
      data.frame(
        sample_id = s,
        patient_id = d$patient_id[1],
        tissue = d$tissue[1],
        n_raw = nrow(d),
        n_primary_qc = sum(d$pass_primary_qc),
        retention = mean(d$pass_primary_qc)
      )
    }))
    write.csv(qc_summary, p(DIR$tables, "02_qc_retention_by_sample.csv"), row.names = FALSE)
    
    qc_patient <- aggregate(
      cbind(n_raw, n_primary_qc) ~ patient_id,
      qc_summary, sum
    )
    qc_patient$retention <- qc_patient$n_primary_qc / qc_patient$n_raw
    write.csv(
      qc_patient,
      p(DIR$tables, "02b_qc_retention_by_patient.csv"),
      row.names = FALSE
    )
    qc_tissue <- aggregate(
      cbind(n_raw, n_primary_qc) ~ tissue,
      qc_summary, sum
    )
    qc_tissue$retention <- qc_tissue$n_primary_qc / qc_tissue$n_raw
    write.csv(
      qc_tissue,
      p(DIR$tables, "02c_qc_retention_by_tissue.csv"),
      row.names = FALSE
    )
    
    con <- gzfile(p(DIR$tables, "03_cell_qc_metadata.csv.gz"), "wt")
    write.csv(meta, con, row.names = FALSE)
    close(con)
    qc_step4_finished_at <- Sys.time()
    cat("[01_import_qc] STEP 4/5 DONE | QC 指标与元数据已完成 | completed time: ",
        format(qc_step4_finished_at, tz = "Asia/Shanghai"),
        " | actual elapsed: ", format_elapsed(qc_step4_started_at, qc_step4_finished_at),
        " seconds\n", sep = "")
    
    stage_step(
      "01_import_qc", 5L, 5L,
      "写入 QC 汇总、筛选通过 QC 的细胞并写入断点缓存",
      {
        save_step_cache(obj, out_all)
        obj_qc <- subset(obj, cells = meta$barcode[meta$pass_primary_qc])
        save_step_cache(obj_qc, out_qc)
        invisible(NULL)
      }, "1–5 min")
    
    write_session_info("01_import_qc")
    cat("Raw cells:", ncol(obj), "\n")
    cat("Primary QC cells:", ncol(obj_qc), "\n")
    obj_qc
  })
}


############################
# Stage 2：15 文库分别 scDblFinder
############################

stage2_doublet <- function(force = FALSE) {
  stage_log("02_samplewise_scdblfinder", {
    check_packages(c("Seurat","SingleCellExperiment","scDblFinder","BiocParallel"))
    
    input <- p(DIR$cache_input, "02_GSE178318_primary_qc.rds")
    assert_file(input, "Stage1 QC对象")
    if (!cache_has_current_version(input)) {
      stop("Stage 2 输入 QC 缓存缺失版本标记或版本不匹配；请先重新运行 Stage 1。")
    }
    out <- p(DIR$cache_input, "03_GSE178318_qc_singlets.rds")
    if (!force && cache_has_current_version(out)) {
      message("[02_samplewise_scdblfinder] CACHE HIT：直接返回 singlet 对象。")
      return(load_step_cache(out))
    }
    
    obj <- load_step_cache(input)
    samples <- sort(unique(obj$sample_id))
    stopifnot(length(samples) == 15L)
    message("Stage 2 starting sample-wise scDblFinder for ", length(samples), " libraries.")
    
    cache_dir <- p(DIR$cache_input, "doublet_by_sample")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    
    calls <- vector("list", length(samples))
    for (j in seq_along(samples)) {
      s <- samples[j]
      cf <- p(cache_dir, paste0(s, "_scDblFinder.rds"))
      if (!force && cache_has_current_version(cf)) {
        calls[[j]] <- load_step_cache(cf)
        cat(sprintf("[02_samplewise_scdblfinder] STEP %02d/15 CACHE HIT | %s | current time: %s\n",
                    j, s, format(Sys.time(), tz = "Asia/Shanghai")))
        next
      }
      
      cells <- colnames(obj)[obj$sample_id == s]
      sub <- obj[, cells]
      sce <- Seurat::as.SingleCellExperiment(sub, assay = "RNA")
      run_seed <- PARAM$seed + j
      set.seed(run_seed)
      sce <- stage_step(
        "02_samplewise_scdblfinder", j, length(samples),
        paste0("scDblFinder：", s, "（", length(cells), " QC-passing cells）"),
        {
          scDblFinder::scDblFinder(
            sce,
            dbr.per1k = PARAM$scdblfinder_dbr_per1k,
            BPPARAM = BiocParallel::SerialParam()
          )
        }, "1–10 min"
      )
      z <- data.frame(
        barcode = colnames(sce),
        sample_id = s,
        patient_id = as.character(sub$patient_id),
        tissue = as.character(sub$tissue),
        scDblFinder.score = SummarizedExperiment::colData(sce)$scDblFinder.score,
        scDblFinder.class = as.character(
          SummarizedExperiment::colData(sce)$scDblFinder.class
        ),
        seed = run_seed,
        stringsAsFactors = FALSE
      )
      save_step_cache(z, cf)
      calls[[j]] <- z
      rm(sub, sce); invisible(gc())
    }
    
    calls <- do.call(rbind, calls)
    stopifnot(nrow(calls) == ncol(obj), !anyDuplicated(calls$barcode))
    idx <- match(colnames(obj), calls$barcode)
    stopifnot(!anyNA(idx))
    obj$scDblFinder.score <- calls$scDblFinder.score[idx]
    obj$scDblFinder.class <- calls$scDblFinder.class[idx]
    obj$pass_singlet <- obj$scDblFinder.class == "singlet"
    
    con <- gzfile(p(DIR$tables, "04_scdblfinder_cell_calls.csv.gz"), "wt")
    write.csv(calls, con, row.names = FALSE); close(con)
    
    sm <- do.call(rbind, lapply(samples, function(s) {
      d <- calls[calls$sample_id == s, , drop = FALSE]
      data.frame(
        sample_id = s,
        patient_id = d$patient_id[1],
        tissue = d$tissue[1],
        n_qc = nrow(d),
        expected_doublet_rate = (nrow(d) / 1000) * PARAM$scdblfinder_dbr_per1k,
        expected_doublets = nrow(d) * (nrow(d) / 1000) * PARAM$scdblfinder_dbr_per1k,
        n_singlet = sum(d$scDblFinder.class == "singlet"),
        n_doublet = sum(d$scDblFinder.class != "singlet"),
        singlet_rate = mean(d$scDblFinder.class == "singlet"),
        dbr_per1k = PARAM$scdblfinder_dbr_per1k,
        dbr_method = PARAM$scdblfinder_dbr_method
      )
    }))
    write.csv(sm, p(DIR$tables, "05_scdblfinder_by_sample.csv"), row.names = FALSE)
    
    db_patient <- aggregate(
      cbind(n_qc, n_singlet, n_doublet) ~ patient_id,
      sm, sum
    )
    db_patient$singlet_rate <- db_patient$n_singlet / db_patient$n_qc
    write.csv(
      db_patient,
      p(DIR$tables, "05b_scdblfinder_by_patient.csv"),
      row.names = FALSE
    )
    db_tissue <- aggregate(
      cbind(n_qc, n_singlet, n_doublet) ~ tissue,
      sm, sum
    )
    db_tissue$singlet_rate <- db_tissue$n_singlet / db_tissue$n_qc
    write.csv(
      db_tissue,
      p(DIR$tables, "05c_scdblfinder_by_tissue.csv"),
      row.names = FALSE
    )
    
    singlets <- colnames(obj)[obj$pass_singlet]
    obj_s <- subset(obj, cells = singlets)
    singlet_audit <- data.frame(
      n_qc_cells = ncol(obj), n_singlets = ncol(obj_s),
      n_doublets = ncol(obj) - ncol(obj_s),
      legacy_reference_singlets = 108353L,
      count_is_descriptive = TRUE,
      note = "New rerun does not force the legacy singlet count; the fixed QC and sample-wise doublet rule determines retention."
    )
    write.csv(singlet_audit, p(DIR$audit, "02_singlet_count_audit.csv"), row.names = FALSE)
    save_step_cache(obj_s, out)
    write_session_info("02_scdblfinder")
    
    cat("QC cells:", ncol(obj), "\n")
    cat("Singlets:", ncol(obj_s), "\n")
    obj_s
  })
}


############################
# Stage 3：全局 atlas（冻结主参数 + 新版注释审计）
############################

stage3_global_atlas <- function(force = FALSE) {
  stage_log("03_global_atlas_and_sensitivities", {
    future::plan("sequential")
    check_packages(c("Seurat","harmony"))
    
    input <- p(DIR$cache_input, "03_GSE178318_qc_singlets.rds")
    assert_file(input)
    out <- p(DIR$cache_annotation, "04_GSE178318_global_patient_harmony.rds")
    out_sens <- p(DIR$cache_annotation, "05_GSE178318_global_sensitivity.rds")
    if (!force && cache_has_current_version(out) &&
        cache_has_current_version(out_sens)) {
      message("[03_global_atlas_and_sensitivities] CACHE HIT：直接返回全局敏感性对象。")
      return(load_step_cache(out_sens))
    }
    
    obj <- load_step_cache(input)
    if (is.null(obj)) stop("Stage 3 输入缓存缺失或版本不匹配：", input)
    stopifnot(ncol(obj) > 0,
              length(unique(obj$patient_id)) == 6L,
              length(unique(obj$sample_id)) == 15L)
    SeuratObject::DefaultAssay(obj) <- "RNA"
    
    seed <- PARAM$seed + 100L
    set.seed(seed)
    cache_file <- p(DIR$cache_annotation,
                    "stage3_01_normalized.rds")
    
    if(cache_has_current_version(cache_file) && !force){
      message("[03_global_atlas_and_sensitivities] STEP 1/6 CACHE HIT | NormalizeData")
      obj <- load_step_cache(cache_file)
      
    }else{
      
      obj <- stage_step(
        "03_global_atlas_and_sensitivities", 1L, 6L,
        "NormalizeData：全局 LogNormalize",
        {
          Seurat::NormalizeData(
            obj,
            normalization.method="LogNormalize",
            scale.factor=10000,
            verbose=FALSE
          )
        }, "3–15 min"
      )
      
      save_step_cache(obj, cache_file)
      
    }
    cache_file <- p(DIR$cache_annotation,
                    "stage3_02_hvg.rds")
    
    if(cache_has_current_version(cache_file) && !force){
      message("[03_global_atlas_and_sensitivities] STEP 2/6 CACHE HIT | FindVariableFeatures")
      obj <- load_step_cache(cache_file)
      
    }else{
      
      obj <- stage_step(
        "03_global_atlas_and_sensitivities", 2L, 6L,
        "FindVariableFeatures：选择全局 HVG",
        {
          Seurat::FindVariableFeatures(
            obj,
            selection.method="vst",
            nfeatures=PARAM$global_hvg,
            verbose=FALSE
          )
        }, "3–15 min"
      )
      
      save_step_cache(obj, cache_file)
      
    }
    cache_file <- p(DIR$cache_annotation,
                    "stage3_03_pca.rds")
    
    if(cache_has_current_version(cache_file) && !force){
      message("[03_global_atlas_and_sensitivities] STEP 3/6 CACHE HIT | ScaleData + PCA")
      obj <- load_step_cache(cache_file)
      
    }else{
      
      obj <- stage_step(
        "03_global_atlas_and_sensitivities", 3L, 6L,
        "ScaleData + RunPCA：建立全局 PCA 表示",
        {
          ans <- Seurat::ScaleData(
            obj,
            features=Seurat::VariableFeatures(obj),
            verbose=FALSE
          )
          Seurat::RunPCA(
            ans,
            npcs=PARAM$global_pcs,
            verbose=FALSE
          )
        }, "5–25 min"
      )
      
      
      save_step_cache(obj, cache_file)
      
    }
    
    # Primary representation: Harmony by patient
    cache_file <- p(DIR$cache_annotation,
                    "stage3_04_harmony_patient.rds")
    
    if(cache_has_current_version(cache_file) && !force){
      message("[03_global_atlas_and_sensitivities] STEP 4/6 CACHE HIT | Harmony(patient)")
      obj <- load_step_cache(cache_file)
      
    }else{
      
      obj <- stage_step(
        "03_global_atlas_and_sensitivities", 4L, 6L,
        "RunHarmony：按 patient_id 校正批次",
        {
          future::plan("sequential")
          harmony::RunHarmony(
            obj,
            group.by.vars="patient_id",
            reduction="pca",
            dims.use=1:PARAM$global_pcs,
            reduction.save="harmony.patient",
            verbose=TRUE
          )
        }, "5–30 min"
      )
      
      save_step_cache(obj, cache_file)
      
    }
    cache_file <- p(DIR$cache_annotation,
                    "stage3_05_umap_neighbors.rds")
    
    if(cache_has_current_version(cache_file) && !force){
      message("[03_global_atlas_and_sensitivities] STEP 5/6 CACHE HIT | 邻居图 + UMAP")
      obj <- load_step_cache(cache_file)
      
    }else{
      
      obj <- stage_step(
        "03_global_atlas_and_sensitivities", 5L, 6L,
        "FindNeighbors + RunUMAP：建立全局图与嵌入",
        {
          ans <- Seurat::FindNeighbors(
            obj, reduction = "harmony.patient",
            dims = PARAM$global_neighbor_dims,
            graph.name = c("patient_nn","patient_snn"), verbose = FALSE
          )
          Seurat::RunUMAP(
            ans, reduction = "harmony.patient",
            dims = PARAM$global_neighbor_dims,
            reduction.name = "umap.patient", reduction.key = "UMAPP_",
            seed.use = seed, verbose = FALSE
          )
        }, "5–30 min"
      )
      
      save_step_cache(obj, cache_file)
      
    }
    
    resolutions <- c(0.2,0.4,0.6,0.8,1.0)
    for (res in resolutions) {
      obj <- Seurat::FindClusters(
        obj, graph.name = "patient_snn", resolution = res,
        algorithm = 1, random.seed = seed,
        cluster.name = paste0("global_patient_res.", res), verbose = FALSE
      )
    }
    obj$global_cluster_primary <- obj[["global_patient_res.0.6"]][, 1]
    SeuratObject::Idents(obj) <- factor(obj$global_cluster_primary)
    
    write.csv(
      data.frame(
        analysis = c("global_primary", "global_sample_sensitivity"),
        normalization = "Seurat LogNormalize scale.factor=10000",
        hvg = PARAM$global_hvg,
        pca_npcs = PARAM$global_pcs,
        harmony_group_by = c("patient_id", "sample_id"),
        harmony_dims = c("1:50", "1:50"),
        graph_dims = c("1:30", "1:30"),
        resolution = PARAM$global_primary_resolution,
        seed = c(PARAM$seed + 100L, PARAM$seed + 200L),
        batch_correction_excludes_biology = TRUE,
        stringsAsFactors = FALSE
      ),
      p(DIR$audit, "03_Harmony_parameter_manifest.csv"), row.names = FALSE
    )
    
    markers <- stage_step(
      "03_global_atlas_and_sensitivities", 6L, 6L,
      "FindAllMarkers：全局 marker 与后续敏感性审计",
      {
        Seurat::FindAllMarkers(
          obj,
          assay="RNA",
          layer="data",
          only.pos=TRUE, min.pct = 0.20,
          logfc.threshold = 0.25, test.use = "wilcox",
          max.cells.per.ident = 5000, random.seed = seed, verbose = FALSE
        )
      }, "5–30 min"
    )
    write.csv(markers, p(DIR$tables, "06_global_primary_positive_markers.csv"),
              row.names = FALSE)
    
    # Primary patient-Harmony representation plus no-batch PCA and
    # sample-Harmony sensitivity representations retained for audit.
    set.seed(PARAM$seed + 200L)
    obj <- Seurat::FindNeighbors(
      obj, reduction = "pca", dims = 1:30,
      graph.name = c("pca_nn","pca_snn"), verbose = FALSE
    )
    obj <- Seurat::FindClusters(
      obj, graph.name = "pca_snn", resolution = 0.6, algorithm = 1,
      random.seed = PARAM$seed + 200L,
      cluster.name = "global_pca_res.0.6", verbose = FALSE
    )
    obj <- Seurat::RunUMAP(
      obj, reduction = "pca", dims = 1:30,
      reduction.name = "umap.pca", reduction.key = "UMAPN_",
      seed.use = PARAM$seed + 200L, verbose = FALSE
    )
    
    obj <- harmony::RunHarmony(
      obj, group.by.vars = "sample_id", reduction.use = "pca",
      dims.use = 1:50, reduction.save = "harmony.sample",
      plot_convergence = FALSE, verbose = TRUE
    )
    obj <- Seurat::FindNeighbors(
      obj, reduction = "harmony.sample", dims = 1:30,
      graph.name = c("sample_nn","sample_snn"), verbose = FALSE
    )
    obj <- Seurat::FindClusters(
      obj, graph.name = "sample_snn", resolution = 0.6, algorithm = 1,
      random.seed = PARAM$seed + 200L,
      cluster.name = "global_sample_res.0.6", verbose = FALSE
    )
    obj <- Seurat::RunUMAP(
      obj, reduction = "harmony.sample", dims = 1:30,
      reduction.name = "umap.sample", reduction.key = "UMAPS_",
      seed.use = PARAM$seed + 200L, verbose = FALSE
    )
    
    d <- obj[[]]
    for (nm in c("patient_vs_pca","patient_vs_sample","pca_vs_sample")) {
      tab <- switch(
        nm,
        patient_vs_pca = table(d$global_patient_res.0.6, d$global_pca_res.0.6),
        patient_vs_sample = table(d$global_patient_res.0.6, d$global_sample_res.0.6),
        pca_vs_sample = table(d$global_pca_res.0.6, d$global_sample_res.0.6)
      )
      write.csv(as.data.frame.matrix(tab),
                p(DIR$tables, paste0("07_cluster_crosstab_", nm, ".csv")))
    }
    
    # Canonical marker audit
    panels <- list(
      immune = c("PTPRC"),
      epithelial = c("EPCAM","KRT8","KRT18","KRT19","KRT20"),
      t_cell = c("CD3D","CD3E","TRAC","IL7R","CD8A"),
      nk = c("NKG7","GNLY","KLRD1","NCR1"),
      b_cell = c("CD79A","MS4A1","CD37","CD74","CD19"),
      plasma = c("JCHAIN","MZB1","SDC1","IGHG1"),
      myeloid = c("LST1","TYROBP","FCER1G","LILRB1","CTSS"),
      mast = c("TPSAB1","TPSB2","KIT","MS4A2"),
      fibroblast = c("COL1A1","COL1A2","COL3A1","DCN","LUM"),
      endothelial = c("PECAM1","VWF","EMCN","KDR","ESAM"),
      erythroid = c("HBA1","HBA2","HBB","ALAS2"),
      cycling = c("MKI67","TOP2A","STMN1","TUBB")
    )
    genes <- unique(unlist(panels))
    present <- intersect(genes, rownames(obj))
    write.csv(
      data.frame(gene = genes, present = genes %in% present),
      p(DIR$tables, "08_canonical_marker_presence.csv"), row.names = FALSE
    )
    cl <- as.character(obj$global_patient_res.0.6)
    clv <- sort(unique(cl))
    expr <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")[
      present, , drop = FALSE
    ]
    avg <- do.call(cbind, lapply(clv, function(z) {
      Matrix::rowMeans(expr[, cl == z, drop = FALSE])
    }))
    colnames(avg) <- clv
    avg_long <- as.data.frame(as.table(as.matrix(avg)), stringsAsFactors = FALSE)
    names(avg_long) <- c("gene","cluster","average_expression")
    avg_long$panel <- vapply(avg_long$gene, function(g) {
      hit <- which(vapply(panels, function(x) g %in% x, logical(1)))
      if (length(hit)) names(panels)[hit[1]] else NA_character_
    }, character(1))
    write.csv(avg_long, p(DIR$tables, "09_canonical_marker_average_expression.csv"),
              row.names = FALSE)
    
    pct <- do.call(cbind, lapply(clv, function(z) {
      Matrix::rowMeans(expr[, cl == z, drop = FALSE] > 0)
    }))
    colnames(pct) <- clv
    pct_long <- as.data.frame(as.table(as.matrix(pct)), stringsAsFactors = FALSE)
    names(pct_long) <- c("gene","cluster","fraction_detected")
    write.csv(pct_long, p(DIR$tables, "10_canonical_marker_detection.csv"),
              row.names = FALSE)
    
    # PC dimension sensitivity: same sample-Harmony representation, graph dims 20/30/40/50
    adjusted_rand <- function(a, b) {
      tab <- table(a,b)
      choose2 <- function(x) sum(x*(x-1)/2)
      n <- sum(tab)
      idx <- choose2(tab)
      rp <- choose2(rowSums(tab)); cp <- choose2(colSums(tab))
      tp <- choose2(n)
      expected <- rp*cp/tp
      maximum <- 0.5*(rp+cp)
      if (maximum == expected) return(1)
      (idx-expected)/(maximum-expected)
    }
    dims_set <- c(20L,30L,40L,50L)
    cluster_cols <- character()
    for (kk in dims_set) {
      nn <- paste0("pc",kk,"_nn")
      snn <- paste0("pc",kk,"_snn")
      ccol <- paste0("cluster_pc",kk,"_res0.6")
      obj <- Seurat::FindNeighbors(
        obj, reduction = "harmony.sample", dims = seq_len(kk),
        graph.name = c(nn,snn), verbose = FALSE
      )
      obj <- Seurat::FindClusters(
        obj, graph.name = snn, resolution = 0.6, algorithm = 1,
        random.seed = PARAM$seed + 600L, cluster.name = ccol, verbose = FALSE
      )
      cluster_cols <- c(cluster_cols, ccol)
    }
    md <- obj[[]]
    pc_summary <- do.call(rbind, lapply(cluster_cols, function(cc) {
      z <- as.character(md[[cc]])
      data.frame(
        partition = cc, n_clusters = length(unique(z)),
        min_cluster_size = min(table(z)),
        median_cluster_size = median(table(z)),
        max_cluster_size = max(table(z))
      )
    }))
    write.csv(pc_summary, p(DIR$tables, "11_pc_dimension_partition_summary.csv"),
              row.names = FALSE)
    ari <- do.call(rbind, combn(cluster_cols, 2, simplify = FALSE, FUN = function(pair) {
      data.frame(
        partition_a = pair[1], partition_b = pair[2],
        adjusted_rand_index = adjusted_rand(md[[pair[1]]], md[[pair[2]]])
      )
    }))
    write.csv(ari, p(DIR$tables, "12_pc_dimension_adjusted_rand_index.csv"),
              row.names = FALSE)
    
    # Harmony preservation audit: marker-module scores are computed from the
    # same normalized RNA assay and summarized under each representation's
    # partition; batch correction is not allowed to redefine expression scores.
    expr_pres <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
    pres_scores <- do.call(cbind, lapply(BROAD_LINEAGE_MODULES, function(gs)
      module_score(expr_pres, gs)))
    colnames(pres_scores) <- paste0("module_", names(BROAD_LINEAGE_MODULES))
    rownames(pres_scores) <- colnames(obj)
    pres_md <- obj[[]]
    pres_md$barcode <- rownames(pres_md)
    pres_md <- cbind(pres_md, as.data.frame(pres_scores[rownames(pres_md), , drop = FALSE]))
    partition_cols <- c(
      "global_patient_res.0.6", "global_pca_res.0.6", "global_sample_res.0.6"
    )
    preservation_rows <- list(); pp <- 0L
    for (pc in partition_cols) {
      if (!pc %in% names(pres_md)) next
      for (mod in colnames(pres_scores)) {
        z <- aggregate(pres_md[[mod]], list(cluster = pres_md[[pc]]), mean, na.rm = TRUE)
        names(z)[2] <- "mean_module_score"
        z$partition <- pc; z$module <- mod
        pp <- pp + 1L; preservation_rows[[pp]] <- z
      }
    }
    write.csv(
      do.call(rbind, preservation_rows),
      p(DIR$tables, "14_harmony_marker_module_preservation.csv"),
      row.names = FALSE
    )
    
    # composition / metadata outputs
    comp <- do.call(rbind, lapply(c("patient_id","sample_id","tissue"), function(v) {
      x <- as.data.frame(
        table(cluster = obj$global_cluster_primary, group = obj[[v]][,1]),
        stringsAsFactors = FALSE
      )
      names(x)[3] <- "n_cells"; x$group_variable <- v; x
    }))
    write.csv(comp, p(DIR$tables, "13_global_primary_cluster_composition.csv"),
              row.names = FALSE)
    
    meta_out <- obj[[]]
    meta_out$barcode <- rownames(meta_out)
    keep <- c(
      "barcode","sample_id","patient_id","tissue","global_cluster_primary",
      paste0("global_patient_res.",resolutions),
      "global_pca_res.0.6","global_sample_res.0.6",cluster_cols
    )
    con <- gzfile(p(DIR$tables, "14_global_cluster_and_sensitivity_metadata.csv.gz"), "wt")
    write.csv(meta_out[, intersect(keep, names(meta_out)), drop = FALSE],
              con, row.names = FALSE)
    close(con)
    
    save_step_cache(obj, out)
    save_step_cache(obj, out_sens)
    future::plan(
      "multisession",
      workers = 4
    )
    write_session_info("03_global_atlas_and_sensitivities")
    cat("Global cells:", ncol(obj), "\n")
    cat("Primary clusters:", length(unique(obj$global_cluster_primary)), "\n")
    obj
  })
}

############################
# Stage 4：marker-module broad-lineage + Myeloid/DC reclustering + fine-state
############################

stage4_myeloid_annotation <- function(force = FALSE) {
  stage_log("04_all_lineage_subclustering_and_fine_maps", {
    check_packages(c("Seurat","harmony"))
    # Seurat may dispatch NormalizeData/ScaleData through future_lapply.
    # The lineage objects can be several GiB; exporting them to multisession
    # workers is unsafe on typical desktop memory and caused worker death.
    # Keep Stage 4 single-process and restore the caller's plan on exit.
    old_future_plan <- future::plan()
    future::plan(future::sequential)
    on.exit(future::plan(old_future_plan), add = TRUE)
    message("[04_all_lineage_subclustering_and_fine_maps] future plan: sequential (large lineage objects)")
    
    input_primary <- p(DIR$cache_annotation, "05_GSE178318_global_sensitivity.rds")
    input_fallback <- p(DIR$cache_annotation, "04_GSE178318_global_patient_harmony.rds")
    input <- if (cache_has_current_version(input_primary)) {
      input_primary
    } else if (cache_has_current_version(input_fallback)) {
      input_fallback
    } else {
      input_primary
    }
    assert_file(input)
    if (!cache_has_current_version(input)) {
      stop("Stage 4 输入缓存不是当前 pipeline 版本；请先重新运行 Stage 3：", input)
    }
    out_my <- p(DIR$cache_annotation, "fine_myeloid_dc_annotated.rds")
    out_t <- p(DIR$cache_annotation, "fine_t_nk_annotated.rds")
    # A previous run may have written lineage objects before failing during
    # final audit-table assembly.  Treat those objects as reusable only when
    # the Stage-4 reviewer-facing audit outputs also exist; otherwise rebuild
    # the summary layer from the current lineage caches.
    stage4_audit_ready <- file.exists(p(DIR$audit, "04_annotation_audit.csv")) &&
      file.exists(p(DIR$tables, "fine_state_annotation_confidence.csv"))
    if (!force && stage4_audit_ready && cache_has_current_version(out_my) &&
        cache_has_current_version(out_t)) {
      message("[04_all_lineage_subclustering_and_fine_maps] CACHE HIT：直接返回髓系/DC 对象。")
      return(load_step_cache(out_my))
    }
    
    obj <- load_step_cache(input)
    if (is.null(obj)) stop("Stage 4 输入缓存缺失或版本不匹配：", input)
    cl <- as.character(obj$global_cluster_primary)
    lineage_step_total <- 1L + 3L * length(LINEAGES_TO_SUBCLUSTER) + 2L
    annotation_audit_rows <- list()
    
    # Primary broad-lineage assignment is recomputed from a prespecified
    # marker-module score on the new rerun.  Legacy cluster labels are only
    # joined below for concordance auditing and never determine inclusion.
    broad <- stage_step(
      "04_all_lineage_subclustering_and_fine_maps", 1L, lineage_step_total,
      "marker-module broad-lineage 评分与初筛",
      {
        cluster_module_scores(obj, "global_cluster_primary", BROAD_LINEAGE_MODULES)
      }, "2–15 min"
    )
    broad_scores <- broad$scores

    # Stress is a state-like overlay rather than a broad lineage.  The old
    # implementation allowed JUN/FOS/HSP modules to compete directly with
    # lineage modules, which incorrectly converted clear T/NK, fibroblast and
    # endothelial clusters into excluded Stress clusters.  Cycling and
    # erythroid remain exclusion modules, but are evaluated separately from
    # the primary lineage winner.
    primary_lineage_names <- intersect(
      LINEAGES_TO_SUBCLUSTER, colnames(broad_scores)
    )
    artifact_names <- intersect(c("Erythroid", "Cycling", "Stress"),
                                colnames(broad_scores))
    if (!length(primary_lineage_names)) {
      stop("Broad lineage module scores缺少主谱系模块。")
    }
    primary_scores <- broad_scores[, primary_lineage_names, drop = FALSE]
    artifact_scores <- broad_scores[, artifact_names, drop = FALSE]
    primary_winner <- apply(primary_scores, 1, function(v) {
      if (!any(is.finite(v)) || max(v, na.rm = TRUE) <= 0) return("Unknown_mixed")
      oo <- order(v, decreasing = TRUE)
      if (length(oo) > 1L &&
          (v[oo[1]] - v[oo[2]]) < PARAM$annotation_module_margin) {
        return("Unknown_mixed")
      }
      names(v)[oo[1]]
    })
    primary_winner_score <- vapply(seq_len(nrow(primary_scores)), function(ii) {
      v <- primary_scores[ii, ]
      if (is.na(primary_winner[ii]) || primary_winner[ii] == "Unknown_mixed")
        return(NA_real_)
      unname(v[primary_winner[ii]])
    }, numeric(1))
    artifact_winner <- if (ncol(artifact_scores)) {
      apply(artifact_scores, 1, function(v) {
        if (!any(is.finite(v))) return(NA_character_)
        names(v)[which.max(v)]
      })
    } else rep(NA_character_, nrow(broad_scores))
    artifact_score <- if (ncol(artifact_scores)) {
      apply(artifact_scores, 1, max, na.rm = TRUE)
    } else rep(NA_real_, nrow(broad_scores))

    # A strong cycling signal is treated as an exclusion only when it is
    # close to the primary lineage score; stress alone is retained as an
    # overlay so that stressed but biologically identifiable lineages remain.
    cycling_score <- if ("Cycling" %in% colnames(broad_scores))
      broad_scores[, "Cycling"] else rep(0, nrow(broad_scores))
    erythroid_score <- if ("Erythroid" %in% colnames(broad_scores))
      broad_scores[, "Erythroid"] else rep(0, nrow(broad_scores))
    cycling_dominant <- is.finite(primary_winner_score) &
      cycling_score >= 1 & cycling_score >= 0.8 * primary_winner_score
    erythroid_dominant <- is.finite(primary_winner_score) &
      erythroid_score >= 1 & erythroid_score >= 0.8 * primary_winner_score

    broad_winner <- primary_winner
    broad_winner[cycling_dominant] <- "Cycling"
    broad_winner[erythroid_dominant] <- "Erythroid"
    broad_winner[!is.finite(primary_winner_score)] <- "Unknown_mixed"
    excluded_broad <- c("Unknown_mixed", "Cycling", "Erythroid")
    broad_status <- ifelse(
      broad_winner %in% excluded_broad,
      "excluded_from_lineage_primary", "retained"
    )
    broad_audit <- annotation_score_audit(
      primary_scores, broad_winner, broad_status,
      level = "broad", lineage = "global"
    )
    broad_audit$primary_winner_score <- primary_winner_score
    broad_audit$artifact_winner <- artifact_winner
    broad_audit$artifact_score <- artifact_score
    broad_audit$stress_overlay <- artifact_winner == "Stress"
    broad_audit$n_cells <- as.integer(table(cl)[broad_audit$cluster])
    annotation_audit_rows[[1L]] <- broad_audit
    map_new <- data.frame(
      cluster = names(broad_winner),
      broad_lineage = unname(broad_winner),
      annotation_status = unname(broad_status),
      winner_module = broad_audit$winner,
      runner_up_module = broad_audit$runner_up,
      top_module_score = broad_audit$winner_score,
      second_module_score = broad_audit$runner_up_score,
      artifact_winner = broad_audit$artifact_winner,
      artifact_score = broad_audit$artifact_score,
      stress_overlay = broad_audit$stress_overlay,
      annotation_margin = broad_audit$margin,
      annotation_margin_definition = "winner_minus_runner_up_on_cluster_mean_log1p_normalized_module_score",
      stringsAsFactors = FALSE
    )
    write.csv(map_new, p(DIR$tables, "04_new_broad_lineage_module_map.csv"), row.names = FALSE)
    write.csv(
      data.frame(
        cluster = rownames(broad_scores), broad_scores,
        stringsAsFactors = FALSE, check.names = FALSE
      ),
      p(DIR$tables, "04_new_broad_lineage_module_scores.csv"), row.names = FALSE
    )
    
    legacy_idx <- match(map_new$cluster, LEGACY_GLOBAL_MAP$cluster)
    concordance <- data.frame(
      cluster = map_new$cluster,
      n_cells = as.integer(table(cl)[map_new$cluster]),
      new_broad_lineage = map_new$broad_lineage,
      new_annotation_status = map_new$annotation_status,
      legacy_broad_lineage = LEGACY_GLOBAL_MAP$broad_lineage[legacy_idx],
      legacy_annotation_status = LEGACY_GLOBAL_MAP$annotation_status[legacy_idx],
      stringsAsFactors = FALSE
    )
    concordance$concordant <- concordance$new_broad_lineage ==
      concordance$legacy_broad_lineage
    write.csv(concordance, p(DIR$audit, "04_annotation_concordance_legacy.csv"),
              row.names = FALSE)
    
    obj$global_broad_lineage <- map_new$broad_lineage[match(cl, map_new$cluster)]
    obj$global_lineage_status <- map_new$annotation_status[match(cl, map_new$cluster)]
    expr_audit <- Seurat::GetAssayData(obj, assay = "RNA", layer = "data")
    obj$myeloid_marker_score <- module_score(expr_audit, MYELOID_MODULE)
    sel_audit <- concordance
    sel_audit$myeloid_mean_score <- vapply(sel_audit$cluster, function(cc) {
      mean(obj$myeloid_marker_score[cl == cc], na.rm = TRUE)
    }, numeric(1))
    sel_audit$myeloid_detection_fraction <- vapply(sel_audit$cluster, function(cc) {
      mean(obj$myeloid_marker_score[cl == cc] > 0, na.rm = TRUE)
    }, numeric(1))
    write.csv(sel_audit, p(DIR$audit, "04_myeloid_selection_audit.csv"), row.names = FALSE)
    
    write.csv(
      data.frame(
        analysis = "Myeloid_DC_lineage_reclustering",
        normalization = "Seurat LogNormalize scale.factor=10000",
        hvg = PARAM$lineage_hvg,
        pca_npcs = PARAM$lineage_pcs,
        harmony_group_by = "sample_id",
        harmony_dims = "1:30",
        graph_dims = "1:20",
        resolutions_tested = "0.4;0.6;1.0",
        primary_resolution = PARAM$lineage_primary_resolution,
        seed = PARAM$seed + 800L,
        batch_correction_excludes_biology = TRUE,
        stringsAsFactors = FALSE
      ),
      p(DIR$audit, "04_lineage_Harmony_parameter_manifest.csv"), row.names = FALSE
    )
    
    # Explicit manuscript-count audit: singlets -> retained analysis cells.
    n_singlet <- ncol(obj)
    n_excluded_primary <- sum(
      obj$global_lineage_status == "excluded_from_lineage_primary",
      na.rm = TRUE
    )
    n_analysis_cells <- n_singlet - n_excluded_primary
    count_audit <- data.frame(
      n_singlet = n_singlet,
      n_excluded_from_lineage_primary = n_excluded_primary,
      n_analysis_cells = n_analysis_cells,
      expected_singlet = 108353L,
      expected_analysis_cells = NA_integer_,
      pass_singlet_count = n_singlet == 108353L,
      note = "Analysis-cell count is descriptive after new marker-based annotation; it is not constrained to the legacy 106,115 value."
    )
    write.csv(count_audit, p(DIR$audit, "04_global_analysis_cell_count_audit.csv"), row.names = FALSE)
    message(
      "New marker-based annotation retained ", n_analysis_cells,
      " analysis cells from ", n_singlet,
      " singlets; legacy 108353/106115 counts are audit references only."
    )
    
    meta0 <- obj[[]]; meta0$barcode <- rownames(meta0)
    write.csv(
      meta0[, c("barcode","sample_id","patient_id","tissue",
                "global_cluster_primary","global_broad_lineage","global_lineage_status")],
      p(DIR$tables, "15_global_broad_lineage_metadata.csv"), row.names = FALSE
    )
    comp0 <- as.data.frame(
      table(lineage = obj$global_broad_lineage,
            tissue = obj$tissue, patient = obj$patient_id),
      stringsAsFactors = FALSE
    )
    write.csv(comp0, p(DIR$tables, "16_global_broad_lineage_composition.csv"),
              row.names = FALSE)
    save_step_cache(obj, p(DIR$cache_annotation, "06_GSE178318_global_broad_lineage.rds"))
    
    all_summary <- list()
    for (lin in LINEAGES_TO_SUBCLUSTER) {
      cells <- colnames(obj)[
        obj$global_broad_lineage == lin &
          obj$global_lineage_status == "retained"
      ]
      cat("Lineage", lin, "cells", length(cells), "\n")
      if (length(cells) < 100L) {
        skipped_from <- 1L + (match(lin, LINEAGES_TO_SUBCLUSTER) - 1L) * 3L + 1L
        skipped_to <- skipped_from + 2L
        cat("[04_all_lineage_subclustering_and_fine_maps] STEPS ",
            skipped_from, "-", skipped_to, "/", lineage_step_total,
            " SKIPPED | ", lin, " lineage 仅 ", length(cells),
            " 个细胞（少于 100） | current time: ",
            format(Sys.time(), tz = "Asia/Shanghai"), "\n", sep = "")
        next
      }
      
      lineage_cache <- p(
        DIR$cache_annotation,
        paste0("lineage_", tolower(lin), "_subclustered.rds")
      )
      if (!force && cache_has_current_version(lineage_cache)) {
        cached_sub <- load_step_cache(lineage_cache)
        if (!is.null(cached_sub) &&
            setequal(colnames(cached_sub), cells)) {
          cache_step_from <- 1L + (match(lin, LINEAGES_TO_SUBCLUSTER) - 1L) * 3L + 1L
          cache_step_to <- cache_step_from + 2L
          cat(
            "[04_all_lineage_subclustering_and_fine_maps] STEPS ",
            cache_step_from, "-", cache_step_to, "/", lineage_step_total,
            " CACHE HIT | ", lin, " lineage subcluster | current time: ",
            format(Sys.time(), tz = "Asia/Shanghai"), "\n", sep = ""
          )
          all_summary[[lin]] <- data.frame(
            lineage = lin, n_cells = ncol(cached_sub),
            n_clusters = length(unique(cached_sub$lineage_cluster_res0.6))
          )
          next
        }
        message(
          "[04_all_lineage_subclustering_and_fine_maps] lineage cache cells "
          , "与当前 broad-lineage 结果不一致，重新计算：", lin
        )
      }
      
      sub <- subset(obj, cells = cells)
      SeuratObject::DefaultAssay(sub) <- "RNA"
      seed <- PARAM$seed + 800L
      set.seed(seed)
      sub <- stage_step(
        "04_all_lineage_subclustering_and_fine_maps",
        1L + (match(lin, LINEAGES_TO_SUBCLUSTER) - 1L) * 3L + 1L,
        lineage_step_total,
        paste0("NormalizeData + HVG + ScaleData：", lin, " lineage"),
        {
          ans <- Seurat::NormalizeData(
            sub, normalization.method = "LogNormalize",
            scale.factor = 10000, verbose = FALSE
          )
          ans <- Seurat::FindVariableFeatures(
            ans, selection.method = "vst", nfeatures = PARAM$lineage_hvg,
            verbose = FALSE
          )
          Seurat::ScaleData(
            ans, features = Seurat::VariableFeatures(ans),
            vars.to.regress = "percent.mt", verbose = FALSE
          )
        }, "2–15 min"
      )
      sub <- stage_step(
        "04_all_lineage_subclustering_and_fine_maps",
        1L + (match(lin, LINEAGES_TO_SUBCLUSTER) - 1L) * 3L + 2L,
        lineage_step_total,
        paste0("RunPCA：", lin, " lineage"),
        {
          Seurat::RunPCA(
            sub, features = Seurat::VariableFeatures(sub), npcs = PARAM$lineage_pcs,
            seed.use = seed, verbose = FALSE
          )
        }, "3–20 min"
      )
      sub <- stage_step(
        "04_all_lineage_subclustering_and_fine_maps",
        1L + (match(lin, LINEAGES_TO_SUBCLUSTER) - 1L) * 3L + 3L,
        lineage_step_total,
        paste0("RunHarmony：", lin, " lineage 批次校正"),
        {
          harmony::RunHarmony(
            sub, group.by.vars = "sample_id", reduction.use = "pca",
            dims.use = PARAM$lineage_neighbor_dims,
            reduction.save = "harmony.lineage",
            plot_convergence = FALSE, verbose = FALSE
          )
        }, "3–20 min"
      )
      graph <- paste0(tolower(lin), "_nn")
      snn <- paste0(tolower(lin), "_snn")
      sub <- Seurat::FindNeighbors(
        sub, reduction = "harmony.lineage", dims = PARAM$lineage_neighbor_dims,
        graph.name = c(graph,snn), verbose = FALSE
      )
      for (res in c(0.4,0.6,1.0)) {
        nm <- paste0("lineage_cluster_res", format(res, nsmall = 1))
        sub <- Seurat::FindClusters(
          sub, graph.name = snn, resolution = res,
          random.seed = seed, cluster.name = nm, verbose = FALSE
        )
      }
      SeuratObject::Idents(sub) <- "lineage_cluster_res0.6"
      markers <- Seurat::FindAllMarkers(
        sub,
        assay = "RNA",
        layer = "data",
        only.pos = TRUE,min.pct = 0.2,
        logfc.threshold = 0.25, test.use = "wilcox",
        max.cells.per.ident = 3000, random.seed = seed, verbose = FALSE
      )
      markers$lineage <- lin
      
      sub_meta <- sub[[]]; sub_meta$barcode <- rownames(sub_meta); sub_meta$lineage <- lin
      write.csv(
        sub_meta,
        p(DIR$tables, paste0("lineage_",tolower(lin),"_metadata.csv")),
        row.names = FALSE
      )
      write.csv(
        markers,
        p(DIR$tables, paste0("lineage_",tolower(lin),"_markers.csv")),
        row.names = FALSE
      )
      comp <- as.data.frame(
        table(cluster = sub$lineage_cluster_res0.6,
              patient = sub$patient_id, tissue = sub$tissue),
        stringsAsFactors = FALSE
      )
      comp$lineage <- lin
      write.csv(
        comp,
        p(DIR$tables, paste0("lineage_",tolower(lin),"_composition.csv")),
        row.names = FALSE
      )
      save_step_cache(
        sub,
        lineage_cache
      )
      all_summary[[lin]] <- data.frame(
        lineage = lin, n_cells = ncol(sub),
        n_clusters = length(unique(sub$lineage_cluster_res0.6))
      )
      rm(sub); invisible(gc())
    }
    summary_df <- do.call(rbind, all_summary)
    write.csv(summary_df, p(DIR$tables, "17_lineage_subcluster_summary.csv"),
              row.names = FALSE)
    
    # Fine-state assignment is also recomputed from the new rerun's marker
    # modules.  Legacy maps are joined only as a concordance reference.
    apply_new_fine_modules <- function(lin, modules, legacy_map) {
      f <- p(DIR$cache_annotation, paste0("lineage_",tolower(lin),"_subclustered.rds"))
      assert_file(f)
      sub <- load_step_cache(f)
      if (is.null(sub)) stop("Stage 4 lineage 缓存缺失或版本不匹配：", f)
      cc <- as.character(sub$lineage_cluster_res0.6)
      sm <- cluster_module_scores(sub, "lineage_cluster_res0.6", modules)
      scores <- sm$scores
      # Stress, cycling and mixed-lineage modules are quality/state overlays,
      # not biological identities.  Excluding them from the primary winner
      # prevents a heat-shock or cycling signal from masking a clear C1Q,
      # SPP1/TREM2, monocyte or DC state.  They are still retained in the
      # score table and used to flag clusters with a strong artifact signal.
      artifact_names <- intersect(
        c("Myeloid_stress", "T_stress", "Myeloid_cycling",
          "Mixed_T_NK_myeloid", "Mixed_TNK_contamination",
          "Mixed_low_confidence"),
        colnames(scores)
      )
      biological_names <- setdiff(colnames(scores), artifact_names)
      if (!length(biological_names)) {
        stop("Fine-state module scores缺少生物学状态模块：", lin)
      }
      biological_scores <- scores[, biological_names, drop = FALSE]
      winners <- apply(biological_scores, 1, function(v) {
        if (!any(is.finite(v)) || max(v, na.rm = TRUE) <= 0) return("Unresolved")
        oo <- order(v, decreasing = TRUE)
        if (length(oo) > 1L &&
            (v[oo[1]] - v[oo[2]]) < PARAM$annotation_module_margin) {
          return("Unresolved")
        }
        names(v)[oo[1]]
      })
      biological_winner_score <- vapply(seq_len(nrow(biological_scores)), function(ii) {
        v <- biological_scores[ii, ]
        if (is.na(winners[ii]) || winners[ii] == "Unresolved") return(NA_real_)
        unname(v[winners[ii]])
      }, numeric(1))
      artifact_winner <- if (length(artifact_names)) {
        apply(scores[, artifact_names, drop = FALSE], 1, function(v) {
          if (!any(is.finite(v))) return(NA_character_)
          names(v)[which.max(v)]
        })
      } else rep(NA_character_, nrow(scores))
      artifact_score <- if (length(artifact_names)) {
        apply(scores[, artifact_names, drop = FALSE], 1, max, na.rm = TRUE)
      } else rep(NA_real_, nrow(scores))
      stress_name <- intersect(c("Myeloid_stress", "T_stress"), colnames(scores))
      cycling_name <- intersect("Myeloid_cycling", colnames(scores))
      mixed_names <- intersect(
        c("Mixed_TNK_contamination", "Mixed_low_confidence", "Mixed_T_NK_myeloid"),
        colnames(scores)
      )
      stress_score <- if (length(stress_name))
        apply(scores[, stress_name, drop = FALSE], 1, max, na.rm = TRUE)
      else rep(0, nrow(scores))
      cycling_score <- if (length(cycling_name))
        scores[, cycling_name[[1L]]] else rep(0, nrow(scores))
      mixed_score <- if (length(mixed_names))
        apply(scores[, mixed_names, drop = FALSE], 1, max, na.rm = TRUE)
      else rep(0, nrow(scores))
      biological_margin <- vapply(seq_len(nrow(biological_scores)), function(ii) {
        v <- biological_scores[ii, ]
        vv <- sort(v[is.finite(v)], decreasing = TRUE)
        if (length(vv) < 2L) return(NA_real_)
        vv[[1L]] - vv[[2L]]
      }, numeric(1))
      cycling_dominant <- is.finite(biological_winner_score) &
        cycling_score >= 1 & cycling_score >= 0.8 * biological_winner_score
      mixed_dominant <- is.finite(biological_winner_score) &
        mixed_score >= 1 & mixed_score >= 0.8 * biological_winner_score
      stress_dominant <- is.finite(biological_winner_score) &
        stress_score >= 1.8 & stress_score >= 1.2 * biological_winner_score &
        is.finite(biological_margin) & biological_margin < 0.2
      if ("Myeloid_cycling" %in% names(modules)) winners[cycling_dominant] <- "Myeloid_cycling"
      if (length(mixed_names)) {
        mixed_winner <- artifact_winner
        winners[mixed_dominant & mixed_winner %in% mixed_names] <-
          mixed_winner[mixed_dominant & mixed_winner %in% mixed_names]
      }
      if (length(stress_name)) winners[stress_dominant] <- stress_name[[1L]]
      # Marker-audit fallback labels are data-driven, not cluster-number driven.
      # For each prespecified quality module, identify the highest-scoring
      # cluster and require that it is near the top of the distribution and
      # above the across-cluster median.  This prevents a changed clustering
      # solution from silently inheriting a legacy "cluster 8/11" override.
      marker_audit_module <- rep(NA_character_, nrow(scores))
      marker_audit_rank <- rep(NA_integer_, nrow(scores))
      marker_audit_score <- rep(NA_real_, nrow(scores))
      marker_audit_pass <- rep(FALSE, nrow(scores))
      if (lin == "Myeloid_DC") {
        marker_audit_specs <- c(
          Mixed_TNK_contamination = "marker_audit_TNK_contamination_module",
          Mixed_low_confidence = "marker_audit_low_confidence_module"
        )
        override_clusters <- character(0)
        override_rows <- list()
        for (mod in names(marker_audit_specs)) {
          if (!mod %in% colnames(scores)) {
            stop("Marker audit module缺失，停止Stage 4人工核对守门：", mod)
          }
          v <- as.numeric(scores[, mod])
          finite <- is.finite(v)
          if (!any(finite)) {
            stop("Marker audit module全为非有限值，停止Stage 4人工核对守门：", mod)
          }
          med <- stats::median(v[finite])
          candidate_i <- which.max(ifelse(finite, v, -Inf))
          candidate <- rownames(scores)[candidate_i]
          rk <- rank(-v[finite], ties.method = "min")
          rank_all <- rep(NA_real_, length(v)); rank_all[finite] <- rk
          candidate_rank <- unname(rank_all[candidate_i])
          candidate_score <- v[candidate_i]
          passed <- is.finite(candidate_rank) &&
            candidate_rank <= PARAM$marker_audit_top_rank &&
            is.finite(candidate_score) && candidate_score > med
          override_rows[[mod]] <- data.frame(
            lineage = lin,
            module = mod,
            module_description = unname(marker_audit_specs[[mod]]),
            candidate_cluster = candidate,
            candidate_score = candidate_score,
            candidate_rank = as.integer(candidate_rank),
            across_cluster_median = med,
            top_rank_limit = PARAM$marker_audit_top_rank,
            passed = passed,
            stringsAsFactors = FALSE
          )
          if (!passed) {
            stop(
              "Marker audit guard failed for ", mod,
              ": candidate cluster=", candidate,
              ", rank=", candidate_rank,
              ", score=", signif(candidate_score, 5),
              ", median=", signif(med, 5),
              ". 请重新人工核对marker模块后再继续；未应用任何cluster覆盖。"
            )
          }
          if (candidate %in% override_clusters) {
            stop("两个marker audit模块指向同一cluster，停止以避免标签冲突：", candidate)
          }
          override_clusters <- c(override_clusters, candidate)
          marker_audit_module[candidate_i] <- mod
          marker_audit_rank[candidate_i] <- as.integer(candidate_rank)
          marker_audit_score[candidate_i] <- candidate_score
          marker_audit_pass[candidate_i] <- TRUE
          winners[candidate_i] <- mod
        }
        write.csv(
          do.call(rbind, override_rows),
          p(DIR$audit, "04_fine_marker_override_audit.csv"),
          row.names = FALSE
        )
      }
      assignment_reason <- rep("highest_biological_module_with_margin", length(winners))
      assignment_reason[winners == "Unresolved"] <- "winner_runner_up_margin_below_threshold_or_no_positive_module"
      assignment_reason[winners == "Myeloid_cycling"] <- "cycling_overlay_dominant"
      assignment_reason[winners %in% c("Myeloid_stress", "T_stress")] <- "stress_overlay_dominant"
      assignment_reason[winners == "Mixed_TNK_contamination"] <- "marker_audit_TNK_contamination_module"
      assignment_reason[winners == "Mixed_low_confidence"] <- "marker_audit_low_confidence_module"
      bad <- if (lin == "Myeloid_DC")
        c("Myeloid_stress", "Myeloid_cycling", "Mixed_T_NK_myeloid",
          "Mixed_TNK_contamination", "Mixed_low_confidence", "Unresolved")
      else c("T_stress", "Unresolved")
      fine_status <- ifelse(
        winners %in% bad, "flagged_or_excluded", "retained"
      )
      fine_audit <- annotation_score_audit(
        biological_scores, winners, fine_status,
        level = "fine", lineage = lin
      )
      fine_audit$biological_winner_score <- biological_winner_score
      fine_audit$biological_margin <- biological_margin
      fine_audit$artifact_winner <- artifact_winner
      fine_audit$artifact_score <- artifact_score
      fine_audit$stress_overlay <- artifact_winner %in% c("Myeloid_stress", "T_stress")
      fine_audit$assignment_reason <- assignment_reason
      fine_audit$marker_audit_module <- marker_audit_module
      fine_audit$marker_audit_rank <- marker_audit_rank
      fine_audit$marker_audit_score <- marker_audit_score
      fine_audit$marker_audit_pass <- marker_audit_pass
      fine_audit$n_cells <- as.integer(table(cc)[fine_audit$cluster])
      annotation_audit_rows[[length(annotation_audit_rows) + 1L]] <<- fine_audit
      fmap <- data.frame(
        cluster = names(winners),
        fine_label = unname(winners),
        annotation_status = fine_status,
        winner_module = fine_audit$assigned_label,
        runner_up_module = fine_audit$runner_up,
        top_module_score = fine_audit$winner_score,
        second_module_score = fine_audit$runner_up_score,
        artifact_winner = fine_audit$artifact_winner,
        artifact_score = fine_audit$artifact_score,
        stress_overlay = fine_audit$stress_overlay,
        assignment_reason = assignment_reason,
        marker_audit_module = marker_audit_module,
        marker_audit_rank = marker_audit_rank,
        marker_audit_score = marker_audit_score,
        marker_audit_pass = marker_audit_pass,
        annotation_margin = fine_audit$margin,
        annotation_margin_definition = "winner_minus_runner_up_on_cluster_mean_log1p_normalized_module_score",
        stringsAsFactors = FALSE
      )
      write.csv(fmap, p(DIR$tables, paste0("04_new_",tolower(lin),"_fine_state_map.csv")), row.names = FALSE)
      write.csv(
        data.frame(cluster = rownames(scores), scores, check.names = FALSE),
        p(DIR$tables, paste0("04_new_",tolower(lin),"_fine_state_module_scores.csv")),
        row.names = FALSE
      )
      legacy_idx <- match(fmap$cluster, legacy_map$cluster)
      write.csv(
        data.frame(
          cluster = fmap$cluster,
          new_fine_label = fmap$fine_label,
          new_annotation_status = fmap$annotation_status,
          legacy_fine_label = legacy_map$fine_label[legacy_idx],
          legacy_annotation_status = legacy_map$annotation_status[legacy_idx],
          stringsAsFactors = FALSE
        ),
        p(DIR$audit, paste0("04_",tolower(lin),"_fine_annotation_concordance_legacy.csv")),
        row.names = FALSE
      )
      j <- match(cc, fmap$cluster)
      sub$fine_label <- fmap$fine_label[j]
      sub$fine_status <- fmap$annotation_status[j]
      md <- sub[[]]; md$barcode <- rownames(md)
      write.csv(
        md[, c("barcode","sample_id","patient_id","tissue",
               "lineage_cluster_res0.6","fine_label","fine_status")],
        p(DIR$tables, paste0("fine_",tolower(lin),"_metadata.csv")),
        row.names = FALSE
      )
      ctab <- as.data.frame(
        table(label = sub$fine_label, patient = sub$patient_id, tissue = sub$tissue),
        stringsAsFactors = FALSE
      )
      sm <- unique(fmap[,c("fine_label","annotation_status")])
      ctab$status <- sm$annotation_status[match(ctab$label, sm$fine_label)]
      write.csv(
        ctab,
        p(DIR$tables, paste0("fine_",tolower(lin),"_composition.csv")),
        row.names = FALSE
      )
      save_step_cache(
        sub,
        p(DIR$cache_annotation, paste0("fine_",tolower(lin),"_annotated.rds"))
      )
      sub
    }
    
    tnk <- stage_step(
      "04_all_lineage_subclustering_and_fine_maps", lineage_step_total - 1L,
      lineage_step_total,
      "T/NK fine-state module annotation",
      { apply_new_fine_modules("T_NK", FINE_TNK_MODULES, LEGACY_TNK_MAP) },
      "2–15 min"
    )
    my <- stage_step(
      "04_all_lineage_subclustering_and_fine_maps", lineage_step_total,
      lineage_step_total,
      "Myeloid/DC fine-state module annotation",
      { apply_new_fine_modules("Myeloid_DC", FINE_MYELOID_MODULES, LEGACY_MYELOID_MAP) },
      "2–20 min"
    )
    # Broad-lineage and fine-state audits intentionally contain different
    # reviewer-facing fields.  Bind by the union of columns and fill fields
    # that are not applicable at a given level, rather than failing at the
    # final Stage 4 summary step because the column sets differ.
    audit_cols <- unique(unlist(lapply(annotation_audit_rows, names),
                               use.names = FALSE))
    annotation_audit <- do.call(rbind, lapply(annotation_audit_rows, function(z) {
      missing_cols <- setdiff(audit_cols, names(z))
      for (nm in missing_cols) {
        exemplar <- NULL
        for (zz in annotation_audit_rows) {
          if (nm %in% names(zz)) {
            exemplar <- zz[[nm]]
            break
          }
        }
        if (is.character(exemplar)) {
          z[[nm]] <- rep(NA_character_, nrow(z))
        } else if (is.logical(exemplar)) {
          z[[nm]] <- rep(NA, nrow(z))
        } else {
          z[[nm]] <- rep(NA_real_, nrow(z))
        }
      }
      z[, audit_cols, drop = FALSE]
    }))
    write.csv(
      annotation_audit,
      p(DIR$audit, "04_annotation_audit.csv"),
      row.names = FALSE
    )
    audit_group <- interaction(
      annotation_audit$level, annotation_audit$lineage, drop = TRUE,
      sep = "::"
    )
    annotation_audit_summary <- do.call(rbind, lapply(
      split(seq_len(nrow(annotation_audit)), audit_group),
      function(ii) {
        z <- annotation_audit[ii, , drop = FALSE]
        data.frame(
          level = z$level[1L], lineage = z$lineage[1L],
          n_clusters = nrow(z),
          n_unresolved = sum(z$unresolved, na.rm = TRUE),
          unresolved_fraction = mean(z$unresolved, na.rm = TRUE),
          n_flagged_or_excluded = sum(z$annotation_status != "retained", na.rm = TRUE),
          flagged_or_excluded_fraction = mean(z$annotation_status != "retained", na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    ))
    write.csv(
      annotation_audit_summary,
      p(DIR$audit, "04_annotation_audit_summary.csv"),
      row.names = FALSE
    )
    # Compact, reviewer-facing confidence table for every newly assigned
    # fine-state cluster.  Cells remain the aggregation unit; legacy labels
    # are not used in this table or in the assignment rule.
    fine_confidence <- annotation_audit[annotation_audit$level == "fine", , drop = FALSE]
    fine_confidence <- data.frame(
      lineage = fine_confidence$lineage,
      cluster = fine_confidence$cluster,
      state = fine_confidence$assigned_label,
      cells = fine_confidence$n_cells,
      winner = fine_confidence$winner,
      runner_up = fine_confidence$runner_up,
      margin = fine_confidence$margin,
      margin_threshold = fine_confidence$margin_threshold,
      status = fine_confidence$annotation_status,
      artifact_winner = fine_confidence$artifact_winner,
      artifact_score = fine_confidence$artifact_score,
      assignment_reason = fine_confidence$assignment_reason,
      stringsAsFactors = FALSE
    )
    write.csv(
      fine_confidence,
      p(DIR$tables, "fine_state_annotation_confidence.csv"),
      row.names = FALSE
    )
    write_session_info("04_all_lineage_subclustering_and_fine_maps")
    cat("T/NK:", ncol(tnk), " Myeloid/DC:", ncol(my), "\n")
    my
  })
}

############################
# Stage 5：冻结 1500-HVG NMF 输入矩阵
############################

stage5_prepare_nmf_matrix <- function(force = FALSE) {
  stage_log("05_core_patient_statistics_tcell_and_NMF_input", {
    check_packages(c("Seurat","Matrix"))
    
    my_file <- p(DIR$cache_annotation, "fine_myeloid_dc_annotated.rds")
    t_file <- p(DIR$cache_annotation, "fine_t_nk_annotated.rds")
    assert_file(my_file); assert_file(t_file)
    if (!cache_has_current_version(my_file) ||
        !cache_has_current_version(t_file)) {
      stop("Stage 5 输入注释缓存不是当前 pipeline 版本；请先重新运行 Stage 4。")
    }
    out <- p(DIR$cache_nmf_matrix, "GSE178318_myeloid_NMF_input.rds")
    if (!force && cache_has_current_version(out) &&
        file.exists(p(DIR$tables, "18_patient_paired_state_abundance.csv")) &&
        file.exists(p(DIR$tables, "22_tcell_program_means_patient_tissue.csv"))) {
      message("[05_core_patient_statistics_tcell_and_NMF_input] CACHE HIT：直接返回 NMF 输入矩阵。")
      return(load_step_cache(out))
    }
    
    my <- load_step_cache(my_file)
    tnk <- load_step_cache(t_file)
    
    # ---- 原文章 patient-level state abundance ----
    paired_summary <- function(obj, lineage_name, labels_keep) {
      d <- obj[[]]
      d <- d[d$fine_status == "retained" &
               d$tissue %in% c("CRC","LM"), , drop = FALSE]
      d$one <- 1L
      patients <- sort(unique(d$patient_id))
      grid <- expand.grid(
        patient_id = patients, tissue = c("CRC","LM"),
        fine_label = labels_keep, stringsAsFactors = FALSE
      )
      counts <- aggregate(one ~ patient_id+tissue+fine_label, d, sum)
      den <- aggregate(one ~ patient_id+tissue, d, sum)
      x <- merge(grid, counts, all.x = TRUE)
      x$one[is.na(x$one)] <- 0L
      names(x)[names(x)=="one"] <- "n_state"
      x <- merge(x, den, by = c("patient_id","tissue"), all.x = TRUE)
      names(x)[names(x)=="one"] <- "n_lineage"
      x$proportion <- x$n_state/x$n_lineage
      
      sm <- do.call(rbind, lapply(labels_keep, function(label) {
        z <- x[x$fine_label == label, c("patient_id","tissue","proportion","n_state")]
        lm <- z[z$tissue=="LM",]; cr <- z[z$tissue=="CRC",]
        lm <- lm[match(patients,lm$patient_id),]
        cr <- cr[match(patients,cr$patient_id),]
        delta <- lm$proportion-cr$proportion
        dd <- delta[is.finite(delta)]
        nz <- dd[dd != 0]
        p_sign <- if (!length(nz)) 1 else
          stats::binom.test(sum(nz>0),length(nz),0.5)$p.value
        # Exploratory abundance audit; not the primary Figure 4A inference.
        p_wil <- if (length(dd)<2 || all(dd==0)) 1 else
          suppressWarnings(stats::wilcox.test(
            dd, mu=0, exact=FALSE, correct=FALSE
          )$p.value)
        data.frame(
          lineage=lineage_name,fine_label=label,
          n_patients=length(dd),n_nonzero_delta=length(nz),
          median_delta_LM_minus_CRC=median(dd),mean_delta_LM_minus_CRC=mean(dd),
          exact_sign_p=p_sign,paired_wilcoxon_p=p_wil,
          min_LM_cells=min(lm$n_state[is.finite(delta)],na.rm=TRUE),
          min_CRC_cells=min(cr$n_state[is.finite(delta)],na.rm=TRUE)
        )
      }))
      list(summary=sm,proportions=x)
    }
    
    t_labels <- sort(unique(tnk$fine_label[tnk$fine_status=="retained"]))
    m_labels <- sort(unique(my$fine_label[my$fine_status=="retained"]))
    tr <- paired_summary(tnk,"T_NK",t_labels)
    mr <- paired_summary(my,"Myeloid_DC",m_labels)
    st <- rbind(tr$summary,mr$summary)
    st$q_within_all_states <- p.adjust(st$paired_wilcoxon_p,"BH")
    write.csv(st,p(DIR$tables,"18_patient_paired_state_abundance.csv"),row.names=FALSE)
    write.csv(rbind(tr$proportions,mr$proportions),
              p(DIR$tables,"19_patient_state_proportions_long.csv"),row.names=FALSE)
    
    # ---- B7/AP patient-state audit ----
    my$B7_score <- score_genes_z(my,c("CD80","CD86"))
    my$AP_score <- score_genes_z(
      my,c("HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DPA1","HLA-DPB1","CD74","CIITA")
    )
    md <- my[[]]; md$barcode <- rownames(md)
    md$B7_score <- my$B7_score; md$AP_score <- my$AP_score; md$one <- 1L
    lm <- md[md$tissue=="LM" & md$fine_status=="retained",]
    prog <- aggregate(
      cbind(B7_score,AP_score,one) ~ patient_id+fine_label,lm,mean
    )
    counts <- aggregate(one ~ patient_id+fine_label,lm,sum)
    prog$n_cells <- counts$one[match(
      paste(prog$patient_id,prog$fine_label),
      paste(counts$patient_id,counts$fine_label)
    )]
    prog$one <- NULL
    write.csv(prog,p(DIR$tables,"20_LM_patient_state_AP_B7_scores.csv"),row.names=FALSE)
    target <- lm[lm$fine_label=="Mac_SPP1_MARCO_LPL",]
    if (nrow(target)) {
      pm <- aggregate(cbind(B7_score,AP_score)~patient_id,target,mean)
      rho <- if(nrow(pm)>=3) safe_cor(pm$B7_score,pm$AP_score) else NA_real_
      write.csv(data.frame(n_patients=nrow(pm),spearman_rho=rho),
                p(DIR$tables,"21_LM_SPP1_MARCO_LPL_AP_B7_patient_correlation.csv"),
                row.names=FALSE)
    }
    
    # ---- T-cell coexistence audit retained from the original scientific aim ----
    tnk$checkpoint_score <- score_genes_z(
      tnk,c("CXCL13","PDCD1","CTLA4","TOX","LAG3","HAVCR2")
    )
    tnk$cytotoxic_score <- score_genes_z(tnk,c("GZMB","PRF1","IFNG","NKG7","GZMH"))
    tnk$gzmK_score <- score_genes_z(tnk,c("GZMK","IL7R","LTB","MAL"))
    td <- tnk[[]]; td$one <- 1L
    td$comparison_group <- ifelse(
      td$fine_label %in% c("CD8_ZNF683_checkpoint_cytotoxic","CD8_CXCL13_checkpoint"),
      "checkpoint",
      ifelse(td$fine_label=="CD8_GZMK_memory","GZMK_reference",NA)
    )
    td <- td[!is.na(td$comparison_group) & td$tissue %in% c("CRC","LM"),]
    means <- aggregate(
      cbind(checkpoint_score,cytotoxic_score,gzmK_score,one) ~
        patient_id+tissue+comparison_group,
      td,mean
    )
    counts2 <- aggregate(one~patient_id+tissue+comparison_group,td,sum)
    means$n_cells <- counts2$one[match(
      paste(means$patient_id,means$tissue,means$comparison_group),
      paste(counts2$patient_id,counts2$tissue,counts2$comparison_group)
    )]
    means$one <- NULL
    write.csv(means,p(DIR$tables,"22_tcell_program_means_patient_tissue.csv"),
              row.names=FALSE)
    
    paired_t <- function(tissue_name) {
      z <- means[means$tissue==tissue_name,]
      a <- z[z$comparison_group=="checkpoint",]
      b <- z[z$comparison_group=="GZMK_reference",]
      pp <- intersect(a$patient_id,b$patient_id)
      a <- a[match(pp,a$patient_id),]; b <- b[match(pp,b$patient_id),]
      dd <- a$cytotoxic_score-b$cytotoxic_score
      nz <- dd[dd!=0]
      data.frame(
        tissue=tissue_name,n_patients=length(dd),
        median_cytotoxic_delta=if(length(dd))median(dd)else NA_real_,
        exact_sign_p=if(!length(nz))1 else stats::binom.test(sum(nz>0),length(nz),.5)$p.value,
        paired_wilcoxon_p=if(length(dd)<2)1 else
          suppressWarnings(stats::wilcox.test(dd,mu=0,exact=FALSE,correct=FALSE)$p.value)
      )
    }
    write.csv(
      rbind(paired_t("CRC"),paired_t("LM")),
      p(DIR$tables,"23_tcell_checkpoint_vs_gzmk_cytotoxic_patient_paired.csv"),
      row.names=FALSE
    )
    expr_t <- Seurat::GetAssayData(tnk,assay="RNA",layer="data")
    mask <- tnk$fine_label %in% c(
      "CD8_ZNF683_checkpoint_cytotoxic","CD8_CXCL13_checkpoint"
    ) & tnk$tissue=="LM"
    genes_t <- intersect(
      c("CXCL13","PDCD1","LAG3","HAVCR2","GZMB","PRF1","IFNG"),
      rownames(expr_t)
    )
    det <- as.matrix(expr_t[genes_t,mask,drop=FALSE])>0
    cyto <- intersect(c("GZMB","PRF1","IFNG"),genes_t)
    co <- data.frame(
      n_checkpoint_LM=sum(mask),
      n_genes_present=length(genes_t),
      pct_checkpoint_with_any_cytotoxic=
        if(length(cyto)&&ncol(det))mean(colSums(det[cyto,,drop=FALSE])>0)*100 else NA_real_,
      pct_checkpoint_with_GZMB=if("GZMB"%in%genes_t)mean(det["GZMB",])*100 else NA_real_,
      pct_checkpoint_with_PRF1=if("PRF1"%in%genes_t)mean(det["PRF1",])*100 else NA_real_,
      pct_checkpoint_with_IFNG=if("IFNG"%in%genes_t)mean(det["IFNG",])*100 else NA_real_
    )
    write.csv(co,p(DIR$tables,"24_tcell_checkpoint_cytotoxic_codetection_LM.csv"),
              row.names=FALSE)
    
    # ---- 1500 HVG NMF input，保持原文章 ----
    # [FIX] 原脚本此处误用了未定义的对象 `sub`（Stage 4 循环变量的残留），
    # 会导致 Stage 5 在构建 NMF 矩阵之前就以 "object 'sub' not found" 崩溃。
    # 这里应操作的对象是刚刚读入并计算过 T-cell/B7/AP 审计的髓系对象 `my`。
    SeuratObject::DefaultAssay(my) <- "RNA"
    genes <- Seurat::VariableFeatures(my)
    if (length(genes) < PARAM$nmf_hvg) {
      stop("Myeloid/DC VariableFeatures少于1500，不能按原文章继续。")
    }
    candidate <- genes[seq_len(PARAM$nmf_hvg)]
    X <- Seurat::GetAssayData(my,assay="RNA",layer="data")[candidate,,drop=FALSE]
    X <- as.matrix(X)
    X[!is.finite(X)] <- 0; X[X<0] <- 0
    nonzero <- rowSums(X)>0
    X <- X[nonzero,,drop=FALSE]
    genes <- candidate[nonzero]
    
    mm <- my[[]]; mm$barcode <- rownames(mm)
    nmf_input <- list(
      X=X,genes=genes,candidate_hvg_1500=candidate,
      cell_barcodes=colnames(my),metadata=mm,
      method=PARAM$nmf_method,ranks=PARAM$nmf_ranks,
      runs_per_rank=PARAM$nmf_runs_per_rank,
      pipeline_version=PIPELINE_VERSION
    )
    save_step_cache(nmf_input, out)
    write.csv(data.frame(rank=seq_along(genes),gene=genes),
              p(DIR$tables,"25_NMF_input_genes.csv"),row.names=FALSE)
    write.csv(
      data.frame(
        n_candidate_hvg=length(candidate),n_nonzero_genes=nrow(X),
        n_cells=ncol(X),n_patients=length(unique(mm$patient_id)),
        n_samples=length(unique(mm$sample_id)),
        n_fine_states=length(unique(mm$fine_label))
      ),
      p(DIR$audit,"05_NMF_input_summary.csv"),row.names=FALSE
    )
    write_session_info("05_core_patient_statistics_tcell_and_NMF_input")
    cat("NMF matrix:",nrow(X),"genes x",ncol(X),"cells\n")
    nmf_input
  })
}

############################
# Stage 6：全矩阵 rank 4–8 × 30 独立 Brunet fits
############################

# A deliberately small worker function for Windows-safe PSOCK execution.
# Only the 1,500 x n_cells NMF matrix is exported to workers; the full Seurat
# object and the stage closure are never serialized.  Each task receives its
# own explicit seed so parallel execution does not change the analysis design.
nmf_worker_single_fit <- function(task, Xmat, method, max_iter) {
  set.seed(as.integer(task$seed))
  tm <- system.time({
    fit <- NMF::nmf(
      Xmat, rank = as.integer(task$rank), method = method, nrun = 1L,
      maxIter = as.integer(max_iter), .options = ""
    )
  })
  list(
    rank = as.integer(task$rank),
    run = as.integer(task$run),
    seed = as.integer(task$seed),
    W = NMF::basis(fit),
    H = NMF::coef(fit),
    residual = tryCatch(as.numeric(NMF::residuals(fit))[1L],
                        error = function(e) NA_real_),
    n_iterations = tryCatch(as.integer(NMF::niter(fit)),
                            error = function(e) NA_integer_),
    elapsed_seconds = unname(tm[["elapsed"]])
  )
}

# PSOCK wrapper: the large input matrix is exported once per worker rather
# than serialized again as an argument for every task.  This is an execution
# and memory-safety change only; task seeds and NMF estimands are unchanged.
nmf_worker_single_fit_psock <- function(task, method, max_iter) {
  nmf_worker_single_fit(
    task, Xmat = .CRLM_NMF_X, method = method, max_iter = max_iter
  )
}

stage6_run_multiseed_nmf <- function(force = FALSE) {
  stage_log("06_multiseed_nmf", {
    check_packages(c("NMF","data.table"))
    # NMF's internal seed dispatch uses the package's S4 generic via
    # getGeneric("seed").  requireNamespace() alone does not attach/register
    # that generic in a fresh R session, which otherwise yields the opaque
    # "'what' must be a function or character string" error.  Attach NMF only
    # for this stage; all model calls remain explicitly namespaced below.
    if (!"package:NMF" %in% search()) {
      suppressPackageStartupMessages(library(NMF))
    }
    # NMF 0.28 performs its own internal function lookup.  Keep the future
    # backend sequential.  Stage 6 defaults to one fit at a time because
    # serializing this large sparse matrix to Windows PSOCK workers can fail;
    # explicit parallelism remains available through CRLM_NMF_WORKERS for
    # controlled environments.
    old_future_plan <- future::plan()
    future::plan(future::sequential)
    on.exit(future::plan(old_future_plan), add = TRUE)
    nmf_workers <- suppressWarnings(as.integer(
      Sys.getenv("CRLM_NMF_WORKERS", unset = "1")
    ))
    if (!is.finite(nmf_workers) || nmf_workers < 1L) nmf_workers <- 1L
    nmf_workers <- min(nmf_workers, parallel::detectCores(logical = TRUE))
    message("[06_multiseed_nmf] future plan: sequential; NMF PSOCK workers: ",
            nmf_workers)
    
    input <- p(DIR$cache_nmf_matrix, "GSE178318_myeloid_NMF_input.rds")
    assert_file(input)
    dat <- load_step_cache(input)
    if (is.null(dat)) stop("Stage 6 NMF 输入缓存缺失或版本不匹配：", input)
    X <- dat$X
    .CRLM_NMF_X <- X

    cl <- NULL
    if (nmf_workers > 1L) {
      cl <- parallel::makeCluster(nmf_workers, type = "PSOCK")
      # A worker can disappear on Windows when a large NMF object exhausts
      # transient memory.  Guard cleanup so a failed PSOCK launch does not
      # mask the original error or abort a resumable run.
      on.exit(
        if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE),
        add = TRUE
      )
      parallel::clusterEvalQ(
        cl,
        suppressPackageStartupMessages(library(NMF))
      )
      # Export only X and the small worker function.  In particular, do not
      # export `dat`, the Seurat object, or this stage's parent environment.
      parallel::clusterExport(
        cl,
        varlist = c(".CRLM_NMF_X", "nmf_worker_single_fit",
                    "nmf_worker_single_fit_psock"),
        envir = environment()
      )
    }
    
    run_log <- data.table::data.table()
    total_fits <- length(PARAM$nmf_ranks) * PARAM$nmf_runs_per_rank
    progress_stride <- 5L
    progress_file <- p(DIR$logs, "06_NMF_progress_current.txt")
    write_nmf_progress <- function(status, rank = NA_integer_, run = NA_integer_,
                                    completed = NA_integer_, started_at = NA,
                                    elapsed_seconds = NA_real_, message = "") {
      writeLines(c(
        paste0("status=", status),
        paste0("rank=", rank),
        paste0("run=", run),
        paste0("completed_fits=", completed),
        paste0("total_fits=", total_fits),
        paste0("started_at=", if (length(started_at) && !is.na(started_at))
          format(started_at, tz = "Asia/Shanghai") else ""),
        paste0("elapsed_seconds=", elapsed_seconds),
        paste0("updated_at=", format(Sys.time(), tz = "Asia/Shanghai")),
        paste0("message=", message)
      ), progress_file, useBytes = TRUE)
    }
    stage6_completed <- FALSE
    write_nmf_progress("started", completed = 0L,
                       message = "NMF input loaded; waiting for first fit")
    on.exit(
      if (!stage6_completed)
        write_nmf_progress("exited", message = "Stage 6 process exited before completion"),
      add = TRUE
    )
    
    stage_started_at <- Sys.time()
    completed_count <- 0L
    if (!force) {
      for (k0 in PARAM$nmf_ranks) {
        rank_dir0 <- p(DIR$cache_nmf_fits, paste0("rank_", k0))
        fs0 <- p(rank_dir0, sprintf("rank%d_run%02d.rds", k0,
                                    seq_len(PARAM$nmf_runs_per_rank)))
        completed_count <- completed_count + sum(vapply(
          fs0, cache_has_current_version, logical(1)
        ))
      }
    }

    for (k in PARAM$nmf_ranks) {
      rank_dir <- p(DIR$cache_nmf_fits, paste0("rank_", k))
      dir.create(rank_dir, recursive = TRUE, showWarnings = FALSE)
      runs <- seq_len(PARAM$nmf_runs_per_rank)
      files <- p(rank_dir, sprintf("rank%d_run%02d.rds", k, runs))
      keep <- if (force) rep(TRUE, length(runs)) else
        !vapply(files, cache_has_current_version, logical(1))
      pending <- runs[keep]
      if (!length(pending)) {
        message("[06_multiseed_nmf] rank ", k,
                ": all runs available in current/compatible cache")
        next
      }

      for (batch_start in seq(1L, length(pending), by = nmf_workers)) {
        batch_runs <- pending[batch_start:min(batch_start + nmf_workers - 1L,
                                              length(pending))]
        tasks <- lapply(batch_runs, function(run) list(
          rank = k,
          run = run,
          seed = PARAM$seed + 600000L + k * 1000L + run
        ))
        write_nmf_progress(
          "running", rank = k, run = min(batch_runs),
          completed = completed_count, started_at = stage_started_at,
          message = paste0("parallel batch runs ",
                           paste(sprintf("%02d", batch_runs), collapse = ","),
                           "; workers=", nmf_workers)
        )
        message("[06_multiseed_nmf] rank ", k, " batch ",
                paste(sprintf("run%02d", batch_runs), collapse = ", "),
                " / ", total_fits, " total fits; workers=", nmf_workers)
        batch_results <- if (nmf_workers > 1L) {
          tryCatch(
            parallel::parLapplyLB(
              cl, tasks, nmf_worker_single_fit_psock,
              method = PARAM$nmf_method,
              max_iter = PARAM$nmf_max_iter
            ),
            error = function(e) {
              # Resumable fallback: keep the exact task list and seeds, stop
              # the failed cluster, and finish this batch serially.  This is
              # an execution-level recovery only; it does not alter the NMF
              # design or any downstream estimand.
              message(
                "[06_multiseed_nmf] PSOCK batch failed (",
                conditionMessage(e),
                "); falling back to sequential execution for this and later batches."
              )
              try(parallel::stopCluster(cl), silent = TRUE)
              cl <<- NULL
              nmf_workers <<- 1L
              write_nmf_progress(
                "running", rank = k, run = min(batch_runs),
                completed = completed_count, started_at = stage_started_at,
                message = paste0(
                  "sequential fallback after PSOCK failure: ",
                  conditionMessage(e)
                )
              )
              lapply(
                tasks, nmf_worker_single_fit, Xmat = X,
                method = PARAM$nmf_method, max_iter = PARAM$nmf_max_iter
              )
            }
          )
        } else {
          lapply(tasks, nmf_worker_single_fit, Xmat = X,
                 method = PARAM$nmf_method, max_iter = PARAM$nmf_max_iter)
        }

        for (j in seq_along(batch_results)) {
          res <- batch_results[[j]]
          run <- res$run
          of <- p(rank_dir, sprintf("rank%d_run%02d.rds", k, run))
          W <- res$W
          H <- res$H
          rownames(W) <- rownames(X)
          colnames(W) <- paste0("factor_", seq_len(k))
          rownames(H) <- colnames(W)
          colnames(H) <- colnames(X)

          normed <- normalize_nmf_scale(W, H)
          W <- normed$W
          H <- normed$H
          rss <- compute_rss(X, W, H)
          n_iterations <- res$n_iterations
          converged <- is.finite(n_iterations) &&
            n_iterations < PARAM$nmf_max_iter
          convergence_status <- ifelse(
            !is.finite(n_iterations),
            "unknown",
            ifelse(
              n_iterations < PARAM$nmf_max_iter,
              "criterion_met_before_cap",
              "reached_cap_stop_reason_unavailable"
            )
          )

          save_step_cache(
            list(
              rank = k, run = run, seed = res$seed,
              method = PARAM$nmf_method, W = W, H = H,
              rss = rss, residual = res$residual,
              n_iterations = n_iterations, max_iter = PARAM$nmf_max_iter,
              converged = converged,
              convergence_status = convergence_status,
              elapsed_seconds = res$elapsed_seconds,
              genes = rownames(X), cells = colnames(X),
              n_genes = nrow(X), n_cells = ncol(X),
              R_version = R.version.string,
              NMF_version = as.character(packageVersion("NMF")),
              pipeline_version = PIPELINE_VERSION
            ),
            of
          )

          run_log <- data.table::rbindlist(list(
            run_log,
            data.table::data.table(
              rank = k, run = run, seed = res$seed,
              rss = rss, residual = res$residual,
              n_iterations = n_iterations, max_iter = PARAM$nmf_max_iter,
              converged = converged,
              convergence_status = convergence_status,
              elapsed_seconds = res$elapsed_seconds,
              n_genes = nrow(X), n_cells = ncol(X),
              R_version = R.version.string,
              NMF_version = as.character(packageVersion("NMF")),
              parallel_workers = nmf_workers,
              completed_at = format(Sys.time(), tz = "Asia/Shanghai")
            )
          ), fill = TRUE)
          data.table::fwrite(
            run_log,
            p(DIR$tables, "12_NMF_multiseed_run_log_current_session.csv")
          )
          completed_count <- completed_count + 1L
          mean_elapsed <- mean(run_log$elapsed_seconds, na.rm = TRUE)
          # Fits in a batch run concurrently; divide the serial estimate by
          # the number of active workers so the progress report reflects the
          # actual wall-clock throughput rather than CPU-time throughput.
          remaining_hours <- if (is.finite(mean_elapsed))
            (total_fits - completed_count) * mean_elapsed /
              (3600 * nmf_workers) else NA_real_
          write_nmf_progress(
            "completed_fit", rank = k, run = run,
            completed = completed_count, started_at = stage_started_at,
            elapsed_seconds = res$elapsed_seconds,
            message = paste0("iterations=", n_iterations,
                             "; converged_before_max_iter=", converged,
                             "; workers=", nmf_workers,
                             "; estimated_remaining_hours=",
                             if (is.finite(remaining_hours))
                               signif(remaining_hours, 3) else "NA")
          )
          message("[06_multiseed_nmf] completed rank ", k, " run ",
                  sprintf("%02d", run), " | ", completed_count, "/",
                  total_fits, " | elapsed=", round(res$elapsed_seconds, 1),
                  " s | converged=", converged)
          rm(W, H, normed, res); invisible(gc())
        }
      }
    }
    
    # 完成度审计
    rows <- list()
    for (k in PARAM$nmf_ranks) {
      rank_dir <- p(DIR$cache_nmf_fits, paste0("rank_", k))
      fs <- p(rank_dir, sprintf("rank%d_run%02d.rds", k, seq_len(PARAM$nmf_runs_per_rank)))
      fs_current <- vapply(fs, cache_has_current_version, logical(1))
      rows[[as.character(k)]] <- data.frame(
        rank = k,
        requested = PARAM$nmf_runs_per_rank,
        completed = sum(fs_current),
        complete = all(fs_current)
      )
    }
    status <- do.call(rbind, rows)
    write.csv(status, p(DIR$audit, "06_NMF_fit_completion.csv"), row.names = FALSE)
    if (!all(status$complete)) {
      stop("仍有NMF fit未完成。可直接重新运行Stage 6断点续跑。")
    }
    manifest_rows <- list(); mm <- 0L
    for (k in PARAM$nmf_ranks) {
      rank_dir <- p(DIR$cache_nmf_fits, paste0("rank_", k))
      for (run in seq_len(PARAM$nmf_runs_per_rank)) {
        fit_file <- p(rank_dir, sprintf("rank%d_run%02d.rds", k, run))
        f <- load_step_cache(fit_file)
        if (is.null(f)) stop("NMF fit 缓存缺失或版本不匹配：", fit_file)
        mm <- mm + 1L
        manifest_rows[[mm]] <- data.frame(
          rank = f$rank, run = f$run, seed = f$seed,
          method = f$method, rss = f$rss, residual = f$residual,
          n_iterations = f$n_iterations, max_iter = f$max_iter,
          converged = f$converged,
          convergence_status = if (!is.null(f$convergence_status))
            f$convergence_status else if (isTRUE(f$converged))
              "criterion_met_before_cap" else
              "reached_cap_stop_reason_unavailable",
          elapsed_seconds = f$elapsed_seconds,
          n_genes = f$n_genes, n_cells = f$n_cells,
          R_version = f$R_version, NMF_version = f$NMF_version,
          pipeline_version = f$pipeline_version,
          fit_file = p(rank_dir, sprintf("rank%d_run%02d.rds", k, run)),
          stringsAsFactors = FALSE
        )
      }
    }
    write.csv(do.call(rbind, manifest_rows),
              p(DIR$tables, "NMF_run_manifest.csv"), row.names = FALSE)
    write_session_info("06_multiseed_nmf")
    stage6_completed <- TRUE
    write_nmf_progress("completed", completed = total_fits,
                       message = "All requested NMF fits completed")
    invisible(status)
  })
}


############################
# Stage 7：LSAP 因子对齐、稳定性、rank重新选择
############################

stage7_nmf_stability_rank_selection <- function(force = FALSE) {
  stage_log("07_fullmatrix_NMF_consensus_stability_rank_selection", {
    check_packages(c("clue","data.table","ggplot2"))
    
    final_out <- p(DIR$cache_nmf_alignment, "rank_selection.rds")
    if (!force && cache_has_current_version(final_out)) {
      message("[07_fullmatrix_NMF_consensus_stability_rank_selection] CACHE HIT：直接返回 rank 选择结果。")
      return(load_step_cache(final_out))
    }
    
    nmf_input_file <- p(DIR$cache_nmf_matrix, "GSE178318_myeloid_NMF_input.rds")
    dat <- load_step_cache(nmf_input_file)
    if (is.null(dat)) stop("Stage 7 NMF 输入缓存缺失或版本不匹配：", nmf_input_file)
    X <- dat$X
    
    pair_rows <- list()
    factor_rows <- list()
    rank_rows <- list()
    convergence_rows <- list()
    alignment_rows <- list()
    consensus_diag_rows <- list()
    consensus_list <- list()
    
    # Full-matrix stability is primary.  Consensus matrices are retained on a
    # fixed, deterministic 999-cell balanced audit subset so that the
    # consensus object and cophenetic calculation remain inspectable without
    # allocating a 12,412 x 12,412 dense matrix for every rank.
    rank_step_total <- length(PARAM$nmf_ranks) + 1L
    cons_idx <- stage_step(
      "07_fullmatrix_NMF_consensus_stability_rank_selection", 1L,
      rank_step_total,
      "选择固定的 999-cell consensus audit 子集",
      {
        select_consensus_cells(dat$metadata, target = 999L, seed = 257319L)
      }, "<5 min"
    )
    write.csv(
      data.frame(
        index = cons_idx,
        barcode = dat$cell_barcodes[cons_idx],
        fine_label = if ("fine_label" %in% names(dat$metadata))
          as.character(dat$metadata$fine_label[cons_idx]) else "all",
        seed = 257319L
      ),
      p(DIR$tables, "NMF_consensus_audit_cells.csv"), row.names = FALSE
    )
    
    for (k in PARAM$nmf_ranks) {
      stage_step(
        "07_fullmatrix_NMF_consensus_stability_rank_selection",
        match(k, PARAM$nmf_ranks) + 1L, rank_step_total,
        paste0("rank ", k, "：因子对齐、稳定性与 consensus 计算"),
        {
          cat("\n===== stability rank",k,"=====\n")
          rank_dir <- p(DIR$cache_nmf_fits,paste0("rank_",k))
          fits <- lapply(seq_len(PARAM$nmf_runs_per_rank),function(run){
            f <- p(rank_dir,sprintf("rank%d_run%02d.rds",k,run))
            assert_file(f)
            fit <- load_step_cache(f)
            if (is.null(fit)) stop("NMF fit 缓存缺失或版本不匹配：", f)
            fit
          })
          # NMF 0.28 does not provide a method-independent convergence flag
          # that can be safely serialized and compared across fits.  The
          # historical `converged` field is therefore an iteration-ceiling
          # proxy, not a package-certified convergence result.  Block only
          # malformed fits with missing/non-finite iteration counts; fits that
          # reach maxIter remain eligible for the prespecified alignment,
          # stability and RSS diagnostics, with an explicit audit message.
          bad_runs <- which(!vapply(
            fits,
            function(z) is.list(z) && length(z$n_iterations) == 1L &&
              is.finite(z$n_iterations),
            logical(1)
          ))
          if (length(bad_runs)) {
            stop(
              "Stage 7发现缺少有效迭代记录的NMF fits：rank ", k,
              " runs ", paste(bad_runs, collapse = ", "),
              "。禁止将无法审计的fit用于稳定性和rank选择。"
            )
          }
          ceiling_runs <- which(vapply(
            fits,
            function(z) z$n_iterations >= z$max_iter,
            logical(1)
          ))
          if (length(ceiling_runs)) {
            message(
              "[07_NMF] rank ", k, ": ", length(ceiling_runs),
              "/", length(fits),
              " fits reached maxIter (", fits[[1]]$max_iter,
              "); iteration-ceiling status is retained for audit and is not",
              " treated as an automatic fit failure."
            )
          }
          
          # 1) 全pair matching，用于run-level稳定性和medoid选择
          pr <- list(); qn <- 0L
          pair_mean <- list(); pmn <- 0L
          for(i in 1:(length(fits)-1L)){
            for(j in (i+1L):length(fits)){
              mt <- match_factors_lsap(fits[[i]]$W,fits[[j]]$W)
              tmp <- mt
              tmp$rank <- k; tmp$run_i <- i; tmp$run_j <- j
              tmp$usage_spearman <- vapply(seq_len(nrow(tmp)),function(z){
                safe_cor(
                  fits[[i]]$H[tmp$factor_ref[z],],
                  fits[[j]]$H[tmp$factor_other[z],]
                )
              },numeric(1))
              tmp$top15_overlap_fraction <- vapply(seq_len(nrow(tmp)),function(z){
                top_gene_overlap(
                  fits[[i]]$W,tmp$factor_ref[z],
                  fits[[j]]$W,tmp$factor_other[z],15L
                )
              },numeric(1))
              qn <- qn+1L; pr[[qn]] <- tmp
              pmn <- pmn+1L
              pair_mean[[pmn]] <- data.frame(
                run_i=i,run_j=j,
                mean_basis_spearman=mean(tmp$basis_spearman,na.rm=TRUE),
                mean_usage_spearman=mean(tmp$usage_spearman,na.rm=TRUE)
              )
            }
          }
          pr <- do.call(rbind,pr)
          pm <- do.call(rbind,pair_mean)
          pair_rows[[as.character(k)]] <- pr
          
          run_stab_basis <- vapply(seq_len(length(fits)),function(r){
            z <- pm[pm$run_i==r | pm$run_j==r,]
            if(!nrow(z)) return(NA_real_)
            mean(z$mean_basis_spearman,na.rm=TRUE)
          },numeric(1))
          run_stab_usage <- vapply(seq_len(length(fits)),function(r){
            z <- pm[pm$run_i==r | pm$run_j==r,]
            if(!nrow(z)) return(NA_real_)
            mean(z$mean_usage_spearman,na.rm=TRUE)
          },numeric(1))
          # Pairwise basis Spearman is the primary run-level metric; larger values
          # indicate greater stability. Usage stability is used only to break a
          # primary-metric tie, or as a rank-wide fallback when every basis metric
          # is unavailable. This avoids comparing the two metrics as if they were
          # on an interchangeable scale.
          basis_candidates <- which(is.finite(run_stab_basis))
          if(length(basis_candidates)){
            best_stab <- max(run_stab_basis[basis_candidates])
            medoid_candidates <- basis_candidates[
              abs(run_stab_basis[basis_candidates] - best_stab) <=
                sqrt(.Machine$double.eps)
            ]
            medoid_selection_metric <- "mean_pairwise_basis_spearman"
            if(length(medoid_candidates) > 1L){
              usage_tie <- run_stab_usage[medoid_candidates]
              if(any(is.finite(usage_tie))){
                best_usage <- max(usage_tie[is.finite(usage_tie)])
                medoid_candidates <- medoid_candidates[
                  is.finite(usage_tie) &
                    abs(usage_tie - best_usage) <= sqrt(.Machine$double.eps)
                ]
                medoid_selection_metric <-
                  "mean_pairwise_basis_spearman_then_usage_spearman"
              }
            }
            run_stab <- run_stab_basis
          } else {
            usage_candidates <- which(is.finite(run_stab_usage))
            if(!length(usage_candidates)){
              stop("NMF run-level stability metrics are all unavailable; cannot choose a reference run.")
            }
            best_stab <- max(run_stab_usage[usage_candidates])
            medoid_candidates <- usage_candidates[
              abs(run_stab_usage[usage_candidates] - best_stab) <=
                sqrt(.Machine$double.eps)
            ]
            medoid_selection_metric <- "mean_pairwise_usage_spearman_fallback"
            run_stab <- run_stab_usage
          }
          if(!length(medoid_candidates)){
            stop("NMF run-level stability metrics are all unavailable; cannot choose a reference run.")
          }
          # Deterministic tie rule: retain the smallest run index among tied
          # highest-stability runs, rather than defaulting to the first run.
          medoid_run <- min(medoid_candidates)
          med <- fits[[medoid_run]]
          
          # 2) 每个run对齐到medoid，同一factor编号才能形成真正的共识摘要
          nrun <- length(fits)
          W_arr <- array(
            NA_real_,
            dim=c(nrow(med$W),k,nrun),
            dimnames=list(rownames(med$W),paste0("factor_",1:k),paste0("run_",1:nrun))
          )
          H_arr <- array(
            NA_real_,
            dim=c(k,ncol(med$H),nrun),
            dimnames=list(paste0("factor_",1:k),colnames(med$H),paste0("run_",1:nrun))
          )
          medoid_match <- vector("list",nrun)
          for(r in seq_len(nrun)){
            if(r==medoid_run){
              mt <- data.frame(
                factor_ref=1:k,factor_other=1:k,
                basis_spearman=rep(1,k)
              )
            } else {
              mt <- match_factors_lsap(med$W,fits[[r]]$W)
            }
            ord <- mt$factor_other[match(1:k,mt$factor_ref)]
            W_arr[,,r] <- fits[[r]]$W[,ord,drop=FALSE]
            H_arr[,,r] <- fits[[r]]$H[ord,,drop=FALSE]
            medoid_match[[r]] <- data.frame(
              rank=k,run=r,medoid_run=medoid_run,
              factor=1:k,matched_factor_in_run=ord,
              loading_spearman=mt$basis_spearman[match(1:k,mt$factor_ref)]
            )
          }
          med_match <- do.call(rbind,medoid_match)
          alignment_rows[[as.character(k)]] <- med_match
          
          # 3) factor-level recurrence/usage/top15
          fr <- lapply(1:k,function(f){
            rho_w <- vapply(seq_len(nrun),function(r){
              safe_cor(med$W[,f],W_arr[,f,r])
            },numeric(1))
            rho_h <- vapply(seq_len(nrun),function(r){
              safe_cor(med$H[f,],H_arr[f,,r])
            },numeric(1))
            ov <- vapply(seq_len(nrun),function(r){
              if(r==medoid_run) return(1)
              g1 <- names(sort(med$W[,f],decreasing=TRUE))[1:15]
              g2 <- names(sort(W_arr[,f,r],decreasing=TRUE))[1:15]
              length(intersect(g1,g2))/15
            },numeric(1))
            data.frame(
              rank=k,factor=f,medoid_run=medoid_run,
              recurrence_fraction=mean(rho_w>=PARAM$nmf_loading_spearman_threshold,na.rm=TRUE),
              usage_stability_fraction=mean(rho_h>=PARAM$nmf_usage_spearman_threshold,na.rm=TRUE),
              median_loading_spearman=median(rho_w,na.rm=TRUE),
              q05_loading_spearman=quantile(rho_w,.05,na.rm=TRUE),
              median_usage_spearman=median(rho_h,na.rm=TRUE),
              median_top15_overlap_fraction=median(ov,na.rm=TRUE),
              stable_factor=
                mean(rho_w>=PARAM$nmf_loading_spearman_threshold,na.rm=TRUE) >=
                PARAM$factor_recurrence_threshold &&
                median(rho_w,na.rm=TRUE) >= PARAM$nmf_loading_spearman_threshold &&
                median(ov,na.rm=TRUE) >= PARAM$top15_overlap_threshold
            )
          })
          fr <- do.call(rbind,fr)
          factor_rows[[as.character(k)]] <- fr
          
          # 4) consensus：对尺度已统一、factor已对齐的30个解逐元素取median
          Wc <- apply(W_arr,c(1,2),median,na.rm=TRUE)
          Hc <- apply(H_arr,c(1,2),median,na.rm=TRUE)
          rownames(Wc) <- rownames(med$W); colnames(Wc) <- paste0("factor_",1:k)
          rownames(Hc) <- colnames(Wc); colnames(Hc) <- colnames(med$H)
          normed <- normalize_nmf_scale(Wc,Hc)
          Wc <- normed$W; Hc <- normed$H
          consensus_rss <- compute_rss(X,Wc,Hc)
          
          # Consensus matrix and cophenetic coefficient on the fixed balanced
          # audit subset.  This is an auxiliary rank diagnostic, not a hidden
          # replacement for the full-matrix recurrence metrics.
          hard <- apply(
            H_arr[, cons_idx, , drop = FALSE], c(2, 3), which.max
          )
          if (is.null(dim(hard))) hard <- matrix(hard, ncol = nrun)
          C <- matrix(0, nrow = length(cons_idx), ncol = length(cons_idx))
          for (rr in seq_len(nrun)) {
            C <- C + outer(hard[, rr], hard[, rr], FUN = "==")
          }
          C <- C / nrun
          diag(C) <- 1
          consensus_dist <- stats::as.dist(1 - C)
          hc <- stats::hclust(consensus_dist, method = "average")
          cophenetic_dist <- stats::cophenetic(hc)
          cophenetic_r <- safe_cor(
            as.numeric(consensus_dist),
            as.numeric(cophenetic_dist),
            method = "pearson"
          )
          consensus_file <- p(
            DIR$cache_nmf_alignment,
            paste0("consensus_matrix_rank", k, ".rds")
          )
          consensus_pdf <- p(DIR$figures, paste0("NMF_rank", k, "_consensus.pdf"))
          save_step_cache(
            list(
              rank = k, matrix = C, cells = dat$cell_barcodes[cons_idx],
              n_runs = nrun, subset_size = length(cons_idx), seed = 257319L,
              cophenetic_coefficient = cophenetic_r,
              method = "hard-assignment consensus on balanced audit subset"
            ),
            consensus_file
          )
          # Human-readable consensus diagnostic.  The matrix is computed on the
          # fixed balanced 999-cell audit subset; this PDF is diagnostic only and
          # is not used to choose P1 or to inspect external outcomes.
          grDevices::pdf(consensus_pdf, width = 6, height = 6, useDingbats = FALSE)
          graphics::par(mar = c(4, 4, 3, 1))
          graphics::image(
            x = seq_len(nrow(C)), y = seq_len(ncol(C)),
            z = C[nrow(C):1, , drop = FALSE],
            col = grDevices::colorRampPalette(c("white", "#C6DBEF", "#2171B5"))(101),
            axes = FALSE, xlab = "Audit cells", ylab = "Audit cells",
            main = paste0("NMF rank ", k, " consensus matrix")
          )
          graphics::box()
          grDevices::dev.off()
          consensus_diag_rows[[as.character(k)]] <- data.frame(
              rank = k, subset_size = length(cons_idx), n_runs = nrun,
              cophenetic_coefficient = cophenetic_r,
              cophenetic_input = "1-consensus_matrix",
              consensus_distance_definition = "D = 1 - C",
              consensus_matrix_file = consensus_file,
              consensus_pdf = consensus_pdf,
            stringsAsFactors = FALSE
          )
          
          cons <- list(
            rank=k,medoid_run=medoid_run,W=Wc,H=Hc,
            factor_stability=fr,medoid_matching=med_match,
            consensus_rss=consensus_rss,
            method=paste0(PARAM$nmf_method,"_aligned_median_consensus"),
            n_runs=nrun,pipeline_version=PIPELINE_VERSION
          )
          cf <- p(DIR$cache_nmf_alignment,paste0("rank",k,"_consensus.rds"))
          save_step_cache(cons,cf)
          consensus_list[[as.character(k)]] <- cons
          
          # 5) 10→20→30 共识收敛：固定同一medoid alignment，不重新挑factor
          finalW <- Wc
          for(nn in c(10L,20L,30L)){
            Wn <- apply(W_arr[,,seq_len(nn),drop=FALSE],c(1,2),median,na.rm=TRUE)
            Hn <- apply(H_arr[,,seq_len(nn),drop=FALSE],c(1,2),median,na.rm=TRUE)
            nrm <- normalize_nmf_scale(Wn,Hn); Wn <- nrm$W
            rhos <- vapply(1:k,function(f) safe_cor(Wn[,f],finalW[,f]),numeric(1))
            convergence_rows[[length(convergence_rows)+1L]] <- data.frame(
              rank=k,first_n_runs=nn,
              median_consensus_loading_rho_vs_30=median(rhos,na.rm=TRUE),
              min_consensus_loading_rho_vs_30=min(rhos,na.rm=TRUE),
              n_factors=k
            )
          }
          
          rss_runs <- vapply(fits,function(z)z$rss,numeric(1))
          rank_rows[[as.character(k)]] <- data.frame(
            rank=k,n_runs=nrun,medoid_run=medoid_run,
            median_pairwise_basis_spearman=median(pr$basis_spearman,na.rm=TRUE),
            q05_pairwise_basis_spearman=quantile(pr$basis_spearman,.05,na.rm=TRUE),
            median_pairwise_usage_spearman=median(pr$usage_spearman,na.rm=TRUE),
            median_pairwise_top15_overlap=median(pr$top15_overlap_fraction,na.rm=TRUE),
            n_stable_factors=sum(fr$stable_factor),
            stable_factor_fraction=mean(fr$stable_factor),
            median_factor_recurrence=median(fr$recurrence_fraction),
            median_factor_usage_stability=median(fr$usage_stability_fraction),
            median_factor_loading_spearman=median(fr$median_loading_spearman),
            median_factor_top15_overlap=median(fr$median_top15_overlap_fraction),
            median_run_rss=median(rss_runs,na.rm=TRUE),
            consensus_rss=consensus_rss,
            medoid_mean_pairwise_basis=run_stab_basis[medoid_run],
            medoid_mean_pairwise_usage=run_stab_usage[medoid_run],
            medoid_selection_stability=run_stab[medoid_run],
            medoid_selection_metric=medoid_selection_metric
          )
          rm(fits,W_arr,H_arr); invisible(gc())
        }, "5–30 min"
      )
    }
    
    pair_all <- do.call(rbind,pair_rows)
    factor_all <- do.call(rbind,factor_rows)
    rank_summary <- do.call(rbind,rank_rows)
    convergence <- do.call(rbind,convergence_rows)
    alignment_all <- do.call(rbind, alignment_rows)
    consensus_diag <- do.call(rbind, consensus_diag_rows)
    rank_summary <- rank_summary[order(rank_summary$rank),]
    rank_summary <- merge(rank_summary, consensus_diag, by = "rank", all.x = TRUE,
                          sort = FALSE)
    rank_summary <- rank_summary[order(rank_summary$rank),]
    rank_summary$rss_improvement_from_previous <- c(
      NA,
      (head(rank_summary$consensus_rss,-1)-tail(rank_summary$consensus_rss,-1))/
        head(rank_summary$consensus_rss,-1)
    )
    
    # rank选择完全机械：不看PT-LM、不看外部、不看SPP1/生物解释
    elbow_rank <- compute_elbow_rank(rank_summary$rank,rank_summary$consensus_rss)
    eligible <- rank_summary$rank[
      rank_summary$stable_factor_fraction >=
        PARAM$eligible_rank_stable_factor_fraction &
        rank_summary$median_factor_loading_spearman >=
        PARAM$nmf_loading_spearman_threshold &
        rank_summary$median_factor_top15_overlap >=
        PARAM$top15_overlap_threshold
    ]
    if(length(eligible)){
      dist <- abs(eligible-elbow_rank)
      cand <- eligible[dist==min(dist)]
      selected_rank <- min(cand)  # 同距离时取较小rank，避免过度分解
      selection_reason <- paste0(
        "mechanical: eligible stability ranks; choose nearest to consensus-RSS elbow; ",
        "tie -> smaller rank"
      )
    } else {
      # 无rank达硬阈值时，优先stable-factor fraction，其次loading稳定，再较小rank
      q <- rank_summary[order(
        -rank_summary$stable_factor_fraction,
        -rank_summary$median_factor_loading_spearman,
        rank_summary$rank
      ),]
      selected_rank <- q$rank[1]
      selection_reason <- paste0(
        "no rank met hard eligibility; fallback maximum stable-factor fraction, ",
        "then loading stability, then smaller rank"
      )
      warning(selection_reason)
    }
    
    # 跨相邻rank program persistence（审计，不参与rank选择）
    cr <- list(); cc <- 0L
    for(a in seq_len(length(PARAM$nmf_ranks)-1L)){
      k1 <- PARAM$nmf_ranks[a]; k2 <- PARAM$nmf_ranks[a+1L]
      W1 <- consensus_list[[as.character(k1)]]$W
      W2 <- consensus_list[[as.character(k2)]]$W
      S <- spearman_matrix(W1,W2)
      for(f in seq_len(ncol(W1))){
        cc <- cc+1L
        cr[[cc]] <- data.frame(
          rank_from=k1,factor_from=f,rank_to=k2,
          best_factor_to=which.max(S[f,]),
          best_loading_spearman=max(S[f,])
        )
      }
    }
    cross_rank <- do.call(rbind,cr)
    
    write.csv(pair_all,p(DIR$tables,"26_NMF_pairwise_factor_stability.csv"),row.names=FALSE)
    write.csv(alignment_all,p(DIR$tables,"25_NMF_factor_alignment.csv"),row.names=FALSE)
    write.csv(factor_all,p(DIR$tables,"27_NMF_factor_recurrence_stability.csv"),row.names=FALSE)
    write.csv(rank_summary,p(DIR$tables,"28_NMF_rank_summary.csv"),row.names=FALSE)
    write.csv(convergence,p(DIR$tables,"29_NMF_10_20_30_consensus_convergence.csv"),row.names=FALSE)
    write.csv(consensus_diag,p(DIR$tables,"31b_NMF_cophenetic_coefficients.csv"),row.names=FALSE)
    write.csv(consensus_diag,p(DIR$tables,"cophenetic_rank_summary.csv"),row.names=FALSE)
    write.csv(cross_rank,p(DIR$tables,"30_NMF_cross_rank_persistence.csv"),row.names=FALSE)
    write.csv(
      data.frame(
        selected_rank=selected_rank,
        rss_elbow_rank=elbow_rank,
        nmf_loading_spearman_threshold=PARAM$nmf_loading_spearman_threshold,
        nmf_usage_spearman_threshold=PARAM$nmf_usage_spearman_threshold,
        recurrence_threshold=PARAM$factor_recurrence_threshold,
        top15_overlap_threshold=PARAM$top15_overlap_threshold,
        eligible_rank_stable_factor_fraction=
          PARAM$eligible_rank_stable_factor_fraction,
        eligible_ranks=paste(eligible,collapse=";"),
        selection_reason=selection_reason
      ),
      p(DIR$tables,"31_NMF_selected_rank.csv"),row.names=FALSE
    )
    
    # 两张独立诊断图（避免复合图）
    g1 <- ggplot2::ggplot(
      rank_summary,ggplot2::aes(rank,stable_factor_fraction)
    ) +
      ggplot2::geom_line() + ggplot2::geom_point(size=2) +
      ggplot2::geom_hline(
        yintercept=PARAM$eligible_rank_stable_factor_fraction,linetype="dashed"
      ) +
      ggplot2::geom_vline(xintercept=selected_rank,linetype="dotted") +
      ggplot2::scale_x_continuous(breaks=PARAM$nmf_ranks) +
      ggplot2::labs(x="Rank",y="Stable factor fraction",
                    title="Full-matrix NMF stability") +
      ggplot2::theme_classic()
    ggplot2::ggsave(
      p(DIR$figures,"NMF_rank_stable_factor_fraction.pdf"),
      g1,width=6,height=4
    )
    g2 <- ggplot2::ggplot(
      rank_summary,ggplot2::aes(rank,consensus_rss)
    ) +
      ggplot2::geom_line()+ggplot2::geom_point(size=2) +
      ggplot2::geom_vline(xintercept=elbow_rank,linetype="dashed") +
      ggplot2::geom_vline(xintercept=selected_rank,linetype="dotted") +
      ggplot2::scale_x_continuous(breaks=PARAM$nmf_ranks) +
      ggplot2::labs(x="Rank",y="Consensus reconstruction RSS",
                    title="NMF reconstruction curve") +
      ggplot2::theme_classic()
    ggplot2::ggsave(
      p(DIR$figures,"NMF_rank_consensus_RSS.pdf"),
      g2,width=6,height=4
    )
    
    ans <- list(
      pairwise=pair_all,factor_stability=factor_all,
      rank_summary=rank_summary,convergence=convergence,
      cross_rank=cross_rank,selected_rank=selected_rank,
      elbow_rank=elbow_rank,eligible_ranks=eligible,
      selection_reason=selection_reason,
      consensus_files=setNames(
        p(DIR$cache_nmf_alignment,paste0("rank",PARAM$nmf_ranks,"_consensus.rds")),
        PARAM$nmf_ranks
      )
    )
    save_step_cache(ans, final_out)
    write_session_info("07_fullmatrix_NMF_consensus_stability_rank_selection")
    cat("Selected rank:",selected_rank,"\n")
    cat("RSS elbow:",elbow_rank,"\n")
    cat("Eligible:",paste(eligible,collapse=","),"\n")
    ans
  })
}

############################
# Stage 8：P1 correspondence audit + 新版P1冻结 + discovery usage
############################

stage8_freeze_p1 <- function(force = FALSE) {
  stage_log("08_post_rank_P1_correspondence_and_freeze", {
    check_packages(c("Seurat","data.table","ggplot2"))
    
    sel_file <- p(DIR$cache_nmf_alignment,"rank_selection.rds")
    obj_file <- p(DIR$cache_annotation,"fine_myeloid_dc_annotated.rds")
    assert_file(sel_file); assert_file(obj_file)
    out <- p(DIR$cache_projection,"frozen_P1_v4.rds")
    if(!force && cache_has_current_version(out)){
      message("[08_post_rank_P1_correspondence_and_freeze] CACHE HIT：直接返回冻结 P1。")
      return(must_load_step_cache(out))
    }
    
    sel <- must_load_step_cache(sel_file, "Stage 7 rank selection")
    k <- sel$selected_rank
    cons_file <- p(DIR$cache_nmf_alignment,paste0("rank",k,"_consensus.rds"))
    assert_file(cons_file)
    cons <- must_load_step_cache(cons_file, "Stage 7 consensus")
    W <- cons$W; H <- cons$H
    
    # P1选择只使用发现队列内预先冻结的anchor module；
    # 不读取外部队列、空间结果、患者分解或legacy P1来决定P1。
    fr <- cons$factor_stability
    stable_factors <- fr$factor[fr$stable_factor]
    if(!length(stable_factors)){
      stop("选定rank没有满足稳定性阈值的factor，不能冻结P1。")
    }
    anchor_genes <- intersect(P1_ANCHOR_GENES, rownames(W))
    if (length(intersect(P1_CORE_ANCHOR_GENES, anchor_genes)) < 2L ||
        length(anchor_genes) < 4L) {
      stop(
        "P1 anchor module可用基因不足：需要两个core anchors及至少4个模块基因；实际为 ",
        paste(anchor_genes, collapse = ",")
      )
    }
    anchor_mat <- W[anchor_genes, stable_factors, drop=FALSE]
    anchor_z <- t(scale(t(anchor_mat)))
    anchor_z[!is.finite(anchor_z)] <- 0
    anchor_score <- rep(NA_real_, ncol(W))
    anchor_score[stable_factors] <- colMeans(anchor_z)
    stable_anchor_scores <- anchor_score[stable_factors]
    score_order <- order(stable_anchor_scores, decreasing = TRUE)
    best_anchor <- stable_anchor_scores[score_order[1L]]
    runner_up_anchor <- if (length(score_order) > 1L) {
      stable_anchor_scores[score_order[2L]]
    } else {
      NA_real_
    }
    anchor_margin <- best_anchor - runner_up_anchor
    p1_selection_low_confidence <- is.finite(anchor_margin) &&
      anchor_margin < PARAM$p1_anchor_min_margin
    p1_selection_status <- if (p1_selection_low_confidence) {
      "P1_selection_low_confidence"
    } else if (is.finite(anchor_margin)) {
      "P1_selection_above_min_margin"
    } else {
      "P1_selection_no_runner_up"
    }
    cand <- stable_factors[abs(anchor_score[stable_factors]-best_anchor) < 1e-12]
    if(length(cand) > 1L){
      cand <- cand[which.max(fr$median_loading_spearman[match(cand, fr$factor)])]
    }
    p1_factor <- cand[1]
    matching_source <- "discovery_prespecified_anchor_module"
    legacy_corv <- rep(NA_real_, ncol(W))
    legacy_mass <- rep(NA_real_, ncol(W))
    g_legacy <- intersect(LEGACY_P1_TOP15$gene, rownames(W))
    if(length(g_legacy) >= 10L){
      ref_legacy <- LEGACY_P1_TOP15$weight_norm[match(g_legacy,LEGACY_P1_TOP15$gene)]
      legacy_corv <- apply(W[g_legacy,,drop=FALSE],2,function(v)safe_cor(v,ref_legacy))
      legacy_mass <- colSums(W[g_legacy,,drop=FALSE])
    }
    
    audit <- data.frame(
      factor=seq_len(ncol(W)),
      stable_factor=seq_len(ncol(W)) %in% stable_factors,
      anchor_genes=paste(anchor_genes, collapse=";"),
      anchor_score=anchor_score,
      selected_anchor_score=best_anchor,
      runner_up_anchor_factor=if (length(score_order) > 1L)
        stable_factors[score_order[2L]] else NA_integer_,
      runner_up_anchor_score=runner_up_anchor,
      anchor_margin=anchor_margin,
      p1_anchor_min_margin=PARAM$p1_anchor_min_margin,
      p1_selection_status=p1_selection_status,
      legacy_spearman=legacy_corv,
      legacy_top15_loading_mass=legacy_mass,
      selected_as_P1=seq_len(ncol(W))==p1_factor,
      matching_source=matching_source
    )
    cat(
      "[08_post_rank_P1_correspondence_and_freeze] P1 anchor winner factor ",
      p1_factor, " | winner score: ", signif(best_anchor, 5),
      " | runner-up score: ", signif(runner_up_anchor, 5),
      " | margin: ", signif(anchor_margin, 5),
      " | threshold: ", PARAM$p1_anchor_min_margin,
      " | status: ", p1_selection_status, "\n", sep = ""
    )
    write.csv(audit,p(DIR$tables,"32_P1_correspondence_after_rank_selection.csv"),
              row.names=FALSE)
    write.csv(audit,p(DIR$tables,"P1_selection_audit.csv"),row.names=FALSE)
    
    # P1移到第1列，其他program维持consensus factor顺序。
    ord <- c(p1_factor,setdiff(seq_len(k),p1_factor))
    W2 <- W[,ord,drop=FALSE]
    H2 <- H[ord,,drop=FALSE]
    colnames(W2) <- paste0("P",seq_len(k))
    rownames(H2) <- colnames(W2)
    
    top50 <- do.call(rbind,lapply(seq_len(k),function(j){
      oo <- order(W2[,j],decreasing=TRUE)[seq_len(min(50L,nrow(W2)))]
      data.frame(
        program=paste0("P",j),rank=seq_along(oo),
        gene=rownames(W2)[oo],weight=W2[oo,j]
      )
    }))
    top50$weight_norm_within_top15 <- NA_real_
    for(pp in unique(top50$program)){
      ii <- top50$program==pp & top50$rank<=15
      top50$weight_norm_within_top15[ii] <- top50$weight[ii]/sum(top50$weight[ii])
    }
    write.csv(top50,p(DIR$tables,"33_NMF_program_top50_genes.csv"),row.names=FALSE)
    
    p1_top15 <- top50[top50$program=="P1" & top50$rank<=15,]
    p1_top15$weight_norm <- p1_top15$weight/sum(p1_top15$weight)
    write.csv(p1_top15,p(DIR$tables,"34_frozen_P1_top15_weights.csv"),row.names=FALSE)
    write.csv(p1_top15,p(DIR$tables,"frozen_P1_v4_weights.csv"),row.names=FALSE)
    
    all_top15 <- top50[top50$rank<=15,]
    write.csv(all_top15,p(DIR$tables,"35_frozen_all_program_top15_weights.csv"),
              row.names=FALSE)
    
    # discovery usage
    obj <- must_load_step_cache(obj_file, "Stage 4 myeloid annotation")
    meta <- obj[[]]; meta$barcode <- rownames(meta)
    usage <- data.frame(barcode=colnames(H2),t(H2),check.names=FALSE)
    usage <- merge(
      meta[,c("barcode","sample_id","patient_id","tissue","fine_label","fine_status")],
      usage,by="barcode",sort=FALSE
    )
    con <- gzfile(p(DIR$tables,"36_discovery_program_usage.csv.gz"),"wt")
    write.csv(usage,con,row.names=FALSE); close(con)
    patient_usage <- aggregate(
      usage[,paste0("P",seq_len(k)),drop=FALSE],
      usage[,c("patient_id","sample_id","tissue")],mean
    )
    write.csv(patient_usage,p(DIR$tables,"37_discovery_program_patient_usage.csv"),
              row.names=FALSE)
    
    # frozen top15 projection fidelity against consensus P1 H
    SeuratObject::DefaultAssay(obj) <- "RNA"
    norm <- Seurat::GetAssayData(obj,assay="RNA",layer="data")
    wdf <- data.frame(gene=p1_top15$gene,weight=p1_top15$weight_norm)
    proj <- weighted_score(norm,wdf)
    h1 <- as.numeric(H2["P1",colnames(obj)])
    cell <- data.frame(
      barcode=colnames(obj),patient=as.character(obj$patient_id),
      sample=as.character(obj$sample_id),tissue=as.character(obj$tissue),
      fine_state=as.character(obj$fine_label),
      consensus_nmf_P1_usage=h1,frozen_top15_weighted_score=proj
    )
    cell_rho <- safe_cor(cell$consensus_nmf_P1_usage,cell$frozen_top15_weighted_score)
    agg <- aggregate(
      cbind(consensus_nmf_P1_usage,frozen_top15_weighted_score)~
        patient+fine_state,
      cell,mean
    )
    agg_rho <- safe_cor(
      agg$consensus_nmf_P1_usage,agg$frozen_top15_weighted_score
    )
    within <- do.call(rbind,lapply(sort(unique(cell$patient)),function(pt){
      z <- cell[cell$patient==pt,]
      data.frame(
        analysis="within_patient",excluded_patient=pt,n_cells=nrow(z),
        rho_spearman=safe_cor(
          z$consensus_nmf_P1_usage,z$frozen_top15_weighted_score
        )
      )
    }))
    loo <- do.call(rbind,lapply(sort(unique(cell$patient)),function(pt){
      z <- cell[cell$patient!=pt,]
      data.frame(
        analysis="leave_one_patient_out",excluded_patient=pt,n_cells=nrow(z),
        rho_spearman=safe_cor(
          z$consensus_nmf_P1_usage,z$frozen_top15_weighted_score
        )
      )
    }))
    fullrow <- data.frame(
      analysis="full",excluded_patient="none",n_cells=nrow(cell),
      rho_spearman=cell_rho
    )
    fid_sens <- rbind(fullrow,within,loo)
    con <- gzfile(p(DIR$tables,"38_P1_projection_fidelity_cell_level.csv.gz"),"wt")
    write.csv(cell,con,row.names=FALSE); close(con)
    # Export a direct Figure-1 source table: the same discovery myeloid cells,
    # frozen-state labels, P1 activity and the patient-Harmony UMAP whenever
    # that reduction is retained in the annotated object.
    red_name <- if ("umap.patient" %in% Seurat::Reductions(obj)) {
      "umap.patient"
    } else if ("umap" %in% Seurat::Reductions(obj)) {
      "umap"
    } else {
      NA_character_
    }
    if (is.na(red_name)) {
      write.csv(
        data.frame(reduction_available = FALSE, note = "No UMAP reduction retained in fine_myeloid_dc_annotated.rds"),
        p(DIR$audit, "08_discovery_UMAP_export_status.csv"), row.names = FALSE
      )
    } else {
      emb <- as.data.frame(Seurat::Embeddings(obj, red_name)[colnames(obj), 1:2, drop = FALSE])
      names(emb) <- c("UMAP_1", "UMAP_2")
      umap_out <- cbind(
        data.frame(
          barcode = colnames(obj),
          patient = as.character(obj$patient_id),
          sample = as.character(obj$sample_id),
          tissue = as.character(obj$tissue),
          fine_state = as.character(obj$fine_label),
          fixed_state = as.character(obj$fine_status),
          P1_usage = h1,
          P1_top15_score = proj,
          stringsAsFactors = FALSE
        ),
        emb
      )
      con <- gzfile(p(DIR$tables, "36b_GSE178318_myeloid_UMAP_P1.csv.gz"), "wt")
      write.csv(umap_out, con, row.names = FALSE); close(con)
      write.csv(
        data.frame(reduction = red_name, n_cells = nrow(umap_out),
                   file = p(DIR$tables, "36b_GSE178318_myeloid_UMAP_P1.csv.gz")),
        p(DIR$audit, "08_discovery_UMAP_export_status.csv"), row.names = FALSE
      )
    }
    write.csv(agg,p(DIR$tables,"39_P1_projection_fidelity_patient_fine_state.csv"),
              row.names=FALSE)
    write.csv(fid_sens,p(DIR$tables,"40_P1_projection_fidelity_LOPO.csv"),
              row.names=FALSE)
    write.csv(
      data.frame(
        level=c("cell","patient_fine_state","within_patient","LOPO"),
        n_units=c(nrow(cell),nrow(agg),nrow(within),nrow(loo)),
        spearman_rho=c(cell_rho,agg_rho,NA,NA),
        range_min=c(NA,NA,min(within$rho_spearman,na.rm=TRUE),min(loo$rho_spearman,na.rm=TRUE)),
        range_median=c(NA,NA,median(within$rho_spearman,na.rm=TRUE),median(loo$rho_spearman,na.rm=TRUE)),
        range_max=c(NA,NA,max(within$rho_spearman,na.rm=TRUE),max(loo$rho_spearman,na.rm=TRUE))
      ),
      p(DIR$tables,"41_P1_projection_fidelity_summary.csv"),row.names=FALSE
    )
    
    # discovery NMF input gene expression/detection stats：matched-null只用发现队列属性
    xnmf <- must_load_step_cache(
      p(DIR$cache_nmf_matrix,"GSE178318_myeloid_NMF_input.rds"),
      "Stage 5 NMF input"
    )$X
    gene_stats <- data.frame(
      gene=rownames(xnmf),
      mean_log_expression=rowMeans(xnmf),
      detection_frequency=rowMeans(xnmf>0)
    )
    write.csv(
      gene_stats,
      p(DIR$tables,"42_GSE178318_NMF_input_gene_expression_detection_stats.csv"),
      row.names=FALSE
    )
    
    # 30个selected-rank runs中，P1 correspondence与final consensus P1的审计
    run_audit <- list()
    for(r in seq_len(PARAM$nmf_runs_per_rank)){
      fit <- must_load_step_cache(p(
        DIR$cache_nmf_fits,paste0("rank_",k),
        sprintf("rank%d_run%02d.rds",k,r)
      ), "Stage 6 NMF fit")
      g <- intersect(LEGACY_P1_TOP15$gene,rownames(fit$W))
      ref <- LEGACY_P1_TOP15$weight_norm[match(g,LEGACY_P1_TOP15$gene)]
      rr <- apply(fit$W[g,,drop=FALSE],2,function(v)safe_cor(v,ref))
      ff <- which.max(rr)
      finalrho <- safe_cor(fit$W[,ff],W2[,"P1"])
      top_run <- names(sort(fit$W[,ff],decreasing=TRUE))[1:15]
      top_final <- p1_top15$gene
      run_audit[[r]] <- data.frame(
        run=r,matched_factor=ff,legacy_top15_spearman=max(rr,na.rm=TRUE),
        loading_rho_vs_final_consensus_P1=finalrho,
        top15_overlap_vs_final=length(intersect(top_run,top_final))/15
      )
    }
    run_audit <- do.call(rbind,run_audit)
    write.csv(run_audit,p(DIR$tables,"43_selected_rank_P1_30run_audit.csv"),
              row.names=FALSE)
    
    convergence_p1 <- do.call(rbind,lapply(c(10L,20L,30L),function(nn){
      z <- run_audit[seq_len(nn),]
      data.frame(
        first_n_runs=nn,
        median_loading_rho_vs_final=median(z$loading_rho_vs_final_consensus_P1),
        median_top15_overlap_vs_final=median(z$top15_overlap_vs_final),
        median_legacy_correspondence=median(z$legacy_top15_spearman)
      )
    }))
    write.csv(convergence_p1,p(DIR$tables,"44_P1_10_20_30_run_convergence.csv"),
              row.names=FALSE)
    
    # fidelity scatter
    gp <- ggplot2::ggplot(
      cell,ggplot2::aes(frozen_top15_weighted_score,consensus_nmf_P1_usage)
    ) +
      ggplot2::geom_point(size=.35,alpha=.2) +
      ggplot2::labs(
        x="Frozen top-15 weighted P1 score",
        y="Consensus NMF P1 usage",
        title="Discovery-cohort P1 projection fidelity"
      ) + ggplot2::theme_classic()
    ggplot2::ggsave(
      p(DIR$figures,"P1_projection_fidelity.pdf"),gp,width=6,height=5
    )
    
    frozen <- list(
      selected_rank=k,
      P1_original_consensus_factor=p1_factor,
      factor_order=ord,
      matching_source=matching_source,
      p1_anchor_genes=anchor_genes,
      p1_core_anchor_genes=intersect(P1_CORE_ANCHOR_GENES, anchor_genes),
      p1_auxiliary_anchor_genes=intersect(P1_AUXILIARY_ANCHOR_GENES, anchor_genes),
      p1_anchor_score=best_anchor,
      p1_anchor_runner_up_score=runner_up_anchor,
      p1_anchor_margin=anchor_margin,
      p1_anchor_min_margin=PARAM$p1_anchor_min_margin,
      p1_selection_status=p1_selection_status,
      p1_selection_low_confidence=p1_selection_low_confidence,
      p1_selection_rule="highest mean within-factor z-scored loading of prespecified SPP1/APOC1 core plus APOE/GPNMB/MARCO/LPL/CST3/LGALS3 auxiliary module among stable factors; tie -> higher loading stability",
      W=W2,H=H2,
      all_program_top50=top50,
      all_program_top15=all_top15,
      P1_top15=p1_top15,
      legacy_correspondence=audit,
      P1_30run_audit=run_audit,
      cell_projection_fidelity=cell_rho,
      patient_fine_state_projection_fidelity=agg_rho,
      consensus_method=cons$method,
      pipeline_version=PIPELINE_VERSION
    )
    save_step_cache(frozen,out)
    # 保存新版冻结对象；后续阶段统一读取 frozen_P1_v4.rds
    save_step_cache(frozen,p(DIR$cache_projection,"frozen_P1_v4.rds"))
    jsonlite::write_json(
      list(
        pipeline_version = PIPELINE_VERSION,
        selected_rank = k,
        nmf_ranks = PARAM$nmf_ranks,
        runs_per_rank = PARAM$nmf_runs_per_rank,
        p1_anchor_genes = anchor_genes,
        p1_core_anchor_genes = intersect(P1_CORE_ANCHOR_GENES, anchor_genes),
        p1_auxiliary_anchor_genes = intersect(P1_AUXILIARY_ANCHOR_GENES, anchor_genes),
        p1_anchor_score = best_anchor,
        p1_anchor_runner_up_score = runner_up_anchor,
        p1_anchor_margin = anchor_margin,
        p1_anchor_min_margin = PARAM$p1_anchor_min_margin,
        p1_selection_status = p1_selection_status,
        p1_selection_low_confidence = p1_selection_low_confidence,
        p1_selection_rule = frozen$p1_selection_rule,
        marker_excluded_genes_provisional = P1_MARKER_EXCLUDED_GENES,
        marker_exclusion_confirmed = P1_MARKER_EXCLUSION_CONFIRMED,
        frozen_top15_file = p(DIR$tables, "34_frozen_P1_top15_weights.csv"),
        frozen_rds = out,
        external_data_used_for_selection = FALSE,
        spatial_data_used_for_selection = FALSE,
        patient_level_results_used_for_selection = FALSE
      ),
      p(DIR$audit, "P1_freeze_manifest.json"),
      auto_unbox = TRUE, pretty = TRUE
    )
    
    write_session_info("08_post_rank_P1_correspondence_and_freeze")
    cat("Selected rank:",k,"\n")
    cat("P1 consensus factor before reorder:",p1_factor,"\n")
    cat("Cell projection fidelity rho:",round(cell_rho,4),"\n")
    frozen
  })
}

############################
# Stage 9：GSE245552 冻结投影（不重新NMF/聚类）
############################

parse_gse245552_soft <- function(soft_file) {
  assert_file(soft_file)
  con <- gzfile(soft_file, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  out <- list()
  cur <- NULL
  tissue <- NULL
  flush_one <- function() {
    if (!is.null(cur) && !is.null(tissue)) out[[cur]] <<- tissue
  }
  for (line in lines) {
    if (startsWith(line, "!Sample_title")) {
      flush_one()
      cur <- trimws(strsplit(line, "=", fixed = TRUE)[[1]][2])
      cur <- trimws(strsplit(cur, ",", fixed = TRUE)[[1]][1])
      tissue <- NULL
    } else if (startsWith(line, "!Sample_characteristics_ch1") && !is.null(cur)) {
      val <- trimws(strsplit(line, "=", fixed = TRUE)[[1]][2])
      if (grepl("^tissue:", val, ignore.case = TRUE)) {
        tissue <- gsub(" ", "_", trimws(sub("^[^:]+:", "", val)))
      }
    }
  }
  flush_one()
  unlist(out)
}

stage9_gse245552_projection <- function(force = FALSE) {
  stage_log("09_GSE245552_frozen_projection_all_programs", {
    frozen_file <- p(DIR$cache_projection,"frozen_P1_v4.rds")
    assert_file(frozen_file)
    frozen <- must_load_step_cache(frozen_file, "frozen P1")
    
    raw_dir <- p(GSE245552_DIR,"raw")
    soft_file <- p(GSE245552_DIR,"GSE245552_family.soft.gz")
    if(!dir.exists(raw_dir) || !file.exists(soft_file)){
      stop(
        "GSE245552尚未就绪。\n期望：\n",raw_dir,"\n",soft_file,
        "\n下载完成后重新运行Stage 9。"
      )
    }
    
    out <- p(DIR$cache_projection,"GSE245552_frozen_projection_all.csv")
    if(!force && file.exists(out)){
      message("[09_GSE245552_frozen_projection_all_programs] CACHE HIT：读取冻结投影结果。")
      return(read.csv(out,stringsAsFactors=FALSE))
    }
    
    p1 <- frozen$P1_top15
    p1w <- data.frame(gene=p1$gene,weight=p1$weight_norm)
    p1_weight_map <- setNames(p1w$weight,p1w$gene)
    
    # Full P1 and the single-gene sensitivity are available immediately.
    # The multi-gene marker-excluded set is added only after explicit review.
    p1_sets <- make_p1_score_sets(p1w$gene)
    
    all15 <- frozen$all_program_top15
    program_weights <- lapply(unique(all15$program),function(pp){
      z <- all15[all15$program==pp & all15$rank<=15,]
      data.frame(gene=z$gene,weight=z$weight/sum(z$weight))
    })
    names(program_weights) <- unique(all15$program)
    
    tissue_map <- parse_gse245552_soft(soft_file)
    allowed <- c(
      "primary_tumor","liver_metastasis",
      "colon_adjacent_tissue","liver_adjacent_tissue"
    )
    
    state_score_fun <- function(norm,modules){
      sm <- sapply(modules,function(gs){
        g <- intersect(gs,rownames(norm))
        if(!length(g)) return(rep(0,ncol(norm)))
        ww <- vapply(g,function(x){
          if(x %in% names(p1_weight_map)) p1_weight_map[[x]] else 1
        },numeric(1))
        ww <- ww/sum(ww)
        as.numeric(Matrix::t(norm[g,,drop=FALSE])%*%ww)
      })
      if(!is.matrix(sm)) sm <- matrix(sm,ncol=length(modules))
      colnames(sm) <- names(modules)
      sm
    }
    
    mtxs <- sort(list.files(raw_dir,pattern="_matrix\\.mtx\\.gz$",full.names=TRUE))
    if(!length(mtxs)) stop("GSE245552 raw目录没有 *_matrix.mtx.gz。")
    
    rows <- list(); rr <- 0L
    coverage <- list(); cv <- 0L
    for(mf in mtxs){
      stem <- sub("_matrix\\.mtx\\.gz$","",basename(mf))
      bf <- p(raw_dir,paste0(stem,"_barcodes.tsv.gz"))
      ff <- p(raw_dir,paste0(stem,"_features.tsv.gz"))
      X <- read_10x_triplet(mf,bf,ff)
      norm <- normalize_log1p_cp10k_sparse(X)
      
      # Frozen myeloid inclusion
      my <- module_score(norm,MYELOID_MODULE)
      ep <- module_score(norm,EPITHELIAL_MODULE)
      lyz <- if("LYZ"%in%rownames(X))as.numeric(X["LYZ",]>0)>0 else rep(FALSE,ncol(X))
      inc <- my>ep & lyz
      
      # 两套state map：主冻结模块；P1 top15完全移除后的overlap sensitivity
      modules_main <- STATE_MODULES
      modules_no_overlap <- lapply(STATE_MODULES,setdiff,y=p1w$gene)
      sm1 <- state_score_fun(norm,modules_main)
      sm2 <- state_score_fun(norm,modules_no_overlap)
      lab1 <- colnames(sm1)[max.col(sm1,ties.method="first")]
      lab2 <- colnames(sm2)[max.col(sm2,ties.method="first")]
      lab1[apply(sm1,1,max)<=0] <- "Unassigned"
      lab2[apply(sm2,1,max)<=0] <- "Unassigned"
      
      m <- regexpr("(Pt[0-9]+_[TZ][0-9]+)$",stem,perl=TRUE)
      sample_id <- if(m[1]>0)regmatches(stem,m) else stem
      patient <- sub("_.*$","",sample_id)
      tissue <- unname(tissue_map[sample_id])
      if(is.na(tissue)||!nzchar(tissue)) stop("缺少GSE245552 tissue metadata：",sample_id)
      if(!tissue%in%allowed) stop("未预期GSE245552 tissue：",tissue)
      
      # score vectors一次算完
      p1_scores <- lapply(p1_sets,function(gs){
        ww <- p1w[p1w$gene%in%gs,,drop=FALSE]
        weighted_score(norm,ww)
      })
      prog_scores <- lapply(program_weights,function(ww)weighted_score(norm,ww))
      
      for(mapver in c("original_frozen_modules","P1_top15_excluded_modules")){
        labels <- if(mapver=="original_frozen_modules")lab1 else lab2
        for(st in names(STATE_MODULES)){
          mask <- inc & labels==st
          if(!any(mask)) next
          # P1主/敏感性
          for(ver in names(p1_scores)){
            rr <- rr+1L
            sv <- p1_scores[[ver]]
            rows[[rr]] <- data.frame(
              sample=sample_id,patient=patient,tissue=tissue,state=st,
              state_mapping_version=mapver,
              projection_normalization=PARAM$projection_normalization,
              score_version=ver,score_family="P1_sensitivity",
              n_cells=sum(mask),
              mean_score=mean(sv[mask]),median_score=median(sv[mask])
            )
          }
          # P1-Pk program benchmark
          for(pp in names(prog_scores)){
            rr <- rr+1L
            sv <- prog_scores[[pp]]
            rows[[rr]] <- data.frame(
              sample=sample_id,patient=patient,tissue=tissue,state=st,
              state_mapping_version=mapver,
              projection_normalization=PARAM$projection_normalization,
              score_version=pp,score_family="all_programs",
              n_cells=sum(mask),
              mean_score=mean(sv[mask]),median_score=median(sv[mask])
            )
          }
        }
      }
      
      # gene coverage：按sample记录
      for(nm in names(p1_sets)){
        cv <- cv+1L
        gs <- p1_sets[[nm]]
        coverage[[cv]] <- data.frame(
          sample=sample_id,projection_normalization=PARAM$projection_normalization,
          score_version=nm,
          genes_expected=length(gs),
          genes_present=sum(gs%in%rownames(norm)),
          coverage_fraction=mean(gs%in%rownames(norm))
        )
      }
      for(pp in names(program_weights)){
        cv <- cv+1L
        gs <- program_weights[[pp]]$gene
        coverage[[cv]] <- data.frame(
          sample=sample_id,projection_normalization=PARAM$projection_normalization,
          score_version=pp,
          genes_expected=length(gs),
          genes_present=sum(gs%in%rownames(norm)),
          coverage_fraction=mean(gs%in%rownames(norm))
        )
      }
      rm(X,norm); invisible(gc())
    }
    
    ans <- do.call(rbind,rows)
    cov <- do.call(rbind,coverage)
    if(any(cov$coverage_fraction<0.80)){
      warning("GSE245552至少一个sample×score的gene coverage <80%；见coverage表。")
    }
    if (!isTRUE(P1_MARKER_EXCLUSION_CONFIRMED)) {
      write.csv(
        data.frame(
          marker_exclusion_confirmed = FALSE,
          note = "P1_full and P1_no_SPP1 were projected; the multi-gene marker-excluded score is pending post-Stage-8 review.",
          stringsAsFactors = FALSE
        ),
        p(DIR$audit, "09_marker_exclusion_pending.csv"), row.names = FALSE
      )
    }
    write.csv(ans,out,row.names=FALSE)
    write.csv(cov,p(DIR$tables,"45_GSE245552_projection_gene_coverage.csv"),
              row.names=FALSE)
    
    # 兼容文章旧命名：只输出主state map的P1三版本
    legacy_like <- ans[
      ans$state_mapping_version=="original_frozen_modules" &
        ans$score_family=="P1_sensitivity",,
      drop=FALSE
    ]
    write.csv(
      legacy_like,
      p(DIR$tables,"46_GSE245552_frozen_replication_by_state.csv"),
      row.names=FALSE
    )
    allprog <- ans[
      ans$state_mapping_version=="original_frozen_modules" &
        ans$score_family=="all_programs",,
      drop=FALSE
    ]
    write.csv(
      allprog,
      p(DIR$tables,"47_GSE245552_all_programs_by_state.csv"),
      row.names=FALSE
    )
    
    write_session_info("09_GSE245552_frozen_projection_all_programs")
    ans
  })
}

############################
# Stage 10：GSE245552 患者级统计 + composition/activity
############################

stage10_gse245552_statistics <- function(force = FALSE) {
  stage_log("10_GSE245552_full_patient_level_statistics", {
    input <- p(DIR$cache_projection,"GSE245552_frozen_projection_all.csv")
    assert_file(input)
    raw <- read.csv(input,stringsAsFactors=FALSE)
    
    PRIMARY_STATES <- c("Mono","P3_C1Q_MHCII","Other","P2_inflammatory")
    ALL_STATES <- c(PRIMARY_STATES,"P1_state")
    PRIMARY_SCORES <- intersect(
      c("P1_full", "P1_marker_excluded"),
      unique(raw$score_version)
    )
    if (!"P1_full" %in% PRIMARY_SCORES) {
      stop("GSE245552投影结果缺少P1_full；请先完成Stage 9。")
    }
    TISSUES <- c("primary_tumor","liver_metastasis")
    
    # 聚合：重复sample row在patient x tissue x state内先等权平均；cell数求和。
    aggregate_external <- function(d){
      key <- c("score_version","patient","tissue","state")
      m <- aggregate(mean_score~score_version+patient+tissue+state,d,mean,na.rm=TRUE)
      med <- aggregate(median_score~score_version+patient+tissue+state,d,mean,na.rm=TRUE)
      nc <- aggregate(n_cells~score_version+patient+tissue+state,d,sum,na.rm=TRUE)
      out <- Reduce(function(x,y)merge(x,y,by=key,all=TRUE),list(m,med,nc))
      out
    }
    
    primary_raw <- raw[
      raw$state_mapping_version=="original_frozen_modules" &
        raw$score_version%in%PRIMARY_SCORES &
        raw$tissue%in%TISSUES,,
      drop=FALSE
    ]
    agg <- aggregate_external(primary_raw)
    write.csv(agg,p(DIR$tables,"48_GSE245552_patient_tissue_state_aggregated.csv"),
              row.names=FALSE)
    
    # -------- A. exploratory state-specific: two-sided paired Wilcoxon + BH --------
    # 预设的 Wilcoxon 规则：能 exact 则 exact，否则使用渐近法。
    state_rows <- list(); sr <- 0L
    for(ver in PRIMARY_SCORES){
      for(st in ALL_STATES){
        z <- agg[agg$score_version==ver & agg$state==st,]
        piv <- reshape(
          z[,c("patient","tissue","mean_score")],
          idvar="patient",timevar="tissue",direction="wide"
        )
        need <- paste0("mean_score.",TISSUES)
        if(!all(need%in%names(piv))) next
        piv <- piv[complete.cases(piv[,need]),]
        delta <- piv[[paste0("mean_score.liver_metastasis")]]-
          piv[[paste0("mean_score.primary_tumor")]]
        if(!length(delta)) next
        set.seed(20260719L+sr)
        boots <- replicate(PARAM$n_boot,median(sample(delta,length(delta),replace=TRUE)))
        ci <- quantile(boots,c(.025,.975),names=FALSE)
        wt2 <- wilcox_signed_rank_auto(delta,"two.sided")
        wtg <- wilcox_signed_rank_auto(delta,"greater")
        sr <- sr+1L
        state_rows[[sr]] <- data.frame(
          score_version=ver,state=st,n_paired_patients=length(delta),
          n_positive=sum(delta>0),mean_delta=mean(delta),median_delta=median(delta),
          median_ci_low=ci[1],median_ci_high=ci[2],
          wilcoxon_stat_two_sided=wt2$statistic,
          wilcoxon_p_two_sided=wt2$p.value,
          wilcoxon_exact_used_two_sided=wt2$exact_used,
          wilcoxon_stat_greater=wtg$statistic,
          wilcoxon_p_greater=wtg$p.value,
          wilcoxon_exact_used_greater=wtg$exact_used
        )
      }
    }
    state_stats <- do.call(rbind,state_rows)
    state_stats$BH_q_primary_two_sided <- NA_real_
    state_stats$BH_q_primary_greater <- NA_real_
    for(ver in PRIMARY_SCORES){
      ii <- state_stats$score_version==ver & state_stats$state%in%PRIMARY_STATES
      state_stats$BH_q_primary_two_sided[ii] <- p.adjust(
        state_stats$wilcoxon_p_two_sided[ii],"BH"
      )
      state_stats$BH_q_primary_greater[ii] <- p.adjust(
        state_stats$wilcoxon_p_greater[ii],"BH"
      )
    }
    write.csv(state_stats,p(DIR$tables,"49_GSE245552_state_specific_primary.csv"),
              row.names=FALSE)
    
    # -------- B. 每端≥1/≥5/≥10 cells state-specific sensitivity --------
    threshold_rows <- list(); tr <- 0L
    for(th in c(1L,5L,10L)){
      tmp_rows <- list(); tt <- 0L
      for(ver in PRIMARY_SCORES){
        for(st in PRIMARY_STATES){
          z <- agg[agg$score_version==ver & agg$state==st & agg$n_cells>=th,]
          piv <- reshape(
            z[,c("patient","tissue","mean_score","n_cells")],
            idvar="patient",timevar="tissue",direction="wide"
          )
          if(!all(paste0("mean_score.",TISSUES)%in%names(piv))) next
          piv <- piv[complete.cases(piv[,paste0("mean_score.",TISSUES)]),]
          delta <- piv$mean_score.liver_metastasis-piv$mean_score.primary_tumor
          if(!length(delta)) next
          wt <- wilcox_signed_rank_auto(delta,"two.sided")
          tt <- tt+1L
          tmp_rows[[tt]] <- data.frame(
            threshold=th,score_version=ver,state=st,n_patients=length(delta),
            n_positive=sum(delta>0),median_delta=median(delta),
            wilcoxon_p_two_sided=wt$p.value,exact_used=wt$exact_used
          )
        }
      }
      tmp <- if(length(tmp_rows))do.call(rbind,tmp_rows)else data.frame()
      if(nrow(tmp)){
        for(ver in PRIMARY_SCORES){
          ii <- tmp$score_version==ver
          tmp$BH_q <- NA_real_
          tmp$BH_q[ii] <- p.adjust(tmp$wilcoxon_p_two_sided[ii],"BH")
        }
        tr <- tr+1L; threshold_rows[[tr]] <- tmp
      }
    }
    threshold_state <- do.call(rbind,threshold_rows)
    write.csv(
      threshold_state,
      p(DIR$tables,"50_GSE245552_state_specific_cell_threshold_sensitivity.csv"),
      row.names=FALSE
    )
    
    # -------- C. 四状态对称composition/activity --------
    decomp_one <- function(d,ver,patient,threshold=1L){
      z <- d[d$score_version==ver & d$patient==patient &
               d$state%in%PRIMARY_STATES,]
      # 所有四状态必须在两端达到threshold
      for(st in PRIMARY_STATES){
        w <- z[z$state==st & z$n_cells>=threshold,]
        if(!all(TISSUES%in%w$tissue)) return(NULL)
        if(any(!is.finite(w$mean_score[match(TISSUES,w$tissue)]))) return(NULL)
      }
      muPT <- vapply(PRIMARY_STATES,function(st){
        z$mean_score[z$state==st & z$tissue=="primary_tumor"][1]
      },numeric(1))
      muLM <- vapply(PRIMARY_STATES,function(st){
        z$mean_score[z$state==st & z$tissue=="liver_metastasis"][1]
      },numeric(1))
      nPT <- vapply(PRIMARY_STATES,function(st){
        z$n_cells[z$state==st & z$tissue=="primary_tumor"][1]
      },numeric(1))
      nLM <- vapply(PRIMARY_STATES,function(st){
        z$n_cells[z$state==st & z$tissue=="liver_metastasis"][1]
      },numeric(1))
      pPT <- nPT/sum(nPT); pLM <- nLM/sum(nLM)
      total <- sum(pLM*muLM)-sum(pPT*muPT)
      comp <- sum((pLM-pPT)*(muLM+muPT)/2)
      activity <- sum((muLM-muPT)*(pLM+pPT)/2)
      err <- total-comp-activity
      if(abs(err)>1e-10)stop("四状态分解closure失败。")
      data.frame(
        score_version=ver,patient=patient,threshold=threshold,
        overall_delta=total,state_composition=comp,
        within_state_activity=activity,
        activity_minus_composition=activity-comp,
        closure_error=err
      )
    }
    
    drows <- list(); dr <- 0L
    for(ver in PRIMARY_SCORES){
      for(th in c(1L,5L,10L)){
        for(pt in sort(unique(agg$patient))){
          z <- decomp_one(agg,ver,pt,th)
          if(!is.null(z)){dr<-dr+1L;drows[[dr]]<-z}
        }
      }
    }
    decomp <- do.call(rbind,drows)
    write.csv(decomp,p(DIR$tables,"51_GSE245552_composition_activity_by_patient.csv"),
              row.names=FALSE)
    
    ds_rows <- list(); ds <- 0L
    for(ver in PRIMARY_SCORES){
      for(th in c(1L,5L,10L)){
        g <- decomp[decomp$score_version==ver & decomp$threshold==th,]
        if(!nrow(g))next
        for(metric in c(
          "overall_delta","state_composition","within_state_activity",
          "activity_minus_composition"
        )){
          v <- g[[metric]]
          ci <- bootstrap_median_ci(
            v,PARAM$n_boot,seed=20260815L+ds
          )
          ds<-ds+1L
          # LOPO：仅报告中位值范围，不新增显著性检验
          loo_med <- if(length(v)>1L)vapply(seq_along(v),function(i)median(v[-i]),numeric(1)) else NA_real_
          ds_rows[[ds]] <- data.frame(
            score_version=ver,threshold=th,component=metric,
            n_patients=length(v),median=median(v),
            ci_low=ci[1],ci_high=ci[2],
            n_positive=sum(v>0),
            exact_sign_p_greater=exact_sign_p_greater(v),
            LOPO_median_min=if(all(is.na(loo_med)))NA_real_ else min(loo_med,na.rm=TRUE),
            LOPO_median_max=if(all(is.na(loo_med)))NA_real_ else max(loo_med,na.rm=TRUE)
          )
        }
      }
    }
    decomp_summary <- do.call(rbind,ds_rows)
    write.csv(
      decomp_summary,
      p(DIR$tables,"52_GSE245552_composition_activity_summary_threshold_LOPO.csv"),
      row.names=FALSE
    )
    
    # -------- D. 五状态 P1-state-inclusive；结构零显式保留 --------
    five_rows <- list(); fr <- 0L
    for(ver in PRIMARY_SCORES){
      base_patients <- unique(
        decomp$patient[decomp$score_version==ver & decomp$threshold==1L]
      )
      for(pt in sort(base_patients)){
        z <- agg[agg$score_version==ver & agg$patient==pt,]
        # four non-P1
        muPT4 <- vapply(PRIMARY_STATES,function(st)
          z$mean_score[z$state==st & z$tissue=="primary_tumor"][1],numeric(1))
        muLM4 <- vapply(PRIMARY_STATES,function(st)
          z$mean_score[z$state==st & z$tissue=="liver_metastasis"][1],numeric(1))
        nPT4 <- vapply(PRIMARY_STATES,function(st)
          z$n_cells[z$state==st & z$tissue=="primary_tumor"][1],numeric(1))
        nLM4 <- vapply(PRIMARY_STATES,function(st)
          z$n_cells[z$state==st & z$tissue=="liver_metastasis"][1],numeric(1))
        
        getone <- function(tissue,col){
          q <- z[z$state=="P1_state" & z$tissue==tissue,col]
          if(!length(q))return(NA_real_)
          q[1]
        }
        nPPT <- getone("primary_tumor","n_cells"); if(!is.finite(nPPT))nPPT<-0
        nPLM <- getone("liver_metastasis","n_cells"); if(!is.finite(nPLM))nPLM<-0
        muPPT <- getone("primary_tumor","mean_score")
        muPLM <- getone("liver_metastasis","mean_score")
        totNPT <- sum(nPT4)+nPPT; totNLM <- sum(nLM4)+nPLM
        pPT4 <- nPT4/totNPT; pLM4 <- nLM4/totNLM
        pPPT <- nPPT/totNPT; pPLM <- nPLM/totNLM
        comp <- sum((pLM4-pPT4)*(muLM4+muPT4)/2)
        act <- sum((muLM4-muPT4)*(pLM4+pPT4)/2)
        termPT <- if(nPPT>0 && is.finite(muPPT))pPPT*muPPT else 0
        termLM <- if(nPLM>0 && is.finite(muPLM))pPLM*muPLM else 0
        p1con <- termLM-termPT
        overall <- sum(pLM4*muLM4)+termLM - (sum(pPT4*muPT4)+termPT)
        err <- overall-comp-act-p1con
        if(abs(err)>1e-10)stop("五状态closure失败。")
        presence <- if(nPPT>0 && nPLM>0)"both_endpoints_positive" else
          if(nPPT>0)"primary_tumour_only" else
            if(nPLM>0)"liver_metastasis_only" else "neither_endpoint"
        fr<-fr+1L
        five_rows[[fr]] <- data.frame(
          score_version=ver,patient=pt,overall_delta=overall,
          non_P1_composition=comp,non_P1_activity=act,
          P1_state_contribution=p1con,
          activity_minus_P1_state=act-p1con,
          P1_state_PT_cells=nPPT,P1_state_LM_cells=nPLM,
          P1_state_presence=presence,closure_error=err
        )
      }
    }
    five <- do.call(rbind,five_rows)
    write.csv(five,p(DIR$tables,"53_GSE245552_five_state_by_patient.csv"),
              row.names=FALSE)
    
    fsum <- list(); fs <- 0L
    for(ver in PRIMARY_SCORES){
      for(subset_name in c("all_complete_four_state","both_endpoints_positive")){
        g <- five[five$score_version==ver,]
        if(subset_name=="both_endpoints_positive"){
          g <- g[g$P1_state_presence=="both_endpoints_positive",]
        }
        if(!nrow(g))next
        for(metric in c(
          "overall_delta","non_P1_composition","non_P1_activity",
          "P1_state_contribution","activity_minus_P1_state"
        )){
          v <- g[[metric]]
          ci <- bootstrap_median_ci(v,PARAM$n_boot,seed=20260826L+fs)
          fs<-fs+1L
          fsum[[fs]] <- data.frame(
            score_version=ver,subset=subset_name,component=metric,
            n_patients=length(v),median=median(v),ci_low=ci[1],ci_high=ci[2],
            n_positive=sum(v>0),exact_sign_p_greater=exact_sign_p_greater(v)
          )
        }
      }
    }
    five_summary <- do.call(rbind,fsum)
    write.csv(
      five_summary,
      p(DIR$tables,"54_GSE245552_five_state_presence_sensitivity_summary.csv"),
      row.names=FALSE
    )
    
    # -------- E. unified paired model：patient-cluster bootstrap + shared sign flip --------
    unified_input <- list(); ui <- 0L
    for(ver in PRIMARY_SCORES){
      z <- agg[agg$score_version==ver & agg$state%in%PRIMARY_STATES,]
      for(pt in unique(z$patient)){
        for(st in PRIMARY_STATES){
          w <- z[z$patient==pt & z$state==st,]
          if(!all(TISSUES%in%w$tissue))next
          a <- w$mean_score[w$tissue=="primary_tumor"][1]
          b <- w$mean_score[w$tissue=="liver_metastasis"][1]
          if(!is.finite(a)||!is.finite(b))next
          ui<-ui+1L
          unified_input[[ui]] <- data.frame(
            score_version=ver,patient=pt,state=st,
            primary_tumor=a,liver_metastasis=b,delta=b-a
          )
        }
      }
    }
    unified_input <- do.call(rbind,unified_input)
    write.csv(unified_input,p(DIR$tables,"55_GSE245552_unified_model_input.csv"),
              row.names=FALSE)
    
    unified_summary <- list(); us <- 0L
    unified_state <- list(); uc <- 0L
    for(ver in PRIMARY_SCORES){
      z <- unified_input[unified_input$score_version==ver,]
      pats <- sort(unique(z$patient))
      mat <- matrix(
        NA_real_,nrow=length(pats),ncol=length(PRIMARY_STATES),
        dimnames=list(pats,PRIMARY_STATES)
      )
      for(i in seq_len(nrow(z)))mat[z$patient[i],z$state[i]]<-z$delta[i]
      direct <- colMeans(mat,na.rm=TRUE)
      
      set.seed(20260813L+match(ver,PRIMARY_SCORES))
      B <- PARAM$n_boot
      boots <- matrix(NA_real_,B,ncol(mat))
      for(b in seq_len(B)){
        ix <- sample(seq_len(nrow(mat)),nrow(mat),replace=TRUE)
        boots[b,] <- colMeans(mat[ix,,drop=FALSE],na.rm=TRUE)
      }
      P <- PARAM$n_signflip
      perms <- matrix(NA_real_,P,ncol(mat))
      for(b in seq_len(P)){
        sg <- sample(c(-1,1),nrow(mat),replace=TRUE)
        perms[b,] <- colMeans(mat*sg,na.rm=TRUE)
      }
      two_p <- function(null,obs){
        # The permutation denominator is the number of finite null draws;
        # failed/NA draws must not dilute or inflate the empirical P value.
        n_finite <- sum(is.finite(null))
        if (n_finite < 1L) return(NA_real_)
        up <- plus_one_p(sum(null>=obs,na.rm=TRUE), n_finite)
        lo <- plus_one_p(sum(null<=obs,na.rm=TRUE), n_finite)
        min(1,2*min(up,lo))
      }
      for(j in seq_along(PRIMARY_STATES)){
        uc<-uc+1L
        unified_state[[uc]] <- data.frame(
          score_version=ver,state=PRIMARY_STATES[j],
          estimate=direct[j],
          bootstrap_ci_low=quantile(boots[,j],.025,na.rm=TRUE),
          bootstrap_ci_high=quantile(boots[,j],.975,na.rm=TRUE),
          sign_flip_p_two_sided=two_p(perms[,j],direct[j]),
          n_patient_state_pairs=sum(is.finite(mat[,j])),
          n_patients=length(pats)
        )
      }
      marg <- mean(direct)
      bm <- rowMeans(boots)
      pm <- rowMeans(perms)
      qobs <- mean((direct-marg)^2)
      qperm <- rowMeans((perms-rowMeans(perms))^2)
      us<-us+1L
      unified_summary[[us]] <- data.frame(
        score_version=ver,
        estimand="equal-weight marginal LM_minus_PT across four prespecified states",
        estimate=marg,
        bootstrap_ci_low=quantile(bm,.025),
        bootstrap_ci_high=quantile(bm,.975),
        sign_flip_p_two_sided=two_p(pm,marg),
        state_heterogeneity_Q=qobs,
        state_heterogeneity_permutation_p=
          plus_one_p(sum(qperm>=qobs),length(qperm)),
        n_patient_state_pairs=sum(is.finite(mat)),
        n_patients=length(pats),seed=20260813L,
        n_bootstrap=B,n_permutations=P
      )
    }
    unified_state <- do.call(rbind,unified_state)
    unified_summary <- do.call(rbind,unified_summary)
    unified_state$BH_q_within_score_version <- NA_real_
    for(ver in PRIMARY_SCORES){
      ii <- unified_state$score_version==ver
      unified_state$BH_q_within_score_version[ii] <-
        p.adjust(unified_state$sign_flip_p_two_sided[ii],"BH")
    }
    write.csv(unified_state,p(DIR$tables,"56_GSE245552_unified_state_effects.csv"),
              row.names=FALSE)
    write.csv(unified_summary,p(DIR$tables,"57_GSE245552_unified_summary.csv"),
              row.names=FALSE)
    
    # -------- F. P1-Pk external specificity benchmark --------
    prog_raw <- raw[
      raw$state_mapping_version=="original_frozen_modules" &
        raw$score_family=="all_programs" &
        raw$tissue%in%TISSUES & raw$state%in%PRIMARY_STATES,,
      drop=FALSE
    ]
    prog_agg <- aggregate_external(prog_raw)
    prog_rows <- list(); pr <- 0L
    prog_delta_rows <- list(); pdr <- 0L
    for(pp in sort(unique(prog_agg$score_version))){
      z <- prog_agg[prog_agg$score_version==pp,]
      pvals <- c()
      state_pairs <- c()
      for(pt in sort(unique(z$patient))){
        vals <- c()
        for(st in PRIMARY_STATES){
          w <- z[z$patient==pt & z$state==st,]
          if(all(TISSUES%in%w$tissue)){
            a <- w$mean_score[w$tissue=="primary_tumor"][1]
            b <- w$mean_score[w$tissue=="liver_metastasis"][1]
            if(is.finite(a)&&is.finite(b)){
              vals <- c(vals,b-a)
              pdr<-pdr+1L
              prog_delta_rows[[pdr]] <- data.frame(
                program=pp,patient=pt,state=st,delta_LM_minus_PT=b-a
              )
            }
          }
        }
        if(length(vals))pvals[pt] <- mean(vals)
      }
      if(length(pvals)){
        ci <- bootstrap_median_ci(pvals,PARAM$n_boot,seed=20260815L+pr)
        pr<-pr+1L
        dg <- do.call(rbind,prog_delta_rows)
        zdel <- dg[dg$program==pp,]
        prog_rows[[pr]] <- data.frame(
          program=pp,n_patients=length(pvals),
          median_effect=median(pvals),ci_low=ci[1],ci_high=ci[2],
          n_positive_patients=sum(pvals>0),
          exact_sign_p_greater=exact_sign_p_greater(pvals),
          positive_state_pair_fraction=mean(zdel$delta_LM_minus_PT>0)
        )
      }
    }
    prog_summary <- do.call(rbind,prog_rows)
    prog_deltas <- do.call(rbind,prog_delta_rows)
    write.csv(prog_summary,p(DIR$tables,"58_GSE245552_P1_Pk_specificity_summary.csv"),
              row.names=FALSE)
    write.csv(prog_deltas,p(DIR$tables,"59_GSE245552_P1_Pk_state_pair_deltas.csv"),
              row.names=FALSE)
    
    # -------- G. P1-gene/state-module overlap sensitivity --------
    ovraw <- raw[
      raw$state_mapping_version=="P1_top15_excluded_modules" &
        raw$score_version%in%PRIMARY_SCORES &
        raw$tissue%in%TISSUES,,
      drop=FALSE
    ]
    ovagg <- aggregate_external(ovraw)
    ovrows <- list(); ov <- 0L
    for(ver in PRIMARY_SCORES){
      for(st in PRIMARY_STATES){
        z <- ovagg[ovagg$score_version==ver & ovagg$state==st,]
        piv <- reshape(z[,c("patient","tissue","mean_score")],
                       idvar="patient",timevar="tissue",direction="wide")
        if(!all(paste0("mean_score.",TISSUES)%in%names(piv)))next
        piv <- piv[complete.cases(piv[,paste0("mean_score.",TISSUES)]),]
        dlt <- piv$mean_score.liver_metastasis-piv$mean_score.primary_tumor
        if(!length(dlt))next
        wt <- wilcox_signed_rank_auto(dlt,"two.sided")
        ov<-ov+1L
        ovrows[[ov]] <- data.frame(
          score_version=ver,state=st,n_patients=length(dlt),
          median_delta=median(dlt),n_positive=sum(dlt>0),
          wilcoxon_p_two_sided=wt$p.value
        )
      }
    }
    overlap_summary <- if(length(ovrows))do.call(rbind,ovrows)else data.frame()
    if(nrow(overlap_summary)){
      overlap_summary$BH_q <- NA_real_
      for(ver in PRIMARY_SCORES){
        ii<-overlap_summary$score_version==ver
        overlap_summary$BH_q[ii] <-
          p.adjust(overlap_summary$wilcoxon_p_two_sided[ii],"BH")
      }
    }
    write.csv(
      overlap_summary,
      p(DIR$tables,"60_GSE245552_P1_state_module_overlap_sensitivity.csv"),
      row.names=FALSE
    )
    
    write_session_info("10_GSE245552_full_patient_level_statistics")
    list(
      state_specific=state_stats,
      threshold_state=threshold_state,
      decomposition=decomp,
      decomposition_summary=decomp_summary,
      five_state=five,
      five_state_summary=five_summary,
      unified_state=unified_state,
      unified_summary=unified_summary,
      program_specificity=prog_summary,
      overlap_sensitivity=overlap_summary
    )
  })
}

############################
# Stage 11：GSE231559 冻结投影（方向可迁移性支持）
############################

parse_gse231559_soft <- function(soft_file) {
  assert_file(soft_file)
  con <- gzfile(soft_file,"rt"); on.exit(close(con))
  lines <- readLines(con,warn=FALSE)
  rows <- list(); rr <- 0L
  cur <- NULL; gsm <- NULL; tissue <- NULL; source <- NULL
  flush <- function(){
    if(!is.null(cur)){
      rr <<- rr+1L
      rows[[rr]] <<- data.frame(
        sample_id=cur,gsm=if(is.null(gsm))NA_character_ else gsm,
        source_name=if(is.null(source))NA_character_ else source,
        tissue=if(is.null(tissue))NA_character_ else tissue,
        stringsAsFactors=FALSE
      )
    }
  }
  for(line in lines){
    if(startsWith(line,"^SAMPLE")){
      flush(); cur<-NULL;gsm<-NULL;tissue<-NULL;source<-NULL
    } else if(startsWith(line,"!Sample_title")){
      cur <- trimws(strsplit(line,"=",fixed=TRUE)[[1]][2])
    } else if(startsWith(line,"!Sample_geo_accession") && !is.null(cur)){
      gsm <- trimws(strsplit(line,"=",fixed=TRUE)[[1]][2])
    } else if(startsWith(line,"!Sample_source_name_ch1") && !is.null(cur)){
      source <- trimws(strsplit(line,"=",fixed=TRUE)[[1]][2])
    } else if(startsWith(line,"!Sample_characteristics_ch1") && !is.null(cur)){
      v <- trimws(strsplit(line,"=",fixed=TRUE)[[1]][2])
      if(grepl("^tissue:",v,ignore.case=TRUE)){
        tissue <- trimws(sub("^[^:]+:","",v))
      }
    }
  }
  flush()
  meta <- do.call(rbind,rows)
  meta$patient <- sub("([NT][0-9]*)$","",meta$sample_id)
  meta$side <- sub("^.*?([NT][0-9]*)$","\\1",meta$sample_id)
  meta$type <- substr(meta$side,1,1)
  classify <- function(v){
    key <- gsub("[^a-z]+"," ",tolower(as.character(v)))
    if(grepl("liver",key) && grepl("normal|adjacent",key))
      return("normal_liver_reference")
    if(grepl("colon",key) && grepl("normal|adjacent",key))
      return("normal_colon_reference")
    if(grepl("liver",key) && grepl("tumor|tumour|metast",key))
      return("liver_tumor")
    if(grepl("colon",key) && grepl("tumor|tumour|primary",key))
      return("colon_tumor")
    NA_character_
  }
  meta$reference_class <- vapply(meta$tissue,classify,character(1))
  meta
}

stage11_gse231559_projection <- function(force = FALSE) {
  stage_log("11_GSE231559_directional_portability", {
    frozen_file <- p(DIR$cache_projection,"frozen_P1_v4.rds")
    assert_file(frozen_file)
    frozen <- must_load_step_cache(frozen_file, "frozen P1")
    
    raw_dir <- p(GSE231559_DIR,"raw")
    soft_file <- p(GSE231559_DIR,"GSE231559_family.soft.gz")
    if(!dir.exists(raw_dir)||!file.exists(soft_file)){
      stop("GSE231559尚未就绪：\n",raw_dir,"\n",soft_file)
    }
    
    out <- p(DIR$cache_projection,"GSE231559_frozen_projection.csv")
    if(!force && file.exists(out)){
      message("[11_GSE231559_directional_portability] CACHE HIT：读取冻结投影结果。")
      return(read.csv(out,stringsAsFactors=FALSE))
    }
    
    md <- parse_gse231559_soft(soft_file)
    if(any(is.na(md$reference_class))){
      write.csv(md,p(DIR$audit,"11_GSE231559_unclassified_metadata.csv"),row.names=FALSE)
      stop("GSE231559存在无法按官方tissue字段分类的样本。")
    }
    write.csv(md,p(DIR$tables,"61_GSE231559_reference_metadata.csv"),row.names=FALSE)
    
    p1 <- frozen$P1_top15
    w <- data.frame(gene=p1$gene,weight=p1$weight_norm)
    sets <- make_p1_score_sets(w$gene)
    sets$P1_no_top5 <- setdiff(w$gene, head(w$gene, 5))
    
    mtxs <- sort(list.files(raw_dir,pattern="_matrix\\.mtx\\.gz$",full.names=TRUE))
    if(!length(mtxs))stop("GSE231559 raw目录没有 *_matrix.mtx.gz。")
    rows <- list(); rr <- 0L
    
    for(mf in mtxs){
      stem <- sub("_matrix\\.mtx\\.gz$","",basename(mf))
      bf <- p(raw_dir,paste0(stem,"_barcodes.tsv.gz"))
      ff <- p(raw_dir,paste0(stem,"_features.tsv.gz"))
      X <- read_10x_triplet(mf,bf,ff)
      norm <- normalize_log1p_cp10k_sparse(X)
      
      # 该队列旧流程：state/lineage module均等权；P1权重只用于projection
      my <- module_score(norm,MYELOID_MODULE)
      ep <- module_score(norm,EPITHELIAL_MODULE)
      lyz <- if("LYZ"%in%rownames(X))as.numeric(X["LYZ",]>0)>0 else rep(FALSE,ncol(X))
      inc <- my>ep & lyz
      sm <- sapply(STATE_MODULES,function(gs)module_score(norm,gs))
      if(!is.matrix(sm))sm<-matrix(sm,ncol=length(STATE_MODULES))
      colnames(sm)<-names(STATE_MODULES)
      labels <- colnames(sm)[max.col(sm,ties.method="first")]
      labels[apply(sm,1,max)<=0] <- "Unassigned"
      
      accession <- strsplit(stem,"_",fixed=TRUE)[[1]][1]
      hit <- md[md$gsm==accession,]
      if(!nrow(hit)){
        # 某些文件stem本身可能是sample title；第二路径匹配
        hit <- md[md$sample_id==stem,]
      }
      if(!nrow(hit)) stop("GSE231559 metadata缺失：",stem)
      sample_id <- hit$sample_id[1]
      patient <- hit$patient[1]
      tissue <- hit$type[1]
      if(!tissue%in%c("N","T")) stop("GSE231559无法解析N/T：",sample_id)
      
      for(ver in names(sets)){
        ww <- w[w$gene%in%sets[[ver]],,drop=FALSE]
        sv <- weighted_score(norm,ww)
        for(st in names(STATE_MODULES)){
          mask <- inc & labels==st
          rr<-rr+1L
          rows[[rr]] <- data.frame(
            sample=stem,sample_id=sample_id,patient=patient,tissue=tissue,
            reference_class=hit$reference_class[1],
            state=st,projection_normalization=PARAM$projection_normalization,
            score_version=ver,n_cells=sum(mask),
            mean_score=if(any(mask))mean(sv[mask])else NA_real_,
            median_score=if(any(mask))median(sv[mask])else NA_real_
          )
        }
      }
      rm(X,norm); invisible(gc())
    }
    
    ans <- do.call(rbind,rows)
    write.csv(ans,out,row.names=FALSE)
    
    # 旧版方向可迁移性：每endpoint state至少20 cells
    z <- ans[ans$n_cells>=20 & ans$tissue%in%c("N","T"),]
    agg <- aggregate(
      mean_score~patient+state+score_version+tissue,z,mean,na.rm=TRUE
    )
    pair_rows <- list(); pr <- 0L
    for(ver in unique(agg$score_version)){
      for(st in unique(agg$state)){
        q <- agg[agg$score_version==ver & agg$state==st,]
        piv <- reshape(q,idvar="patient",timevar="tissue",direction="wide")
        if(!all(c("mean_score.N","mean_score.T")%in%names(piv)))next
        piv <- piv[complete.cases(piv[,c("mean_score.N","mean_score.T")]),]
        if(!nrow(piv))next
        for(i in seq_len(nrow(piv))){
          pr<-pr+1L
          pair_rows[[pr]]<-data.frame(
            score_version=ver,patient=piv$patient[i],state=st,
            delta_T_minus_N=piv$mean_score.T[i]-piv$mean_score.N[i]
          )
        }
      }
    }
    deltas <- if(length(pair_rows))do.call(rbind,pair_rows)else data.frame()
    write.csv(deltas,p(DIR$tables,"62_GSE231559_replication_deltas.csv"),row.names=FALSE)
    
    summary_rows <- list(); ss <- 0L
    for(ver in unique(deltas$score_version)){
      g <- deltas[deltas$score_version==ver,]
      ss<-ss+1L
      summary_rows[[ss]]<-data.frame(
        score_version=ver,n_pairs=nrow(g),
        n_positive=sum(g$delta_T_minus_N>0),
        median_delta=median(g$delta_T_minus_N)
      )
    }
    summary <- do.call(rbind,summary_rows)
    write.csv(summary,p(DIR$tables,"63_GSE231559_replication_summary.csv"),row.names=FALSE)
    
    loo <- list(); ll <- 0L
    for(drop in unique(deltas$patient)){
      q <- deltas[deltas$patient!=drop,]
      for(ver in unique(q$score_version)){
        g <- q[q$score_version==ver,]
        ll<-ll+1L
        loo[[ll]]<-data.frame(
          omitted_patient=drop,score_version=ver,
          median_delta=median(g$delta_T_minus_N,na.rm=TRUE),
          n_positive=sum(g$delta_T_minus_N>0,na.rm=TRUE),
          n_pairs=sum(is.finite(g$delta_T_minus_N))
        )
      }
    }
    loo <- do.call(rbind,loo)
    write.csv(loo,p(DIR$tables,"64_GSE231559_LOPO.csv"),row.names=FALSE)
    
    no_p1 <- deltas[deltas$state!="P1_state",]
    no_p1_sum <- aggregate(
      delta_T_minus_N~score_version,no_p1,
      function(v)c(
        n_pairs=length(v),n_positive=sum(v>0),median_delta=median(v)
      )
    )
    # flatten aggregate matrix column
    vals <- no_p1_sum$delta_T_minus_N
    no_p1_sum$delta_T_minus_N <- NULL
    no_p1_sum$n_pairs <- vals[,"n_pairs"]
    no_p1_sum$n_positive <- vals[,"n_positive"]
    no_p1_sum$median_delta <- vals[,"median_delta"]
    write.csv(
      no_p1_sum,
      p(DIR$tables,"65_GSE231559_exclude_P1_state_summary.csv"),
      row.names=FALSE
    )
    
    write_session_info("11_GSE231559_directional_portability")
    ans
  })
}

############################
# Stage 12：空间 projection + residualization
############################

create_spatial_manifest_template <- function() {
  if (file.exists(SPATIAL_MANIFEST)) return(invisible(FALSE))
  template <- data.frame(
    dataset = c(rep("GSE217414",4), rep("GSE225857",6)),
    section = c("19G081","19G0619","19G0635","19G02977","C1","C2","C3","C4","L1","L2"),
    matrix_type = c(rep("h5",4), rep("mtx",6)),
    matrix_path = "",
    barcodes_path = "",
    features_path = "",
    positions_path = "",
    pixel_size_um = NA_real_,
    primary_graph = paste0("knn", PARAM$spatial_knn_k),
    radius_um_sensitivity = paste(PARAM$spatial_radius_um_sensitivity, collapse = ";"),
    # Legacy pixel-radius column retained only for provenance; it is not used
    # by the primary graph or the new physical-distance sensitivity analysis.
    radius_px = NA_real_,
    stringsAsFactors = FALSE
  )
  write.csv(template, SPATIAL_MANIFEST, row.names = FALSE)
  invisible(TRUE)
}

read_10x_h5_rhdf5 <- function(path) {
  # Seurat::Read10X_h5 currently delegates to hdf5r.  Keep an equivalent,
  # read-only fallback for environments where rhdf5 (but not hdf5r) is
  # available; this does not alter the matrix, normalization, or spatial
  # analysis parameters.
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("读取10x HDF5需要hdf5r或rhdf5；当前环境两者均不可用。")
  }
  info <- rhdf5::h5ls(path, recursive = TRUE)
  root <- if (any(info$group == "/matrix" & info$name == "data")) "/matrix" else ""
  read_dataset <- function(name) {
    rhdf5::h5read(path, paste0(root, "/", name))
  }
  data <- read_dataset("data")
  indices <- read_dataset("indices")
  indptr <- read_dataset("indptr")
  shape <- as.integer(read_dataset("shape"))
  if (length(shape) != 2L || length(indptr) != shape[2L] + 1L) {
    stop("10x HDF5矩阵维度或indptr长度不一致：", path)
  }
  X <- Matrix::sparseMatrix(
    i = as.integer(indices) + 1L,
    p = as.integer(indptr),
    x = as.numeric(data),
    dims = shape,
    repr = "C"
  )
  feature_group <- paste0(root, "/features")
  feature_name_path <- paste0(feature_group, "/name")
  if (any(info$group == feature_group & info$name == "name")) {
    genes <- as.character(rhdf5::h5read(path, feature_name_path))
  } else if (any(info$group == root & info$name == "gene_names")) {
    genes <- as.character(read_dataset("gene_names"))
  } else {
    stop("10x HDF5中未找到features/name或gene_names：", path)
  }
  barcodes <- as.character(read_dataset("barcodes"))
  if (length(genes) != nrow(X) || length(barcodes) != ncol(X)) {
    stop("10x HDF5基因/条形码数量与矩阵维度不一致：", path)
  }
  rownames(X) <- make.unique(genes)
  colnames(X) <- make.unique(barcodes)
  X
}

read_spatial_section <- function(row) {
  if (row$matrix_type == "h5") {
    assert_file(row$matrix_path)
    if (requireNamespace("hdf5r", quietly = TRUE)) {
      X <- Seurat::Read10X_h5(row$matrix_path)
      if (is.list(X)) X <- X[[1]]
      X <- methods::as(X, "dgCMatrix")
    } else {
      X <- read_10x_h5_rhdf5(row$matrix_path)
    }
  } else if (row$matrix_type == "mtx") {
    X <- read_10x_triplet(row$matrix_path, row$barcodes_path, row$features_path)
  } else {
    stop("未知matrix_type：", row$matrix_type)
  }
  assert_file(row$positions_path)
  pos <- read.csv(row$positions_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(pos) < 6) stop("Visium positions文件少于6列：", row$positions_path)
  names(pos)[1:6] <- c("barcode","in_tissue","array_row","array_col","pxl_row","pxl_col")
  common <- intersect(colnames(X), pos$barcode[pos$in_tissue == 1])
  X <- X[, common, drop = FALSE]
  pos <- pos[match(common, pos$barcode), , drop = FALSE]
  list(X = X, pos = pos)
}

stage12_spatial_projection <- function(force = FALSE) {
  stage_log("12_spatial_projection", {
    message("Stage 12: checking frozen P1 and spatial manifest ...")
    frozen_file <- p(DIR$cache_projection, "frozen_P1_v4.rds")
    assert_file(frozen_file)
    create_spatial_manifest_template()
    man <- read.csv(SPATIAL_MANIFEST, stringsAsFactors = FALSE)
    message("Stage 12 manifest: ", SPATIAL_MANIFEST, " (", nrow(man), " sections)")
    
    if (any(!nzchar(man$matrix_path)) || any(!nzchar(man$positions_path))) {
      stop(
        "已生成/检测到 spatial_manifest.csv，但路径尚未填写完整：\n",
        SPATIAL_MANIFEST,
        "\n填好两个空间队列的实际文件路径后重新运行Stage 12。"
      )
    }
    
    out <- p(DIR$cache_spatial, "spatial_P1_residuals.csv.gz")
    if (!force && file.exists(out)) {
      message("[12_spatial_projection] CACHE HIT：读取空间 P1 residual 结果。")
      con <- gzfile(out, "rt"); z <- read.csv(con); close(con); return(z)
    }
    
    frozen <- must_load_step_cache(frozen_file, "frozen P1")
    w <- frozen$P1_top15[, c("gene","weight_norm")]
    names(w)[2] <- "weight"
    
    rows <- list()
    for (i in seq_len(nrow(man))) {
      cat("Spatial:", man$dataset[i], man$section[i], "\n")
      z <- read_spatial_section(man[i, ])
      X <- z$X; pos <- z$pos
      norm <- normalize_log1p_cp10k_sparse(X)
      
      score <- weighted_score(norm, w)
      my <- module_score(norm, SPATIAL_COVARIATE_MODULES$myeloid)
      ep <- module_score(norm, SPATIAL_COVARIATE_MODULES$epithelial)
      lv <- module_score(norm, SPATIAL_COVARIATE_MODULES$liver)
      lib <- Matrix::colSums(X)
      
      df <- data.frame(
        dataset = man$dataset[i],
        section = man$section[i],
        barcode = colnames(X),
        pxl_row = pos$pxl_row,
        pxl_col = pos$pxl_col,
        radius_px = man$radius_px[i],
        pixel_size_um = if ("pixel_size_um" %in% names(man))
          suppressWarnings(as.numeric(man$pixel_size_um[i])) else NA_real_,
        projection_normalization = PARAM$projection_normalization,
        P1_score = score,
        myeloid_module = my,
        epithelial_module = ep,
        liver_module = lv,
        log1p_total_UMI = log1p(lib)
      )
      fit <- lm(
        P1_score ~ myeloid_module + epithelial_module + liver_module + log1p_total_UMI,
        data = df
      )
      df$P1_residual <- residuals(fit)
      rows[[i]] <- df
      rm(X,norm); invisible(gc())
    }
    ans <- do.call(rbind, rows)
    con <- gzfile(out, "wt"); write.csv(ans, con, row.names = FALSE); close(con)
    write_session_info("12_spatial_projection")
    ans
  })
}


############################
# Stage 13：空间 Moran I + local High-High（999 permutations）
############################

make_radius_neighbors <- function(coords, radius) {
  check_packages("dbscan")
  fr <- dbscan::frNN(coords, eps = radius, sort = FALSE)
  lapply(seq_along(fr$id), function(i) setdiff(fr$id[[i]], i))
}

make_knn_neighbors <- function(coords, k = 6L) {
  check_packages("FNN")
  nn <- FNN::get.knn(coords, k = k)$nn.index
  lapply(seq_len(nrow(nn)), function(i) nn[i, ])
}

symmetrize_neighbors <- function(nb) {
  n <- length(nb)
  out <- lapply(nb, unique)
  for (i in seq_len(n)) {
    for (j in out[[i]]) {
      if (j >= 1 && j <= n && !i %in% out[[j]]) out[[j]] <- c(out[[j]], i)
    }
  }
  lapply(out, function(x) sort(unique(x[x >= 1 & x <= n])))
}

component_labels <- function(nb, selected) {
  n <- length(nb)
  selected <- as.logical(selected)
  lab <- rep(NA_integer_, n)
  cid <- 0L
  for (start in which(selected)) {
    if (!is.na(lab[start])) next
    cid <- cid + 1L
    queue <- start
    lab[start] <- cid
    while (length(queue)) {
      v <- queue[1]
      queue <- queue[-1]
      nei <- nb[[v]]
      nei <- nei[selected[nei]]
      new <- nei[is.na(lab[nei])]
      if (length(new)) {
        lab[new] <- cid
        queue <- c(queue, new)
      }
    }
  }
  lab
}

moran_I_from_neighbors <- function(x, nb) {
  x <- as.numeric(scale(x))
  n <- length(x)
  s0 <- sum(lengths(nb))
  if (s0 == 0) return(NA_real_)
  num <- 0
  for (i in seq_len(n)) {
    if (length(nb[[i]])) num <- num + sum(x[i] * x[nb[[i]]])
  }
  (n / s0) * num / sum(x^2)
}

spatial_permutation_test <- function(x, nb, B, seed) {
  obs <- moran_I_from_neighbors(x, nb)
  set.seed(seed)
  perm <- replicate(B, moran_I_from_neighbors(sample(x), nb))
  pval <- (1 + sum(perm >= obs, na.rm = TRUE)) / (1 + B)
  c(I = obs, p = pval)
}

local_high_high <- function(x, nb, B, seed) {
  z <- as.numeric(scale(x))
  lag <- sapply(seq_along(nb), function(i) {
    if (!length(nb[[i]])) return(NA_real_)
    mean(z[nb[[i]]])
  })
  obs <- z * lag
  set.seed(seed)
  # Conditional local permutation: keep the focal spot value z_i fixed and
  # randomize only the values drawn for its neighbor set.  This avoids the
  # less specific global permutation in which focal and neighbor values are
  # shuffled together.
  exceed <- integer(length(z))
  for (b in seq_len(B)) {
    for (i in seq_along(nb)) {
      if (!length(nb[[i]]) || !is.finite(obs[i])) next
      pool <- z[-i]
      # A large physical-radius graph can contain more neighbors than the
      # available pool; without replacement is used whenever possible, while
      # replacement keeps the conditional null defined for that edge case.
      draw <- sample(
        pool, length(nb[[i]]),
        replace = length(nb[[i]]) > length(pool)
      )
      null_i <- z[i] * mean(draw)
      if (is.finite(null_i) && null_i >= obs[i]) exceed[i] <- exceed[i] + 1L
    }
  }
  pval <- vapply(seq_along(obs), function(i) {
    if (!is.finite(obs[i])) return(NA_real_)
    (1 + exceed[i]) / (1 + B)
  }, numeric(1))
  q <- p.adjust(pval, method = "BH")
  data.frame(
    z = z, spatial_lag = lag, local_I = obs,
    local_p = pval, local_q = q,
    high_high = is.finite(q) & z > 0 & lag > 0 & q < 0.05,
    permutation_scheme = "conditional_neighbor_randomization_fixed_focal",
    permutations = B
  )
}

stage13_spatial_statistics <- function(force = FALSE) {
  stage_log("13_spatial_statistics", {
    input <- p(DIR$cache_spatial, "spatial_P1_residuals.csv.gz")
    assert_file(input)
    con <- gzfile(input, "rt"); d <- read.csv(con); close(con)
    
    if (!"pixel_size_um" %in% names(d) ||
        any(!is.finite(as.numeric(d$pixel_size_um)) | as.numeric(d$pixel_size_um) <= 0)) {
      stop(
        "空间敏感性分析需要每张切片的 pixel_size_um。请在 spatial_manifest.csv " ,
        "中填写真实物理像素尺寸后重新运行Stage 12和Stage 13。"
      )
    }
    
    global_rows <- list()
    local_rows <- list()
    gg <- 0L; ll <- 0L
    
    section_cols <- c("dataset", "section", "radius_px")
    if ("pixel_size_um" %in% names(d)) section_cols <- c(section_cols, "pixel_size_um")
    sections <- unique(d[, section_cols, drop = FALSE])
    for (i in seq_len(nrow(sections))) {
      ds <- sections$dataset[i]
      sec <- sections$section[i]
      z <- d[d$dataset == ds & d$section == sec, , drop = FALSE]
      coords <- as.matrix(z[, c("pxl_col","pxl_row")])
      message(
        "[13_spatial_statistics] section ", i, "/", nrow(sections),
        " ", ds, "/", sec, " (", nrow(z), " spots)"
      )
      
      # Primary graph is kNN(k=6), independent of platform-specific pixel
      # spacing.  Physical-distance graphs are sensitivity analyses only and
      # require the section-specific Visium pixel size in the manifest.
      graphs <- list(
        knn6 = symmetrize_neighbors(make_knn_neighbors(coords, PARAM$spatial_knn_k))
      )
      px_um <- if ("pixel_size_um" %in% names(sections))
        suppressWarnings(as.numeric(sections$pixel_size_um[i])) else NA_real_
      if (is.finite(px_um) && px_um > 0) {
        for (radius_um in PARAM$spatial_radius_um_sensitivity) {
          radius_px <- radius_um / px_um
          graphs[[paste0("radius_", radius_um, "um")]] <-
            symmetrize_neighbors(make_radius_neighbors(coords, radius_px))
        }
      }
      
      for (gname in names(graphs)) {
        nb <- graphs[[gname]]
        seed <- PARAM$seed + 1300000L + i*100L + match(gname, names(graphs))
        gt <- spatial_permutation_test(
          z$P1_residual, nb, PARAM$n_spatial_perm, seed
        )
        gg <- gg + 1L
        global_rows[[gg]] <- data.frame(
          dataset = ds, section = sec, graph = gname,
          radius_px = if (grepl("^radius_", gname))
            as.numeric(sub("^radius_([0-9.]+)um$", "\\1", gname)) / px_um else NA_real_,
          radius_um = if (grepl("^radius_", gname))
            as.numeric(sub("^radius_([0-9.]+)um$", "\\1", gname)) else NA_real_,
          pixel_size_um = px_um,
          k = if (gname=="knn6") PARAM$spatial_knn_k else NA,
          n_spots = nrow(z),
          Moran_I = unname(gt["I"]),
          empirical_p = unname(gt["p"]),
          permutations = PARAM$n_spatial_perm
        )
        
        loc <- local_high_high(
          z$P1_residual, nb, PARAM$n_spatial_perm, seed + 50L
        )
        lab <- component_labels(nb, loc$high_high)
        comp_size <- rep(0L, length(lab))
        if (any(!is.na(lab))) {
          tt <- table(lab)
          comp_size[!is.na(lab)] <- as.integer(tt[as.character(lab[!is.na(lab)])])
        }
        loc$HH_component <- lab
        loc$HH_component_size <- comp_size
        loc$robust_HH_ge3 <- loc$high_high & comp_size >= 3L
        loc$sensitivity_HH_ge2 <- loc$high_high & comp_size >= 2L
        loc$sensitivity_HH_ge5 <- loc$high_high & comp_size >= 5L
        
        ll <- ll + 1L
        local_rows[[ll]] <- cbind(
          z[, c("dataset","section","barcode","pxl_row","pxl_col")],
          graph = gname,
          loc
        )
        message(
          "[13_spatial_statistics] completed graph ", gname,
          " for ", ds, "/", sec
        )
      }
    }
    
    global <- do.call(rbind, global_rows)
    global$BH_q <- ave(
      global$empirical_p,
      interaction(global$dataset, global$graph, drop = TRUE),
      FUN = function(x) p.adjust(x, "BH")
    )
    local <- do.call(rbind, local_rows)
    
    write.csv(global, p(DIR$tables, "31_spatial_global_Moran_999.csv"),
              row.names = FALSE)
    con <- gzfile(p(DIR$tables, "32_spatial_local_high_high_999.csv.gz"), "wt")
    write.csv(local, con, row.names = FALSE); close(con)
    
    burden <- aggregate(
      cbind(high_high, robust_HH_ge3, sensitivity_HH_ge2, sensitivity_HH_ge5) ~
        dataset + section + graph,
      local,
      FUN = sum
    )
    total_spots <- aggregate(
      barcode ~ dataset + section + graph, local, length
    )
    names(total_spots)[4] <- "n_spots"
    burden <- merge(burden, total_spots, by = c("dataset","section","graph"))
    burden$robust_HH_ge3_fraction <- burden$robust_HH_ge3 / burden$n_spots
    write.csv(
      burden,
      p(DIR$tables, "33_spatial_hotspot_burden_999.csv"),
      row.names = FALSE
    )
    
    write_session_info("13_spatial_statistics")
    list(global = global, local = local, burden = burden)
  })
}



############################
# Stage 14：发现队列补充审计
# - P1–Pk跨状态承载
# - P1–Pk patient×fine-state within-state contrasts
# - GSE178318 14 fine-state composition/activity
# - full-loading preranked enrichment（若msigdbr/fgsea已安装）
############################

stage14_discovery_supplementary <- function(force = FALSE) {
  stage_log("14_discovery_supplementary_robustness", {
    frozen_file <- p(DIR$cache_projection,"frozen_P1_v4.rds")
    my_file <- p(DIR$cache_annotation,"fine_myeloid_dc_annotated.rds")
    assert_file(frozen_file); assert_file(my_file)
    frozen <- must_load_step_cache(frozen_file, "frozen P1")
    my <- must_load_step_cache(my_file, "fine myeloid annotation")
    
    usage <- data.frame(
      barcode=colnames(frozen$H),
      t(frozen$H),check.names=FALSE
    )
    md <- my[[]]; md$barcode <- rownames(md)
    usage <- merge(
      md[,c("barcode","sample_id","patient_id","tissue","fine_label","fine_status")],
      usage,by="barcode",sort=FALSE
    )
    programs <- paste0("P",seq_len(frozen$selected_rank))
    
    # P1-Pk carriage matrix：fine-state平均program usage
    carriage <- aggregate(
      usage[,programs,drop=FALSE],
      list(fine_label=usage$fine_label),
      mean,na.rm=TRUE
    )
    write.csv(
      carriage,
      p(DIR$tables,"66_discovery_fine_state_P1_Pk_carriage_matrix.csv"),
      row.names=FALSE
    )
    
    # patient×state within-state PT-LM contrasts for every program
    z <- aggregate(
      usage[,programs,drop=FALSE],
      usage[,c("patient_id","tissue","fine_label")],
      mean,na.rm=TRUE
    )
    contrast_rows <- list(); cr <- 0L
    for(pp in programs){
      for(pt in unique(z$patient_id)){
        for(st in unique(z$fine_label)){
          q <- z[z$patient_id==pt & z$fine_label==st,]
          if(all(c("CRC","LM")%in%q$tissue)){
            # Single-column data-frame extraction drops to a vector by default;
            # use an explicit first element rather than a two-dimensional index.
            a <- q[q$tissue=="CRC", pp, drop=TRUE][1]
            b <- q[q$tissue=="LM", pp, drop=TRUE][1]
            if(is.finite(a)&&is.finite(b)){
              cr<-cr+1L
              contrast_rows[[cr]]<-data.frame(
                program=pp,patient_id=pt,fine_label=st,
                delta_within_state=b-a
              )
            }
          }
        }
      }
    }
    contrasts <- do.call(rbind,contrast_rows)
    write.csv(
      contrasts,
      p(DIR$tables,"67_discovery_NMF_within_state_patient_contrasts.csv"),
      row.names=FALSE
    )
    csum <- aggregate(
      delta_within_state~program+fine_label,
      contrasts,
      function(v)c(
        n_pairs=length(v),n_positive=sum(v>0),
        median_delta=median(v),mean_delta=mean(v)
      )
    )
    vv <- csum$delta_within_state
    csum$delta_within_state <- NULL
    csum$n_pairs <- vv[,"n_pairs"]; csum$n_positive <- vv[,"n_positive"]
    csum$median_delta <- vv[,"median_delta"]; csum$mean_delta <- vv[,"mean_delta"]
    write.csv(
      csum,
      p(DIR$tables,"68_discovery_NMF_within_state_summary.csv"),
      row.names=FALSE
    )
    
    # 14 fine-state decomposition on consensus NMF P1 usage.
    # Denominator is the entire 14-state myeloid/DC compartment for each patient×tissue.
    p1col <- "P1"
    cell <- usage[,c("patient_id","tissue","fine_label",p1col)]
    names(cell)[4] <- "score"
    patient_rows <- list(); pr <- 0L
    for(pt in sort(unique(cell$patient_id))){
      q <- cell[cell$patient_id==pt & cell$tissue%in%c("CRC","LM"),]
      if(!all(c("CRC","LM")%in%q$tissue))next
      states <- sort(unique(q$fine_label))
      # observed p/mu for every side
      get <- function(tissue,st){
        v <- q$score[q$tissue==tissue & q$fine_label==st]
        c(n=length(v),mu=if(length(v))mean(v)else NA_real_)
      }
      nC <- muC <- nL <- muL <- numeric(length(states))
      for(i in seq_along(states)){
        a<-get("CRC",states[i]);b<-get("LM",states[i])
        nC[i]<-a["n"];muC[i]<-a["mu"];nL[i]<-b["n"];muL[i]<-b["mu"]
      }
      pC <- nC/sum(nC);pL <- nL/sum(nL)
      common <- nC>0 & nL>0
      oneC <- nC>0 & nL==0
      oneL <- nL>0 & nC==0
      
      comp <- sum((pL[common]-pC[common])*(muL[common]+muC[common])/2)
      act <- sum((muL[common]-muC[common])*(pL[common]+pC[common])/2)
      one_contrib <- 0
      if(any(oneL))one_contrib <- one_contrib + sum(pL[oneL]*muL[oneL])
      if(any(oneC))one_contrib <- one_contrib - sum(pC[oneC]*muC[oneC])
      total <- sum(pL[nL>0]*muL[nL>0])-sum(pC[nC>0]*muC[nC>0])
      err <- total-comp-act-one_contrib
      if(abs(err)>1e-10)stop("GSE178318 14-state closure失败。")
      pr<-pr+1L
      patient_rows[[pr]] <- data.frame(
        patient=pt,overall_fine_state_delta=total,
        fine_state_composition=comp,
        within_fine_state_activity=act,
        one_sided_state_contribution=one_contrib,
        activity_minus_composition=act-comp,
        n_common_states=sum(common),
        n_one_sided_states=sum(oneC|oneL),
        closure_error=err
      )
    }
    fine_decomp <- do.call(rbind,patient_rows)
    write.csv(
      fine_decomp,
      p(DIR$tables,"69_GSE178318_14fine_state_decomposition_by_patient.csv"),
      row.names=FALSE
    )
    sum_rows <- list(); sr <- 0L
    for(metric in c(
      "overall_fine_state_delta","fine_state_composition",
      "within_fine_state_activity","one_sided_state_contribution",
      "activity_minus_composition"
    )){
      v <- fine_decomp[[metric]]
      ci <- bootstrap_median_ci(v,PARAM$n_boot,seed=20260818L+sr)
      sr<-sr+1L
      sum_rows[[sr]] <- data.frame(
        component=metric,n_patients=length(v),median=median(v),
        ci_low=ci[1],ci_high=ci[2],n_positive=sum(v>0),
        exact_sign_p_greater=exact_sign_p_greater(v)
      )
    }
    fine_sum <- do.call(rbind,sum_rows)
    write.csv(
      fine_sum,
      p(DIR$tables,"70_GSE178318_14fine_state_decomposition_summary.csv"),
      row.names=FALSE
    )
    
    # Full loading enrichment，严格保持“全loading排序，不用top15截断”
    enrichment_status <- data.frame(
      attempted=FALSE,completed=FALSE,n_loading_genes=nrow(frozen$W),
      note="msigdbr/fgsea not checked"
    )
    if(requireNamespace("msigdbr",quietly=TRUE) &&
       requireNamespace("fgsea",quietly=TRUE) &&
       requireNamespace("BiocParallel",quietly=TRUE)){
      enrichment_status$attempted <- TRUE
      loading <- frozen$W[,"P1"]
      names(loading) <- rownames(frozen$W)
      loading <- sort(loading[is.finite(loading)],decreasing=TRUE)
      run_collection <- function(collection,subcollection,label){
        ms <- msigdbr::msigdbr(
          species="Homo sapiens",collection=collection,subcollection=subcollection
        )
        pathways <- split(ms$gene_symbol,ms$gs_name)
        pathways <- lapply(pathways,unique)
        pathways <- pathways[
          vapply(pathways,function(x)sum(x%in%names(loading))>=10L,logical(1))
        ]
        set.seed(20260815L)
        res <- fgsea::fgseaSimple(
          pathways=pathways,stats=loading,nperm=1000L,
          minSize=10L,maxSize=500L,scoreType="pos",
          BPPARAM=BiocParallel::SerialParam()
        )
        res <- as.data.frame(res)
        res$collection <- label
        res$leadingEdge <- vapply(res$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
        res
      }
      go <- run_collection("C5","GO:BP","GO_BP")
      re <- run_collection("C2","CP:REACTOME","Reactome")
      enrich <- rbind(go,re)
      write.csv(
        enrich,
        p(DIR$tables,"71_P1_full_loading_preranked_enrichment_all_terms.csv"),
        row.names=FALSE
      )
      target <- enrich[
        grepl(
          "lysosom|phagocyt|endocyt|lipid transport|lipid localization|cholesterol|lipoprotein|vesicle",
          tolower(enrich$pathway)
        ),,
        drop=FALSE
      ]
      write.csv(
        target,
        p(DIR$tables,"72_P1_full_loading_targeted_enrichment_terms.csv"),
        row.names=FALSE
      )
      enrichment_status$completed <- TRUE
      enrichment_status$note <- paste0(
        "GO_BP tested=",sum(enrich$collection=="GO_BP"),
        "; Reactome tested=",sum(enrich$collection=="Reactome"),
        "; FDR<0.05=",sum(enrich$padj<0.05,na.rm=TRUE)
      )
    } else {
      enrichment_status$note <-
        "Skipped only because msigdbr/fgsea/BiocParallel is not installed; no substitute analysis used."
      warning(enrichment_status$note)
    }
    write.csv(
      enrichment_status,
      p(DIR$audit,"14_full_loading_enrichment_status.csv"),
      row.names=FALSE
    )
    write_session_info("14_discovery_supplementary_robustness")
    list(
      carriage=carriage,within_state=contrasts,
      fine_state_decomposition=fine_decomp,
      fine_state_summary=fine_sum,
      enrichment_status=enrichment_status
    )
  })
}


############################
# Stage 14b：P1技术复杂度与HC3协变量审计
############################

stage14b_complexity_hc3 <- function(force = FALSE) {
  stage_log("14b_complexity_HC3_audit", {
    check_packages(c("Seurat", "data.table", "sandwich"))
    frozen_file <- p(DIR$cache_projection, "frozen_P1_v4.rds")
    my_file <- p(DIR$cache_annotation, "fine_myeloid_dc_annotated.rds")
    assert_file(frozen_file); assert_file(my_file)
    
    # The stage writes one file per model plus a combined audit table.  The
    # former cache check looked for a filename that is never produced
    # (100_HC3_full_covariate_model.csv), so Stage 14b was needlessly rerun
    # on every invocation.  Treat the actual complete set of outputs as the
    # cache unit and preserve the existing filenames for backwards
    # compatibility.
    model_files <- c(
      p(DIR$tables, "100_HC3_full_covariate_model_loo_P1.csv"),
      p(DIR$tables, "100_HC3_complexity_only_model.csv"),
      p(DIR$tables, "100_HC3_models_all_terms.csv")
    )
    if (!force && all(file.exists(model_files))) {
      message("[14b_complexity_HC3_audit] CACHE HIT：读取已有 HC3 模型结果。")
      return(read.csv(model_files[[3]], stringsAsFactors = FALSE))
    }
    
    frozen <- must_load_step_cache(frozen_file, "frozen P1")
    my <- must_load_step_cache(my_file, "fine myeloid annotation")
    md <- my[[]]
    md$barcode <- rownames(md)
    H <- frozen$H
    usage <- data.frame(barcode = colnames(H), t(H), check.names = FALSE)
    usage <- merge(
      md[, intersect(c("barcode", "patient_id", "tissue", "fine_label",
                       "nCount_RNA", "nFeature_RNA", "percent.mt"), names(md)),
         drop = FALSE],
      usage, by = "barcode", sort = FALSE
    )
    required <- c("P1", paste0("P", 2:6), "patient_id", "tissue",
                  "fine_label", "nCount_RNA", "nFeature_RNA", "percent.mt")
    if (!all(required %in% names(usage))) {
      stop("HC3审计缺少字段：", paste(setdiff(required, names(usage)), collapse = ","))
    }
    usage$log_nCount <- log1p(usage$nCount_RNA)
    usage$log_nFeature <- log1p(usage$nFeature_RNA)
    usage$total_nmf_loo_p1 <- rowSums(usage[, paste0("P", 2:6), drop = FALSE])
    
    corr_vars <- c("log_nCount", "log_nFeature", "percent.mt", "total_nmf_loo_p1")
    corr <- do.call(rbind, lapply(corr_vars, function(v) {
      ok <- is.finite(usage$P1) & is.finite(usage[[v]])
      data.frame(
        variable = v,
        spearman_rho = suppressWarnings(cor(usage$P1[ok], usage[[v]][ok], method = "spearman")),
        n = sum(ok),
        note = "descriptive only; cells are not independent biological replicates",
        stringsAsFactors = FALSE
      )
    }))
    write.csv(corr, p(DIR$tables, "98_P1_complexity_cell_correlations.csv"), row.names = FALSE)
    
    z <- aggregate(
      cbind(P1, log_nCount, log_nFeature, percent.mt, total_nmf_loo_p1) ~
        patient_id + tissue + fine_label,
      usage, mean, na.action = na.omit
    )
    n_cells <- aggregate(
      barcode ~ patient_id + tissue + fine_label,
      usage, length
    )
    names(n_cells)[4] <- "n_cells"
    z <- merge(z, n_cells, by = c("patient_id", "tissue", "fine_label"), sort = FALSE)
    write.csv(z, p(DIR$tables, "99_P1_complexity_patient_state.csv"), row.names = FALSE)
    # The aggregate-unit count is determined by the current marker-defined
    # patient x tissue x fine-state support.  Derive the expected count from
    # the same complete-case rows used by aggregate(), rather than carrying a
    # legacy hard-coded value (the old annotation happened to yield 185).
    vars_needed <- c(
      "P1", "log_nCount", "log_nFeature", "percent.mt",
      "total_nmf_loo_p1", "patient_id", "tissue", "fine_label"
    )
    u_complete <- usage[
      complete.cases(usage[, vars_needed, drop = FALSE]),
      c("patient_id", "tissue", "fine_label"),
      drop = FALSE
    ]
    expected_units <- nrow(unique(u_complete))
    aggregate_audit <- data.frame(
      expected_complete_case_units = expected_units,
      observed_aggregate_units = nrow(z),
      aggregation_unit = "patient x tissue x fine_label",
      status = if (expected_units == nrow(z)) "PASS" else "FAIL",
      stringsAsFactors = FALSE
    )
    write.csv(
      aggregate_audit,
      p(DIR$audit, "14b_HC3_aggregate_unit_audit.csv"),
      row.names = FALSE
    )
    if (expected_units != nrow(z)) {
      stop(
        "HC3 aggregate-unit mismatch: expected ", expected_units,
        ", obtained ", nrow(z),
        ". See 14b_HC3_aggregate_unit_audit.csv."
      )
    }
    
    fit_one <- function(formula_text, model_name) {
      fit <- stats::lm(stats::as.formula(formula_text), data = z)
      V <- sandwich::vcovHC(fit, type = "HC3")
      est <- stats::coef(fit)
      se <- sqrt(diag(V))
      stat <- est / se
      out <- data.frame(
        model = model_name,
        term = names(est),
        estimate = unname(est),
        robust_se_HC3 = unname(se),
        ci_low = unname(est - stats::qnorm(0.975) * se),
        ci_high = unname(est + stats::qnorm(0.975) * se),
        p_value = unname(2 * stats::pnorm(-abs(stat))),
        n_aggregate_units = nrow(z),
        formula = formula_text,
        stringsAsFactors = FALSE
      )
      write.csv(out, p(DIR$tables, paste0("100_", model_name, ".csv")), row.names = FALSE)
      out
    }
    
    full_formula <- "P1 ~ tissue + log_nCount + log_nFeature + percent.mt + total_nmf_loo_p1 + patient_id + fine_label"
    complexity_formula <- "P1 ~ tissue + log_nCount + log_nFeature + percent.mt + patient_id + fine_label"
    full <- fit_one(full_formula, "HC3_full_covariate_model_loo_P1")
    simple <- fit_one(complexity_formula, "HC3_complexity_only_model")
    write.csv(rbind(full, simple), p(DIR$tables, "100_HC3_models_all_terms.csv"), row.names = FALSE)
    write_session_info("14b_complexity_HC3_audit")
    list(cell_correlations = corr, aggregate = z, full_model = full, complexity_only = simple)
  })
}

############################
# Stage 15：5,000 expression/detection-matched random-signature null
############################

stage15_matched_random_signature_null <- function(force = FALSE) {
  stage_log("15_matched_random_signature_null_5000", {
    frozen_file <- p(DIR$cache_projection,"frozen_P1_v4.rds")
    stats_file <- p(DIR$tables,"42_GSE178318_NMF_input_gene_expression_detection_stats.csv")
    raw_dir <- p(GSE245552_DIR,"raw")
    # GEO archives in some local bundles are unpacked under `extracted`
    # rather than `raw`.  Prefer the canonical directory but accept the
    # equivalent extracted layout without copying or altering the source data.
    if (!dir.exists(raw_dir) && dir.exists(p(GSE245552_DIR, "extracted"))) {
      raw_dir <- p(GSE245552_DIR, "extracted")
      message("[15_matched_random_signature_null_5000] 使用GSE245552 extracted目录：", raw_dir)
    }
    soft_file <- p(GSE245552_DIR,"GSE245552_family.soft.gz")
    assert_file(frozen_file);assert_file(stats_file);assert_file(soft_file)
    if(!dir.exists(raw_dir))stop("GSE245552 raw目录不存在。")
    
    out_summary <- p(DIR$tables,"73_matched_random_signature_null_results_5000.csv")
    balance_file <- p(DIR$tables,"74b_matched_random_signature_balance_audit.csv")
    sensitivity_summary_file <- p(DIR$tables,"73b_matched_random_signature_null_sensitivity_5000.csv")
    sensitivity_balance_file <- p(DIR$tables,"74c_matched_random_signature_sensitivity_balance_5000.csv")
    sensitivity_definition_file <- p(DIR$tables,"74d_matched_random_signature_sensitivity_definitions.csv")
    if(!force && file.exists(out_summary) && file.exists(balance_file) &&
       file.exists(sensitivity_summary_file) && file.exists(sensitivity_balance_file) &&
       file.exists(sensitivity_definition_file)){
      message("[15_matched_random_signature_null_5000] CACHE HIT：读取 matched-null 结果及匹配平衡审计。")
      return(read.csv(out_summary,stringsAsFactors=FALSE))
    }
    if(!force && file.exists(out_summary) && !file.exists(balance_file)){
      message("[15_matched_random_signature_null_5000] 旧版结果缺少匹配平衡审计，将重新计算。")
    }
    
    frozen <- must_load_step_cache(frozen_file, "frozen P1")
    p1 <- frozen$P1_top15[order(frozen$P1_top15$rank),]
    p1_genes <- p1$gene
    p1_weights <- p1$weight_norm
    stats <- read.csv(stats_file,stringsAsFactors=FALSE)
    
    # 所有GSE245552 feature matrix的基因交集
    ffs <- sort(list.files(raw_dir,pattern="_features\\.tsv\\.gz$",full.names=TRUE))
    if(!length(ffs))stop("GSE245552没有features文件。")
    feat_sets <- lapply(ffs,function(ff){
      q<-read.delim(gzfile(ff),header=FALSE,stringsAsFactors=FALSE)
      unique(as.character(q[[min(2,ncol(q))]]))
    })
    ext_genes <- Reduce(intersect,feat_sets)
    
    stats <- stats[stats$gene%in%ext_genes,]
    stats$q_mean <- rank(stats$mean_log_expression,ties.method="average")/nrow(stats)
    stats$q_detect <- rank(stats$detection_frequency,ties.method="average")/nrow(stats)
    ps <- stats[match(p1_genes,stats$gene),]
    if(anyNA(ps$gene))stop("P1 top15至少一个基因缺少discovery expression/detection stats。")
    cand <- stats[!stats$gene%in%p1_genes,]
    if(nrow(cand)<PARAM$matched_null_knn*2)stop("matched-null候选基因不足。")
    
    # Match on the original expression/detection scales first, and use the
    # quantile distance only to rank candidates within the admissible pool.
    # This avoids the previous behaviour in which a broad 50-gene pool was
    # sampled probabilistically even when its expression level was plainly
    # incompatible with the target gene.
    expr_caliper <- as.numeric(PARAM$matched_null_expr_caliper)
    detect_caliper <- as.numeric(PARAM$matched_null_detection_caliper)
    if (!is.finite(expr_caliper) || expr_caliper <= 0 ||
        !is.finite(detect_caliper) || detect_caliper <= 0) {
      stop("matched-null caliper必须为正数。")
    }
    matching_info <- lapply(seq_along(p1_genes), function(j) {
      expr_gap <- abs(cand$mean_log_expression - ps$mean_log_expression[j])
      detect_gap <- abs(cand$detection_frequency - ps$detection_frequency[j])
      qdist <- sqrt(
        (cand$q_mean - ps$q_mean[j])^2 +
          (cand$q_detect - ps$q_detect[j])^2
      )
      admissible <- which(
        is.finite(expr_gap) & is.finite(detect_gap) &
          is.finite(qdist) & expr_gap <= expr_caliper &
          detect_gap <= detect_caliper
      )
      list(
        expr_gap = expr_gap,
        detect_gap = detect_gap,
        qdist = qdist,
        admissible = admissible
      )
    })

    set.seed(20260815L)
    sig_rows <- vector("list",PARAM$matched_null_n*length(p1_genes))
    sr <- 0L
    for(sid in seq_len(PARAM$matched_null_n)){
      order_target <- sample(seq_along(p1_genes))
      used <- character()
      chosen <- rep(NA_character_,length(p1_genes))
      chosen_status <- rep(NA_character_, length(p1_genes))
      chosen_pool_size <- integer(length(p1_genes))
      chosen_caliper_pool_size <- integer(length(p1_genes))
      chosen_qdist <- rep(NA_real_, length(p1_genes))
      for(j in order_target){
        mi <- matching_info[[j]]
        admissible <- mi$admissible
        available_admissible <- admissible[!cand$gene[admissible] %in% used]
        if (length(available_admissible)) {
          ids <- available_admissible[order(mi$qdist[available_admissible])]
          status <- "within_caliper"
        } else {
          # If the caliper is empty (or all admissible genes are already
          # used), retain the nearest unused candidate and label the draw.
          # It remains in the null so the requested signature count is fixed,
          # while the balance audit makes this limitation explicit.
          available_all <- which(!cand$gene %in% used & is.finite(mi$qdist))
          if (!length(available_all))
            stop("matched-null无法为目标基因分配未使用候选基因：", p1_genes[j])
          ids <- available_all[order(mi$qdist[available_all])]
          status <- if (length(admissible))
            "caliper_blocked_by_uniqueness" else "extreme_target_no_caliper"
        }
        pool_size <- length(ids)
        ids <- ids[seq_len(min(PARAM$matched_null_knn, length(ids)))]
        dd <- mi$qdist[ids]
        temp <- max(median(dd),1e-3)
        prob <- exp(-(dd-min(dd))/temp);prob<-prob/sum(prob)
        # `sample(x, 1)` treats a length-one numeric x as the integer range
        # 1:x; use an index explicitly so a one-candidate pool remains valid
        # and its probability vector cannot be misinterpreted.
        pick <- ids[sample.int(length(ids), 1L, prob = prob)]
        chosen[j] <- cand$gene[pick]
        chosen_status[j] <- status
        chosen_pool_size[j] <- pool_size
        chosen_caliper_pool_size[j] <- length(available_admissible)
        chosen_qdist[j] <- mi$qdist[pick]
        used <- c(used,chosen[j])
      }
      for(j in seq_along(p1_genes)){
        sr<-sr+1L
        cc <- cand[cand$gene==chosen[j],][1,]
        sig_rows[[sr]] <- data.frame(
          signature_id=sid,target_rank=j,target_gene=p1_genes[j],
          matched_gene=chosen[j],target_weight=p1_weights[j],
          # Keep the target-specific weight attached to the matched feature;
          # randomizing weights would change two null properties at once.
          assigned_weight=p1_weights[j],
          target_mean_log_expression=ps$mean_log_expression[j],
          target_detection_frequency=ps$detection_frequency[j],
          matched_mean_log_expression=cc$mean_log_expression,
          matched_detection_frequency=cc$detection_frequency,
          abs_expression_gap=abs(cc$mean_log_expression-ps$mean_log_expression[j]),
          abs_detection_gap=abs(cc$detection_frequency-ps$detection_frequency[j]),
          quantile_distance=chosen_qdist[j],
          matching_status=chosen_status[j],
          candidate_pool_size=chosen_pool_size[j],
          caliper_pool_size=chosen_caliper_pool_size[j],
          expression_caliper=expr_caliper,
          detection_caliper=detect_caliper,
          stringsAsFactors=FALSE
        )
      }
    }
    signatures <- do.call(rbind,sig_rows)
    write.csv(
      signatures,
      p(DIR$tables,"74_matched_random_signatures_5000.csv"),
      row.names=FALSE
    )
    balance <- do.call(rbind, lapply(split(signatures, signatures$target_gene),
      function(q) {
        data.frame(
          target_gene = q$target_gene[[1]],
          n_draws = nrow(q),
          n_within_caliper = sum(q$matching_status == "within_caliper"),
          fraction_within_caliper = mean(q$matching_status == "within_caliper"),
          n_caliper_blocked = sum(q$matching_status == "caliper_blocked_by_uniqueness"),
          n_extreme_target = sum(q$matching_status == "extreme_target_no_caliper"),
          median_abs_expression_gap = median(q$abs_expression_gap, na.rm=TRUE),
          median_abs_detection_gap = median(q$abs_detection_gap, na.rm=TRUE),
          median_quantile_distance = median(q$quantile_distance, na.rm=TRUE),
          target_mean_log_expression = q$target_mean_log_expression[[1]],
          target_detection_frequency = q$target_detection_frequency[[1]],
          expression_caliper = expr_caliper,
          detection_caliper = detect_caliper,
          stringsAsFactors=FALSE
        )
      }))
    write.csv(
      balance,
      p(DIR$tables,"74b_matched_random_signature_balance_audit.csv"),
      row.names=FALSE
    )
    
    # 再沿用同一GSE245552 state assignment，提取候选基因的patient×sample×state均值
    tissue_map <- parse_gse245552_soft(soft_file)
    p1_weight_map <- setNames(p1_weights,p1_genes)
    state_score_fun <- function(norm){
      sm <- sapply(STATE_MODULES,function(gs){
        g <- intersect(gs,rownames(norm))
        if(!length(g))return(rep(0,ncol(norm)))
        ww <- vapply(g,function(x)if(x%in%names(p1_weight_map))p1_weight_map[[x]]else 1,
                     numeric(1))
        ww<-ww/sum(ww)
        as.numeric(Matrix::t(norm[g,,drop=FALSE])%*%ww)
      })
      if(!is.matrix(sm))sm<-matrix(sm,ncol=length(STATE_MODULES))
      colnames(sm)<-names(STATE_MODULES);sm
    }
    candidate_genes <- unique(c(p1_genes,signatures$matched_gene))
    mtxs <- sort(list.files(raw_dir,pattern="_matrix\\.mtx\\.gz$",full.names=TRUE))
    state_rows <- list(); rr <- 0L
    for(mf in mtxs){
      stem<-sub("_matrix\\.mtx\\.gz$","",basename(mf))
      X<-read_10x_triplet(
        mf,p(raw_dir,paste0(stem,"_barcodes.tsv.gz")),
        p(raw_dir,paste0(stem,"_features.tsv.gz"))
      )
      norm<-normalize_log1p_cp10k_sparse(X)
      my<-module_score(norm,MYELOID_MODULE);ep<-module_score(norm,EPITHELIAL_MODULE)
      lyz<-if("LYZ"%in%rownames(X))as.numeric(X["LYZ",]>0)>0 else rep(FALSE,ncol(X))
      inc<-my>ep&lyz
      sm<-state_score_fun(norm)
      lab<-colnames(sm)[max.col(sm,ties.method="first")]
      lab[apply(sm,1,max)<=0]<-"Unassigned"
      
      mm<-regexpr("(Pt[0-9]+_[TZ][0-9]+)$",stem,perl=TRUE)
      sample_id<-if(mm[1]>0)regmatches(stem,mm)else stem
      patient<-sub("_.*$","",sample_id)
      tissue<-unname(tissue_map[sample_id])
      g<-intersect(candidate_genes,rownames(norm))
      for(st in c("Mono","P3_C1Q_MHCII","Other","P2_inflammatory")){
        mask<-inc&lab==st
        if(!any(mask))next
        vals<-Matrix::rowMeans(norm[g,mask,drop=FALSE])
        rr<-rr+1L
        row<-data.frame(
          sample=sample_id,patient=patient,tissue=tissue,state=st,
          n_cells=sum(mask),stringsAsFactors=FALSE
        )
        for(j in seq_along(g))row[[g[j]]]<-vals[j]
        state_rows[[rr]]<-row
      }
      rm(X,norm);invisible(gc())
    }
    state_df<-do.call(rbind,state_rows)
    # 缺失gene列由0补齐是不允许的；所有signature gene必须在全部external matrices交集中
    needed<-candidate_genes
    missing_cols<-setdiff(needed,names(state_df))
    if(length(missing_cols))stop("matched-null state mean表缺少基因：",paste(missing_cols,collapse=","))
    
    score_signature <- function(genes,weights){
      ww<-weights/sum(weights)
      as.numeric(as.matrix(state_df[,genes,drop=FALSE])%*%ww)
    }
    make_deltas <- function(scores){
      tmp<-state_df[,c("patient","tissue","state")]
      tmp$score<-scores
      rows<-list();ii<-0L
      for(pt in unique(tmp$patient)){
        for(st in c("Mono","P3_C1Q_MHCII","Other","P2_inflammatory")){
          q<-tmp[tmp$patient==pt&tmp$state==st,]
          if(all(c("primary_tumor","liver_metastasis")%in%q$tissue)){
            a<-mean(q$score[q$tissue=="primary_tumor"])
            b<-mean(q$score[q$tissue=="liver_metastasis"])
            if(is.finite(a)&&is.finite(b)){
              ii<-ii+1L
              rows[[ii]]<-data.frame(patient=pt,state=st,delta=b-a)
            }
          }
        }
      }
      if(length(rows))do.call(rbind,rows)else data.frame()
    }
    
    p1_delta<-make_deltas(score_signature(p1_genes,p1_weights))
    p1_pat<-aggregate(delta~patient,p1_delta,mean)
    p1_effect<-median(p1_pat$delta)
    p1_pat_frac<-mean(p1_pat$delta>0)
    p1_state_frac<-mean(p1_delta$delta>0)
    
    null_rows<-vector("list",PARAM$matched_null_n)
    for(sid in seq_len(PARAM$matched_null_n)){
      q<-signatures[signatures$signature_id==sid,]
      dd<-make_deltas(score_signature(q$matched_gene,q$assigned_weight))
      pd<-aggregate(delta~patient,dd,mean)
      null_rows[[sid]]<-data.frame(
        signature_id=sid,
        median_patient_equal_weight_effect=median(pd$delta),
        mean_patient_equal_weight_effect=mean(pd$delta),
        n_patients=nrow(pd),n_positive_patients=sum(pd$delta>0),
        patient_positive_fraction=mean(pd$delta>0),
        n_state_pairs=nrow(dd),n_positive_state_pairs=sum(dd$delta>0),
        state_positive_fraction=mean(dd$delta>0)
      )
      if(sid%%500==0)cat("matched-null processed",sid,"\n")
    }
    null<-do.call(rbind,null_rows)
    write.csv(
      null,p(DIR$tables,"75_matched_random_signature_summary_5000.csv"),
      row.names=FALSE
    )
    
    eff<-null$median_patient_equal_weight_effect
    pf<-null$patient_positive_fraction
    sf<-null$state_positive_fraction
    result<-data.frame(
      n_random_signatures=nrow(null),
      matching_method="raw_scale_caliper_then_quantile_topK_without_replacement",
      assigned_weight_policy="target_specific_weight_preserved",
      matching_knn=PARAM$matched_null_knn,
      expression_caliper=expr_caliper,
      detection_caliper=detect_caliper,
      P1_median_patient_equal_weight_effect=p1_effect,
      P1_patient_positive_fraction=p1_pat_frac,
      P1_state_positive_fraction=p1_state_frac,
      P1_n_patients=nrow(p1_pat),P1_n_state_pairs=nrow(p1_delta),
      random_effect_ge_P1=sum(eff>=p1_effect),
      random_patient_fraction_ge_P1=sum(pf>=p1_pat_frac),
      random_state_fraction_ge_P1=sum(sf>=p1_state_frac),
      random_joint_effect_state_ge_P1=sum(eff>=p1_effect & sf>=p1_state_frac),
      random_joint_effect_patient_ge_P1=sum(eff>=p1_effect & pf>=p1_pat_frac)
    )
    result$empirical_p_effect <- plus_one_p(result$random_effect_ge_P1,nrow(null))
    result$empirical_p_patient_fraction <-
      plus_one_p(result$random_patient_fraction_ge_P1,nrow(null))
    result$empirical_p_state_fraction <-
      plus_one_p(result$random_state_fraction_ge_P1,nrow(null))
    result$empirical_p_joint_effect_state <-
      plus_one_p(result$random_joint_effect_state_ge_P1,nrow(null))
    result$empirical_p_joint_effect_patient <-
      plus_one_p(result$random_joint_effect_patient_ge_P1,nrow(null))
    write.csv(result,out_summary,row.names=FALSE)

    # Target-subset matched-null sensitivities.  These are independent
    # 5,000-signature draws generated with the same frozen matching rule,
    # while removing (i) genes with no in-caliper support at all and (ii)
    # genes whose in-caliper support is below the prespecified 0.8 threshold.
    # The full-P1 null above remains the primary result; these analyses audit
    # sensitivity to limited common support rather than redefine P1.
    frac_lookup <- setNames(balance$fraction_within_caliper, balance$target_gene)
    sensitivity_sets <- list(
      P1_no_extreme = which(is.finite(frac_lookup[p1_genes]) & frac_lookup[p1_genes] > 0),
      P1_well_matchable = which(is.finite(frac_lookup[p1_genes]) & frac_lookup[p1_genes] >= 0.8)
    )

    generate_matched_subset <- function(target_idx, seed) {
      target_idx <- as.integer(target_idx)
      n_target <- length(target_idx)
      if (n_target < 2L) return(data.frame())
      set.seed(seed)
      rows2 <- vector("list", PARAM$matched_null_n * n_target)
      rr2 <- 0L
      for (sid2 in seq_len(PARAM$matched_null_n)) {
        order_local <- sample(seq_len(n_target))
        used2 <- character()
        for (j_local in order_local) {
          j <- target_idx[j_local]
          mi <- matching_info[[j]]
          admissible <- mi$admissible
          available_admissible <- admissible[!cand$gene[admissible] %in% used2]
          if (length(available_admissible)) {
            ids2 <- available_admissible[order(mi$qdist[available_admissible])]
            status2 <- "within_caliper"
          } else {
            available_all <- which(!cand$gene %in% used2 & is.finite(mi$qdist))
            if (!length(available_all))
              stop("matched-null敏感性无法为目标基因分配候选基因：", p1_genes[j])
            ids2 <- available_all[order(mi$qdist[available_all])]
            status2 <- if (length(admissible))
              "caliper_blocked_by_uniqueness" else "extreme_target_no_caliper"
          }
          pool_size2 <- length(ids2)
          ids2 <- ids2[seq_len(min(PARAM$matched_null_knn, length(ids2)))]
          dd2 <- mi$qdist[ids2]
          temp2 <- max(median(dd2), 1e-3)
          prob2 <- exp(-(dd2 - min(dd2)) / temp2)
          prob2 <- prob2 / sum(prob2)
          pick2 <- ids2[sample.int(length(ids2), 1L, prob = prob2)]
          chosen2 <- cand$gene[pick2]
          rr2 <- rr2 + 1L
          rows2[[rr2]] <- data.frame(
            sensitivity_signature = TRUE,
            signature_id = sid2,
            target_rank = j_local,
            target_gene = p1_genes[j],
            matched_gene = chosen2,
            target_weight = p1_weights[j],
            assigned_weight = p1_weights[j],
            target_mean_log_expression = ps$mean_log_expression[j],
            target_detection_frequency = ps$detection_frequency[j],
            matched_mean_log_expression = cand$mean_log_expression[pick2],
            matched_detection_frequency = cand$detection_frequency[pick2],
            abs_expression_gap = abs(cand$mean_log_expression[pick2] - ps$mean_log_expression[j]),
            abs_detection_gap = abs(cand$detection_frequency[pick2] - ps$detection_frequency[j]),
            quantile_distance = mi$qdist[pick2],
            matching_status = status2,
            candidate_pool_size = pool_size2,
            caliper_pool_size = length(available_admissible),
            expression_caliper = expr_caliper,
            detection_caliper = detect_caliper,
            stringsAsFactors = FALSE
          )
          used2 <- c(used2, chosen2)
        }
      }
      do.call(rbind, rows2[seq_len(rr2)])
    }

    sensitivity_summary_rows <- list()
    sensitivity_balance_rows <- list()
    sensitivity_definition_rows <- list()
    for (sn in names(sensitivity_sets)) {
      idx2 <- sensitivity_sets[[sn]]
      genes2 <- p1_genes[idx2]
      sensitivity_definition_rows[[sn]] <- data.frame(
        sensitivity = sn,
        rule = if (sn == "P1_no_extreme")
          "Retain P1 genes with fraction_within_caliper > 0" else
          "Retain P1 genes with fraction_within_caliper >= 0.8",
        n_target_genes = length(idx2),
        target_genes = paste(genes2, collapse = ";"),
        draw_method = "Independent 5,000 matched signatures; target-specific weights preserved",
        stringsAsFactors = FALSE
      )
      if (length(idx2) < 2L) next
      sig2 <- generate_matched_subset(idx2, 20260815L + match(sn, names(sensitivity_sets)) * 1000L)
      bal2 <- do.call(rbind, lapply(split(sig2, sig2$target_gene), function(q2) {
        data.frame(
          sensitivity = sn,
          target_gene = q2$target_gene[[1]],
          n_draws = nrow(q2),
          n_within_caliper = sum(q2$matching_status == "within_caliper"),
          fraction_within_caliper = mean(q2$matching_status == "within_caliper"),
          n_caliper_blocked = sum(q2$matching_status == "caliper_blocked_by_uniqueness"),
          n_extreme_target = sum(q2$matching_status == "extreme_target_no_caliper"),
          median_abs_expression_gap = median(q2$abs_expression_gap, na.rm = TRUE),
          median_abs_detection_gap = median(q2$abs_detection_gap, na.rm = TRUE),
          median_quantile_distance = median(q2$quantile_distance, na.rm = TRUE),
          expression_caliper = expr_caliper,
          detection_caliper = detect_caliper,
          stringsAsFactors = FALSE
        )
      }))
      sensitivity_balance_rows[[sn]] <- bal2
      p1d2 <- make_deltas(score_signature(genes2, p1_weights[idx2]))
      p1pat2 <- aggregate(delta ~ patient, p1d2, mean)
      p1eff2 <- median(p1pat2$delta)
      p1pf2 <- mean(p1pat2$delta > 0)
      p1sf2 <- mean(p1d2$delta > 0)
      null2_rows <- vector("list", PARAM$matched_null_n)
      for (sid2 in seq_len(PARAM$matched_null_n)) {
        q2 <- sig2[sig2$signature_id == sid2, , drop = FALSE]
        dd3 <- make_deltas(score_signature(q2$matched_gene, q2$assigned_weight))
        pd3 <- aggregate(delta ~ patient, dd3, mean)
        null2_rows[[sid2]] <- data.frame(
          signature_id = sid2,
          median_patient_equal_weight_effect = median(pd3$delta),
          patient_positive_fraction = mean(pd3$delta > 0),
          state_positive_fraction = mean(dd3$delta > 0),
          stringsAsFactors = FALSE
        )
      }
      null2 <- do.call(rbind, null2_rows)
      ne2 <- sum(null2$median_patient_equal_weight_effect >= p1eff2)
      np2 <- sum(null2$patient_positive_fraction >= p1pf2)
      ns2 <- sum(null2$state_positive_fraction >= p1sf2)
      sensitivity_summary_rows[[sn]] <- data.frame(
        sensitivity = sn,
        n_target_genes = length(idx2),
        target_genes = paste(genes2, collapse = ";"),
        P1_median_patient_equal_weight_effect = p1eff2,
        P1_patient_positive_fraction = p1pf2,
        P1_state_positive_fraction = p1sf2,
        n_patients = nrow(p1pat2),
        n_state_pairs = nrow(p1d2),
        random_effect_ge_P1 = ne2,
        random_patient_fraction_ge_P1 = np2,
        random_state_fraction_ge_P1 = ns2,
        empirical_p_effect = plus_one_p(ne2, nrow(null2)),
        empirical_p_patient_fraction = plus_one_p(np2, nrow(null2)),
        empirical_p_state_fraction = plus_one_p(ns2, nrow(null2)),
        draw_method = "independent matched signatures; target-specific weights preserved",
        stringsAsFactors = FALSE
      )
      cat("matched-null sensitivity processed ", sn, " (", length(idx2), " genes)\n", sep = "")
    }
    write.csv(do.call(rbind, sensitivity_definition_rows), sensitivity_definition_file, row.names = FALSE)
    write.csv(do.call(rbind, sensitivity_summary_rows), sensitivity_summary_file, row.names = FALSE)
    if (length(sensitivity_balance_rows)) {
      write.csv(do.call(rbind, sensitivity_balance_rows), sensitivity_balance_file, row.names = FALSE)
    } else {
      write.csv(data.frame(), sensitivity_balance_file, row.names = FALSE)
    }
    write_session_info("15_matched_random_signature_null_5000")
    result
  })
}


############################
# Stage 16：非恶性肝脏 / 肝病环境边界
############################

stage16_nonmalignant_liver_boundary <- function(force = FALSE) {
  stage_log("16_nonmalignant_liver_boundary", {
    frozen_file <- p(DIR$cache_projection,"frozen_P1_v4.rds")
    assert_file(frozen_file)
    frozen<-must_load_step_cache(frozen_file, "frozen P1")
    p1<-frozen$P1_top15
    w<-data.frame(gene=p1$gene,weight=p1$weight_norm)
    no_core4 <- if (isTRUE(P1_MARKER_EXCLUSION_CONFIRMED))
      setdiff(w$gene, P1_MARKER_EXCLUDED_GENES) else character()
    
    results<-list()
    
    # GSE115469：作者处理表达矩阵+CellClusterType
    data_file<-p(GSE115469_DIR,"raw","GSE115469_Data.csv.gz")
    meta_file<-p(GSE115469_DIR,"raw","GSE115469_CellClusterType.txt.gz")
    if(file.exists(data_file)&&file.exists(meta_file)){
      expr<-read.csv(gzfile(data_file),check.names=FALSE)
      names(expr)[1]<-"gene";rownames(expr)<-expr$gene;expr$gene<-NULL
      mm<-read.delim(gzfile(meta_file),stringsAsFactors=FALSE)
      mm<-mm[mm$CellName%in%colnames(expr),]
      lib<-colSums(expr)
      norm<-log1p(sweep(expr,2,pmax(lib,1),"/")*1e4)
      score_mat<-function(gs){
        g<-intersect(gs,rownames(norm))
        if (!length(g)) return(rep(NA_real_, ncol(norm)))
        ww<-w$weight[match(g,w$gene)]
        ww<-ww/sum(ww);as.numeric(t(norm[g,,drop=FALSE])%*%ww)
      }
      sc1<-score_mat(w$gene);sc2<-score_mat(no_core4)
      names(sc1)<-colnames(norm);names(sc2)<-colnames(norm)
      mm$P1_full<-sc1[mm$CellName];mm$P1_no_core4<-sc2[mm$CellName]
      mm$myeloid<-grepl("Macrophage|Monocyte",mm$CellType,ignore.case=TRUE)
      z<-mm[mm$myeloid,]
      write.csv(z,p(DIR$tables,"76_GSE115469_frozen_projection_cells.csv"),row.names=FALSE)
      # When marker exclusion is still pending, P1_no_core4 is intentionally
      # all NA.  Do not pass that column to aggregate(), whose default
      # na.omit would otherwise discard every myeloid cell.
      if (length(no_core4)) {
        sm<-aggregate(
          cbind(P1_full,P1_no_core4)~Sample+CellType,z,mean
        )
      } else {
        sm<-aggregate(
          P1_full~Sample+CellType,z,mean
        )
        sm$P1_no_core4 <- NA_real_
      }
      nc<-aggregate(CellName~Sample+CellType,z,length);names(nc)[3]<-"n_cells"
      sm<-merge(sm,nc,by=c("Sample","CellType"))
      sm$environment<-"healthy_liver"
      write.csv(sm,p(DIR$tables,"77_GSE115469_frozen_projection_summary.csv"),row.names=FALSE)
      results$GSE115469<-sm
    } else {
      warning("GSE115469 raw files absent; Stage16 skips this dataset.")
    }
    
    # 通用10x肝病投影
    run_liver_10x<-function(root,study){
      rawdir<-p(root,"raw")
      mtxs<-if(dir.exists(rawdir))sort(list.files(rawdir,pattern="_matrix\\.mtx\\.gz$",full.names=TRUE))else character()
      rows<-list();rr<-0L
      fixed235<-c(
        GSM7494563_LF004="MAFLD",
        GSM7494564_LF005="MAFLD",
        GSM7494565_LF008="MAFLD"
      )
      for(mf in mtxs){
        stem<-sub("_matrix\\.mtx\\.gz$","",basename(mf))
        if(study=="GSE136103"&&!grepl("cd45\\+",stem,ignore.case=TRUE))next
        bf<-p(rawdir,paste0(stem,"_barcodes.tsv.gz"))
        ff<-p(rawdir,paste0(stem,"_features.tsv.gz"))
        if(!file.exists(ff))ff<-p(rawdir,paste0(stem,"_genes.tsv.gz"))
        if(!file.exists(bf)||!file.exists(ff))next
        X<-read_10x_triplet(mf,bf,ff);norm<-normalize_log1p_cp10k_sparse(X)
        my<-module_score(norm,MYELOID_MODULE);ep<-module_score(norm,EPITHELIAL_MODULE)
        lyz<-if("LYZ"%in%rownames(X))as.numeric(X["LYZ",]>0)>0 else rep(FALSE,ncol(X))
        inc<-my>ep&lyz
        if(study=="GSE235079"){
          if(!stem%in%names(fixed235))stop("未映射GSE235079 sample：",stem)
          env<-fixed235[[stem]]
        } else {
          env<-if(grepl("healthy",stem,ignore.case=TRUE))"healthy" else
            if(grepl("cirrh",stem,ignore.case=TRUE))"cirrhotic" else "unknown"
        }
        s1<-weighted_score(norm,w)
        s2<-if (length(no_core4)) weighted_score(norm,w[w$gene%in%no_core4,,drop=FALSE]) else rep(NA_real_, ncol(norm))
        rr<-rr+1L
        rows[[rr]]<-data.frame(
          study=study,sample=stem,environment=env,n_myeloid=sum(inc),
          mean_P1=if(any(inc))mean(s1[inc])else NA_real_,
          median_P1=if(any(inc))median(s1[inc])else NA_real_,
          mean_no_core4=if(any(inc))mean(s2[inc])else NA_real_,
          q90_P1=if(any(inc))quantile(s1[inc],.9)else NA_real_
        )
        rm(X,norm);invisible(gc())
      }
      if(length(rows))do.call(rbind,rows)else data.frame()
    }
    
    g136<-run_liver_10x(GSE136103_DIR,"GSE136103")
    g235<-run_liver_10x(GSE235079_DIR,"GSE235079")
    liver_disease<-rbind(g136,g235)
    if(nrow(liver_disease)){
      write.csv(
        liver_disease,
        p(DIR$tables,"78_liver_disease_frozen_projection_summary.csv"),
        row.names=FALSE
      )
      results$liver_disease<-liver_disease
    }
    
    # GSE178318 CRC/LM同一冻结score环境摘要
    my<-must_load_step_cache(
      p(DIR$cache_annotation,"fine_myeloid_dc_annotated.rds"),
      "fine myeloid annotation"
    )
    norm<-Seurat::GetAssayData(my,assay="RNA",layer="data")
    s1<-weighted_score(norm,w)
    s2<-if (length(no_core4)) weighted_score(norm,w[w$gene%in%no_core4,,drop=FALSE]) else rep(NA_real_, ncol(norm))
    dm<-data.frame(
      tissue=as.character(my$tissue),P1_full=s1,P1_no_core4=s2
    )
    dm<-dm[dm$tissue%in%c("CRC","LM"),]
    if (length(no_core4)) {
      cr<-aggregate(cbind(P1_full,P1_no_core4)~tissue,dm,mean)
      names(cr)[2:3]<-c("mean_P1","mean_no_core4")
    } else {
      cr<-aggregate(P1_full~tissue,dm,mean)
      names(cr)[2]<-"mean_P1"
      cr$mean_no_core4 <- NA_real_
    }
    write.csv(cr,p(DIR$tables,"79_GSE178318_environment_summary.csv"),row.names=FALSE)
    results$GSE178318<-cr
    
    # 环境梯度表：明确非恶性是sample-level，GSE178318是pooled tissue summary
    env_rows<-list();er<-0L
    # GSE115469 is a sample-level non-malignant liver reference.  The source
    # table is summarized by CellType, so collapse it back to one weighted
    # value per biological sample before adding it to the environment table.
    # This prevents the two macrophage labels from being treated as separate
    # biological replicates and keeps the cross-dataset interpretation
    # descriptive only.
    if (!is.null(results$GSE115469) && nrow(results$GSE115469)) {
      g115 <- results$GSE115469
      g115$n_cells <- as.numeric(g115$n_cells)
      g115$P1_full <- as.numeric(g115$P1_full)
      g115$P1_no_core4 <- as.numeric(g115$P1_no_core4)
      g115_sample <- do.call(rbind, lapply(split(g115, g115$Sample), function(q) {
        ok_full <- is.finite(q$P1_full) & is.finite(q$n_cells) & q$n_cells > 0
        ok_excl <- is.finite(q$P1_no_core4) & is.finite(q$n_cells) & q$n_cells > 0
        data.frame(
          sample = q$Sample[[1]],
          environment = "healthy_liver",
          mean_P1 = if (any(ok_full))
            stats::weighted.mean(q$P1_full[ok_full], q$n_cells[ok_full]) else NA_real_,
          mean_no_core4 = if (any(ok_excl))
            stats::weighted.mean(q$P1_no_core4[ok_excl], q$n_cells[ok_excl]) else NA_real_,
          stringsAsFactors = FALSE
        )
      }))
      g115_sample <- g115_sample[is.finite(g115_sample$mean_P1), , drop = FALSE]
      if (nrow(g115_sample)) {
        er <- er + 1L
        env_rows[[er]] <- data.frame(
          study = "GSE115469",
          environment = "healthy_liver",
          n_samples = nrow(g115_sample),
          aggregation_unit = "sample",
          median_mean_P1 = stats::median(g115_sample$mean_P1, na.rm = TRUE),
          iqr_low = stats::quantile(g115_sample$mean_P1, .25, na.rm = TRUE),
          iqr_high = stats::quantile(g115_sample$mean_P1, .75, na.rm = TRUE),
          median_no_core4 = if (any(is.finite(g115_sample$mean_no_core4)))
            stats::median(g115_sample$mean_no_core4, na.rm = TRUE) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
    if(nrow(liver_disease)){
      ld<-liver_disease[liver_disease$n_myeloid>=100,]
      for(en in unique(ld$environment)){
        q<-ld[ld$environment==en,]
        er<-er+1L
        env_rows[[er]]<-data.frame(
          study="nonmalignant_liver",environment=en,n_samples=nrow(q),
          aggregation_unit="sample",
          median_mean_P1=median(q$mean_P1),
          iqr_low=quantile(q$mean_P1,.25),iqr_high=quantile(q$mean_P1,.75),
          median_no_core4=median(q$mean_no_core4)
        )
      }
    }
    for(i in seq_len(nrow(cr))){
      er<-er+1L
      env_rows[[er]]<-data.frame(
        study="GSE178318",
        environment=if(cr$tissue[i]=="LM")"CRLM"else"CRC",
        n_samples=NA_integer_,aggregation_unit="pooled_tissue_summary",
        median_mean_P1=cr$mean_P1[i],
        iqr_low=NA_real_,iqr_high=NA_real_,
        median_no_core4=cr$mean_no_core4[i]
      )
    }
    env<-do.call(rbind,env_rows)
    write.csv(env,p(DIR$tables,"80_environment_gradient_summary.csv"),row.names=FALSE)
    
    # GSE168933：当前代码包无可执行的原始重建脚本，不能凭空发明。
    # 若用户提供预冻结reference-marker map，则只做透明的可重建入口。
    ref_status<-data.frame(
      dataset="GSE168933",
      executable_source_found_in_current_bundle=FALSE,
      reference_marker_file=GSE168933_REFERENCE_MARKER_FILE,
      marker_file_available=file.exists(GSE168933_REFERENCE_MARKER_FILE),
      action=if(file.exists(GSE168933_REFERENCE_MARKER_FILE))
        "reference marker file available; dedicated mapping can be added without using legacy result"
      else
        "not reconstructed: source method absent; no legacy result imported"
    )
    write.csv(ref_status,p(DIR$audit,"16_GSE168933_rebuild_status.csv"),row.names=FALSE)
    
    write_session_info("16_nonmalignant_liver_boundary")
    list(environment=env,GSE168933_status=ref_status)
  })
}


############################
# Stage 17：空间补充审计
# - P1-Pk residual correlations
# - GSE225857 marker-excluded P1
# - scRNA reference-projected state adjustment（探索性，不称deconvolution）
############################

stage17_spatial_supplementary <- function(force = FALSE) {
  stage_log("17_spatial_supplementary_audits", {
    check_packages(c("dbscan","FNN"))
    frozen<-must_load_step_cache(
      p(DIR$cache_projection,"frozen_P1_v4.rds"), "frozen P1"
    )
    create_spatial_manifest_template()
    man<-read.csv(SPATIAL_MANIFEST,stringsAsFactors=FALSE)
    if(any(!nzchar(man$matrix_path))||any(!nzchar(man$positions_path))){
      stop("spatial_manifest尚未填写完整。")
    }
    
    all15<-frozen$all_program_top15
    prog_weights<-lapply(unique(all15$program),function(pp){
      z<-all15[all15$program==pp&all15$rank<=15,]
      data.frame(gene=z$gene,weight=z$weight/sum(z$weight))
    })
    names(prog_weights)<-unique(all15$program)
    p1w<-prog_weights$P1
    p1_no_core4 <- if (isTRUE(P1_MARKER_EXCLUSION_CONFIRMED))
      p1w[!p1w$gene%in%P1_MARKER_EXCLUDED_GENES,] else NULL
    
    resid_rows<-list();rr<-0L
    marker_rows<-list();mr<-0L
    
    for(i in seq_len(nrow(man))){
      z<-read_spatial_section(man[i,])
      X<-z$X;pos<-z$pos;norm<-normalize_log1p_cp10k_sparse(X)
      mymod<-module_score(norm,SPATIAL_COVARIATE_MODULES$myeloid)
      epmod<-module_score(norm,SPATIAL_COVARIATE_MODULES$epithelial)
      livmod<-module_score(norm,SPATIAL_COVARIATE_MODULES$liver)
      lib<-Matrix::colSums(X)
      covdf<-data.frame(
        myeloid=mymod,epithelial=epmod,liver=livmod,logUMI=log1p(lib)
      )
      rmat<-matrix(
        NA_real_,nrow=ncol(X),ncol=length(prog_weights),
        dimnames=list(colnames(X),names(prog_weights))
      )
      for(pp in names(prog_weights)){
        sc<-weighted_score(norm,prog_weights[[pp]])
        fit<-lm(sc~myeloid+epithelial+liver+logUMI,data=covdf)
        rmat[,pp]<-residuals(fit)
      }
      for(pp in setdiff(colnames(rmat),"P1")){
        rr<-rr+1L
        resid_rows[[rr]]<-data.frame(
          dataset=man$dataset[i],section=man$section[i],
          comparison=paste0("P1-",pp),
          residual_spearman=safe_cor(rmat[,"P1"],rmat[,pp]),
          n_spots=nrow(rmat)
        )
      }
      
      # marker-excluded spatial sensitivity只在GSE225857
      if(man$dataset[i]=="GSE225857" && !is.null(p1_no_core4)){
        sc<-weighted_score(norm,p1_no_core4)
        fit<-lm(sc~myeloid+epithelial+liver+logUMI,data=covdf)
        re<-residuals(fit)
        coords<-as.matrix(pos[,c("pxl_col","pxl_row")])
        nb<-symmetrize_neighbors(make_knn_neighbors(coords, PARAM$spatial_knn_k))
        tst<-spatial_permutation_test(
          re,nb,PARAM$n_spatial_perm,
          PARAM$seed+1700000L+i
        )
        mr<-mr+1L
        marker_rows[[mr]]<-data.frame(
          dataset=man$dataset[i],section=man$section[i],
          score_version="P1_marker_excluded",
          excluded_genes=paste(P1_MARKER_EXCLUDED_GENES, collapse=";"),
          n_excluded_from_frozen_top15=length(intersect(P1_MARKER_EXCLUDED_GENES, p1w$gene)),
          graph=paste0("knn", PARAM$spatial_knn_k),
          Moran_I=unname(tst["I"]),empirical_p=unname(tst["p"]),
          permutations=PARAM$n_spatial_perm
        )
      }
      rm(X,norm);invisible(gc())
    }
    rc<-do.call(rbind,resid_rows)
    write.csv(rc,p(DIR$tables,"81_spatial_P1_Pk_residual_correlations.csv"),
              row.names=FALSE)
    me<-if(length(marker_rows))do.call(rbind,marker_rows)else data.frame()
    if(nrow(me))me$BH_q<-p.adjust(me$empirical_p,"BH")
    write.csv(
      me,p(DIR$tables,"82_GSE225857_marker_excluded_spatial_sensitivity.csv"),
      row.names=FALSE
    )
    if (is.null(p1_no_core4)) {
      write.csv(
        data.frame(
          marker_exclusion_confirmed = FALSE,
          note = "Marker-excluded spatial sensitivity was not run pending post-Stage-8 P1 top-15 review.",
          stringsAsFactors = FALSE
        ),
        p(DIR$audit, "17_marker_exclusion_pending.csv"), row.names = FALSE
      )
    } else if (file.exists(p(DIR$audit, "17_marker_exclusion_pending.csv"))) {
      # Remove the stale pre-confirmation sentinel so Stage 18 reports the
      # completed sensitivity rather than the old pending state.
      file.remove(p(DIR$audit, "17_marker_exclusion_pending.csv"))
    }
    
    # exploratory state-reference adjustment
    state_adjust_status<-data.frame(
      attempted=FALSE,completed=FALSE,n_sections=0,
      note="nnls package not checked"
    )
    if(requireNamespace("nnls",quietly=TRUE)){
      state_adjust_status$attempted<-TRUE
      my<-must_load_step_cache(
        p(DIR$cache_annotation,"fine_myeloid_dc_annotated.rds"),
        "fine myeloid annotation"
      )
      expr<-Seurat::GetAssayData(my,assay="RNA",layer="data")
      ref_modules<-list(
        `P1-state`=STATE_MODULES$P1_state,
        Inflammatory=STATE_MODULES$P2_inflammatory,
        `C1Q/MHC-II`=STATE_MODULES$P3_C1Q_MHCII,
        Mono=STATE_MODULES$Mono,
        Other=STATE_MODULES$Other
      )
      priority<-names(ref_modules)
      score<-sapply(ref_modules,function(gs){
        g<-intersect(gs,rownames(expr))
        if(length(g))Matrix::colMeans(expr[g,,drop=FALSE])else rep(0,ncol(expr))
      })
      labels<-priority[max.col(score,ties.method="first")]
      labels[apply(score,1,max)<=0]<-NA_character_
      deconv_genes<-setdiff(unique(unlist(ref_modules)),p1w$gene)
      deconv_genes<-intersect(deconv_genes,rownames(expr))
      ref<-sapply(priority,function(st){
        Matrix::rowMeans(expr[deconv_genes,which(labels==st),drop=FALSE])
      })
      rownames(ref)<-deconv_genes;colnames(ref)<-priority
      write.csv(
        data.frame(gene=rownames(ref),ref,check.names=FALSE),
        p(DIR$tables,"83_spatial_reference_state_profiles_P1_excluded.csv"),
        row.names=FALSE
      )
      
      state_sum<-list();ss<-0L
      for(i in seq_len(nrow(man))){
        z<-read_spatial_section(man[i,]);X<-z$X;pos<-z$pos
        norm<-normalize_log1p_cp10k_sparse(X)
        g<-intersect(rownames(ref),rownames(norm))
        if(length(g)<10L){
          warning("Section ",man$section[i]," state-reference genes <10; skipped.")
          next
        }
        A<-ref[g,,drop=FALSE]
        scale<-pmax(rowMeans(A),1e-4)
        As<-A/scale
        Y<-t(as.matrix(norm[g,,drop=FALSE]))/
          matrix(scale,nrow=ncol(norm),ncol=length(scale),byrow=TRUE)
        co<-matrix(0,nrow=nrow(Y),ncol=ncol(As))
        r2<-rep(NA_real_,nrow(Y))
        for(j in seq_len(nrow(Y))){
          fit<-nnls::nnls(As,Y[j,])
          co[j,]<-coef(fit)
          rn<-fit$deviance
          den<-sum((Y[j,]-mean(Y[j,]))^2)
          r2[j]<-if(den>0)1-rn/den else NA_real_
        }
        total<-rowSums(co)
        frac<-co/pmax(total,1e-12)
        colnames(frac)<-priority
        
        p1score<-weighted_score(norm,p1w)
        covdf<-data.frame(
          f1=frac[,1],f2=frac[,2],f3=frac[,3],f4=frac[,4],
          inferred_total=total,
          myeloid=module_score(norm,SPATIAL_COVARIATE_MODULES$myeloid),
          epithelial=module_score(norm,SPATIAL_COVARIATE_MODULES$epithelial),
          liver=module_score(norm,SPATIAL_COVARIATE_MODULES$liver),
          logUMI=log1p(Matrix::colSums(X))
        )
        # Other omitted as reference fraction
        fit<-lm(
          p1score~f1+f2+f3+f4+inferred_total+myeloid+epithelial+liver+logUMI,
          data=covdf
        )
        resid<-residuals(fit)
        coords<-as.matrix(pos[,c("pxl_col","pxl_row")])
        nb<-symmetrize_neighbors(make_knn_neighbors(coords, PARAM$spatial_knn_k))
        tst<-spatial_permutation_test(
          resid,nb,PARAM$n_spatial_perm,7403L+i
        )
        loc<-local_high_high(resid,nb,PARAM$n_spatial_perm,7503L+i)
        ss<-ss+1L
        state_sum[[ss]]<-data.frame(
          dataset=man$dataset[i],section=man$section[i],
          n_spots=length(resid),graph=paste0("knn", PARAM$spatial_knn_k),
          n_reference_genes=length(g),
          mean_projection_r2=mean(r2,na.rm=TRUE),
          median_projection_r2=median(r2,na.rm=TRUE),
          moran_I_state_adjusted=unname(tst["I"]),
          moran_p_state_adjusted=unname(tst["p"]),
          n_HH_spots=sum(loc$high_high),
          HH_fraction=mean(loc$high_high)
        )
        rm(X,norm);invisible(gc())
      }
      if(length(state_sum)){
        ssum<-do.call(rbind,state_sum)
        ssum$BH_q<-ave(
          ssum$moran_p_state_adjusted,ssum$dataset,
          FUN=function(x)p.adjust(x,"BH")
        )
        write.csv(
          ssum,p(DIR$tables,"84_spatial_state_adjusted_Moran_LISA_summary.csv"),
          row.names=FALSE
        )
        state_adjust_status$completed<-TRUE
        state_adjust_status$n_sections<-nrow(ssum)
        state_adjust_status$note <-
          "Exploratory NNLS reference projection; not probabilistic deconvolution."
      }
    } else {
      state_adjust_status$note <-
        "Skipped only because nnls package is not installed; no substitute deconvolution used."
      warning(state_adjust_status$note)
    }
    write.csv(
      state_adjust_status,p(DIR$audit,"17_spatial_state_adjustment_status.csv"),
      row.names=FALSE
    )
    write_session_info("17_spatial_supplementary_audits")
    list(residual_correlations=rc,marker_excluded=me,
         state_adjustment_status=state_adjust_status)
  })
}


############################
# Stage 18：最终审计 + 当前四张主图对应关系
############################

stage18_final_audit <- function(force = FALSE) {
  stage_log("18_final_analysis_audit", {
    essential<-c(
      p(DIR$cache_input,"03_GSE178318_qc_singlets.rds"),
      p(DIR$cache_annotation,"fine_myeloid_dc_annotated.rds"),
      p(DIR$cache_annotation,"fine_t_nk_annotated.rds"),
      p(DIR$cache_nmf_matrix,"GSE178318_myeloid_NMF_input.rds"),
      p(DIR$cache_nmf_alignment,"rank_selection.rds"),
      p(DIR$cache_projection,"frozen_P1_v4.rds")
    )
    audit<-data.frame(
      item=basename(essential),path=essential,exists=file.exists(essential),
      sha256=vapply(essential,sha256_file,character(1))
    )
    write.csv(audit,p(DIR$audit,"18_core_file_audit.csv"),row.names=FALSE)
    if(!all(audit$exists))stop("核心分析文件不完整，见18_core_file_audit.csv。")
    
    sel<-must_load_step_cache(
      p(DIR$cache_nmf_alignment,"rank_selection.rds"),
      "Stage 7 rank selection"
    )
    fit_completion_file<-p(DIR$audit,"06_NMF_fit_completion.csv")
    fit_completion<-if(file.exists(fit_completion_file))
      read.csv(fit_completion_file)else data.frame()
    my<-must_load_step_cache(
      p(DIR$cache_annotation,"fine_myeloid_dc_annotated.rds"),
      "fine myeloid annotation"
    )
    tnk<-must_load_step_cache(
      p(DIR$cache_annotation,"fine_t_nk_annotated.rds"),
      "fine T/NK annotation"
    )
    
    closure_files<-c(
      p(DIR$tables,"51_GSE245552_composition_activity_by_patient.csv"),
      p(DIR$tables,"53_GSE245552_five_state_by_patient.csv"),
      p(DIR$tables,"69_GSE178318_14fine_state_decomposition_by_patient.csv")
    )
    closure_rows<-list()
    for(f in closure_files){
      if(file.exists(f)){
        q<-read.csv(f)
        cc<-grep("closure_error",names(q),value=TRUE)
        closure_rows[[length(closure_rows)+1L]]<-data.frame(
          file=basename(f),n_rows=nrow(q),
          max_abs_closure_error=if(length(cc))max(abs(q[[cc[1]]]),na.rm=TRUE)else NA_real_
        )
      }
    }
    closure<-if(length(closure_rows))do.call(rbind,closure_rows)else data.frame()
    write.csv(closure,p(DIR$audit,"18_decomposition_closure_audit.csv"),row.names=FALSE)
    
    final<-data.frame(
      item=c(
        "pipeline_version","GSE178318_myeloid_cells","GSE178318_TNK_cells",
        "NMF_requested_fits","NMF_completed_fits","selected_rank",
        "rss_elbow_rank","rank_selection_used_external_outcome",
        "rank_selection_used_biological_interpretation"
      ),
      value=c(
        PIPELINE_VERSION,ncol(my),ncol(tnk),
        length(PARAM$nmf_ranks)*PARAM$nmf_runs_per_rank,
        if(nrow(fit_completion))sum(fit_completion$completed)else NA,
        sel$selected_rank,sel$elbow_rank,
        "FALSE","FALSE"
      )
    )
    write.csv(final,p(DIR$audit,"18_final_methodological_audit.csv"),row.names=FALSE)

    # Extended closure checks for analyses whose interpretation depends on a
    # secondary audit rather than on the core cache files alone.  In
    # particular, do not label the matched-null or marker-exclusion analyses
    # as complete merely because an older summary table is present.
    hc3_units_file <- p(DIR$audit, "14b_HC3_aggregate_unit_audit.csv")
    hc3_unit_audit <- if (file.exists(hc3_units_file))
      read.csv(hc3_units_file, stringsAsFactors = FALSE) else data.frame()
    hc3_units <- if (nrow(hc3_unit_audit))
      suppressWarnings(as.numeric(hc3_unit_audit$observed_aggregate_units[[1]]))
    else if (file.exists(p(DIR$tables, "99_P1_complexity_patient_state.csv")))
      nrow(read.csv(p(DIR$tables, "99_P1_complexity_patient_state.csv")))
    else NA_real_
    marker_pending_file <- p(DIR$audit, "17_marker_exclusion_pending.csv")
    marker_pending <- if (file.exists(marker_pending_file))
      read.csv(marker_pending_file, stringsAsFactors = FALSE) else data.frame()
    # Stage 15 can produce a valid 5,000-signature null while still having
    # target genes for which no candidate meets the prespecified calipers.
    # Keep existence and quality as separate audit fields so that a populated
    # table is never silently presented as a fully balanced matched null.
    stage15_balance_file <- p(DIR$tables, "74b_matched_random_signature_balance_audit.csv")
    stage15_balance <- if (file.exists(stage15_balance_file))
      read.csv(stage15_balance_file, stringsAsFactors = FALSE) else data.frame()
    stage15_min_within <- if (nrow(stage15_balance) &&
                              "fraction_within_caliper" %in% names(stage15_balance)) {
      vv <- suppressWarnings(as.numeric(stage15_balance$fraction_within_caliper))
      if (any(is.finite(vv))) min(vv[is.finite(vv)]) else NA_real_
    } else NA_real_
    stage15_extreme_genes <- if (nrow(stage15_balance) &&
                                 "n_extreme_target" %in% names(stage15_balance)) {
      ee <- suppressWarnings(as.numeric(stage15_balance$n_extreme_target))
      sum(is.finite(ee) & ee > 0)
    } else NA_real_
    stage15_balance_status <- if (!file.exists(stage15_balance_file)) {
      "MISSING"
    } else if (!nrow(stage15_balance)) {
      "EMPTY"
    } else if (!is.finite(stage15_min_within)) {
      "INCOMPLETE"
    } else if (stage15_min_within >= (1 - 1e-12) &&
               is.finite(stage15_extreme_genes) && stage15_extreme_genes == 0) {
      "PASS_ALL_TARGETS_WITHIN_CALIPER"
    } else {
      "BALANCE_LIMITED_EXTREME_OR_OUTSIDE_CALIPER"
    }
    extended <- data.frame(
      item = c(
        "stage14b_aggregate_units",
        "stage14b_hc3_aggregate_audit_exists",
        "stage14b_hc3_aggregate_audit_status",
        "stage14b_model_files_complete",
        "stage15_summary_exists",
        "stage15_balance_audit_exists",
        "stage15_matching_status",
        "stage15_min_fraction_within_caliper",
        "stage15_genes_with_extreme_target",
        "stage15_balance_status",
        "stage15_sensitivity_summary_exists",
        "stage15_sensitivity_balance_exists",
        "stage15_sensitivity_definitions_exists",
        "stage17_marker_exclusion_confirmed",
        "stage17_marker_excluded_output_exists",
        "stage17_state_adjusted_sections"
      ),
      value = c(
        hc3_units,
        file.exists(hc3_units_file),
        if (nrow(hc3_unit_audit)) hc3_unit_audit$status[[1]] else "MISSING",
        all(file.exists(c(
          p(DIR$tables, "100_HC3_full_covariate_model_loo_P1.csv"),
          p(DIR$tables, "100_HC3_complexity_only_model.csv"),
          p(DIR$tables, "100_HC3_models_all_terms.csv")
        ))),
        file.exists(p(DIR$tables, "73_matched_random_signature_null_results_5000.csv")),
        file.exists(stage15_balance_file),
        if (file.exists(stage15_balance_file))
          "strict_caliper_audit_available" else "legacy_summary_without_balance_audit",
        stage15_min_within,
        stage15_extreme_genes,
        stage15_balance_status,
        file.exists(p(DIR$tables, "73b_matched_random_signature_null_sensitivity_5000.csv")),
        file.exists(p(DIR$tables, "74c_matched_random_signature_sensitivity_balance_5000.csv")),
        file.exists(p(DIR$tables, "74d_matched_random_signature_sensitivity_definitions.csv")),
        if (nrow(marker_pending)) marker_pending$marker_exclusion_confirmed[[1]]
          else P1_MARKER_EXCLUSION_CONFIRMED,
        file.exists(p(DIR$tables, "82_GSE225857_marker_excluded_spatial_sensitivity.csv")) &&
          file.info(p(DIR$tables, "82_GSE225857_marker_excluded_spatial_sensitivity.csv"))$size > 10,
        if (file.exists(p(DIR$tables, "84_spatial_state_adjusted_Moran_LISA_summary.csv")))
          nrow(read.csv(p(DIR$tables, "84_spatial_state_adjusted_Moran_LISA_summary.csv"))) else 0L
      ),
      stringsAsFactors = FALSE
    )
    write.csv(extended, p(DIR$audit, "18_extended_closure_audit.csv"), row.names=FALSE)
    
    figure_map<-data.frame(
      figure=c("Figure 1","Figure 2","Figure 3","Figure 4"),
      scientific_role=c(
        "Discovery: cluster abundance versus cross-state P1 program activity",
        "Independent GSE245552 validation and composition-activity decomposition",
        "Non-malignant liver / liver-disease contextual boundary",
        "Residual spatial organization in GSE217414 and GSE225857"
      ),
      stages=c(
        "1-8,14","9-10,15","11,16","12-13,17"
      ),
      main_outputs=c(
        "NMF rank/stability; frozen P1; 14-state carriage",
        "state-specific paired effects; four/five-state decomposition; P1-Pk/null",
        "GSE231559 portability; GSE115469/GSE136103/GSE235079 context",
        "global Moran I; local High-High; P1-Pk residual; state-adjusted sensitivity"
      )
    )
    write.csv(figure_map,p(DIR$audit,"18_Figure1_4_stage_mapping.csv"),row.names=FALSE)
    
    # 自动生成一次所有results文件的SHA256 manifest
    rf<-list.files(p(BASE,"results"),recursive=TRUE,full.names=TRUE)
    rf<-rf[file.info(rf)$isdir%in%FALSE]
    base_norm <- normalizePath(BASE,winslash="/",mustWork=FALSE)
    rf_norm <- normalizePath(rf,winslash="/",mustWork=FALSE)
    prefix_n <- nchar(base_norm) + 2L
    rel <- ifelse(
      startsWith(rf_norm,paste0(base_norm,"/")),
      substring(rf_norm,prefix_n),
      rf_norm
    )
    manifest<-data.frame(
      relative_path=rel,
      size_bytes=file.info(rf)$size,
      sha256=vapply(rf,sha256_file,character(1))
    )
    write.csv(manifest,p(DIR$audit,"18_results_SHA256_manifest.csv"),row.names=FALSE)
    write_session_info("18_final_analysis_audit")
    list(final=final,figure_map=figure_map,closure=closure)
  })
}

############################
# 调度器：只运行你指定的Stage
############################

stage0_preflight <- function(force = FALSE) {
  files <- c(
    GSE178318_MATRIX = GSE178318_MATRIX,
    GSE178318_BARCODES = GSE178318_BARCODES,
    GSE178318_FEATURES = GSE178318_FEATURES
  )
  tab <- data.frame(
    input = names(files), path = unname(files),
    exists = file.exists(unname(files)),
    size_bytes = ifelse(file.exists(unname(files)), file.info(unname(files))$size, NA_real_),
    sha256 = vapply(unname(files), sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(tab, p(DIR$audit, "00_stage0_input_preflight.csv"), row.names = FALSE)
  message("Stage 0 preflight: ", sum(tab$exists), "/", nrow(tab), " GSE178318 input files found.")
  print(tab[, c("input", "path", "exists", "size_bytes")], row.names = FALSE)
  invisible(tab)
}

assert_freeze_state <- function(stage) {
  s <- as.character(stage)
  rank_ready <- cache_has_current_version(
    p(DIR$cache_nmf_alignment, "rank_selection.rds")
  )
  p1_ready <- cache_has_current_version(
    p(DIR$cache_projection, "frozen_P1_v4.rds")
  )
  audit <- data.frame(
    requested_stage = s,
    rank_selection_available = rank_ready,
    frozen_P1_available = p1_ready,
    marker_exclusion_confirmed = P1_MARKER_EXCLUSION_CONFIRMED,
    stringsAsFactors = FALSE
  )
  write.csv(audit, p(DIR$audit, "freeze_state_current.csv"), row.names = FALSE)
  if (s %in% c("8", "9", "10", "11", "12", "13", "14", "14b", "15", "16", "17", "18") &&
      !rank_ready) {
    stop("Stage ", s, " 需要先完成Stage 7 rank selection。")
  }
  if (s %in% c("9", "10", "11", "12", "13", "14", "14b", "15", "16", "17", "18") &&
      !p1_ready) {
    stop("Stage ", s, " 需要先完成Stage 8 P1 freeze。")
  }
  invisible(TRUE)
}

run_stage <- function(stage, force = FALSE) {
  message("Starting Stage ", as.character(stage), " ...")
  assert_freeze_state(stage)
  f <- switch(
    as.character(stage),
    "0" = stage0_preflight,
    "1" = stage1_import_qc,
    "2" = stage2_doublet,
    "3" = stage3_global_atlas,
    "4" = stage4_myeloid_annotation,
    "5" = stage5_prepare_nmf_matrix,
    "6" = stage6_run_multiseed_nmf,
    "7" = stage7_nmf_stability_rank_selection,
    "8" = stage8_freeze_p1,
    "9" = stage9_gse245552_projection,
    "10" = stage10_gse245552_statistics,
    "11" = stage11_gse231559_projection,
    "12" = stage12_spatial_projection,
    "13" = stage13_spatial_statistics,
    "14" = stage14_discovery_supplementary,
    "14b" = stage14b_complexity_hc3,
    "15" = stage15_matched_random_signature_null,
    "16" = stage16_nonmalignant_liver_boundary,
    "17" = stage17_spatial_supplementary,
    "18" = stage18_final_audit,
    stop("Stage必须为0–18，或14b。Stage 0仅做输入预检，不执行分析。")
  )
  f(force = force)
}

run_many <- function(stages, force = FALSE) {
  ans <- vector("list", length(stages))
  names(ans) <- paste0("stage", stages)
  for (i in seq_along(stages)) {
    cat("\n\n######## RUN STAGE", stages[i], "########\n")
    ans[[i]] <- run_stage(stages[i], force = force)
  }
  invisible(ans)
}


################################################################################
# 推荐首次运行顺序
#
# source("scripts/analysis/CRLM_methodologically_upgraded_rerun_v1.R")
#
# 发现队列 / standardized rerun：
# run_stage(1)      # raw + frozen QC
# run_stage(2)      # 15 libraries scDblFinder
# run_stage(3)      # global atlas + patient/sample/no-batch + PC sensitivity
# run_stage(4)      # lineage reclustering + marker-module annotations (legacy maps only for concordance)
# run_stage(5)      # patient-level abundance/T-cell audit + 1500-HVG NMF matrix
#
# 方法学修复核心：
# run_stage(6)      # 最耗时：rank4–8 × 30 full-matrix Brunet fits；断点续跑
# run_stage(7)      # LSAP + factor recurrence + consensus + mechanical rank selection
# run_stage(8)      # rank确定后按发现队列anchor module重新选择并冻结P1_v4
#
# 外部scRNA：
# run_stage(9)      # GSE245552 frozen projection + all-program/overlap sensitivity
# run_stage(10)     # exploratory Wilcoxon/BH + 4/5-state decomposition + threshold/LOPO/unified
# run_stage(11)     # GSE231559 directional portability
#
# 空间主分析：
# run_stage(12)     # GSE217414 + GSE225857 frozen spatial projection/residualization
# run_stage(13)     # global Moran + k=6 sensitivity + local High-High
################################################################################
