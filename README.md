# CRLM R code release

This release contains the R analysis pipeline, the five main-figure scripts,
and the supplementary-figure and sensitivity scripts for the frozen
manuscript. It supports two distinct workflows. Choose the workflow that
matches your goal; rendering figures does not require rerunning the full
analysis.

| Goal | Workflow | Required inputs | Output |
|---|---|---|---|
| Reproduce manuscript figures from the frozen results | A | This release, the frozen results package, and the primary rerun output used by Figure 1 | Figure files and figure-specific source-data audits |
| Rerun analysis from GEO inputs | B | This release and the required external GEO files | Recomputed analysis tables, caches, and audit outputs |

Workflow A reproduces figures from the manuscript's frozen inputs. Workflow B
recomputes analysis outputs. Rerun outputs are not automatically substituted
for the frozen tables, figures, or manuscript values; audit them before making
any such replacement.

## Workflow A — reproduce the manuscript figures

### Required folders

1. The unzipped code-release directory (this directory).
2. The final frozen result package, including `tables/` and
   `Figure2_source_data_addendum_20260912/`.
3. The primary rerun output directory, including `results/tables/`. Figure 1
   reads discovery tables from this location; some supplementary scripts also
   require its cache. Figures 2–5 primarily use the frozen result package.
4. A separate writable directory for newly rendered figures.

Raw GEO matrices are not needed when all required frozen and rerun artifacts
are already available. Supplementary scripts can have additional inputs; see
“Supplementary and sensitivity analyses” and each script's input checks.

### Configure paths

Set the paths explicitly before rendering. These figure-workflow variables are
read by `config.R`:

```r
Sys.setenv(
  CRLM_CODE_RELEASE_ROOT = "<path to this unzipped release>",
  CRLM_RESULTS_ROOT      = "<primary rerun output folder>",
  CRLM_FREEZE_DIR        = "<final frozen result package>",
  CRLM_OUTPUT_ROOT       = "<new writable figure-output folder>"
)

source(file.path(Sys.getenv("CRLM_CODE_RELEASE_ROOT"), "config.R"))
```

`CRLM_RESULTS_ROOT` must contain `results/tables/`. `CRLM_FREEZE_DIR` must
contain the final `tables/` and the Figure 2 source-data addendum. Keep
`CRLM_OUTPUT_ROOT` outside both input folders so new exports cannot overwrite
the frozen source data. `CRLM_PRIMARY_RERUN_ROOT` defaults to
`CRLM_RESULTS_ROOT`; set it separately only if your supplementary inputs are
stored elsewhere.

### Render the five main figures

Run these scripts from R or RStudio after configuring paths:

```r
source(file.path(CRLM_CODE_RELEASE_ROOT, "scripts/figures/Figure1_CRLM_P1_actual.R"))
source(file.path(CRLM_CODE_RELEASE_ROOT, "scripts/figures/Figure2_external_scRNA_portability_v7.R"))
source(file.path(CRLM_CODE_RELEASE_ROOT, "scripts/figures/Figure3_CRLM_composition_activity.R"))
source(file.path(CRLM_CODE_RELEASE_ROOT, "scripts/figures/Figure4_CRLM_cross_state_SPP1.R"))
source(file.path(CRLM_CODE_RELEASE_ROOT, "scripts/figures/Figure5_spatial_organization_R.R"))
```

| Figure | Script | Primary inputs |
|---|---|---|
| 1 | `scripts/figures/Figure1_CRLM_P1_actual.R` | Primary rerun tables and frozen results |
| 2 | `scripts/figures/Figure2_external_scRNA_portability_v7.R` | Frozen package and Figure 2 addendum |
| 3 | `scripts/figures/Figure3_CRLM_composition_activity.R` | Frozen result package |
| 4 | `scripts/figures/Figure4_CRLM_cross_state_SPP1.R` | Frozen result package |
| 5 | `scripts/figures/Figure5_spatial_organization_R.R` | Frozen spatial addendum |

The figure scripts are plotting-only: they read source tables and do not
refit the scientific models. Each writes figure exports and a figure-specific
source-data audit beneath `CRLM_OUTPUT_ROOT`. If a required input is absent,
the script stops and reports the unresolved path rather than substituting a
different file.

## Workflow B — rerun the analysis from GEO inputs

Use this workflow only when you intend to recompute the analysis. The code
release does not contain FASTQ files or GEO raw matrices.

### Configure the rerun

The analysis pipeline uses a separate set of path variables; these are not
the figure-workflow variables above:

| Variable | Purpose | Default |
|---|---|---|
| `CRLM_RERUN_DIR` | Project root for the rerun | Directory containing the pipeline script |
| `CRLM_OUTPUT_DIR` | New analysis-output directory | `<CRLM_RERUN_DIR>/outputs/<pipeline version>` |
| `CRLM_RAW_DIR` | Parent directory containing GEO study folders | `<CRLM_RERUN_DIR>/data/GEO_raw` |
| `CRLM_GSE178318_DIR` | GSE178318 input folder | `<CRLM_RAW_DIR>/GSE178318` |
| `CRLM_GSE245552_DIR` … `CRLM_GSE168933_DIR` | Optional per-study overrides | Corresponding folder beneath `CRLM_RAW_DIR` |

Example:

```r
Sys.setenv(
  CRLM_CODE_RELEASE_ROOT = "<path to this unzipped release>",
  CRLM_RERUN_DIR = "<project/input root>",
  CRLM_OUTPUT_DIR = "<new analysis-output directory>",
  CRLM_RAW_DIR = "<directory containing GEO study folders>"
)

source(file.path(CRLM_CODE_RELEASE_ROOT,
                 "scripts/analysis/CRLM_methodologically_upgraded_rerun_v1.R"))
run_stage(0)  # input/path preflight only; no analysis stage is run
```

The GSE168933 mapping stage additionally expects a frozen reference-marker
map CSV based on that study's reference labels; the pipeline does not infer
this mapping from P1. Stage 0 checks required inputs and paths. Resolve any
reported missing files before starting analysis.

After preflight passes, run stages individually in the recommended order at
the end of the pipeline script:

| Stage(s) | Analysis block |
|---|---|
| 0 | Input/path preflight only |
| 1–2 | Discovery import/QC and per-library doublet detection |
| 3–5 | Global atlas, lineage subclustering/state maps, patient summaries and NMF input |
| 6–8 | Multi-seed NMF, rank/stability selection, P1 correspondence and freeze |
| 9–11 | External single-cell projection and directional portability |
| 12–13 | Spatial projection/residualization, global Moran and local High–High statistics |
| 14, 14b–18 | Supplementary robustness, complexity, matched-null, non-malignant, spatial and final audits; run when required |

Stage 6 (rank 4–8, repeated full-matrix NMF fits) is the most computationally
intensive and supports checkpoint recovery. The pipeline does not run all
stages when sourced. Do not launch every stage blindly: later audits can
require additional inputs and should be run only when their outputs are
needed.

## Supplementary and sensitivity analyses

- `scripts/analysis/CRLM_figure4A_state_specific_lmm.R` reproduces the
  prespecified state-specific LM−PT contrasts from the frozen patient ×
  tissue × state aggregate table, using `lmerTest`, `emmeans`, Satterthwaite
  inference, and BH correction. It writes a new table and does not modify the
  frozen results.
- `scripts/sensitivity/CRLM_5gene_SPP1TAM_sensitivity.R` is an alternative
  five-gene comparator branch. It requires the primary rerun output, frozen
  results, and raw GSE245552 input. It does not replace the frozen eight-gene
  comparator or main-figure data.
- `scripts/supplementary/build_supplementary_table_S13.R` formats the
  completed five-gene sensitivity summary as Supplementary Table S13.
- The NMF-depth scripts reproduce the prespecified 1,000- and rank-6
  2,000-iteration sensitivity audits; they do not overwrite the primary NMF
  fits.
- Supplementary Figures S1–S6 are in `scripts/supplementary/`. Some require
  rerun caches in addition to the frozen result package; check each script's
  input-file checks before running it.

The frozen SPP1-comparator residual analyses used 164 complete patient ×
tissue × state aggregate units. For each comparator, `P1_full ~ comparator +
state` was fitted at the aggregate-unit level; tissue was retained as the
LM–PT contrast label rather than included as a residual covariate. Residuals
were paired within patient and state and summarized at the patient level.
Released summaries use percentile bootstrap intervals (2.5th and 97.5th
percentiles of 10,000 patient resamples) and 10,000 patient-level Monte Carlo
sign-flip resamples with a plus-one empirical P-value correction.

## Provenance and naming

The release includes `code_release_manifest.csv`, pipeline and parameter
freeze records, session information, and `CRLM_CODE_RELEASE_SCHEMA_MAP.csv`.
Some frozen CSV columns retain historical identifiers such as
`SPP1_TAM_published_2025` to preserve existing SHA256 records. The
reader-facing term is “literature-derived SPP1-TAM marker score”; the mapping
is documented in the schema map. Do not rename frozen source files without
regenerating and re-auditing their manifests.

R package versions and session information are archived under `provenance/`.
Raw GEO data and other external inputs are not included in this release. See
`provenance/FREEZE_README.txt` for frozen-input and provenance details.
