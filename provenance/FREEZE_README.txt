Current frozen release provenance

This code package supports the synchronized manuscript and Supplementary
Information package dated 2026-09-18. Scientific analysis and figure
generation are R-only. The scripts read frozen source tables unless a script
explicitly documents a reproducible sensitivity branch.

The primary state-specific inference uses a REML linear mixed-effects model
fit with lmerTest::lmer and five state-specific LM-PT contrasts extracted with
emmeans::emmeans, followed by Satterthwaite inference and BH correction. The
Wilcoxon branch is exploratory/sensitivity only and is not the primary Figure
4A workflow.

The patient-level residual comparator model was fitted on 164 complete
patient × tissue × state aggregate units as `P1_full ~ comparator + state`.
Tissue was retained as the LM-PT contrast label rather than included as a
residual covariate because it defines the contrast of interest. Residuals
were paired within patient and state and summarized at the patient level.
Residual intervals use the 2.5th and 97.5th percentiles of 10,000 patient
bootstrap resamples; directionality uses 10,000 patient-level Monte Carlo
sign-flip resamples with a plus-one empirical P-value correction.

The literature-derived SPP1-TAM marker score is the primary eight-gene
comparator; the five-gene macrophage-specific definition is an alternative
sensitivity analysis. The prespecified factor-selection audit is distributed
as `provenance/P1_selection_audit.csv`.

The discovery input is a GEO-supplied 10x-format count matrix, not FASTQ or
raw sequencing data. The public manuscript and Supplementary Information
files are maintained outside this code directory. Their exact file hashes are
recorded in the corresponding package manifests.
