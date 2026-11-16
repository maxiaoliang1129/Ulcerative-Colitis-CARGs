# Final model specification — UC-Dx-2gene

## Model

A logistic regression classifier for ulcerative colitis (UC) using two genes:

```
logit P(UC) = -3.900392 + 0.417981 × CCL11 + 0.409387 × CLDN1
```

| Term | Coefficient | Source |
|---|---:|---|
| Intercept | -3.900392 | fitted GLM on training cohort |
| CCL11 | 0.417981 | fitted GLM |
| CLDN1 | 0.409387 | fitted GLM |

* Training cohort: GSE193677, 1,320 biopsies (865 UC / 455 control), 660 patients.
* Decision threshold: 0.639643 (Youden optimum in the training cohort).
* Machine-readable copy: `model/final_model_spec.json`.
* Reference implementation: `R/00_config.R::predict_uc()` / `R/05_model.R::apply_frozen_model()`.

## Note on the manuscript equation

The equation in the original notes (`公式.txt`, archived in `docs/legacy_scripts/`)
states `-7.735 + 0.416 × CCL11 + 0.199 × CLDN1`; those are the **LASSO**
coefficients of the larger model, not the coefficients of the fitted two-gene
GLM. This repository uses the fitted GLM, which is the only specification
consistent with the reported Youden cut-off of 0.6396. The authors should
confirm the equation in the manuscript text.

## Cut-off transfer

The 0.6396 cut-off is training-derived. Validation cohorts are on different
expression scales; `05_validation.R` therefore also reports logistic
recalibration (slope/intercept, Brier before/after) so that cohort-specific
thresholds can be reported instead of transferring the training cut-off blindly.
