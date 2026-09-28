# basics-chloride-analysis

Analysis code for:

> Zampieri FG, Maia IS, Machado FR, Cavalcanti AB. **Chloride load and serum chloride in critically ill adults: a randomization-anchored analysis of the BaSICS trial.** *Annals of Intensive Care* (under review).

This is a secondary analysis of BaSICS (Balanced Solutions in Intensive Care Study; ClinicalTrials.gov NCT02875873), a 2 × 2 factorial randomized trial of Plasma-Lyte 148 versus 0.9% saline and slow versus fast fluid-challenge infusion rate in 75 Brazilian ICUs ([Zampieri et al., JAMA 2021](https://doi.org/10.1001/jama.2021.11684)).

## Data availability

**The patient-level data are not included in this repository.** The BaSICS dataset is not publicly available. It is governed by the BaSICS/BRICNet data-sharing policy and may be made available on reasonable request to the corresponding author, subject to appropriate ethics and data-transfer approvals. The code expects the file `raw_data/Dados_BASICS_MainPaper.RData`, which contains the data frame `dados`.

All data files, derived data, model fits, figures and outputs are excluded by `.gitignore`.

## Requirements

R 4.6 with `brms` (Stan backend via `rstan`), `posterior`, and `lme4`; `lubridate` must be installed so the date columns in the data load correctly. The exact versions used are in [`docs/session_info.txt`](docs/session_info.txt). All Bayesian models use fixed random seeds.

## Layout

| Path | Contents |
|---|---|
| `R/project_paths.R` | Resolves the project root; run scripts from the repository root or set `BASICS_REVISION_ROOT`. |
| `R/exposure_construction.R` | Post-randomization study-fluid chloride exposure (days 1-3), including the single documented data-entry correction. |
| `R/basics_chloride_analysis.R` | Analytic cohort, key numbers, exposure validation; writes the derived-data cache. |
| `R/basics_chloride_brms.R` | Bayesian hierarchical longitudinal model (`arm × day × eGFR`). |
| `R/analysis1_sensitivity_lme4.R` | Frequentist (lme4) robustness checks of the longitudinal model: participants measured on all three days, and baseline creatinine in place of eGFR. |
| `R/basics_chloride_varpart.R` | Bayesian variance partition of serum chloride. |
| `R/reviewer_sensitivity_16_17.R` | Variance-partition sensitivity analyses (kidney replacement therapy exclusion; change from baseline). |
| `R/reviewer2_analyses.R` | Descriptive intervals, day-2 rationale counts, decomposition of the measured-baseline component, exploratory correlates of patient-specific intercepts, and the randomized-arm sensitivity analysis. |
| `R/table1_final.R` | Table 1. |
| `R/basics_chloride_figures.R` | Figures 1-4, randomized normalization, and inverse-probability-weighted within-saline association. |
| `reproduce_current_submission.R` | Reruns cohort, Table 1 and figures from saved fits and checks the key reported numbers. |

## Running

With the data in `raw_data/`, from the repository root:

```bash
Rscript R/basics_chloride_analysis.R      # cohort and derived-data cache
Rscript R/basics_chloride_brms.R          # longitudinal model (slow)
Rscript R/analysis1_sensitivity_lme4.R    # longitudinal-model robustness checks
Rscript R/basics_chloride_varpart.R       # variance partition (slow)
Rscript R/reviewer_sensitivity_16_17.R    # variance-partition sensitivity analyses (slow)
Rscript R/reviewer2_analyses.R            # additional analyses (slow: one refit)
Rscript reproduce_current_submission.R    # Table 1, Figures 1-4, and checks
```

Outputs are written to `reproduction/`, `figures/`, `derived_data/` and `model_fits/`.

## Licence

MIT; see [LICENSE](LICENSE). If you use this code, please cite the article and the archived software release (see [CITATION.cff](CITATION.cff)).
