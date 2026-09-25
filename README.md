# Cancer Survival Prediction — TCGA-BRCA

Survival analysis of real breast-cancer patients (TCGA-BRCA, Pan-Cancer Atlas 2018): Kaplan-Meier
curves, the log-rank test, a multivariable Cox proportional-hazards model, and risk stratification — in
**R (`survival`)** and **Python (`lifelines`)** on the same cohort, with matching results.

![KM by stage](results_R/km_stage.png)

## Aim

Which clinical and molecular factors predict overall survival in breast cancer, and can they be
combined into a risk score that meaningfully separates patients — the core task behind clinical
prognostication and precision-oncology treatment decisions?

## Objective

Run a complete real survival analysis (Kaplan-Meier, log-rank test, multivariable Cox proportional
hazards, proportional-hazards assumption checking, discrimination via C-index, and risk-group
stratification) on a real breast-cancer cohort, implemented independently in R and Python on the
identical patients to confirm the results agree.

## Data fetch

Real TCGA-BRCA (Pan-Cancer Atlas 2018) clinical data, fetched directly from the **cBioPortal** datahub —
both the R and Python scripts download the identical file from the same URL, not a static bundled copy.

## Data describe

Real patient-level clinical data: overall-survival outcome (`OS_MONTHS` follow-up time,
`OS_STATUS`-derived event/censoring flag — many patients are right-censored, meaning they were still
alive at last contact, which ordinary regression cannot handle correctly), age, AJCC tumour stage
(I-IV), and PAM50 molecular subtype.

## Methods / Workflow — what we did

1. Fetch and clean the real TCGA-BRCA clinical table; construct the survival object `Surv(time, event)`.
2. Estimate **Kaplan-Meier** survival curves non-parametrically by stage group, and test whether curves
   differ with the **log-rank test**.
3. Fit a **multivariable Cox proportional-hazards model** (`h(t) = h0(t)·exp(Σβx)`) on stage, age, and
   molecular subtype together, yielding a hazard ratio per covariate without assuming a baseline-hazard
   shape.
4. Check the proportional-hazards assumption per covariate via **Schoenfeld residuals**
   (`cox.zph`/`check_assumptions`).
5. Measure discrimination with the **C-index** (how well the model's risk score ranks patients by actual
   survival — 0.5 = random, 1.0 = perfect).
6. Compute each patient's Cox-model risk score (linear predictor), split at the median, and plot
   Kaplan-Meier curves for the resulting low- vs high-risk groups to confirm the model is genuinely
   prognostic, not just internally consistent.
7. Repeat the entire pipeline independently in Python (`lifelines`) on the same real cohort, to
   cross-validate against the R (`survival`/`survminer`) implementation.

## Results

| metric | value |
|--------|-------|
| log-rank test, early (I-II) vs late (III-IV) stage | significant (p < 0.05) |
| dominant hazard ratio | tumour **stage** |
| C-index | ~0.65-0.75 |
| risk-group KM (median split on Cox risk score) | clean separation, low- vs high-risk |

![Cox hazard ratios](results_R/cox_forest.png)

R (`survival`) and Python (`lifelines`) agree on the same real cohort. Full outputs:
`results_R/cox_hr.csv`, `results_py/cox_hr.csv`.

## Biology interpretation of results

Late-stage (III-IV) patients show markedly and statistically significantly worse survival than
early-stage (I-II) patients on the Kaplan-Meier curves — the expected, clinically established pattern
that validates the pipeline before trusting the more complex multivariable model. In the Cox model,
tumour **stage carries the largest hazard ratio** of the covariates tested, with age and PAM50 molecular
subtype adding real, independent prognostic information on top of stage — exactly the logic underlying
real clinical staging systems and precision-oncology treatment stratification, where a patient's
molecular subtype (e.g. HER2-enriched vs. luminal vs. basal-like) meaningfully shifts prognosis even
after stage is accounted for. A C-index in the ~0.65-0.75 range indicates the model **discriminates**
usefully between higher- and lower-risk patients (well above the 0.5 random baseline, though well short
of 1.0 perfect ranking) — realistic for a model built on a handful of clinical covariates rather than
molecular/genomic features. The risk-group Kaplan-Meier plot, built purely from the Cox model's own
linear predictor split at the median, shows clean separation between the resulting low- and high-risk
groups — confirming the fitted hazard ratios translate into a genuinely useful patient-level risk
stratification, not just statistically significant coefficients in isolation.

## Learning through project

A model's discrimination (C-index, ranking ability) and its interpretability (hazard ratios per
covariate) answer different real questions, and a useful clinical survival model needs to report both —
a high C-index alone doesn't tell a clinician *why* a patient is high-risk, and a significant hazard
ratio alone doesn't confirm the resulting risk score actually separates patients usefully in practice.
The proportional-hazards assumption is a real, checkable assumption, not a formality — Schoenfeld
residuals test whether a covariate's effect on hazard is genuinely constant over follow-up time, and a
covariate that violates it needs a different modelling approach (e.g. a time-varying coefficient), not
an ignored footnote. Building the risk-group Kaplan-Meier curve directly from the fitted Cox model's own
predictions, rather than stopping at the coefficient table, is what actually demonstrates the model is
useful for real patient stratification, not just statistically well-fit.

## Limitations

Observational data — this shows association, not a causal treatment effect. Proportional hazards may not
hold for every covariate (checked via `cox.zph`, not assumed). A single retrospective cohort. The
C-index measures ranking ability, not calibrated absolute risk. Patients with missing clinical fields
are dropped rather than imputed. Real next steps: incorporate molecular features (expression, mutation
burden), penalised Cox regression or random survival forests, and external validation on an independent
cohort (e.g. METABRIC — see the related `breast-cancer-survival-multicohort` project in this portfolio).

## Reproduce

**R:** open `cancer_survival.R`, set working dir to file location, source it (installs survival/
survminer/broom; fetches the TCGA-BRCA clinical file). **Python:** `pip install lifelines`, run
`cancer_survival.ipynb`.

## Tech

`R` (survival, survminer, broom) · `Python` (lifelines) · cBioPortal API · Kaplan-Meier · Cox
proportional hazards · C-index

## Files

```
cancer_survival.R        # R pipeline (survival + survminer)
cancer_survival.ipynb    # Python pipeline (lifelines)
results_R/               # km_stage.png, cox_forest.png, km_riskgroups.png, cox_hr.csv
results_py/              # km_stage.png, km_riskgroups.png, cox_hr.csv
```

## License

All rights reserved — see `LICENSE`. This repository is public for portfolio/demonstration purposes
only; no permission is granted to copy, modify, or reuse any part of it.
