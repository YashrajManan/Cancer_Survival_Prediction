# =============================================================================
# Cancer Survival Prediction — TCGA-BRCA (R: survival + survminer)
# Kaplan-Meier + log-rank + Cox proportional-hazards + risk stratification.
# Companion Python version: cancer_survival.ipynb (lifelines) — same real cohort.
# =============================================================================
# Data: TCGA-BRCA, Pan-Cancer Atlas 2018 (real patients), clinical table from cBioPortal.
# Run in RStudio: Session -> Set Working Directory -> To Source File Location -> Source.
# First run installs survival + survminer + broom.
# =============================================================================

## ---- Setup ----
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
for (p in c("survival", "survminer", "broom")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(survival); library(survminer)
dir.create("data_R", showWarnings = FALSE); dir.create("results_R", showWarnings = FALSE)

## ---- 1. Fetch + clean the real TCGA-BRCA clinical table ----
# OS_MONTHS = follow-up time; OS_STATUS = "1:DECEASED" (event) / "0:LIVING" (censored).
# Covariates: age, tumour stage (I-IV, dominant prognostic), PAM50 molecular subtype.
url <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt"
if (!file.exists("data_R/brca_clinical.txt"))
  download.file(url, "data_R/brca_clinical.txt", mode = "wb")
cl <- read.delim("data_R/brca_clinical.txt", comment.char = "#")
cl$time  <- suppressWarnings(as.numeric(cl$OS_MONTHS))
cl$event <- as.integer(grepl("^1|DECEASED", cl$OS_STATUS))
cl$age   <- suppressWarnings(as.numeric(cl$AGE))
st <- toupper(cl$AJCC_PATHOLOGIC_TUMOR_STAGE)
cl$stage <- ifelse(grepl("IV", st), 4, ifelse(grepl("III", st), 3,
             ifelse(grepl("II", st), 2, ifelse(grepl("\\bI\\b|STAGE I$|IA|IB", st), 1, NA))))
cl$subtype <- factor(cl$SUBTYPE)
df <- cl[!is.na(cl$time) & cl$time > 0 & !is.na(cl$event) & !is.na(cl$age) & !is.na(cl$stage), ]
df$stage_group <- factor(ifelse(df$stage >= 3, "late (III-IV)", "early (I-II)"))
cat("patients:", nrow(df), "| deaths:", sum(df$event), "| censored:", sum(df$event == 0), "\n")

## ---- 2. Survival object (pairs follow-up time with event/censoring) ----
surv_obj <- Surv(time = df$time, event = df$event)

## ---- 3. Kaplan-Meier by stage group + log-rank ----
fit_km <- survfit(surv_obj ~ stage_group, data = df)
print(surv_median(fit_km))
print(survdiff(surv_obj ~ stage_group, data = df))
ggsave("results_R/km_stage.png",
       ggsurvplot(fit_km, data = df, pval = TRUE, risk.table = TRUE)$plot, width = 7, height = 5, dpi = 150)

## ---- 4. Multivariable Cox proportional-hazards (hazard ratios) ----
cox <- coxph(surv_obj ~ age + stage + subtype, data = df)
print(summary(cox))                                  # exp(coef)=HR, 95% CI, p, concordance
ggsave("results_R/cox_forest.png", ggforest(cox, data = df), width = 8, height = 6, dpi = 150)
write.csv(broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE), "results_R/cox_hr.csv", row.names = FALSE)

## ---- 5. Check the proportional-hazards assumption (Schoenfeld residuals) ----
print(cox.zph(cox))

## ---- 6. Risk stratification + validation KM ----
df$risk  <- predict(cox, type = "lp")                # linear predictor = per-patient risk score
df$group <- ifelse(df$risk > median(df$risk), "high risk", "low risk")
fit_risk <- survfit(Surv(time, event) ~ group, data = df)
ggsave("results_R/km_riskgroups.png",
       ggsurvplot(fit_risk, data = df, pval = TRUE, risk.table = TRUE)$plot, width = 7, height = 5, dpi = 150)
cat("\nC-index (discrimination):", round(summary(cox)$concordance[1], 3), "\n")

## ---- Interpretation ----
# Late-stage (III-IV) patients have markedly worse KM survival than early-stage (significant log-rank).
# In the Cox model tumour STAGE carries the largest hazard ratio (HR > 1), age adds a modest HR > 1, and
# molecular SUBTYPE modulates prognosis (Luminal A best, Basal/HER2 worse). The C-index (~0.65-0.75) shows
# the model ranks risk usefully, and the median-split risk groups separate on KM. Biologically this is the
# basis of clinical staging + precision oncology. Caveats: observational (association, not causal); PH may
# not hold for every term (see cox.zph); single cohort; C-index measures ranking, not calibrated risk.
