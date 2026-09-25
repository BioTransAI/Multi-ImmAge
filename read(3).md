# Multimodal Immune Age in UK Biobank

*Analysis workflow · Hematology, metabolomics, and proteomics*

Three modality-specific age predictors are combined into an integrated immune-age clock, followed by analyses of cross-modal concordance, phenotypes, disease, and mortality.

**Analysis script:** [multimodal_immune_age_analysis.R](multimodal_immune_age_analysis.R)

## Workflow

| Module | Analysis | Main operations |
| :---: | :--- | :--- |
| **1** | Model development | Use an 80/20 split and five-fold outer/inner cross-validation; select modality-specific learners and the fusion method, evaluate the test set, and generate full-cohort out-of-fold predictions. |
| **2** | Model comparison | Compare held-out MAE using 2,000 paired bootstrap samples; calculate modality contributions and prediction correlations. |
| **3** | Cross-modal concordance | Derive standardized age residuals, acceleration patterns, discordance, and directional concordance. |
| **4** | Phenotype associations | Fit separate adjusted OLS models for four clocks, using standardized outcomes, HC1 robust standard errors, and global Bonferroni correction. |
| **5** | Disease and mortality | Fit logistic, cause-specific Cox, and negative-binomial models; retain proportional-hazards checks and regional follow-up dates. |
| **6** | Disease overlap | Apply Bonferroni correction across planned incident endpoints and four clocks; summarize shared, integrated-only, and recovered associations. |
| **7** | Mortality risk and prediction | Estimate spline dose responses and standardized absolute risk; use 200 paired bootstrap samples for out-of-bag AUC and Uno's C. |

## Measures and adjustment

**Age acceleration** is the standardized residual from regressing predicted age on chronological age, using pooled out-of-fold predictions. Four clocks are assessed: integrated, hematologic, metabolomic, and proteomic.

**Covariates:** age, sex, assessment-centre region, ethnicity, Townsend deprivation index, smoking, alcohol, physical activity, and BMI. The supplied formulas, seeds, thresholds, and model parameters are retained.

## Use

1. Replace the `<...>` placeholders with input and output locations.
2. Provide preprocessed features, modality feature lists, chronological age, covariates, numeric phenotypes, blood-collection dates, and hospital/death records. Participant identifiers, row alignment, and factor coding must match the input checks.
3. Run `multimodal_immune_age_analysis.R` from top to bottom. All seven modules, model-fitting functions, and clinical analyses are included in this file.

The script saves model and clinical RDS results and clinical CSV outputs; mortality bootstrap results are cached. Additional summaries are available within their respective modules. Export steps absent from the supplied code are left as comments.

<details>
<summary><strong>Software requirements</strong></summary>

R ≥ 4.1; Linux for the existing fork-based parallel workflow.

Core packages: `data.table`, `mlr3verse`, `mlr3learners`, `mlr3extralearners`, `callr`, `survival`, `MASS`, `timeROC`, and `survAUC`, plus the packages required by the configured learners. Package-availability checks are included.

</details>
