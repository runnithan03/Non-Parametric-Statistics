# Non-Parametric Statistics

This repository contains coursework and projects for the Nonparametric Statistics IV module at Durham University, across both Michaelmas and Epiphany terms. It includes hands-on exploration of resampling techniques, Bayesian methods, kernel smoothing, and local polynomial regression, all implemented in **R**.

---

## 📁 Michaelmas Term

### 🧪 Practicals

- **Practical 1:** Classical bootstrap with `boot` using summary statistics (e.g. medians)
- **Practical 2:** Bayesian bootstrap with `bayesboot` (Cars93 dataset)
- **Practical 3:** Density estimation with Dirichlet Process Mixture Models using `BNPmix` (Iris dataset)
- **Practical 4:** Regression with DPMMs (Abalone dataset from `AppliedPredictiveModeling`)

### 📊 Mini Project: Market Volatility

- **Topic:** *Bayesian and Bootstrap Methods for VIX Volatility Modelling*
- **Goals:**
  - Compare classical and Bayesian bootstrap confidence intervals
  - Cluster VIX log-returns with a Dirichlet Process Mixture Model
  - Predict VIX returns using XGBoost and extracted clusters
- **Tools:** `boot`, `bayesboot`, `BNPmix`, `xgboost`

### 📚 Lecture Notes & Files
- `Michaelmas Notes.pdf`: Lecture notes on EDF, DPMMs, bootstrap theory
- `Mini_Project_Report.pdf`: VIX project write-up
- `Raul_Unnithan_Code.R`: Full code for bootstrap, Bayesian inference, and XGBoost
- `Marked Report.pdf`: Feedback and grade

---

## 📁 Epiphany Term

### 🧪 Practicals

- **Practical 5:** Histogram estimation and bin-width selection (Davis dataset, Scott’s and ROT methods)
- **Practical 6:** Kernel Density Estimation (KDE) and bandwidth selection (LakeHuron dataset)
- **Practical 7:** Multivariate KDE using `ks` package (Soils dataset up to 5D)

### 📊 Mini Project: Local Linear Regression CIs

- **Topic:** *Analytical and Simulation-Based Confidence Intervals for Local Linear Regression*
- **Dataset:** Onion yield vs planting density (South Australia)
- **Goals:**
  - Derive pointwise 95% confidence intervals using:
    - An analytical plug-in method based on asymptotic bias/variance
    - A simulation-based (bootstrap) method using residual resampling
  - Compare both approaches visually and via average CI width
- **Key Equations:**
  - LL estimator: `m̂(x)`
  - Bias: (1/2) · h² · m″(x) · μ₂
  - Variance: σ² · ν₀ / (n · h · f(x))
- **Tools:** `locpol`, `SemiPar`, `density()`, custom KDE and LL functions

### 📚 Lecture & Report Files
- `Nonparametric_Statistics_IV.pdf`: Lecture slides on density estimation, KDE, CI theory
- `Nonparametric Statistics - Epiphany term 2025 - Full.pdf`: Complete lecture notes
- `Report.pdf`: Full write-up of LL CI comparison
- `Report Code.R`: Implementation of LL regression with analytic & bootstrap CIs
- `Guidelines.pdf`: Project instructions and marking rubric

---

## 📦 Key R Packages Used

```r
install.packages(c(
  "boot", "bayesboot", "BNPmix", "ggplot2", "reshape2", "dplyr",
  "xgboost", "caret", "locpol", "SemiPar", "ks", "MASS"
))
