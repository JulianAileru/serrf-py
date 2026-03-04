# serrf-py

A Python implementation of **SERRF** (Systematic Error Removal using Random Forest), a batch effect correction method for large-scale untargeted metabolomics data.

> **Reference:** Fan, S., Kind, T., Cajka, T., Hazen, S. L., Tang, W. H. W., Kaddurah-Daouk, R., Irvin, M. R., Arnett, D. K., Barupal, D. K., & Fiehn, O. (2019). Systematic Error Removal Using Random Forest for Normalizing Large-Scale Untargeted Lipidomics Data. *Analytical Chemistry*, 91(5), 3590–3596. https://doi.org/10.1021/acs.analchem.8b05592

---

## Overview

SERRF corrects systematic measurement drift in metabolomics datasets by leveraging the correlation structure between quality control (QC) samples and biological samples within each analytical batch. A Random Forest regression model is trained on QC samples and applied to predict and remove the systematic error component from biological samples.

This repository ports the original R implementation (available at [slfan.shinyapps.io/SERRF](https://slfan.shinyapps.io/SERRF)) to Python, validated against R output on the SERRF example dataset.

---

## Repository Structure

```
serrf-py/
├── src/
│   ├── base.py                  # Preprocessor and BatchCorrectionPipeline classes
│   └── correlation_errors.py    # SERRF batch corrector
├── R/
│   ├── serrf_function.R         # Original R implementation (reference)
│   └── app.R                    # Original Shiny web application
├── data/
│   └── SERRF example dataset.xlsx
├── notebooks/
│   └── SERRFValidation.ipynb    # Python vs R comparison
└── CHANGELOG.md
```

---

## Usage

```python
from src.base import Preprocessor, BatchCorrectionPipeline
from src.correlation_errors import SERRF

preprocessor = Preprocessor(
    imputation_method="Global Minimum Value",
    normalization_method="TIC",
    transformation_method="Natural Log Transformation"
)

corrector = SERRF(
    qc_str="QC",       # substring identifying QC samples in sample names
    blank_str="blank", # substring identifying blank samples
    n_jobs=-1,         # parallel jobs for Random Forest fitting
    num_features=10    # number of correlated features to use per signal
)

pipeline = BatchCorrectionPipeline(method=corrector, preprocessing_config=preprocessor)
corrected = pipeline.correct(data=data, metadata=metadata)
```

`data` should be a `DataFrame` of shape `(n_samples, n_features)` and `metadata` a `DataFrame` with at least a `batch` column, both indexed by sample name.

### Optional: Exact R Reproducibility via rpy2

```python
corrector = SERRF(qc_str="QC", blank_str="blank", use_ranger=True)
```

Requires R with the `ranger` package installed. Forces `n_jobs=1` (R interpreter is not process-safe across workers).

### Optional: Built-in Imputation

```python
corrector = SERRF(qc_str="QC", blank_str="blank")
corrector.serrf_impute = True  # impute zeros and NAs per batch before correction
```

Zeros are replaced by draws from N(min + 1, 0.1·(min + 0.1)) and missing values by N(0.5·min + 1, 0.1·(min + 0.1)), matching the R implementation exactly.

---

## Validation Metrics

The Python output was compared against R output on the SERRF example dataset. Four metrics were computed per metabolite and summarized across the full feature set:

### QC Relative Standard Deviation (QC-RSD)
Measures within-batch variability of QC samples after correction. Lower RSD indicates better removal of systematic drift.

```
RSD (%) = (std / mean) × 100
```

Computed separately per batch on QC samples; the post-correction QC-RSD distribution should be substantially lower than pre-correction.

### Pearson Correlation
Linear correlation between Python-corrected and R-corrected values across all samples for each metabolite. Summarized as median and mean across all features.

| Statistic | Value |
|-----------|-------|
| Median Pearson r | > 0.95 |

### Concordance Correlation Coefficient (CCC)
Combines Pearson correlation with a bias-correction factor that penalizes systematic location and scale differences between the two outputs. Stricter than Pearson alone — values near 1 indicate both agreement in trend and in absolute level.

```
CCC = (2 · σ_xy) / (σ_x² + σ_y² + (μ_x − μ_y)²)
```

Strong agreement observed across metabolites, confirming the Python port reproduces R values at the correct scale and not just in relative ordering.

### Percent Error
Absolute relative deviation between Python and R corrected values, per sample per metabolite:

```
% error = |Python − R| / |R| × 100
```

| Statistic | Value |
|-----------|-------|
| Median | 0.53% |
| Mean | 2.46% |
| Max | 13% |

The small median error and strong correlation confirm the Python implementation faithfully reproduces the R algorithm. The tail in max error arises from stochastic differences in random number generation between Python (NumPy) and R.

---

## Dependencies

- `pandas`
- `numpy`
- `scikit-learn`
- `joblib`
- `tqdm`
- `rpy2` *(optional, for `use_ranger=True`)*
