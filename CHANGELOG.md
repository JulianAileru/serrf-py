# Changelog

## [2026-03-31] Logging Improvements

### `src/base.py`

- **Console + file logging:** Module-level setup now outputs to stdout by default. The hardcoded `FileHandler` has been removed in favour of a `log_file` parameter on `BatchCorrectionPipeline`.
- **`BatchCorrectionPipeline`:** Added `log_file` parameter to `__init__`. When provided, a `FileHandler` is added to the root logger so all pipeline, preprocessor, and correction logs are written to the named file.
- **`BatchCorrectionPipeline.correct()`:** Added start log (method name, sample count, feature count), stage transition logs, and a completion log with output shape.
- **`BatchCorrectionPipeline.original_sample_names()`:** Added log describing the operation.
- **`BatchCorrectionPipeline.return_blanks()`:** Added log reporting the number of blank samples re-appended.
- **`Preprocessor.adjust_data_labels()`:** Added `logger` parameter; logs the QC/blank identifier strings and resulting per-label sample counts.
- **`Preprocessor.apply()`:** Removed stray `print(self.imputation_method)` call.
- **`Preprocessor.MeanNorm()`:** Fixed typo in log message (`"Normalziation"` → `"Normalization"`).
- **Bug fix — `impute_median_value()` / `impute_mean_value()`:** Both methods were missing `return` statements, causing them to silently return `None`.

### `src/correlation_errors.py`

- **`SERRF.correct()`:** Removed all duplicate `print()` calls and the unused `root_logger` / `logfile` variables. Replaced `tqdm` progress bar (per-signal updates) with structured log messages at each stage: label adjustment, optional imputation, per-batch correlation, feature selection, model fitting start/complete, cross-batch normalization, and final fix.
- **`SERRF.adjust_data_labels()`:** Added logs for the QC/blank strings used and the resulting per-label sample counts.
- **Removed `tqdm` import** (no longer used).

---

## [2026-03-04] `impute()` Method — Fixes and Integration

### Bug Fixes

#### 1. Missing `self` / incorrect method binding
**Location:** `SERRF.impute()` line 54

`impute` was defined as `def impute(current_batch):` without a `self` parameter and without `@staticmethod`. When called as `self.impute(batch)`, Python would pass `self` as `current_batch` and `batch` as an unexpected extra argument, raising a `TypeError`.

- **Fix:** Added `@staticmethod` decorator; call site updated to `SERRF.impute(...)`

#### 2. `impute()` called before `self.batches` and `self.all_data` were initialized
**Location:** `correct()` lines 347–354 (original)

The imputation block ran before `adjust_data_labels()` was called, so `self.batches` and `self.all_data` did not yet exist, causing an `AttributeError`.

- **Fix:** Moved the imputation block to after `adjust_data_labels()`, `self.all_data = self.all_data.T`, and `self.batches` are all set.

#### 3. Imputed data immediately overwritten
**Location:** `correct()` — original ordering

Even if imputation ran without error, `self.all_data` was set to the imputed result and then immediately overwritten by `self.adjust_data_labels(data=data, metadata=metadata)`, discarding all imputed values.

- **Fix:** Imputation now operates on `self.all_data` after labels are adjusted, slicing per batch via `self.all_data.columns[self.metadata['batch'] == batch]`.

#### 4. `min_val` excluded zeros for zero imputation (R includes them)
**Location:** `SERRF.impute()` — zero imputation block

Python used `row[~row_zero_mask & ~row_na_mask]` (non-zero AND non-NA), whereas R's `min(all[j,...][!is.na(...)])` only excludes NAs, leaving zeros in the minimum calculation.

- **R behavior:** `mean = min(non_NA_values) + 1` — zeros contribute to the minimum
- **Python (before):** `min_val = row[~zero & ~NA].min()` — zeros excluded; rows with only zeros and NAs would skip imputation entirely
- **Fix:** Changed to `valid_for_zeros = row[~row_na_mask]` so zeros are included

#### 5. NA mean formula off by 0.5
**Location:** `SERRF.impute()` — NA imputation block

- **R behavior:** `mean = 0.5 * min + 1`
- **Python (before):** `loc = 0.5 * (min_val + 1)` = `0.5 * min + 0.5`
- **Fix:** Changed to `loc = 0.5 * min_val_na + 1`

#### 6. Scale formula constant wrong (both zeros and NAs)
**Location:** `SERRF.impute()` — both imputation blocks

- **R behavior:** `sd = 0.1 * (min + 0.1)` = `0.1 * min + 0.01`
- **Python (before):** `scale = (0.1 * min_val) + 0.1` = `0.1 * min + 0.1`
- **Fix:** Changed to `scale = 0.1 * (min_val + 0.1)` for both blocks

#### 7. NA min not recomputed after zero imputation
**Location:** `SERRF.impute()` — NA imputation block

R processes zeros first (modifying in-place), then recomputes `min(non_NA)` for NA imputation — at this point zeros have been replaced by imputed values. Python computed a single `min_val` upfront used for both.

- **Fix:** After zero imputation, `min_val_na` is recomputed from the updated row: `current_batch.loc[signal, ~row_na_mask].min()`

---

## [2026-02-05] SERRF Port Fixes

Comprehensive fixes to align the Python SERRF implementation (`src/correlation_errors.py`) with the original R version (`R/serrf_function.R`).

### Error Fixes

#### 1. y-scaling condition and factor inversion
**Location:** `fit_predict()` lines 269-278

The condition for when to scale the target variable was inverted, and the scaling factor was computed as the reciprocal of R's.

- **R behavior:** Scale when `sd_QC >= sd_Bio`, using `factor = sd_QC / sd_Bio`
- **Python (before):** Scaled when `sd_QC < sd_Bio`, using `factor = sd_Bio / sd_QC`
- **Fix:** Condition now triggers centering-only when `qc_std < sample_std`, and factor is `qc_std / sample_std`

#### 2. `final_fix()` QC block copy-paste error
**Location:** `final_fix()` lines 399-415

The QC sample fix block was reading from `normed_train` but writing to `normed_target`. QC samples were never actually fixed.

- **Fix:** Changed all `normed_target` references in the QC block to `normed_train`

#### 3. `normalize_all_batches()` correction factor guard
**Location:** `normalize_all_batches()` line 374

Python used `max(c, 1)` which only allowed scaling QC up. R uses `ifelse(c > 0, c, 1)` which allows scaling down when `0 < c < 1`.

- **Fix:** Changed to `c if c > 0 else 1`

#### 4. Feature selection off-by-one
**Location:** `top_correlated()` lines 220-237

Python dropped the signal from correlation lists before slicing, while R includes it in the top-l and removes it after computing the intersection. This caused Python to consider one extra non-self feature per list at each iteration.

- **Fix:** Keep signal in lists, subtract `{signal}` from intersection inside the loop

#### 5. Random Forest hyperparameters
**Location:** `fit_predict()` lines 292-294

`min_samples_leaf=1` caused overfitting. R's `ranger` defaults to `min.node.size=5` for regression.

- **Fix:** Set `min_samples_leaf=5`, removed `min_samples_split` (no ranger equivalent), set `random_state=1` to match R's seed intent

#### 6. Order of operations in `fit_predict()`
**Location:** `fit_predict()` lines 302-327

The median adjustments were interleaved with normalization, and the negative Bio fix came after median adjustments. R does normalization first, then negative fix, then median adjustments.

- **Before:** QC norm → QC median adjust → Bio norm → Bio median adjust → negative fix
- **After:** QC norm → Bio norm → negative fix → QC median adjust → Bio median adjust

#### 7. `final_fix()` negative replacement distribution
**Location:** `final_fix()` lines 420-426, 438-444

R generates ONE random `Uniform(0,1)` value and applies it to ALL negative entries. Python was generating a different random value per entry from `Uniform(0.1, 1.0)`.

- **R behavior:** `runif(1) * min(positive_values)` — single draw, applied to all
- **Python (before):** `np.random.uniform(0.1, 1.0, size=n) * min_val` — n draws, wrong range
- **Fix:** `np.random.uniform(0, 1) * min_val` — single draw, correct range
- **Note:** The actual random values will differ between Python and R due to different RNGs, but the logic now matches.

#### 8. `final_fix()` strictly positive minimum value
**Location:** `final_fix()` lines 422, 440

When computing the minimum value to scale for negative replacement, Python was using `>= 0` (non-negative, including zeros) while R uses `> 0` (strictly positive).

- **R behavior:** `min(x[x > 0])` — strictly positive values only
- **Python (before):** `row[~neg_mask]` which is `row[row >= 0]` — included zeros
- **Fix:** Changed to `row[row > 0]` for strictly positive values
- **Impact:** If any value is exactly 0, it would have been included in min() calculation before, potentially returning 0 and causing replacement values to be 0.

### New Features

#### `use_ranger` flag for rpy2 integration
**Location:** `SERRF.__init__()`, `_ranger_fit_predict()`, `fit_predict()`, `correct()`

Added option to call R's `ranger` package directly via rpy2 for exact reproducibility with R output.

```python
# Default: sklearn RandomForestRegressor
SERRF(qc_str='_QC', blank_str='blank')

# Optional: R's ranger via rpy2
SERRF(qc_str='_QC', blank_str='blank', use_ranger=True)
```

- Automatically forces `n_jobs=1` when `use_ranger=True` (R interpreter is not process-safe)
- Requires R with `ranger` package installed
- Uses `seed=1` parameter in ranger call to match R's `set.seed(1)`

### Validation

Validated against R output using the SERRF example dataset:
- **Percent error:** Median 0.53%, Mean 2.46%, Max 13%
- **Pearson correlation:** >0.95 for all metabolites
- **Concordance correlation coefficient:** Strong agreement
- **RSD alignment:** Confirmed when using same RSD function on both outputs
