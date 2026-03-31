from abc import ABC,abstractmethod
import logging
import pandas as pd
import numpy as np 
from joblib import Parallel,delayed
from sklearn.ensemble import RandomForestRegressor
import sys,os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))



## Correlation of Errors
class BatchCorrector(ABC):
    @abstractmethod
    def correct(self,data,metadata):
        pass

class SERRF(BatchCorrector):
    def __init__(self,qc_str,blank_str,n_jobs=-1,num_features=10,use_ranger=False):
        self.logger = logging.getLogger(self.__class__.__name__)
        self.qc_str = qc_str
        self.blank_str = blank_str
        self.n_jobs = n_jobs
        self.serrf_impute = False
        self.num_features=num_features
        self.use_ranger = use_ranger
    def adjust_data_labels(self,data,metadata,rowvar=False):
        self.logger.info(f"Adjusting sample labels (QC='{self.qc_str}', Blank='{self.blank_str}')")
        mask_qc_data = data.index.str.contains(self.qc_str)
        mask_qc_meta = metadata.index.str.contains(self.qc_str)
        mask_blank_data = data.index.str.contains(self.blank_str) | data.index.str.endswith("_BLANK")
        mask_blank_meta = metadata.index.str.contains(self.blank_str) | metadata.index.str.endswith("_BLANK")

        # Apply BLANK labels first
        data.index = data.index.where(~mask_blank_data, data.index + "_BLANK")
        metadata.index = metadata.index.where(~mask_blank_meta, metadata.index + "_BLANK")

        # Apply QC labels next
        data.index = data.index.where(~mask_qc_data, data.index + "_QualityControl")
        metadata.index = metadata.index.where(~mask_qc_meta, metadata.index + "_QualityControl")

        # Apply Biological labels to everything else
        mask_bio_data = ~(mask_qc_data | mask_blank_data)
        mask_bio_meta = ~(mask_qc_meta | mask_blank_meta)
        data.index = data.index.where(~mask_bio_data, data.index + "_Biological")
        metadata.index = metadata.index.where(~mask_bio_meta, metadata.index + "_Biological")
        self.logger.info(
            f"Label counts — Biological: {mask_bio_data.sum()}, "
            f"QC: {mask_qc_data.sum()}, Blank: {mask_blank_data.sum()}"
        )
        return data,metadata

    @staticmethod
    def impute(current_batch):
        signals = current_batch.index.to_list()
        zero_mask = current_batch == 0
        na_mask = current_batch.isna()

        for signal in signals:
            row = current_batch.loc[signal]
            row_zero_mask = zero_mask.loc[signal]
            row_na_mask = na_mask.loc[signal]
            # R uses min of non-NA values (zeros included) for zero imputation
            valid_for_zeros = row[~row_na_mask]

            if valid_for_zeros.empty:
                continue

            min_val = valid_for_zeros.min()

            if row_zero_mask.any():
                zero_impute = np.random.normal(
                    size=row_zero_mask.sum(),
                    loc=(min_val + 1),
                    scale=(0.1 * (min_val + 0.1))
                )
                current_batch.loc[signal, row_zero_mask] = zero_impute

            if row_na_mask.any():
                # R recomputes min after zeros are replaced; re-read the updated row
                min_val_na = current_batch.loc[signal, ~row_na_mask].min()
                na_impute = np.random.normal(
                    size=row_na_mask.sum(),
                    loc=(0.5 * min_val_na + 1),
                    scale=(0.1 * (min_val_na + 0.1))
                )
                current_batch.loc[signal, row_na_mask] = na_impute
        return current_batch
    def compute_correlation(self):
        """
        Current Batch is of shape (n_signals,n_samples)
        """
        train = self.current_batch.loc[:, self.current_batch.columns.str.contains("_QualityControl")]
        target = self.current_batch.loc[:, self.current_batch.columns.str.contains("_Biological")]
        
        # Scale per row (per signal)
        train_scale = train.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        target_scale = target.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
        
        # Compute correlation matrices
        corr_train = train_scale.T.corr(method='spearman')
        corr_target = target_scale.T.corr(method='spearman')
        
        train_order = {}
        target_order = {}
        
        for signal in self.current_batch.index:
            train_order[signal] = corr_train.loc[:, signal].abs().sort_values(ascending=False)
            target_order[signal] = corr_target.loc[:, signal].abs().sort_values(ascending=False)
        
        self.train_order,self.target_order = train_order,target_order
        return None
    def top_correlated(self,n=10):
        """
        Find top correlated features of each signal
        """
        features = {}

        for signal in self.current_batch.index:
            train_lst = self.train_order[signal].index.tolist()
            target_lst = self.target_order[signal].index.tolist()

            l = n
            max_len = min(len(train_lst), len(target_lst))

            while l <= max_len:
                temp_train = set(train_lst[:l])
                temp_target = set(target_lst[:l])
                inter = (temp_train & temp_target) - {signal}

                if len(inter) >= n:
                    features[signal] = list(inter)
                    break
                l += 1
            else:
                features[signal] = list((temp_train & temp_target) - {signal})
        self.features = features
        return None
    @staticmethod
    def _ranger_fit_predict(train_X, train_y, test_X):
        """Fit ranger via rpy2 and return (train_pred, test_pred) as numpy arrays.
        Requires R with ranger installed. Must run with n_jobs=1 (R is not process-safe)."""
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter

        col_names = [f'V{i+1}' for i in range(train_X.shape[1])]

        train_df = train_X.copy()
        train_df.columns = col_names
        train_df.insert(0, 'y', train_y.values if hasattr(train_y, 'values') else train_y)

        test_df = test_X.copy()
        test_df.columns = col_names

        with localconverter(ro.default_converter + pandas2ri.converter):
            ro.r.assign('serrf_train', train_df)
            ro.r.assign('serrf_test', test_df)

        ro.r('library(ranger)')
        ro.r('set.seed(1)')
        ro.r('serrf_model <- ranger(y ~ ., data = serrf_train)')
        train_pred = np.array(ro.r('predict(serrf_model, data = serrf_train)$predictions'))
        test_pred = np.array(ro.r('predict(serrf_model, data = serrf_test)$predictions'))

        return train_pred, test_pred

    @staticmethod
    def fit_predict(all_data,metadata,batch,features,signal,use_ranger=False):
        """
        Fit RF model and predict, following R logic more closely
        """
        current_batch = all_data.loc[:,metadata['batch'] == batch]
        train_id = current_batch.columns[current_batch.columns.str.contains("_QualityControl")].tolist()
        test_id = current_batch.columns[current_batch.columns.str.contains("_Biological")].tolist()
        
        feat = features[signal]
        if len(feat) == 0:
            return current_batch.loc[signal, :].copy()
        
        # Get training data (QC samples)
        train_X = current_batch.loc[feat, train_id].T
        train_y = current_batch.loc[signal, train_id]
        
        # Get test data (biological samples)  
        test_X = current_batch.loc[feat, test_id].T
        test_y = current_batch.loc[signal, test_id]
        
        # Scale training features per column (per feature)
        train_X_scaled = (train_X - train_X.mean(axis=0)) / train_X.std(axis=0)
        test_X_scaled = (test_X - test_X.mean(axis=0)) / test_X.std(axis=0)
        
        # Handle scaling of target variable
        qc_std = train_y.std()
        sample_std = test_y.std()
        
        if qc_std == 0 or sample_std == 0 or np.isnan(qc_std) or np.isnan(sample_std) or qc_std < sample_std:
            # Center only
            train_y_scaled = train_y - train_y.mean()
        else:
            # Scale by factor (qc_std / sample_std >= 1)
            factor = qc_std / sample_std
            if len(train_id) * 2 >= len(test_id):
                train_y_scaled = (train_y - train_y.mean()) / factor
            else:
                train_y_scaled = train_y - train_y.mean()
        
        # Remove NaN features
        valid_cols = ~(train_X_scaled.isna().any() | test_X_scaled.isna().any())
        train_X_scaled = train_X_scaled.loc[:, valid_cols]
        test_X_scaled = test_X_scaled.loc[:, valid_cols]
        
        if train_X_scaled.shape[1] == 0:
            return current_batch.loc[signal, :].copy()
        
        # Fit model and predict
        if use_ranger:
            train_pred, test_pred = SERRF._ranger_fit_predict(train_X_scaled, train_y_scaled, test_X_scaled)
        else:
            rfr = RandomForestRegressor(n_estimators=500,max_features=1/3,
                                        min_samples_leaf=5,bootstrap=True,
                                        oob_score=True, random_state=1,n_jobs=1)
            rfr.fit(train_X_scaled, train_y_scaled)
            train_pred = rfr.predict(train_X_scaled)
            test_pred = rfr.predict(test_X_scaled)
        
        # Initialize normalized values
        norm = current_batch.loc[signal, :].copy()

        # 1. Normalize QC samples (training)
        train_pred_original_scale = train_pred + train_y.mean()
        train_scale_factor = train_pred_original_scale / all_data.loc[signal, all_data.columns.str.contains("_QualityControl")].mean()
        norm[train_id] = norm[train_id] / train_scale_factor

        # 2. Normalize biological samples (test)
        test_pred_original_scale = test_pred + test_y.mean()
        test_pred_centered = test_pred_original_scale - test_pred.mean()
        test_scale_factor = test_pred_centered / all_data.loc[signal, all_data.columns.str.contains("_Biological")].median()
        norm[test_id] = norm[test_id] / test_scale_factor

        # 3. Fix negative values in biological samples
        negative_mask = norm[test_id] < 0
        if negative_mask.any():
            norm.loc[negative_mask.index[negative_mask]] = current_batch.loc[signal, negative_mask.index[negative_mask]]

        # 4. Adjust QC to match global median
        norm[train_id] = norm[train_id] / (norm[train_id].median() / all_data.loc[signal, all_data.columns.str.contains("_QualityControl")].median())

        # 5. Adjust biological samples to match global median
        norm[test_id] = norm[test_id] / (norm[test_id].median() / all_data.loc[signal, all_data.columns.str.contains("_Biological")].median())

        # 6. Handle infinite values
        inf_mask = np.isinf(norm)
        if inf_mask.any():
            norm[inf_mask] = np.random.normal(0, norm[~inf_mask].std() * 0.01, size=inf_mask.sum())
        
        # Outlier detection and correction (following R logic)
        Q1, Q3 = np.percentile(norm, [25, 75])
        IQR = Q3 - Q1
        outliers = (norm < (Q1 - 3*IQR)) | (norm > (Q3 + 3*IQR))
        
        if outliers.any():
            # Alternative calculation for outliers
            attempt = ((current_batch.loc[signal, test_id]) - 
                    (test_pred_original_scale) + 
                    all_data.loc[signal, all_data.columns.str.contains("_Biological")].median())
            
            outlier_test_samples = outliers & norm.index.isin(test_id)
            
            if outlier_test_samples.any() and len(attempt) > 0:
                outlier_values = norm[outliers]
                if outlier_values.mean() > norm.mean():
                    if attempt[outlier_test_samples].mean() < outlier_values.mean():
                        norm[outlier_test_samples] = attempt[outlier_test_samples]
                else:
                    if attempt[outlier_test_samples].mean() > outlier_values.mean():
                        norm[outlier_test_samples] = attempt[outlier_test_samples]
        
        # Final fix for negative values
        final_negative_mask = norm < 0
        if final_negative_mask.any():
            norm[final_negative_mask] = current_batch.loc[signal, final_negative_mask]
        
        return norm
    def normalize_all_batches(self,normalized_data):
        result = normalized_data.copy()
        qc_cols = self.all_data.columns[self.all_data.columns.str.contains("_QualityControl")]
        sample_cols = self.all_data.columns[self.all_data.columns.str.contains("_Biological")]
        for signal in self.all_data.index:
            # Calculate the correction factor
            sample_median_norm = result.loc[signal, sample_cols].median()
            qc_median_all = self.all_data.loc[signal, qc_cols].median()
            sample_median_all = self.all_data.loc[signal, sample_cols].median()
            sample_std_all = self.all_data.loc[signal, sample_cols].std()
            sample_std_norm = result.loc[signal, sample_cols].std()
            qc_median_norm = result.loc[signal, qc_cols].median()
            if (sample_std_all == 0) or (qc_median_norm == 0):
                continue
            else:
                c = (sample_median_norm + ((qc_median_all - sample_median_all) / sample_std_all) * sample_std_norm) / qc_median_norm
                # Apply correction to QC samples
                result.loc[signal, qc_cols] = result.loc[signal, qc_cols] * (c if c > 0 else 1)
        return result
    def final_fix(self,result):
        normed_target = result.loc[:,result.columns[result.columns.str.contains("_Biological")]]
        normed_train = result.loc[:,result.columns[result.columns.str.contains("_QualityControl")]]
        for signal in result.index:
            row = normed_target.loc[signal]

            # Fix NaNs
            na_mask = row.isna()
            valid_vals = row[~na_mask]
            if not valid_vals.empty:
                min_val = valid_vals.min()
                std_ = valid_vals.std()
                std_ = std_ * 0.1 if std_ > 0 else 1e-8  # prevent 0 std
                imputed_vals = np.random.normal(loc=min_val, scale=std_, size=na_mask.sum())
                normed_target.loc[signal, na_mask.index[na_mask]] = imputed_vals
            # Fix negatives (R: runif(1) * min(values > 0) — single draw applied to all)
            neg_mask = normed_target.loc[signal] < 0
            pos_vals = normed_target.loc[signal][normed_target.loc[signal] > 0]
            if not pos_vals.empty and neg_mask.any():
                min_val = pos_vals.min()
                replacement_val = np.random.uniform(0, 1) * min_val
                normed_target.loc[signal, neg_mask.index[neg_mask]] = replacement_val

            row = normed_train.loc[signal]
            # Fix NaNs
            na_mask = row.isna()
            valid_vals = row[~na_mask]
            if not valid_vals.empty:
                min_val = valid_vals.min()
                std_ = valid_vals.std()
                std_ = std_ * 0.1 if std_ > 0 else 1e-8  # prevent 0 std
                imputed_vals = np.random.normal(loc=min_val, scale=std_, size=na_mask.sum())
                normed_train.loc[signal, na_mask.index[na_mask]] = imputed_vals
            # Fix negatives (R: runif(1) * min(values > 0) — single draw applied to all)
            neg_mask = normed_train.loc[signal] < 0
            pos_vals = normed_train.loc[signal][normed_train.loc[signal] > 0]
            if not pos_vals.empty and neg_mask.any():
                min_val = pos_vals.min()
                replacement_val = np.random.uniform(0, 1) * min_val
                normed_train.loc[signal, neg_mask.index[neg_mask]] = replacement_val
        df = pd.concat([normed_train,normed_target],axis=1).T
        return df
        
    def correct(self,data,metadata):
        self.all_data,self.metadata = self.adjust_data_labels(data=data,metadata=metadata)
        self.all_data = self.all_data.T
        self.signals = self.all_data.index.to_list()
        self.batches = self.metadata['batch'].unique()
        self.logger.info(
            f"Starting SERRF correction — {len(self.signals)} signals, "
            f"{len(self.batches)} batch(es)"
        )
        if self.serrf_impute:
            self.logger.info("Imputing missing values per batch (Gaussian sampling)")
            imputed_all = []
            for batch in self.batches:
                batch_cols = self.all_data.columns[self.metadata['batch'] == batch]
                filled_na = SERRF.impute(self.all_data.loc[:, batch_cols])
                imputed_all.append(filled_na)
            self.all_data = pd.concat(imputed_all,axis=1)
            self.logger.info("Imputation complete")
        normalized_batches = []
        for i, batch in enumerate(self.batches, 1):
            self.current_batch = self.all_data.loc[:,self.metadata['batch'] == batch]
            n_qc = self.current_batch.columns.str.contains("_QualityControl").sum()
            n_bio = self.current_batch.columns.str.contains("_Biological").sum()
            self.logger.info(
                f"Batch {batch} [{i}/{len(self.batches)}] — "
                f"{n_qc} QC samples, {n_bio} biological samples"
            )
            self.logger.info(f"Batch {batch} — computing Spearman correlation matrix")
            self.compute_correlation()
            self.logger.info(f"Batch {batch} — selecting top {self.num_features} correlated features per signal")
            self.top_correlated(n=self.num_features)
            self.logger.info(f"Batch {batch} — fitting random forest models ({len(self.signals)} signals)")
            n_jobs = 1 if self.use_ranger else self.n_jobs
            normalized_signals = Parallel(n_jobs=n_jobs)(
                delayed(SERRF.fit_predict)(
                    all_data=self.all_data, batch=batch, features=self.features,
                    metadata=self.metadata, signal=signal, use_ranger=self.use_ranger
                )
                for signal in self.signals
            )
            self.logger.info(f"Batch {batch} — model fitting complete")
            batch_normalized = pd.concat(normalized_signals,axis=1)
            normalized_batches.append(batch_normalized)
        normalized_batches = pd.concat(normalized_batches).T
        self.logger.info("Aligning QC samples across all batches")
        result = self.normalize_all_batches(normalized_data=normalized_batches)
        self.logger.info("Applying final imputation and negative-value fixes")
        norm = self.final_fix(result)
        self.logger.info(f"SERRF correction complete — output shape: {norm.shape}")
        return norm
     