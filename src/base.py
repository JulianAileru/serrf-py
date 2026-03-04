from abc import ABC,abstractmethod
import logging
import sys,os
import numpy as np
import pandas as pd
sys.path.append(os.path.dirname(os.path.dirname(__file__)))


logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(name)s] %(levelname)s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    filename='BatchEffectsPipeline.log'
)


class Preprocessor:
    def __init__(self,transformation_method=None,normalization_method=None,imputation_method=None,logger=None):
        if logger:
            self.logger = logger
        else:
            self.logger = logging.getLogger(self.__class__.__name__)
        self.transformation_method = transformation_method
        self.normalization_method = normalization_method
        self.imputation_method = imputation_method
    @staticmethod   
    def return_data(data,method,logger):
        logger.info(f"No {method} Applied")
        return data 
    def _no_option(self, stage_name):
        return lambda data, logger: Preprocessor.return_data(data, stage_name, logger)
    # Imputation Methods
    @staticmethod
    def impute_gaussian_sampling(data,logger):
        logger.info('Imputing using Gaussian Sampling (per feature)')
        data = data.copy()
        nan_mask = data.isna()
        means = data.mean(axis=0, skipna=True).values[None,:]
        stds = data.std(axis=0, skipna=True, ddof=1).values[None,:]
        samples = np.random.normal(loc=means, scale=stds, size=data.shape)
        data[nan_mask] = samples[nan_mask]
        return data 
    @staticmethod
    def impute_global_minimum(data,logger):
        min_val = data.where(data > 0).min().min()
        logger.info(f'Imputing Minimum Non-Zero Value: {min_val}')
        return data.fillna(min_val * 0.10)
    @staticmethod
    def impute_median_value(data,logger):
        logger.info(f'Imputing Minimum Value: {data.median(axis=0)}')
        data = data.fillna(data.median(axis=0))
    @staticmethod
    def impute_mean_value(data,logger):
        logger.info(f'Imputing Mean Value: {data.mean(axis=0)}')
        data = data.fillna(data.mean(axis=0))
    # Normalization Methods
    @staticmethod
    def TIC(data,logger):
        logger.info("Applying Total Ion Current Normalization")
        scaler = data.sum(axis=1).mean()
        data = data.div(data.sum(axis=1), axis=0)
        data *= scaler
        return data
    @staticmethod
    def MedianNorm(data,logger):
        logger.info("Applying Median Normalization")
        data = (data / data.median(axis=1))
        return data
    
    @staticmethod
    def MeanNorm(data,logger):
        logger.info("Applying Mean Normalziation")
        data = (data / data.mean(axis=1))
        return data
    @staticmethod
    def PQN(data,logger):
        raise NotImplementedError
    @staticmethod
    def InternalStandards(data,logger):
        raise NotImplementedError
    # Transformation Methods
    @staticmethod
    def log2transform(data,logger):
        logger.info('Applying Log2 Transformation')
        data = data.copy()
        data = np.log2(data+1)
        return data
    @staticmethod
    def lntransform(data,logger):
        logger.info('Applying Natural Log Transformation')
        data = data.copy()
        data = np.log1p(data)
        return data
    def apply(self,data):
        self.logger.info("Starting Preprocessing")
        self.logger.info(f"Filling Missing Values with: {self.imputation_method} Method")
        
        imputation_methods = {
            "Global Minimum Value": Preprocessor.impute_global_minimum,
            "Median Value": Preprocessor.impute_median_value,
            "Mean Value": Preprocessor.impute_mean_value,
            "Gaussian Sampling": Preprocessor.impute_gaussian_sampling,
        }
        normalization_methods = {
            "TIC": Preprocessor.TIC,
            "Internal Standards": Preprocessor.InternalStandards,
            "PQN": Preprocessor.PQN,
            'MedianNorm': Preprocessor.MedianNorm,
            'MeanNorm': Preprocessor.MeanNorm
        }
        transformation_methods = {'Log2 Transformation':Preprocessor.log2transform,
                                  'Natural Log Transformation':Preprocessor.lntransform}
        
        impute = imputation_methods.get(self.imputation_method,self._no_option('Imputation'))
        print(self.imputation_method)
        transform = transformation_methods.get(self.transformation_method,self._no_option("Transformation"))
        normalize = normalization_methods.get(self.normalization_method,self._no_option("Normalization"))
        data = impute(data, logger=self.logger)
        data = normalize(data, logger=self.logger)
        data = transform(data, logger=self.logger)
        return data
    
    @staticmethod
    def adjust_data_labels(qc_str,blank_str,_data,_metadata):
        data = _data.copy()
        metadata = _metadata.copy()
        # data.index = data.index.astype(str)
        # metadata.index = metadata.index.astype(str)
        mask_qc_data = data.index.str.contains(qc_str)
        mask_qc_meta = metadata.index.str.contains(qc_str)
        mask_blank_data = data.index.str.contains(blank_str) | data.index.str.endswith("_BLANK")
        mask_blank_meta = metadata.index.str.contains(blank_str) | metadata.index.str.endswith("_BLANK")
        
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
        return data,metadata
    
class BatchCorrectionPipeline:
    def __init__(self, method, preprocessing_config):
        self.Method = method
        self.Preprocessor = preprocessing_config
    def original_sample_names(self,corrected_data,metadata):
        corrected_data.index = corrected_data.index.str.replace(r'_(QualityControl|Biological|BLANK)$',"",regex=True)
        metadata.index = metadata.index.str.replace(r'_(QualityControl|Biological|BLANK)$',"",regex=True)
        return corrected_data,metadata
    def return_blanks(self,corrected_data,original_data):
        original_data.index = original_data.index.str.replace(r'_(QualityControl|Biological|BLANK)$',"",regex=True)                
        idx = original_data.index[~original_data.index.isin(corrected_data.index)]
        corrected_data = pd.concat([corrected_data, original_data.loc[idx, :]])
        return corrected_data
    def correct(self, data, metadata):
        # Apply preprocessing
        processed_data = self.Preprocessor.apply(data)
        # Apply correction
        corrected_data = self.Method.correct(processed_data, metadata)
        corrected_data,metadata = self.original_sample_names(corrected_data,metadata=metadata)
        corrected_data = self.return_blanks(corrected_data,processed_data)
        return corrected_data



# metadata = pd.read_csv("/Users/jaileru/Projects/Metabolomics/batch-effects-dashboard-playground/data/streamlit_sample_metadata.csv").set_index("sample_name")
# data = pd.read_csv('/Users/jaileru/Projects/Metabolomics/batch-effects-dashboard-playground/data/streamlit_sample_data.csv').set_index("sample_name")
# metadata = metadata[metadata.batch <= 2]
# data = data.loc[metadata.index,:]
# pipeline = BatchCorrectionPipeline(
#     method=Combat(qc_str='SP',blank_str='B'),
#     preprocessing_config=Preprocessor(imputation_method="Global Minimum Value",normalization_method='TIC',transformation_method='Natural Log Transformation')
# )

# pipeline.correct(data=data,metadata=metadata)