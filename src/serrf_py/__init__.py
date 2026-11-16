"""serrf-py: Python implementation of SERRF batch effect correction."""

from serrf_py.base import BatchCorrectionPipeline, Preprocessor
from serrf_py.correlation_errors import SERRF, BatchCorrector

__all__ = ["BatchCorrectionPipeline", "BatchCorrector", "Preprocessor", "SERRF"]
__version__ = "0.1.0"
