"""python-tetracorder: a Python port of the Tetracorder 6 spectral identification expert system."""
from importlib.metadata import version

from .group_evaluator import GroupEvaluator

__version__ = version("tetracorderp")
__all__ = ["GroupEvaluator", "__version__"]