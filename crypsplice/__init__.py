"""
CrypSplice 2.0 - A robust computational tool for detecting cryptic alternative splicing events.

This package provides tools for analyzing RNA-seq data to identify differential 
splicing patterns that are not annotated (cryptic) in standard gene models.
"""

__version__ = "2.1.0"  # Updated to match your script version
__author__ = "Michelle Dias"
__email__ = "michelledias10@gmail.com"

# Import main functionality
from .main import main

__all__ = ["main"]