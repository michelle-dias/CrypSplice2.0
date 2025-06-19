#!/usr/bin/env python3
"""
Command Line Interface for CrypSplice 2.1
"""

import sys
import os

def cli_main():
    """
    CLI entry point that calls the main function from main.py
    """
    # Import the main function from the same package
    from .main import main
    
    # Call the main function - it already handles all argument parsing
    main()

if __name__ == "__main__":
    cli_main()