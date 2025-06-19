"""
Setup script for CrypSplice 2.1
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file for the long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding='utf-8')

# Read requirements.txt
requirements_path = this_directory / "requirements.txt"
if requirements_path.exists():
    with open(requirements_path, 'r', encoding='utf-8') as f:
        requirements = [line.strip() for line in f if line.strip() and not line.startswith('#')]
else:
    # Fallback requirements based on your imports
    requirements = [
        "numpy>=1.21.0",
        "pandas>=1.5.0",
        "pybedtools>=0.8.2",
        "scikit-learn>=1.1.0",
        "matplotlib>=3.5.0",
        "pyBigWig>=0.3.18",
        "pysam>=0.19.0",
        "concurrent-futures; python_version<'3.2'"
    ]

setup(
    name="crypsplice",
    version="2.1.0",
    author="Michelle Dias",
    author_email="michelledias10@gmail.com",
    description="A robust computational tool designed to detect and analyze cryptic alternative splicing events from RNA-seq data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/michelle-dias/CrypSplice2.0",
    project_urls={
        "Bug Reports": "https://github.com/michelle-dias/CrypSplice2.0/issues",
        "Source": "https://github.com/michelle-dias/CrypSplice2.0",
        "Documentation": "https://github.com/michelle-dias/CrypSplice2.0#readme",
    },
    packages=find_packages(include=['crypsplice', 'crypsplice.*']),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
            "black>=22.0.0",
            "isort>=5.10.0",
            "flake8>=5.0.0"
        ]
    },
    entry_points={
        "console_scripts": [
            "crypsplice=crypsplice.cli:cli_main",
        ],
    },
    include_package_data=True,
    package_data={
        "crypsplice": [
            "lib/*.py",
            "*.py"
        ],
    },
    zip_safe=False,
    keywords="bioinformatics RNA-seq alternative-splicing cryptic-splicing genomics",
)