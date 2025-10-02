"""
Dependency checking for pseudomancer pipeline.
"""

import shutil
import sys


def check_dependencies():
    """Check if required dependencies (mmseqs2, datasets, getorf) are available."""
    
    dependencies = ["mmseqs", "datasets", "getorf"]
    missing = []
    
    for dep in dependencies:
        if not shutil.which(dep):
            missing.append(dep)
    
    if missing:
        sys.stderr.write(
            f"Missing required dependencies: {', '.join(missing)}\n"
            "Please install the following tools:\n"
            "- mmseqs2: https://github.com/soedinglab/MMseqs2\n"
            "- datasets: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/\n"
            "- getorf (EMBOSS): conda install -c bioconda emboss\n"
        )
        sys.exit(1)
