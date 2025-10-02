# pseudomancer

A command line tool for reconstructing pseudogenes in prokaryotic genomes.

## Installation

```bash
# Clone the repository
git clone https://github.com/Floto-Lab/pseudomancer.git
cd pseudomancer

# Create conda environment with all dependencies
mamba env create -f environment.yml
mamba activate pseudomancer_env
```

**Requirements**: [Mamba](https://mamba.readthedocs.io/) or [Conda](https://docs.conda.io/) package manager

## Usage

```bash
# Activate environment (if not already active)
mamba activate pseudomancer_env

# Run pipeline
python -m pseudomancer --taxon Mycobacterium --genome target_genome.fasta --out_dir results/
```

### Parameters

- `--taxon`: Taxon name for downloading reference proteins from NCBI RefSeq (e.g., 'Mycobacterium' or 'Mycobacterium tuberculosis')
- `--genome`: Target genome FASTA file to search for pseudogenes
- `--out_dir`: Output directory for results
- `--evalue`: E-value threshold (default: 1e-5)
- `--cluster`: Enable protein clustering at 99% identity (default: False)

## How it works

1. Identifies all open reading frames (ORFs) in the target genome using getorf
2. Downloads all complete, annotated genomes for the specified taxon (can be a species or an entire genus) from NCBI RefSeq
3. Extracts and merges protein sequences from all assemblies
4. Optionally: Clusters proteins at 99% identity using mmseqs2 to create a non-redundant dataset
5. Searches the (clustered or unclustered) proteins against your target genome using mmseqs2 (tblastn-like search)
6. Outputs results in tabular format with alignment statistics

## Documentation

```bash
cd docs

# Create the html documentation
make html
```

The HTML will be available in docs/build/html/index.html.

You can also view the latest online documentation here:
https://floto-lab.github.io/pseudomancer/
