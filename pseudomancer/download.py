"""
NCBI data download and protein file handling for pseudomancer pipeline.
"""

import subprocess
import sys
import os
import zipfile


def download_ncbi_assemblies(taxon: str, output_dir: str) -> str:
    """Download NCBI RefSeq assemblies for a given taxon.
    
    Parameters
    ----------
    taxon : str
        Taxon name (e.g., 'Mycobacterium')
    output_dir : str
        Directory to store downloaded assemblies
        
    Returns
    -------
    str
        Path to the merged protein FASTA file
    """
    
    taxon_name = taxon.replace(" ", "_")
    assemblies_dir = os.path.join(output_dir, f"{taxon_name.lower()}_assemblies")
    os.makedirs(assemblies_dir, exist_ok=True)
    
    # Use datasets to download assemblies
    print(f"Downloading {taxon} assemblies from NCBI RefSeq...")
    cmd = [
        "datasets",
        "download",
        "genome",
        "taxon",
        taxon,
        "--annotated",
        "--reference",
        "--assembly-level",
        "complete",
        "--include",
        "genome,protein,gff3,gtf,seq-report",
        "--filename",
        os.path.join(assemblies_dir, f"{taxon_name.lower()}_dataset.zip"),
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to download assemblies: {e.stderr}\n")
        sys.exit(1)
    
    # Extract and merge protein files
    protein_files = []
    zip_path = os.path.join(assemblies_dir, f"{taxon_name.lower()}_dataset.zip")
    
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        zip_ref.extractall(assemblies_dir)
    
    # Find all protein.faa files
    for root, dirs, files in os.walk(assemblies_dir):
        for file in files:
            if file == "protein.faa":
                protein_files.append(os.path.join(root, file))
    
    # Merge protein files
    merged_proteins = os.path.join(output_dir, f"{taxon_name.lower()}_proteins.faa")
    with open(merged_proteins, "w") as outfile:
        for protein_file in protein_files:
            with open(protein_file, "r") as infile:
                outfile.write(infile.read())
    
    print(f"Merged {len(protein_files)} protein files into {merged_proteins}")
    return merged_proteins
