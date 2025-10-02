"""
MMseqs2 search functionality for pseudomancer pipeline.
"""

import subprocess
import os


def run_mmseqs2_search(protein_db: str, genome_file: str, output_dir: str, evalue: float = 1e-5) -> str:
    """Run mmseqs2 search (tblastn-like) against genome.
    
    Parameters
    ----------
    protein_db : str
        Path to clustered protein sequences
    genome_file : str
        Path to target genome FASTA file
    output_dir : str
        Directory for output files
    evalue : float
        E-value threshold
        
    Returns
    -------
    str
        Path to search results file
    """
    
    mmseqs_dir = os.path.join(output_dir, "mmseqs_dbs")
    results_dir = os.path.join(output_dir, "mmseqs_results")
    os.makedirs(mmseqs_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    protein_db_path = os.path.join(mmseqs_dir, "protein_db")
    genome_db_path = os.path.join(mmseqs_dir, "genome_db")

    print("Creating mmseqs2 databases...")
    subprocess.run(["mmseqs", "createdb", protein_db, protein_db_path], check=True)
    subprocess.run(["mmseqs", "createdb", genome_file, genome_db_path], check=True)

    search_result = os.path.join(results_dir, "search_result")
    tmp_dir = os.path.join(mmseqs_dir, "tmp")

    print("Running mmseqs2 search...")
    search_cmd = [
        "mmseqs",
        "search",
        protein_db_path,
        genome_db_path,
        search_result,
        tmp_dir,
        "--translation-table",
        "11",
        "-e",
        str(evalue),
    ]
    subprocess.run(search_cmd, check=True)

    # Convert to readable format
    output_file = os.path.join(results_dir, "search_results.txt")
    convert_cmd = [
        "mmseqs",
        "convertalis",
        protein_db_path,
        genome_db_path,
        search_result,
        output_file,
        "--format-output",
        "query,qstart,qend,qcov,target,tstart,tend,tframe,evalue,pident,alnlen,qlen",
    ]
    subprocess.run(convert_cmd, check=True)

    print(f"Search results saved to {output_file}")
    return output_file
