"""
Protein clustering with mmseqs2 for pseudomancer pipeline.
"""

import subprocess
import sys
import os


def cluster_proteins_mmseqs2(protein_file: str, output_dir: str, identity: float = 0.99) -> str:
    """Cluster proteins using mmseqs2 at specified identity threshold.
    
    Parameters
    ----------
    protein_file : str
        Path to input protein FASTA file
    output_dir : str
        Directory for output files
    identity : float
        Sequence identity threshold (default: 0.99)
        
    Returns
    -------
    str
        Path to representative sequences FASTA file
    """

    base_name = os.path.splitext(os.path.basename(protein_file))[0]
    output_prefix = os.path.join(output_dir, f"{base_name}_clustered")
    
    print(f"Clustering proteins at {identity*100}% identity...")
    
    cmd = [
        "mmseqs",
        "easy-cluster",
        protein_file,
        output_prefix,
        output_dir,
        "--min-seq-id",
        str(identity),
        "--cov-mode",
        "0",
        "-c",
        "0.8",
        "-s",
        "7",
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to cluster proteins: {e.stderr}\n")
        sys.exit(1)
    
    rep_seq_file = f"{output_prefix}_rep_seq.fasta"
    print(f"Representative sequences saved to {rep_seq_file}")
    return rep_seq_file
