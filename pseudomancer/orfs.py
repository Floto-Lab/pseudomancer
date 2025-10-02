"""
ORF identification and GFF conversion for pseudomancer pipeline.
"""

import subprocess
import sys
import os
import re
from typing import Optional, Tuple


def parse_getorf_header(header: str) -> Optional[dict]:
    """Parse getorf FASTA header to extract coordinates.

    Example headers:
    Forward: >NC_002677.1_1 [93 - 272] Mycobacterium leprae TN, complete sequence
    Reverse: >NC_002677.1_53466 [1611 - 1420] (REVERSE SENSE) Mycobacterium leprae TN, complete sequence

    Parameters
    ----------
    header : str
        FASTA header line (without the leading '>')

    Returns
    -------
    Optional[dict]
        Dictionary with keys: seq_id, orf_num, start, end, strand, frame, length
        or None if parsing fails
    """

    match = re.match(r"^(.+)_(\d+)", header)
    if not match:
        return None

    seq_id = match.group(1)
    orf_num = match.group(2)

    coord_match = re.search(r"\[(\d+) - (\d+)\]", header)
    if not coord_match:
        return None

    start = int(coord_match.group(1))
    end = int(coord_match.group(2))

    if "(REVERSE SENSE)" in header:
        strand = "-"
    else:
        strand = "+"

    if strand == "+":
        frame = ((start - 1) % 3) + 1
    else:
        frame = -((end - 1) % 3 + 1)

    return {
        "seq_id": seq_id,
        "orf_num": orf_num,
        "start": start,
        "end": end,
        "strand": strand,
        "frame": frame,
        "length": abs(end - start) + 1,
    }


def convert_fasta_to_gff(fasta_file: str, gff_file: str) -> None:
    """Convert getorf FASTA output to GFF format.
    
    Parameters
    ----------
    fasta_file : str
        Path to getorf FASTA file
    gff_file : str
        Path to output GFF file
    """

    print(f"Converting {fasta_file} to GFF format...")

    gff_lines = ["##gff-version 3"]
    orf_count = 0

    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Parse header
                header = line.strip()[1:]  # Remove '>'
                orf_info = parse_getorf_header(header)
                if not orf_info:
                    print(f"Warning: Could not parse header: {header}")
                    continue

                orf_count += 1

                # Create GFF attributes
                orf_id = f"getorf_orf_{orf_info['orf_num']}"
                attributes = f"ID={orf_id};Name={orf_id};length={orf_info['length']};frame={orf_info['frame']}"

                # Create GFF line
                gff_line = "\t".join(
                    [
                        orf_info["seq_id"],  # seqid
                        "getorf",  # source
                        "ORF",  # type
                        str(orf_info["start"]),  # start
                        str(orf_info["end"]),  # end
                        ".",  # score
                        orf_info["strand"],  # strand
                        ".",  # phase
                        attributes,  # attributes
                    ]
                )

                gff_lines.append(gff_line)

    with open(gff_file, "w") as f:
        f.write("\n".join(gff_lines) + "\n")

    print(f"Converted {orf_count} ORFs to GFF format")
    print(f"GFF output written to {gff_file}")


def identify_orfs(genome_file: str, output_dir: str, minsize: int = 90, table: int = 11) -> Tuple[str, str]:
    """Identify open reading frames (ORFs) in the genome using getorf.

    Parameters
    ----------
    genome_file : str
        Path to genome FASTA file
    output_dir : str
        Directory for output files
    minsize : int
        Minimum ORF size in codons (default: 90)
    table : int
        Genetic code table (default: 11 for bacteria)

    Returns
    -------
    Tuple[str, str]
        Paths to ORF sequences file and GFF file
    """
    
    base_name = os.path.splitext(os.path.basename(genome_file))[0]
    orf_file = os.path.join(output_dir, f"{base_name}_orfs.fasta")
    gff_file = os.path.join(output_dir, f"{base_name}_orfs.gff")

    print(f"Running getorf ORF identification on {genome_file}...")

    cmd = [
        "getorf",
        "-sequence",
        genome_file,
        "-outseq",
        orf_file,
        "-find",
        "2",
        "-minsize",
        str(minsize),
        "-table",
        str(table),
        "-flanking",
        "0",
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"Failed to identify ORFs: {e.stderr}\n")
        sys.exit(1)

    print(f"ORF sequences saved to {orf_file}")

    convert_fasta_to_gff(orf_file, gff_file)

    return orf_file, gff_file
