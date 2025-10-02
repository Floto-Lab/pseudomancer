"""
ORF overlap analysis for pseudomancer pipeline.
"""

import os
import pandas as pd
from typing import List, Dict, Any


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF attributes string into a dictionary.

    Parameters
    ----------
    attr_string : str
        GFF attributes string (9th column)
    
    Returns
    -------
    Dict[str, str]
        Dictionary of attribute key-value pairs
    """

    attributes = {}
    
    # Handle hit_details specially since it contains semicolons
    if "hit_details=" in attr_string:
        before_hit_details, rest = attr_string.split("hit_details=", 1)
        
        # Find where hit_details ends (look for next attribute with =)
        hit_details_end = -1
        for i, char in enumerate(rest):
            if char == ";" and "=" in rest[i+1:i+50]:
                hit_details_end = i
                break
        
        if hit_details_end == -1:
            hit_details_value = rest
            after_hit_details = ""
        else:
            hit_details_value = rest[:hit_details_end]
            after_hit_details = rest[hit_details_end+1:]
        
        # Parse before and after normally
        for part in [before_hit_details, after_hit_details]:
            if part:
                for item in part.split(";"):
                    if "=" in item:
                        key, value = item.split("=", 1)
                        attributes[key] = value
        
        attributes["hit_details"] = hit_details_value
    else:
        # Normal parsing if no hit_details
        for item in attr_string.split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                attributes[key] = value
    
    return attributes


def parse_hit_details(hit_details: str) -> List[Dict[str, Any]]:
    """Parse hit_details string to extract individual hit information.
    
    Parameters
    ----------
    hit_details : str
        hit_details attribute string from GFF
    
    Returns
    -------
    List[Dict[str, Any]]
        List of hits with keys: q_coords, t_start, t_end, frame, strand
    """

    hits = []
    for hit in hit_details.split(";"):
        if not hit.strip():
            continue
        
        # Parse: q146-189:t2913374-2913508:f1:e8.38e-09
        parts = hit.split(":")
        if len(parts) >= 4:
            try:
                q_coords = parts[0].replace("q", "")
                t_coords = parts[1].replace("t", "")
                mmseqs_frame = int(parts[2].replace("f", ""))
                evalue = parts[3]
                if evalue.startswith("e"):
                    evalue = evalue[1:]
                
                t_start, t_end = map(int, t_coords.split("-"))
                
                # Determine strand from frame
                strand = "+" if mmseqs_frame > 0 else "-"
                
                # Convert MMseqs2 frame to biological frame
                if mmseqs_frame > 0:
                    biological_frame = mmseqs_frame + 1
                    if biological_frame > 3:
                        biological_frame = 1
                else:
                    biological_frame = mmseqs_frame + 1
                    if biological_frame == 0:
                        biological_frame = -3
                
                hits.append({
                    "q_coords": q_coords,
                    "t_start": t_start,
                    "t_end": t_end,
                    "frame": biological_frame,
                    "strand": strand,
                    "evalue": evalue,
                })
            except (ValueError, IndexError):
                continue
    return hits


def find_overlapping_orfs(hit: Dict[str, Any], orfs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Find ORFs that overlap with a hit in the same strand and frame.
    
    Parameters
    ----------
    hit : Dict[str, Any]
        Hit information with keys: t_start, t_end, frame, strand
    orfs : List[Dict[str, Any]]
        List of ORFs with keys: start, end, strand, frame

    Returns
    -------
    List[Dict[str, Any]]
        List of overlapping ORFs with overlap details
    """

    overlapping = []
    
    hit_start = min(hit["t_start"], hit["t_end"])
    hit_end = max(hit["t_start"], hit["t_end"])
    
    for orf in orfs:
        # Check strand and frame match
        if orf["strand"] != hit["strand"] or orf["frame"] != hit["frame"]:
            continue
        
        # Check for overlap
        orf_start = min(orf["start"], orf["end"])
        orf_end = max(orf["start"], orf["end"])
        
        overlap_start = max(hit_start, orf_start)
        overlap_end = min(hit_end, orf_end)
        
        if overlap_start <= overlap_end:
            overlap_length = overlap_end - overlap_start + 1
            overlapping.append({
                "orf": orf,
                "overlap_length": overlap_length,
                "overlap_start": overlap_start,
                "overlap_end": overlap_end,
            })
    
    return overlapping


def parse_getorf_gff(getorf_gff: str) -> List[Dict[str, Any]]:
    """Parse getorf GFF output to extract ORF information.

    Parameters
    ----------
    getorf_gff : str
        Path to getorf GFF file
    
    Returns
    -------
    List[Dict[str, Any]]
        List of ORFs with keys: orf_id, seqid, start, end,
    """

    orfs = []
    
    with open(getorf_gff, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            
            seqid = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Only process ORF features
            if feature_type.lower() in ["orf", "cds"]:
                # Extract ORF ID from attributes
                orf_id = None
                for attr in attributes.split(";"):
                    if attr.startswith("ID="):
                        orf_id = attr.split("=", 1)[1]
                        break
                
                if not orf_id:
                    orf_id = f"orf_{seqid}_{start}_{end}_{strand}"
                
                # Calculate reading frame from coordinates and strand
                if strand == "+":
                    frame = ((start - 1) % 3) + 1
                else:
                    frame = -((end - 1) % 3 + 1)
                
                orfs.append({
                    "orf_id": orf_id,
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "frame": frame,
                    "length": end - start + 1,
                    "attributes": attributes,
                })
    
    return orfs


def overlap_with_getorf(pseudogene_gff: str, getorf_gff: str, output_dir: str) -> str:
    """Find overlapping ORFs for pseudogene candidates.
    
    Parameters
    ----------
    pseudogene_gff : str
        Path to pseudogene candidate GFF file
    getorf_gff : str
        Path to getorf GFF file
    output_dir : str
        Directory to store output TSV file

    Returns
    -------
    str
        Path to the output TSV file with overlap results
    """

    print(f"Finding ORF overlaps between {pseudogene_gff} and {getorf_gff}")
    
    # Parse getorf ORFs
    orfs = parse_getorf_gff(getorf_gff)
    print(f"Parsed {len(orfs)} ORFs from getorf output")
    
    # Parse pseudogene clusters
    clusters = []
    with open(pseudogene_gff, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            
            seqid = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = parse_gff_attributes(fields[8])
            
            # Only process pseudogene features with hit_details
            if feature_type == "pseudogene" and "hit_details" in attributes:
                clusters.append({
                    "seqid": seqid,
                    "cluster_start": start,
                    "cluster_end": end,
                    "cluster_strand": strand,
                    "cluster_id": attributes.get("ID", f"cluster_{len(clusters)}"),
                    "cluster_name": attributes.get("Name", "Unknown"),
                    "hit_details": attributes["hit_details"],
                    "attributes": attributes,
                })
    
    print(f"Found {len(clusters)} pseudogene clusters")
    
    # Process each cluster
    results = []
    
    for cluster in clusters:
        # Parse hits from this cluster
        hits = parse_hit_details(cluster["hit_details"])
        
        # For each hit, find overlapping ORFs
        for hit_idx, hit in enumerate(hits):
            hit_id = f"{cluster['cluster_id']}_hit_{hit_idx + 1}"
            
            # Find overlapping ORFs in same strand and frame
            overlapping_orfs = find_overlapping_orfs(hit, orfs)
            
            if overlapping_orfs:
                # If multiple ORFs overlap, choose the longest one
                best_orf = max(overlapping_orfs, key=lambda x: x["orf"]["length"])
                
                results.append({
                    "cluster_id": cluster["cluster_id"],
                    "cluster_name": cluster["cluster_name"],
                    "hit_id": hit_id,
                    "hit_start": hit["t_start"],
                    "hit_end": hit["t_end"],
                    "hit_frame": hit["frame"],
                    "hit_strand": hit["strand"],
                    "hit_evalue": hit["evalue"],
                    "orf_id": best_orf["orf"]["orf_id"],
                    "orf_start": best_orf["orf"]["start"],
                    "orf_end": best_orf["orf"]["end"],
                    "orf_frame": best_orf["orf"]["frame"],
                    "orf_strand": best_orf["orf"]["strand"],
                    "orf_length": best_orf["orf"]["length"],
                    "overlap_length": best_orf["overlap_length"],
                    "overlap_start": best_orf["overlap_start"],
                    "overlap_end": best_orf["overlap_end"],
                    "num_overlapping_orfs": len(overlapping_orfs),
                })
            else:
                # Still record it but with no ORF info
                results.append({
                    "cluster_id": cluster["cluster_id"],
                    "cluster_name": cluster["cluster_name"],
                    "hit_id": hit_id,
                    "hit_start": hit["t_start"],
                    "hit_end": hit["t_end"],
                    "hit_frame": hit["frame"],
                    "hit_strand": hit["strand"],
                    "hit_evalue": hit["evalue"],
                    "orf_id": None,
                    "orf_start": None,
                    "orf_end": None,
                    "orf_frame": None,
                    "orf_strand": None,
                    "orf_length": None,
                    "overlap_length": None,
                    "overlap_start": None,
                    "overlap_end": None,
                    "num_overlapping_orfs": 0,
                })
    
    # Convert to DataFrame and save TSV
    df = pd.DataFrame(results)
    output_file = os.path.join(output_dir, "cluster_orf_mappings.tsv")
    df.to_csv(output_file, sep="\t", index=False)
    
    print(f"Results saved to {output_file}")
    print(f"Total hit-ORF mappings: {len(results)}")
    print(f"Hits with matching ORFs: {len(df[df['orf_id'].notna()])}")
    print(f"Hits without matching ORFs: {len(df[df['orf_id'].isna()])}")
    
    return output_file
