"""
Pseudogene candidate identification and processing for pseudomancer pipeline.
"""

import os
import math
import pandas as pd
from typing import List, Dict, Any, Tuple
from collections import defaultdict


def get_genome_info(genome_file: str) -> Tuple[str, int]:
    """Extract sequence ID and length from genome FASTA file.
    
    Parameters
    ----------
    genome_file : str
        Path to the genome FASTA file.

    Returns
    -------
    Tuple[str, int]
        A tuple containing the sequence ID and the total length of the genome.
    """

    with open(genome_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                # Extract sequence ID (first word after >)
                seq_id = line.strip()[1:].split()[0]
                break
        else:
            raise ValueError(f"No sequences found in {genome_file}")

    # Calculate total length
    total_length = 0
    with open(genome_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                total_length += len(line.strip())

    return seq_id, total_length


def calculate_distance(start1: int, end1: int, start2: int, end2: int) -> int:
    """Calculate distance between two genomic regions.

    Parameters
    ----------
    start1 : int
        Start position of the first region.
    end1 : int
        End position of the first region.
    start2 : int
        Start position of the second region.
    end2 : int
        End position of the second region.

    Returns
    -------
    int
        Distance between the two regions. Returns 0 if they overlap.
    """

    return max(0, max(min(start1, end1), min(start2, end2)) - min(max(start1, end1), max(start2, end2)))


def are_adjacent(start1: int, end1: int, start2: int, end2: int, max_distance: int = 100) -> bool:
    """Check if two genomic regions are adjacent within max_distance.

    Parameters
    ----------
    start1 : int
        Start position of the first region.
    end1 : int
        End position of the first region.
    start2 : int
        Start position of the second region.
    end2 : int
        End position of the second region.

    max_distance : int, optional
        Maximum distance to consider regions as adjacent, by default 100.
    """

    return calculate_distance(start1, end1, start2, end2) <= max_distance


def group_adjacent_hits(hits_df: pd.DataFrame, max_distance: int = 100) -> pd.DataFrame:
    """Group adjacent hits within max_distance.

    Parameters
    ----------
    hits_df : pd.DataFrame
        DataFrame containing hits with 'tstart' and 'tend' columns.
    max_distance : int, optional
        Maximum distance to consider hits as adjacent, by default 100.

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional 'cluster' column indicating cluster membership.
    """

    if len(hits_df) == 0:
        return hits_df

    # Sort by genomic position
    hits_df = hits_df.sort_values(by=["tstart", "tend"])
    hits_df = hits_df.reset_index(drop=True)

    # Initialize clustering
    clusters = []
    hits_df["cluster"] = 0
    current_cluster = 1

    hits_df.loc[0, "cluster"] = current_cluster

    for i in range(1, len(hits_df)):
        assigned = False

        # Check against all existing clusters
        for cluster_id in range(1, current_cluster + 1):
            cluster_hits = hits_df[hits_df["cluster"] == cluster_id]

            # Check if current hit is adjacent to any hit in this cluster
            for j in cluster_hits.index:
                if are_adjacent(
                    hits_df.loc[i, "tstart"],
                    hits_df.loc[i, "tend"],
                    hits_df.loc[j, "tstart"],
                    hits_df.loc[j, "tend"],
                    max_distance,
                ):
                    hits_df.loc[i, "cluster"] = cluster_id
                    assigned = True
                    break
            if assigned:
                break

        # If not assigned to any existing cluster, create new cluster
        if not assigned:
            current_cluster += 1
            hits_df.loc[i, "cluster"] = current_cluster

    return hits_df


def create_cluster_summary(
    cluster_hits: pd.DataFrame, query_name: str, event_type: str, strand_name: str
) -> Dict[str, Any]:
    """Create summary for a cluster of hits.
    
    Parameters
    ----------
    cluster_hits : pd.DataFrame
        DataFrame containing hits in the cluster.
    query_name : str
        Name of the query protein.
    event_type : str
        Type of event ("frameshift" or "nonsense").
    strand_name : str
        Strand name ("forward" or "reverse").

    Returns
    -------
    Dict[str, Any]
        Dictionary summarizing the cluster.
    """

    # Sort hits by genomic position
    cluster_hits = cluster_hits.sort_values(by=["tstart", "tend"])

    # Calculate cluster boundaries
    cluster_start = min(cluster_hits["tstart"].min(), cluster_hits["tend"].min())
    cluster_end = max(cluster_hits["tstart"].max(), cluster_hits["tend"].max())

    # Calculate query boundaries
    query_start = cluster_hits["qstart"].min()
    query_end = cluster_hits["qend"].max()

    # Get reading frames involved
    frames = ",".join(map(str, sorted(cluster_hits["tframe"].unique())))

    # Calculate best e-value and mean percent identity
    best_evalue = cluster_hits["evalue"].min()
    mean_pident = cluster_hits["pident"].mean()

    # Create hit details string
    hit_details = []
    for _, hit in cluster_hits.iterrows():
        hit_detail = (
            f"q{hit['qstart']}-{hit['qend']}:t{hit['tstart']}-{hit['tend']}:f{hit['tframe']}:e{hit['evalue']:.2e}"
        )
        hit_details.append(hit_detail)

    return {
        "query": query_name,
        "target": cluster_hits["target"].iloc[0],
        "cluster_id": cluster_hits["cluster"].iloc[0],
        "event_type": event_type,
        "strand": strand_name,
        "num_hits": len(cluster_hits),
        "frames": frames,
        "query_start": query_start,
        "query_end": query_end,
        "query_span": query_end - query_start + 1,
        "target_start": cluster_start,
        "target_end": cluster_end,
        "target_span": cluster_end - cluster_start + 1,
        "best_evalue": best_evalue,
        "mean_pident": mean_pident,
        "hit_details": ";".join(hit_details),
    }


def process_query(query_data: pd.DataFrame) -> List[Dict[str, Any]]:
    """Process a single query to identify frameshift and nonsense candidates.

    Parameters
    ----------
    query_data : pd.DataFrame
        DataFrame containing hits for a single query.
    
    Returns
    -------
    List[Dict[str, Any]]
        List of candidate summaries.
    """

    query_name = query_data["query"].iloc[0]

    # If only one hit, return empty list
    if len(query_data) == 1:
        return []

    # Separate hits by strand
    forward_hits = query_data[query_data["tframe"] > 0]
    reverse_hits = query_data[query_data["tframe"] < 0]

    candidates = []

    # Process forward strand hits
    if len(forward_hits) > 1:
        forward_clustered = group_adjacent_hits(forward_hits)
        for cluster_id in forward_clustered["cluster"].unique():
            cluster_hits = forward_clustered[forward_clustered["cluster"] == cluster_id]

            if len(cluster_hits) == 1:
                continue  # Skip single-hit clusters

            # Check if cluster contains multiple reading frames
            unique_frames = cluster_hits["tframe"].abs().unique()

            if len(unique_frames) > 1:
                # Multiple reading frames - frameshift candidate
                candidate = create_cluster_summary(cluster_hits, query_name, "frameshift", "forward")
                candidates.append(candidate)
            else:
                # Same reading frame - nonsense candidate
                candidate = create_cluster_summary(cluster_hits, query_name, "nonsense", "forward")
                candidates.append(candidate)

    # Process reverse strand hits
    if len(reverse_hits) > 1:
        reverse_clustered = group_adjacent_hits(reverse_hits)
        for cluster_id in reverse_clustered["cluster"].unique():
            cluster_hits = reverse_clustered[reverse_clustered["cluster"] == cluster_id]

            if len(cluster_hits) == 1:
                continue  # Skip single-hit clusters

            # Check if cluster contains multiple reading frames
            unique_frames = cluster_hits["tframe"].abs().unique()

            if len(unique_frames) > 1:
                # Multiple reading frames - frameshift candidate
                candidate = create_cluster_summary(cluster_hits, query_name, "frameshift", "reverse")
                candidates.append(candidate)
            else:
                # Same reading frame - nonsense candidate
                candidate = create_cluster_summary(cluster_hits, query_name, "nonsense", "reverse")
                candidates.append(candidate)

    return candidates


def cluster_overlapping_candidates(candidates_df: pd.DataFrame, min_overlap: int = 1) -> pd.DataFrame:
    """Cluster overlapping genomic regions from different query proteins.

    Parameters
    ----------
    candidates_df : pd.DataFrame
        DataFrame containing candidate summaries.
    min_overlap : int, optional
        Minimum overlap in base pairs to consider candidates as overlapping, by default 1.

    Returns
    -------
    pd.DataFrame
        DataFrame with clustered candidates.
    """

    if len(candidates_df) == 0:
        return candidates_df

    # Determine strand from frames column
    candidates_df["strand_direction"] = candidates_df["frames"].apply(
        lambda x: "reverse" if any(int(f) < 0 for f in str(x).split(",")) else "forward"
    )

    # Group by target and strand
    grouped = candidates_df.groupby(["target", "strand_direction"])

    merged_results = []

    for (target, strand_direction), group_df in grouped:
        if len(group_df) <= 1:
            # Single candidate - create merged entry
            for _, row in group_df.iterrows():
                merged_results.append(
                    {
                        "merged_id": f"merged_{target}_{strand_direction}_1",
                        "target": target,
                        "strand": strand_direction,
                        "target_start": row["target_start"],
                        "target_end": row["target_end"],
                        "target_span": row["target_span"],
                        "num_contributing_queries": 1,
                        "contributing_queries": row["query"],
                        "all_frames": row["frames"],
                        "best_evalue": row["best_evalue"],
                        "mean_pident": row["mean_pident"],
                        "max_pident": row["mean_pident"],
                        "total_hits": row["num_hits"],
                        "min_query_start": row["query_start"],
                        "max_query_end": row["query_end"],
                        "total_query_span": row["query_span"],
                        "event_type": row["event_type"],
                        "combined_hit_details": row["hit_details"],
                    }
                )
            continue

        # Find overlapping candidates using simple approach
        group_df = group_df.reset_index(drop=True)
        visited = [False] * len(group_df)
        overlap_group = 0

        for i in range(len(group_df)):
            if visited[i]:
                continue

            overlap_group += 1
            current_group = [i]
            visited[i] = True

            # Find all overlapping candidates
            for j in range(i + 1, len(group_df)):
                if visited[j]:
                    continue

                # Check if j overlaps with any in current_group
                for k in current_group:
                    overlap_start = max(group_df.iloc[j]["target_start"], group_df.iloc[k]["target_start"])
                    overlap_end = min(group_df.iloc[j]["target_end"], group_df.iloc[k]["target_end"])

                    if overlap_start <= overlap_end and (overlap_end - overlap_start + 1) >= min_overlap:
                        current_group.append(j)
                        visited[j] = True
                        break

            # Create merged candidate for this overlap group
            group_candidates = group_df.iloc[current_group]

            merged_results.append(
                {
                    "merged_id": f"merged_{target}_{strand_direction}_{overlap_group}",
                    "target": target,
                    "strand": strand_direction,
                    "target_start": group_candidates["target_start"].min(),
                    "target_end": group_candidates["target_end"].max(),
                    "target_span": group_candidates["target_end"].max() - group_candidates["target_start"].min() + 1,
                    "num_contributing_queries": len(group_candidates),
                    "contributing_queries": ",".join(group_candidates["query"].unique()),
                    "all_frames": ",".join(
                        sorted(set(frame for frames in group_candidates["frames"] for frame in str(frames).split(",")))
                    ),
                    "best_evalue": group_candidates["best_evalue"].min(),
                    "mean_pident": group_candidates["mean_pident"].mean(),
                    "max_pident": group_candidates["mean_pident"].max(),
                    "total_hits": group_candidates["num_hits"].sum(),
                    "min_query_start": group_candidates["query_start"].min(),
                    "max_query_end": group_candidates["query_end"].max(),
                    "total_query_span": group_candidates["query_span"].sum(),
                    "event_type": ",".join(group_candidates["event_type"].unique()),
                    "combined_hit_details": ";".join(group_candidates["hit_details"]),
                }
            )

    return pd.DataFrame(merged_results)


def export_candidates_to_gff(frameshift_df: pd.DataFrame, output_file: str, genome_file: str) -> None:
    """Export frameshift candidates to GFF3 format.

    Parameters
    ----------
    frameshift_df : pd.DataFrame
        DataFrame containing frameshift candidate summaries.
    output_file : str
        Path to the output GFF3 file.
    genome_file : str
        Path to the genome FASTA file.
    """

    # Get genome information
    seq_id, genome_length = get_genome_info(genome_file)

    # Create GFF3 header
    gff_lines = ["##gff-version 3", f"##sequence-region {seq_id} 1 {genome_length}"]

    # Process frameshift candidates
    for i, row in frameshift_df.iterrows():
        # Create display name
        gene_display = row["contributing_queries"].split(";")[0].replace(".", "_")

        # Create attributes string
        attributes = [
            f"ID={row['merged_id']}",
            f"Name={gene_display}",
            f"Note=Merged frameshift candidate with {row['num_contributing_queries']} contributing queries, {row['total_hits']} total hits in frames {row['all_frames']}",
            f"merged_id={row['merged_id']}",
            f"num_contributing_queries={row['num_contributing_queries']}",
            f"contributing_queries={row['contributing_queries']}",
            f"total_hits={row['total_hits']}",
            f"all_frames={row['all_frames']}",
            f"total_query_span={row['total_query_span']}",
            f"best_evalue={row['best_evalue']:.2e}",
            f"mean_pident={row['mean_pident']:.1f}",
            f"max_pident={row['max_pident']:.1f}",
            f"event_type={row['event_type']}",
            f"hit_details={row['combined_hit_details']}",
        ]

        # Create GFF line
        # Calculate score, handling edge case where evalue is 0
        if row["best_evalue"] <= 0:
            score = "999.0"  # Use maximum score for perfect matches
        else:
            score = str(round(-math.log10(row["best_evalue"]), 1))

        gff_line = "\t".join(
            [
                row["target"],  # seqid
                "MMseqs2",  # source
                "pseudogene",  # type
                str(int(row["target_start"])),  # start
                str(int(row["target_end"])),  # end
                score,  # score
                "+" if row["strand"] == "forward" else "-",  # strand
                ".",  # phase
                ";".join(attributes),  # attributes
            ]
        )

        gff_lines.append(gff_line)

    # Write to file
    with open(output_file, "w") as f:
        f.write("\n".join(gff_lines) + "\n")

    print(f"GFF3 file written to: {output_file}")
    print(f"Features exported: {len(frameshift_df)} merged frameshift candidates")


def process_mmseqs2_results(results_file: str, output_dir: str, genome_file: str) -> str:
    """Process mmseqs2 results to identify pseudogene candidates and generate GFF.

    Parameters
    ----------
    results_file : str
        Path to the mmseqs2 results TSV file.
    output_dir : str
        Directory to save output files.
    genome_file : str
        Path to the genome FASTA file.

    Returns
    -------
    str
        Path to the generated GFF file.
    """

    print(f"Processing mmseqs2 results from {results_file}")

    # Read mmseqs2 results
    column_names = [
        "query",
        "qstart",
        "qend",
        "qcov",
        "target",
        "tstart",
        "tend",
        "tframe",
        "evalue",
        "pident",
        "alnlen",
        "qlen",
    ]
    mmseqs_results = pd.read_csv(results_file, sep="\t", names=column_names)

    print(f"Loaded {len(mmseqs_results)} mmseqs2 hits")

    # Process each query separately
    all_candidates = []

    for query in mmseqs_results["query"].unique():
        query_data = mmseqs_results[mmseqs_results["query"] == query]
        candidates = process_query(query_data)
        all_candidates.extend(candidates)

    print(f"Found {len(all_candidates)} individual candidates")

    # Convert to DataFrame
    candidates_df = pd.DataFrame(all_candidates)

    if len(candidates_df) == 0:
        print("No candidates found")
        return ""

    # Filter for frameshift candidates only
    frameshift_candidates = candidates_df[candidates_df["event_type"] == "frameshift"]

    print(f"Found {len(frameshift_candidates)} frameshift candidates")

    # Remove duplicates based on target_start
    frameshift_candidates = frameshift_candidates.drop_duplicates(subset=["target_start"])

    # Cluster overlapping candidates
    frameshift_merged = cluster_overlapping_candidates(frameshift_candidates)

    print(f"Merged into {len(frameshift_merged)} clustered candidates")

    # Generate GFF file
    gff_file = os.path.join(output_dir, "pseudogene_candidates_annotated_merged.gff")
    export_candidates_to_gff(frameshift_merged, gff_file, genome_file)

    return gff_file
