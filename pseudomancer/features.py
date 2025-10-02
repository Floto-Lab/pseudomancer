"""
Feature extraction functions for genomic sequences, including motif searches
and promoter identification using PWMs from motif databases.
"""

from Bio import motifs
from Bio.Seq import Seq
from Bio.motifs.matrix import PositionSpecificScoringMatrix
from pathlib import Path
from typing import Dict, List, Optional, Union, Mapping, Tuple


# Caches
_MOTIF_DB_CACHE: Dict[Tuple[str, str], Dict[str, motifs.Motif]] = {}
_PWM_CACHE: Dict[Tuple[Tuple[Tuple[float, ...], ...], float], PositionSpecificScoringMatrix] = {}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        Input DNA sequence.

    Returns
    -------
    str
        Reverse complement of the input sequence.
    """

    return str(Seq(seq).reverse_complement())


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate the Hamming distance between two strings.
    
    Parameters
    ----------
    s1 : str
        First string.
    s2 : str
        Second string.

    Returns
    -------
    int
        Hamming distance (number of differing characters).
    """

    if len(s1) != len(s2):
        raise ValueError("hamming_distance requires strings of equal length.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def search_motifs(
        sequence: str,
        region_start: int,
        region_end: int,
        strand_name: str,
        mode: str,
        motifs_list: List[str],
        threshold: Optional[float] = None,
        pseudocounts: float = 0.1,
        mismatches: int = 0
    ) -> List[Dict]:
    """Generic motif search (exact or PWM) in a defined genomic region.

    Parameters
    ----------
    sequence : str
        Full DNA sequence.
    region_start : int
        Start coordinate of the region (inclusive).
    region_end : int
        End coordinate of the region (exclusive).
    strand_name : str
        "forward" or "reverse".
    mode : str
        "exact" or "pwm". Exact matches or PWM scoring.
    motifs_list : list of str
        Motifs to search or seed for PWM.
    threshold : float, optional
        Minimum PWM score cutoff.
    pseudocounts : float
        For PWM normalization.
    mismatches : int
        Allowed mismatches in exact mode.

    Returns
    -------
    list of dict
        Each dict contains motif, relative_pos, absolute_pos, strand_name, distance, (score if PWM).
    """

    sequence = sequence.upper()
    hits = []

    if strand_name == "forward":
        subseq = sequence[region_start:region_end]
        anchor = region_start
    elif strand_name == "reverse":
        subseq = reverse_complement(sequence[region_start:region_end])
        anchor = region_end  # reverse anchor
    else:
        raise ValueError("strand_name must be 'forward' or 'reverse'")

    # Normalise motif seeds to uppercase to match the uppercased search sequence
    motifs_list = [m.upper() for m in motifs_list]

    if mode == "exact":
        for motif in motifs_list:
            mlen = len(motif)
            for i in range(len(subseq) - mlen + 1):
                window = subseq[i:i + mlen]
                if hamming_distance(window, motif) <= mismatches:
                    if strand_name == "forward":
                        abs_pos = anchor + i
                    else:
                        abs_pos = anchor - (i + mlen)
                    hits.append({
                        "motif": window,
                        "strand_name": strand_name,
                        "relative_pos": i,
                        "absolute_pos": abs_pos,
                    })

    elif mode == "pwm":
        grouped = {}
        for s in motifs_list:
            grouped.setdefault(len(s), []).append(s)

        for length, group in grouped.items():
            inst = [Seq(s) for s in group]
            m = motifs.create(inst)
            pwm = m.counts.normalize(pseudocounts=pseudocounts).log_odds()
            thres = threshold if threshold is not None else 0.8 * pwm.max_score()

            for i in range(len(subseq) - length + 1):
                frag = Seq(subseq[i:i + length])
                score = pwm.calculate(frag)
                if score >= thres:
                    if strand_name == "forward":
                        abs_pos = anchor + i
                    else:
                        abs_pos = anchor - (i + length)
                    hits.append({
                        "motif": str(frag),
                        "strand_name": strand_name,
                        "relative_pos": i,
                        "absolute_pos": abs_pos,
                        "score": score,
                    })
    else:
        raise ValueError("mode must be 'exact' or 'pwm'")

    return hits


def find_shine_dalgarno(
        sequence: str,
        start: int,
        strand_name: str,
        window: int = 20,
        mode: str = "exact",
        sd_sequences: List[str] = ["AGGAGG", "GGAG", "AGGA", "AAGGAG", "GGAGG"],
        threshold: Optional[float] = None,
        pseudocounts: float = 0.1,
        return_all: bool = False,
        mismatches: int = 0
    ) -> Optional[Union[Dict, List[Dict]]]:
    """Wrapper to find Shine-Dalgarno motifs upstream of a start codon.

    Parameters
    ----------
    sequence : str
        Full DNA sequence.
    start : int
        Start codon position (0-based).
    strand_name : str
        "forward" or "reverse".
    window : int
        Search window size upstream of start codon.
    mode : str
        "exact" or "pwm". Exact matches or PWM scoring.
    sd_sequences : list of str
        Known Shine-Dalgarno motifs.
    threshold : float, optional
        Minimum PWM score cutoff.
    pseudocounts : float
        For PWM normalization.
    return_all : bool
        If True, return all hits; else return best hit.
    mismatches : int
        Allowed mismatches in exact mode.

    Returns
    -------
    dict or list of dict or None
        Best hit dict, list of all hits, or None if no hits found.
    """

    if strand_name == "forward":
        region_start = max(0, start - window)
        region_end = start
    else:
        region_start = start
        region_end = min(len(sequence), start + window)

    hits = search_motifs(
        sequence, region_start, region_end, strand_name,
        mode, sd_sequences, threshold, pseudocounts, mismatches
    )

    # compute distance to start codon
    for h in hits:
        if strand_name == "forward":
            h["distance"] = start - h["absolute_pos"]
        else:
            h["distance"] = h["absolute_pos"] - start

    if not hits:
        return None
    if return_all:
        return hits
    if mode == "exact":
        return max(hits, key=lambda h: (len(h["motif"]), -abs(h["distance"])))
    else:
        return max(hits, key=lambda h: h["score"])


def load_pwms(db_path: Union[str, Path], fmt: str = "meme") -> Dict[str, motifs.Motif]:
    """Load promoter motifs (PWMs) from a file.

    Parameters
    ----------
    db_path : str or Path
        Path to motif database file (e.g. Prodoric2 export).
    fmt : str, optional
        File format: "meme" (default), "transfac", "jaspar", etc.

    Returns
    -------
    dict
        Mapping of motif name -> Biopython Motif object.
    """

    with open(db_path) as handle:
        motif_list = list(motifs.parse(handle, fmt))
    
    return {m.name: m for m in motif_list}


def _motif_signature(motif: motifs.Motif) -> Tuple[Tuple[float, ...], ...]:
    """Create a hashable signature of a motif's counts for caching.
    
    Parameters
    ----------
    motif : Bio.motifs.Motif
        Motif object.

    Returns
    -------
    tuple of tuples
        Signature representing the motif's counts.
    """

    return tuple(tuple(float(x) for x in motif.counts[base]) for base in "ACGT")


def _get_cached_log_odds(motif: motifs.Motif, pseudocounts: float) -> PositionSpecificScoringMatrix:
    """Get or compute the log-odds PWM for a motif, with caching.

    Parameters
    ----------
    motif : Bio.motifs.Motif
        Motif object.
    pseudocounts : float
        For PWM normalization.

    Returns
    -------
    Bio.motifs.matrix.PositionSpecificScoringMatrix
        Log-odds PWM.
    """

    key = (_motif_signature(motif), float(pseudocounts))
    pwm = _PWM_CACHE.get(key)

    if pwm is None:
        pwm = motif.counts.normalize(pseudocounts=pseudocounts).log_odds()
        _PWM_CACHE[key] = pwm

    return pwm


def _cache_key_for_db(db_path: Union[str, Path], fmt: str) -> Optional[Tuple[str, str]]:
    """Create a cache key for a motif database path and format.

    Parameters
    ----------
    db_path : str or Path
        Path to motif database file.
    fmt : str
        File format.

    Returns
    -------
    tuple or None
        Cache key (path, fmt) or None if db_path is not a path.
    """

    if not isinstance(db_path, (str, Path)):
        return None
    
    path = Path(db_path).expanduser()
    
    try:
        path = path.resolve()
    except OSError:
        path = path.resolve(strict=False)

    return (str(path), fmt)


def find_promoter(
        sequence: str,
        start: int,
        strand_name: str,
        db_path: Union[str, Path, dict, Mapping[str, motifs.Motif]],
        fmt: str = "json",  # <-- now default is Prodoric2 JSON
        window: int = 120,
        sigma: Optional[str] = None,
        threshold: Optional[float] = None,
        pseudocounts: float = 0.1,
        return_all: bool = False,
    ) -> Optional[Union[Dict, List[Dict]]]:
    """
    Identify bacterial promoters upstream of an ORF using PWMs from a database.
    By default this uses Prodoric2 JSON matrices, which is the recommended format.

    Parameters
    ----------
    sequence : str
        Full DNA sequence.
    start : int
        Start codon position (0-based).
    strand_name : str
        "forward" or "reverse".
    db_path : str, Path, or Mapping
        Path to motif database file, a Prodoric2 JSON dict, or a preloaded motif mapping.
    fmt : str, default "json"
        File format. Options:

        * "json": Prodoric2 JSON (preferred and default).
        * "meme", "transfac", "jaspar": other Biopython-supported formats
          (discouraged, kept for legacy use).
    window : int, default 120
        Search window size upstream of start codon.
    sigma : str, optional
        If set, filter motifs to those associated with this sigma factor (e.g., "sigma70").
    threshold : float, optional
        Minimum PWM score cutoff.
    pseudocounts : float, default 0.1
        For PWM normalization.
    return_all : bool, default False
        If True, return all hits; else return best hit.

    Returns
    -------
    dict or list of dict or None
        Best hit dict, list of all hits, or None if no hits found.
    """

    from pseudomancer.motifdb import load_prodoric_json

    sequence = sequence.upper()
    if strand_name == "forward":
        region_start = max(0, start - window)
        region_end = start
        upstream = sequence[region_start:region_end]
        anchor = region_start
    elif strand_name == "reverse":
        region_start = start
        region_end = min(len(sequence), start + window)
        upstream = reverse_complement(sequence[region_start:region_end])
        anchor = region_end
    else:
        raise ValueError("strand_name must be 'forward' or 'reverse'.")

    # Load motifs depending on format
    motif_dict: Dict[str, motifs.Motif]

    def _is_preloaded(mapping: Mapping) -> bool:
        if not mapping:
            return True
        return all(isinstance(v, motifs.Motif) for v in mapping.values())

    cache_key = _cache_key_for_db(db_path, fmt)

    if isinstance(db_path, Mapping) and _is_preloaded(db_path):
        motif_dict = dict(db_path)
    elif cache_key and cache_key in _MOTIF_DB_CACHE:
        motif_dict = _MOTIF_DB_CACHE[cache_key]
    elif fmt == "json":
        motif_dict = load_prodoric_json(db_path)
        if cache_key:
            _MOTIF_DB_CACHE[cache_key] = motif_dict
    else:
        import warnings
        warnings.warn(
            f"Non-JSON motif format '{fmt}' is discouraged. "
            "Please use Prodoric2 JSON whenever possible.",
            UserWarning,
        )
        motif_dict = load_pwms(db_path, fmt=fmt)
        if cache_key:
            _MOTIF_DB_CACHE[cache_key] = motif_dict

    if sigma:
        motif_dict = {k: v for k, v in motif_dict.items() if sigma in k}
        if not motif_dict:
            raise ValueError(f"No motifs for sigma factor '{sigma}' in {db_path}")

    hits = []
    for name, pwm in motif_dict.items():
        pwm = _get_cached_log_odds(pwm, pseudocounts)
        cutoff = threshold if threshold is not None else 0.8 * pwm.max_score()

        for pos in range(len(upstream) - pwm.length + 1):
            frag = Seq(upstream[pos:pos + pwm.length])
            score = pwm.calculate(frag)
            if score >= cutoff:
                abs_pos = anchor + pos if strand_name == "forward" else anchor - (pos + pwm.length)
                hits.append({
                    "motif_name": name,
                    "motif_seq": str(frag),
                    "strand_name": strand_name,
                    "relative_pos": pos,
                    "absolute_pos": abs_pos,
                    "score": score,
                })

    for h in hits:
        if strand_name == "forward":
            h["distance"] = start - h["absolute_pos"]
        else:
            h["distance"] = h["absolute_pos"] - start

    if not hits:
        return None
    if return_all:
        return hits
    return max(hits, key=lambda h: h["score"])


def find_start_codons(
        sequence: str,
        region_start: int,
        region_end: int,
        strand_name: str,
        codons: Optional[List[str]] = None,
        mismatches: int = 0,
        return_all: bool = False,
    ) -> Optional[Union[Dict, List[Dict]]]:
    """Locate putative start codons within a genomic window.

    Parameters
    ----------
    sequence : str
        Full genomic DNA sequence.
    region_start : int
        Start coordinate (inclusive) of the window to inspect.
    region_end : int
        End coordinate (exclusive) of the window to inspect.
    strand_name : str
        "forward" or "reverse".
    codons : list of str, optional
        Start codon triplets to consider. Defaults to ``["ATG", "GTG", "TTG"]``.
        Earlier entries in the list are treated as higher priority when
        selecting a single best hit.
    mismatches : int, optional
        Maximum number of mismatches permitted when matching codons
        (defaults to 0 for exact matches).
    return_all : bool, optional
        If ``True`` return a list of all detected codons. Otherwise, return the
        highest priority hit.

    Returns
    -------
    dict or list of dict or None
        A single hit (default), a list of hits, or ``None`` if no codons are
        detected. Each hit dictionary contains ``codon``, ``absolute_pos``,
        ``relative_pos`` (offset within the inspected window), ``distance``
        (offset to the window boundary: ``region_start`` for the forward strand
        or ``region_end`` for the reverse strand), and ``strand_name``.
    """

    if region_start < 0 or region_end > len(sequence) or region_start >= region_end:
        raise ValueError("Invalid region boundaries provided to find_start_codons.")

    codon_list = tuple(c.upper() for c in (codons or ["ATG", "GTG", "TTG"]))
    
    hits = search_motifs(
        sequence,
        region_start,
        region_end,
        strand_name,
        mode="exact",
        motifs_list=list(codon_list),
        mismatches=mismatches,
    )

    if not hits:
        return None

    preference = {codon: idx for idx, codon in enumerate(codon_list)}

    for hit in hits:
        hit["codon"] = hit.pop("motif")
        if strand_name == "forward":
            hit["distance"] = hit["absolute_pos"] - region_start
        else:
            hit["distance"] = region_end - (hit["absolute_pos"] + len(hit["codon"]))

    hits.sort(key=lambda h: (preference.get(h["codon"], len(codon_list)), h["relative_pos"]))

    if return_all:
        return hits
    return hits[0]
