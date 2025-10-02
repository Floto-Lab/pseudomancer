"""
Motif database access and management for Pseudomancer.
"""

from __future__ import annotations

import json
import requests

from Bio import motifs
from pathlib import Path
from typing import Dict, Optional, Union, Mapping, Sequence


def fetch_prodoric_motif(
        tf_id: str,
        base_url: str = "https://www.prodoric.de/api",
        timeout: float = 15.0,
        session: Optional[requests.Session] = None,
    ) -> dict:
    """Fetch a PRODORIC2 motif for a transcription factor by ID (JSON only).

    Parameters
    ----------
    tf_id : str
        PRODORIC matrix ID (e.g., 'MX000001').
    base_url : str, optional
        Base URL for the PRODORIC API (default: "https://www.prodoric.de/api").
    timeout : float, optional
        Request timeout in seconds (default: 15).
    session : requests.Session, optional
        Reuse a requests session for performance if calling repeatedly.

    Returns
    -------
    dict
        Motif data as a Python dictionary.
    """

    endpoint = f"{base_url}/matrix/{tf_id}"
    req = session.get if session else requests.get
    resp = req(endpoint, timeout=timeout)
    resp.raise_for_status()

    return resp.json()


def fetch_prodoric_motif_json(tf_id: str, **kwargs) -> dict:
    """Wrapper: fetch motif as JSON (dict).

    Parameters
    ----------
    tf_id : str
        PRODORIC matrix ID (e.g., 'MX000001').
    **kwargs
        Additional arguments passed to `fetch_prodoric_motif`.

    Returns
    -------
    dict
        Motif data as a Python dictionary.
    """

    return fetch_prodoric_motif(tf_id, **kwargs)


def load_prodoric_json(data: Union[dict, list, str, Path]) -> Dict[str, motifs.Motif]:
    """Load one or more PWMs from a PRODORIC2 JSON entry and return Biopython Motifs.

    Parameters
    ----------
    data : dict, list, str, or pathlib.Path
        Either:
        - a single parsed JSON dict (one motif),
        - a list of such dicts (multiple motifs),
        - or a path/str to a JSON file containing one of the above.

    Returns
    -------
    dict
        Mapping ``{motif_name: Bio.motifs.Motif}``.
    """

    # Load from file if path or string path
    if isinstance(data, (str, Path)):
        with open(data) as f:
            data = json.load(f)

    # Ensure list for uniform handling
    if isinstance(data, Mapping):
        entries = [data]
    elif isinstance(data, Sequence) and not isinstance(data, (str, bytes)):
        entries = data
    else:
        raise ValueError("Expected dict, list of dicts, or path to JSON file.")

    result: Dict[str, motifs.Motif] = {}

    for entry in entries:
        pwm = entry.get("pwm")
        if pwm is None:
            raise ValueError(f"Invalid Prodoric JSON: missing 'pwm' field in {entry.get('id', 'unknown')}")

        def get_col(key: str) -> Sequence[float]:
            if key in pwm:
                col = pwm[key]
            elif key.lower() in pwm:
                col = pwm[key.lower()]
            else:
                raise ValueError(f"PWM missing nucleotide column '{key}' in {entry.get('id', 'unknown')}.")
            return [float(x) for x in col]

        A, C, G, T = map(get_col, ["A", "C", "G", "T"])

        lengths = {len(A), len(C), len(G), len(T)}
        if len(lengths) != 1 or next(iter(lengths)) == 0:
            raise ValueError(
                f"PWM columns must all be same non-zero length in {entry.get('id', 'unknown')} "
                f"(found A={len(A)}, C={len(C)}, G={len(G)}, T={len(T)})."
            )

        counts = {"A": A, "C": C, "G": G, "T": T}
        m = motifs.Motif(counts=counts)
        m.name = entry.get("name", entry.get("id", "unknown"))

        result[m.name] = m

    return result
