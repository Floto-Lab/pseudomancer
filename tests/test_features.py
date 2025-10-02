import json
from pathlib import Path

import pytest

import pseudomancer.features as features_mod
import pseudomancer.motifdb as motifdb_mod

from pseudomancer.features import find_shine_dalgarno, find_promoter, find_start_codons
from pseudomancer.motifdb import load_prodoric_json


def _make_prodoric_entry(name: str, pattern: str) -> dict:
    """Construct a minimal Prodoric-style PWM entry for a literal pattern."""

    length = len(pattern)
    counts = {base: [0.0] * length for base in "ACGT"}
    for idx, base in enumerate(pattern):
        counts[base][idx] = 10.0

    return {
        "id": name,
        "name": name,
        "pwm": counts,
    }


def test_find_start_codons_forward_all_and_best() -> None:
    """Detect canonical and alternative start codons on the forward strand."""

    seq = "AAAATGAAAGTGTTG"

    best = find_start_codons(seq, 0, len(seq), strand_name="forward")
    assert best is not None
    assert best["codon"] == "ATG"
    assert best["absolute_pos"] == seq.index("ATG")
    assert best["distance"] == best["absolute_pos"]

    all_hits = find_start_codons(seq, 0, len(seq), strand_name="forward", return_all=True)
    assert isinstance(all_hits, list)
    assert [hit["codon"] for hit in all_hits] == ["ATG", "GTG", "TTG"]


def test_find_start_codons_reverse_detection() -> None:
    """Detect start codons on the reverse strand (forward sequence contains complements)."""

    seq = "CATCACCAA"  # complements of ATG, GTG, TTG on the reverse strand

    best = find_start_codons(seq, 0, len(seq), strand_name="reverse")
    assert best is not None
    assert best["codon"] == "ATG"
    assert best["absolute_pos"] == seq.index("CAT")
    assert best["distance"] == len(seq) - (best["absolute_pos"] + 3)

    all_hits = find_start_codons(seq, 0, len(seq), strand_name="reverse", return_all=True)
    assert [hit["codon"] for hit in all_hits] == ["ATG", "GTG", "TTG"]


def test_find_start_codons_invalid_region() -> None:
    """Invalid window coordinates should raise an informative error."""

    seq = "ATG"
    with pytest.raises(ValueError):
        find_start_codons(seq, 5, 2, strand_name="forward")


def test_find_shine_dalgarno_forward_exact() -> None:
    """Exact mode: canonical Shine-Dalgarno on forward strand."""

    seq = "TTTTTAGGAGGATGAAA"
    start = seq.find("ATG")
    result = find_shine_dalgarno(seq, start, strand_name="forward", window=12, mode="exact")

    assert result is not None
    assert result["motif"] == "AGGAGG"
    assert result["strand_name"] == "forward"
    assert "score" not in result


def test_find_shine_dalgarno_forward_pwm() -> None:
    """PWM mode: canonical Shine-Dalgarno on forward strand."""

    seq = "TTTTTAGGAGGATGAAA"
    start = seq.find("ATG")
    result = find_shine_dalgarno(seq, start, strand_name="forward", window=12, mode="pwm", threshold=0.0)

    assert result is not None
    assert result["motif"] in ["AGGAGG", "GGAG", "AGGA", "AAGGAG", "GGAGG"]
    assert result["strand_name"] == "forward"
    assert "score" in result


def test_find_shine_dalgarno_reverse_exact() -> None:
    """Exact mode: Shine-Dalgarno detection on reverse strand."""

    seq = "TTTTTATGAAACCTCCT"
    start = seq.find("ATG")
    result = find_shine_dalgarno(seq, start, strand_name="reverse", window=12, mode="exact")

    assert result is not None
    assert result["motif"] == "AGGAGG"
    assert result["strand_name"] == "reverse"
    assert result["distance"] == result["absolute_pos"] - start


def test_no_shine_dalgarno_exact() -> None:
    """Exact mode: no Shine-Dalgarno motif present."""

    seq = "TTTTTCCCCCCATGAAA"
    start = seq.find("ATG")
    result = find_shine_dalgarno(seq, start, strand_name="forward", window=12, mode="exact")

    assert result is None


def test_edge_case_near_start_exact() -> None:
    """Exact mode: start codon at beginning (edge case)."""

    seq = "ATGAGGAGGTTTT"
    start = 0
    result = find_shine_dalgarno(seq, start, strand_name="forward", window=20, mode="exact")

    assert result is None


def test_variant_detected_with_pwm_thresholds() -> None:
    """PWM mode: weak variant only detected with permissive threshold."""

    seq = "TTTTTAGGACGATGAAA"
    start = seq.find("ATG")

    # stricter threshold ensures weak site is not called
    strict = find_shine_dalgarno(seq, start, strand_name="forward", window=12, mode="pwm", threshold=6.0)
    loose = find_shine_dalgarno(seq, start, strand_name="forward", window=12, mode="pwm", threshold=-2.0)

    assert strict is None
    assert loose is not None
    assert "score" in loose


def test_custom_sd_list_exact() -> None:
    """Exact mode: custom Shine-Dalgarno motif list."""

    seq = "TTTTTAGGAATGAAA"
    start = seq.find("ATG")

    result = find_shine_dalgarno(
        seq,
        start,
        strand_name="forward",
        window=12,
        mode="exact",
        sd_sequences=["AGGA"],
    )

    assert result is not None
    assert result["motif"] == "AGGA"

def test_forward_absolute_position_and_distance() -> None:
    """Check absolute position and distance calculation on forward strand."""

    seq = "TTTTTAGGAGGATGAAA"
    start = seq.find("ATG")
    result = find_shine_dalgarno(seq, start, strand_name="forward", window=12, mode="exact")

    assert result is not None
    assert result["absolute_pos"] == seq.find("AGGAGG")
    assert result["distance"] == start - result["absolute_pos"]

def test_return_all_hits() -> None:
    """Verify multiple hits are returned when return_all=True."""

    seq = "TTTTTAGGAGGAGGAATGAAA"
    start = seq.find("ATG")
    results = find_shine_dalgarno(seq, start, strand_name="forward", window=15, mode="exact", return_all=True)

    assert isinstance(results, list)
    assert len(results) > 1
    assert all("motif" in r for r in results)

def test_pwm_threshold_respected_strict_vs_loose() -> None:
    """With a weak SD-like site, a loose threshold should detect it,
    while a stricter threshold (just above the detected score) should not.
    """
    
    seq = "TTTTTAGGACGATGAAA"  # weak-ish SD variant upstream of ATG
    start = seq.find("ATG")

    # Loose threshold: should find something and report a score
    loose = find_shine_dalgarno(
        seq, start, strand_name="forward",
        window=12, mode="pwm", threshold=-2.0
    )
    assert loose is not None
    assert "score" in loose

    # Now set a threshold a hair above the found score → should yield no hit
    strict_cutoff = loose["score"] + 1e-6
    strict = find_shine_dalgarno(
        seq, start, strand_name="forward",
        window=12, mode="pwm", threshold=strict_cutoff
    )
    assert strict is None


def test_pwm_threshold_override_default() -> None:
    """A very high explicit threshold should override the default (0.8*max_score)
    and return no hits even when a normal call would find one.
    """

    seq = "TTTTTAGGAGGATGAAA"  # strong SD upstream of ATG
    start = seq.find("ATG")

    # Normal permissive search should find a hit
    normal = find_shine_dalgarno(
        seq, start, strand_name="forward",
        window=12, mode="pwm", threshold=0.0
    )
    assert normal is not None
    assert "score" in normal

    # An absurdly high threshold should force no hits
    none_hit = find_shine_dalgarno(
        seq, start, strand_name="forward",
        window=12, mode="pwm", threshold=1e9
    )
    assert none_hit is None


def test_find_promoter_forward_pwm_hit() -> None:
    """find_promoter detects a hand-crafted PWM on the forward strand."""

    seq = "TATAATTTTTATGAAA"
    start = seq.find("ATG")
    pwm_entry = _make_prodoric_entry("TestProm", "TATAAT")

    result = find_promoter(
        seq,
        start,
        strand_name="forward",
        db_path=pwm_entry,
        window=12,
        threshold=-1.0,
    )

    assert result is not None
    assert result["motif_name"] == "TestProm"
    assert result["motif_seq"] == "TATAAT"
    assert result["strand_name"] == "forward"
    assert result["absolute_pos"] == 0
    assert result["distance"] == start


def test_find_promoter_reverse_coordinates() -> None:
    """Reverse-strand promoter hit reports correct absolute coordinates."""

    seq = "GGGGATGATTATAGGG"
    start = seq.find("ATG")
    pwm_entry = _make_prodoric_entry("ReverseProm", "TATAAT")

    result = find_promoter(
        seq,
        start,
        strand_name="reverse",
        db_path=pwm_entry,
        window=12,
        threshold=-1.0,
    )

    assert result is not None
    assert result["motif_name"] == "ReverseProm"
    assert result["motif_seq"] == "TATAAT"
    assert result["strand_name"] == "reverse"
    assert result["relative_pos"] == 3
    assert result["absolute_pos"] == seq.index("ATTATA")
    assert result["distance"] == result["absolute_pos"] - start


def test_find_promoter_accepts_preloaded_motifs() -> None:
    """Passing a preloaded motif dictionary should avoid re-parsing JSON."""

    pwm_entry = _make_prodoric_entry("CachedProm", "TTGACA")
    motif_dict = load_prodoric_json(pwm_entry)

    seq = "TTGACAACCCATG"
    start = seq.find("ATG")

    result = find_promoter(
        seq,
        start,
        strand_name="forward",
        db_path=motif_dict,
        window=12,
        threshold=-1.0,
    )

    assert result is not None
    assert result["motif_name"] == "CachedProm"
    assert result["motif_seq"] == "TTGACA"
    assert result["distance"] == start - result["absolute_pos"]


def test_find_promoter_caches_loaded_database(tmp_path: Path, monkeypatch: any) -> None:
    """Repeated calls with the same path should hit the cached motifs only once."""

    features_mod._MOTIF_DB_CACHE.clear()
    features_mod._PWM_CACHE.clear()

    pwm_entry = _make_prodoric_entry("CachedFile", "TATAAT")
    json_path = tmp_path / "motifs.json"
    json_path.write_text(json.dumps(pwm_entry))

    call_count = {"n": 0}
    original_loader = motifdb_mod.load_prodoric_json

    def wrapped_loader(data):
        call_count["n"] += 1
        return original_loader(data)

    monkeypatch.setattr(motifdb_mod, "load_prodoric_json", wrapped_loader)

    seq = "TATAATTTTTATGAAA"
    start = seq.find("ATG")

    features_mod.find_promoter(seq, start, "forward", db_path=json_path, window=12, threshold=-1.0)
    features_mod.find_promoter(seq, start, "forward", db_path=json_path, window=12, threshold=-1.0)

    assert call_count["n"] == 1
