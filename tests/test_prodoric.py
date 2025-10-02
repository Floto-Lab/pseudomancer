"""
Integration and unit tests for pseudomancer.motifdb (Prodoric2 provider).

This suite contains:
- Offline unit tests that mock HTTP calls to ensure deterministic parsing of JSON payloads.
- Optional online integration tests (marked `@pytest.mark.online`) that query
  the live Prodoric2 API. These are skipped automatically if unreachable.

Usage
-----
Run only offline tests (default):
    pytest -q tests/test_prodoric.py

Run including online tests:
    pytest -q -m online tests/test_prodoric.py
"""

import json
import pytest
import requests

from typing import Optional
from pathlib import Path
from pseudomancer import motifdb


class DummyResponse:
    """Minimal stand-in for `requests.Response` used by the mock.

    Parameters
    ----------
    json_data : dict, optional
        Parsed JSON payload to be returned by `.json()`.
    status : int, optional
        HTTP status code to simulate (default: 200).

    Notes
    -----
    This mock is deliberately minimal: it only implements
    `.raise_for_status()` and `.json()` to satisfy the
    subset of the interface used by `motifdb.fetch_prodoric_motif`.
    """

    def __init__(self, json_data: dict = None, status: int = 200) -> None:
        self._json = json_data
        self._status = status

    def raise_for_status(self) -> None:
        """Raise `requests.HTTPError` if the mock status indicates failure."""
        if self._status >= 400:
            raise requests.HTTPError(f"HTTP {self._status}")

    def json(self) -> dict:
        """Return parsed JSON payload (may be None)."""
        return self._json


@pytest.fixture
def mock_get(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Patch `requests.get` with deterministic mock responses.

    Returns
    -------
    dict
        A dictionary capturing the last called URL under key `"url"`,
        useful for downstream assertions.

    Behavior
    --------
    - Requests ending with `/matrix/MX000001` return a small but valid JSON
      payload including a dummy PWM block.
    - Any other URL returns a 404-like response.
    """

    called = {}

    def _fake_get(url: str, timeout: Optional[float] = None) -> DummyResponse:
        called["url"] = url
        if url.endswith("/matrix/MX000001"):
            return DummyResponse(json_data={
                "id": "MX000001",
                "name": "DummyTF",
                "pwm": {
                    "A": [1, 0, 0],
                    "C": [0, 1, 0],
                    "G": [0, 0, 1],
                    "T": [0, 0, 0],
                },
            })
        return DummyResponse(status=404)

    monkeypatch.setattr(requests, "get", _fake_get)
    return called


def test_fetch_json(mock_get: dict) -> None:
    """Unit: JSON endpoint returns dict with expected keys and URL suffix.

    Uses the `mock_get` fixture to patch `requests.get`.

    Parameters
    ----------
    mock_get : dict
        Captures the last called URL under key `"url"`.
    """

    result = motifdb.fetch_prodoric_motif_json("MX000001")
    assert isinstance(result, dict)
    assert result["id"] == "MX000001"
    assert "name" in result
    assert mock_get["url"].endswith("/matrix/MX000001")


def test_load_prodoric_json_from_dict(mock_get: dict) -> None:
    """Unit: `load_prodoric_json` builds a Motif from a dict.

    This test fetches a mocked PRODORIC2 JSON payload and
    passes it to `load_prodoric_json`. The resulting motif
    object should have the expected name, length, and counts.

    Parameters
    ----------
    mock_get : dict
        Captures the last called URL under key `"url".
    """

    data = motifdb.fetch_prodoric_motif_json("MX000001")
    motif_dict = motifdb.load_prodoric_json(data)
    
    assert isinstance(motif_dict, dict)
    
    m = next(iter(motif_dict.values()))

    assert m.name == "DummyTF"
    assert m.length == 3
    assert "A" in m.counts and "C" in m.counts


def test_load_prodoric_json_from_file(tmp_path: Path) -> None:
    """Unit: `load_prodoric_json` also works from a JSON file.

    This test writes a minimal PWM to a temporary JSON file,
    loads it, and checks that the Motif is constructed properly.

    Parameters
    ----------
    tmp_path : pathlib.Path
        Temporary directory provided by pytest.
    """

    fake_json = {
        "id": "MX999999",
        "name": "FileTF",
        "pwm": {
            "A": [2, 0],
            "C": [0, 2],
            "G": [0, 0],
            "T": [0, 0],
        },
    }
    json_path = tmp_path / "test.json"
    json_path.write_text(json.dumps(fake_json))

    motif_dict = motifdb.load_prodoric_json(json_path)
    assert "FileTF" in motif_dict
    m = motif_dict["FileTF"]
    assert m.length == 2
    assert m.counts["A"][0] == 2.0


@pytest.mark.online
def test_fetch_real_motif_json() -> None:
    """Integration (online): hit live Prodoric2 API and ensure JSON looks sane.

    This test will be skipped automatically if the Prodoric
    API is not reachable within the timeout window.
    """

    tf_id = "MX000001"
    try:
        data = motifdb.fetch_prodoric_motif_json(tf_id, timeout=10)
    except requests.RequestException as e:
        pytest.skip(f"Prodoric API not reachable: {e}")
    assert isinstance(data, dict)
    assert "id" in data
    assert "name" in data
