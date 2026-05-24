"""Parity gate against R destiny on Guo fixture."""
import json, sys
from pathlib import Path
import numpy as np
import pytest
import yaml

PORT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PORT))
sys.path.insert(0, str(PORT.parent / "omicverse-rebuildr" / "engine"))
from parity_metrics import compute_parity


@pytest.fixture(scope="session")
def manifest():
    return yaml.safe_load((PORT / "data" / "manifest.yaml").read_text())


@pytest.fixture(scope="session")
def outputs():
    ref = PORT / "data" / "reference_output.json"
    cand = PORT / "data" / "candidate_output.json"
    if not (ref.exists() and cand.exists()):
        pytest.skip("R reference / Py candidate outputs not generated")
    return json.loads(ref.read_text()), json.loads(cand.read_text())


def test_eigvec_procrustes(manifest, outputs):
    r, p = outputs
    proc = compute_parity(np.array(r["eigenvectors"]), np.array(p["eigenvectors"]), "embedding")
    spec = next(o for o in manifest["outputs"] if o["name"] == "diffusion_map")
    assert proc >= spec["threshold"], f"Procrustes {proc:.4f} < {spec['threshold']}"


def test_eigvalue_pearson(manifest, outputs):
    r, p = outputs
    ev_r = np.array(r["eigenvalues"]); ev_p = np.array(p["eigenvalues"])
    m = compute_parity(ev_r, ev_p, "ordinal")
    spec = next(o for o in manifest["outputs"] if o["name"] == "eigenvalues")
    assert m >= spec["threshold"], f"eigenvalue Pearson {m:.4f} < {spec['threshold']}"


def test_dpt_pearson(manifest, outputs):
    r, p = outputs
    dpt_r = np.array(r["dpt"]); dpt_p = np.array(p["dpt"])
    mask = np.isfinite(dpt_r) & np.isfinite(dpt_p)
    m = compute_parity(dpt_r[mask], dpt_p[mask], "ordinal")
    spec = next(o for o in manifest["outputs"] if o["name"] == "dpt")
    assert m >= spec["threshold"], f"DPT Pearson {m:.4f} < {spec['threshold']}"
