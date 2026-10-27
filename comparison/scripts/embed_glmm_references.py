#!/usr/bin/env python3
"""Validate GLMM pins and generate portable test-only Scala fixtures."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRANCHE = ROOT / "comparison/glmm-tranche"
FRAMES = TRANCHE / "frames"
SPECS = {
    "cbpp": {"numeric": {"incidence", "size"}, "factors": {"herd": [str(i) for i in range(1, 16)], "period": [str(i) for i in range(1, 5)]}},
    "culcitalogreg": {"numeric": {"predation"}, "factors": {"block": [str(i) for i in range(1, 11)], "ttt": ["none", "crabs", "shrimp", "both"]}},
    "synthetic_poisson_offset": {"numeric": {"count", "x", "exposure"}, "factors": {"group": [f"g{i}" for i in range(1, 7)]}},
    "synthetic_poisson_crossed": {"numeric": {"count", "x", "exposure", "visit"}, "factors": {"site": [f"s{i}" for i in range(1, 13)], "batch": [f"b{i}" for i in range(1, 6)]}},
}

def scala_string(s): return json.dumps(s, ensure_ascii=False)
def vector(values): return "Vector(" + ", ".join(values) + ")"
def numeric(values): return vector([repr(float(v)) for v in values])

def frame_object():
    blocks = []
    for name, spec in SPECS.items():
        with (FRAMES / f"{name}.csv").open(newline="") as f: rows = list(csv.DictReader(f))
        columns = []
        for col in rows[0]:
            values = [r[col] for r in rows]
            if col in spec["numeric"]:
                expr = f'ModelFrame.numeric({numeric(values)})'
            else:
                levels = vector([scala_string(x) for x in spec["factors"][col]])
                expr = f'ModelFrame.factor({vector([scala_string(x) for x in values])}, {levels}).toOption.get'
            columns.append(f'      {scala_string(col)} -> {expr}')
        blocks.append(f'  val {name} = ModelFrame.of(\n' + ",\n".join(columns) + "\n    ).toOption.get")
    return "package mixeff4s.fixtures\n\nimport mixeff4s.data.ModelFrame\n\n/** Generated from comparison/glmm-tranche/frames; factor levels are explicit. */\nobject GlmmFrames:\n" + "\n\n".join(blocks) + "\n"

def finite(x):
    if isinstance(x, (int, float)): return math.isfinite(x)
    if isinstance(x, list): return all(finite(v) for v in x)
    if isinstance(x, dict): return all(finite(v) for v in x.values())
    return True

def main():
    refs = json.loads((TRANCHE / "references.json").read_text())
    assert refs["schema_version"] == 1 and len(refs["results"]) == 8
    for row in refs["results"]:
        assert row["status"] == "ok" and finite(row), row["dataset"]
        assert len(row["conditional_eta"]) == row["n"] == len(row["conditional_mu"])
        assert len(row["theta"]) == len(row["theta_names"])
        assert row["objective_minus2loglik_with_constants"] > 0
    with (FRAMES / "synthetic_poisson_offset.csv").open(newline="") as f:
        offset_rows = list(csv.DictReader(f))
    with (FRAMES / "synthetic_poisson_crossed.csv").open(newline="") as f:
        crossed_rows = list(csv.DictReader(f))
    assert len(offset_rows) == 24 and any(int(r["count"]) == 0 for r in offset_rows)
    assert len(crossed_rows) == 300 and any(int(r["count"]) == 0 for r in crossed_rows)
    provenance = {
        "schema_version": 1,
        "inputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(FRAMES.glob("*.csv"))},
        "upstream_inputs": {
            "cbpp": "/Users/bbuchsbaum/code/rust/mixeff-rs/datasets/cbpp/data.csv",
            "culcitalogreg": "/Users/bbuchsbaum/code/rust/mixeff-rs/datasets/culcitalogreg/data.csv",
        },
        "upstream_input_sha256": {
            "cbpp": "2798bc7704a4bdb8ade4e164d8c3a6397efa7ba5fb25c1ec11f517b4ffb940a7",
            "culcitalogreg": "0147190ecdea59319cf3143cee5e1ca5c2afcc94d88e00a8c88c321079e44138",
        },
        "generator_sha256": {
            "freeze_glmm_references.R": hashlib.sha256((ROOT / "comparison/scripts/freeze_glmm_references.R").read_bytes()).hexdigest(),
            "embed_glmm_references.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "factor_levels": {k: v["factors"] for k, v in SPECS.items()},
        "n_agq_0_policy": "documented_fast_divergence only; not lme4 joint-Laplace parity",
    }
    prior = TRANCHE / "prior-default-pwrss"
    if prior.exists():
        provenance["adverse_prior_default_pwrss"] = {
            "reason": "Default tolPwrss=1e-7 gave cbpp -2logLik 184.05313277509163; retained before tightening to 1e-10.",
            "artifacts_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(prior.glob("*")) if p.is_file()},
        }
    (TRANCHE / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    # Scala/JVM rejects a single modified-UTF8 constant above 65535 bytes.
    # Keep chunks comfortably below that limit while preserving byte-for-byte JSON.
    rendered = json.dumps(refs, indent=2)
    chunks = [rendered[i:i + 12000] for i in range(0, len(rendered), 12000)]
    frozen = "package mixeff4s.glmm\n\n/** Generated by comparison/scripts/freeze_glmm_references.R; no R dependency in tests. */\nprivate[mixeff4s] object GlmmFrozen:\n  val json: String = Vector(\n"
    frozen += ",\n".join(f'    \"\"\"{chunk}\"\"\"' for chunk in chunks)
    frozen += "\n  ).mkString\n"
    (ROOT / "src/test/scala/mixeff4s/glmm/GlmmFrozen.scala").write_text(frozen)
    (ROOT / "src/test/scala/mixeff4s/fixtures/GlmmFrames.scala").write_text(frame_object())
    print(f"validated {len(refs['results'])} finite GLMM references and generated Scala fixtures")

if __name__ == "__main__": main()
