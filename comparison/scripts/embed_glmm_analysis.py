#!/usr/bin/env python3
"""Validate and embed the independent GLMM-analysis JSON as portable Scala."""
import hashlib
import json
import math
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRANCHE = ROOT / "comparison" / "glmm-analysis-tranche"
FRAMES = ROOT / "comparison" / "glmm-tranche" / "frames"
OUT = ROOT / "src" / "test" / "scala" / "mixeff4s" / "glmm" / "GlmmAnalysisFrozen.scala"
REQUIRED = ("cbpp", "culcitalogreg", "synthetic_poisson_offset")

def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def finite(value):
    if isinstance(value, float): return math.isfinite(value)
    if isinstance(value, list): return all(finite(x) for x in value)
    if isinstance(value, dict): return all(finite(x) for x in value.values())
    return True

def q(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)

def add_glmm_analysis_pins(catalog):
    refs = json.loads((TRANCHE / "references.json").read_text())
    integrals = json.loads((TRANCHE / "integrals.json").read_text())
    formulas = {"cbpp": "proportion ~ 1 + period + (1 | herd)", "culcitalogreg": "predation ~ 1 + ttt + (1 | block)",
                "synthetic_poisson_offset": "count ~ 1 + x + (1 | group)"}
    key = lambda row: tuple(row[k] for k in ("dataset", "formula", "family", "link", "estimator", "engine"))
    for name in REQUIRED:
        shift = integrals["integrals"][name]["response_only_saturated_minus2loglik"]
        for nodes in (5, 9, 15):
            fit = refs["fits"][name][f"nAGQ{nodes}"]
            family = "Binomial" if name == "cbpp" else "Bernoulli" if name == "culcitalogreg" else "Poisson"
            pin = dict(dataset=name, formula=formulas[name], family=family, link="Log" if family == "Poisson" else "Logit",
                       estimator=f"joint_agq_{nodes}", engine="lme4", status="ok", n_obs=integrals["integrals"][name]["n"],
                       n_agq=nodes, objective=fit["minus2log_lme4_raw"]+shift,
                       objective_definition="minus_two_loglikelihood_with_response_constants", objective_comparable=True,
                       beta=fit["beta"], coef_names=fit["beta_names"], theta=fit["theta"],
                       tolerances=dict(objective=1e-4, beta=1e-3, theta=1e-3),
                       raw_reference_objective=fit["minus2log_lme4_raw"], response_only_normalization_shift=shift,
                       provenance="comparison/glmm-analysis-tranche/reference-provenance.json")
            existing = [row for row in catalog["results"] if key(row) == key(pin)]
            if existing: assert existing == [pin], f"conflicting AGQ pin {key(pin)}"
            else: catalog["results"].append(pin)
    catalog["source"]["glmm_analysis"] = dict(R=refs["source"]["R"], lme4=refs["source"]["lme4"],
        provenance="comparison/glmm-analysis-tranche/reference-provenance.json",
        normalization="raw high-order lme4 objective plus independently checked saturated-density constants")
    return catalog

def embed_workflow_frames():
    import csv
    def vector(items): return "Vector(" + ", ".join(items) + ")"
    binomial_path = TRANCHE / "frames/binomial_workflow.csv"
    rare_path = TRANCHE / "frames/rare_event_bernoulli.csv"
    if not binomial_path.exists(): return
    rows = list(csv.DictReader(binomial_path.open()))
    code = "package mixeff4s.fixtures\n\nimport mixeff4s.data.ModelFrame\n\n/** Frozen synthetic workflow and rare-event diagnostic frames. */\nobject GlmmAnalysisFrames:\n"
    code += "  val trials: Vector[Int] = " + vector([str(int(r["trials"])) for r in rows]) + "\n"
    code += "  val successes: Vector[Double] = " + vector([str(float(r["successes"])) for r in rows]) + "\n"
    code += "  val response: Vector[Double] = successes.zip(trials).map((k, n) => k / n.toDouble)\n"
    code += '  def binomial: ModelFrame = ModelFrame.of(\n    "y" -> ModelFrame.numeric(response),\n    "x" -> ModelFrame.numeric(' + vector([str(float(r["x"])) for r in rows]) + '),\n    "g" -> ModelFrame.factor(' + vector([q(r["group"]) for r in rows]) + ')\n  ).toOption.get\n'
    if rare_path.exists():
        rare = list(csv.DictReader(rare_path.open()))
        assert len(rare) == 120 and sum(int(r["y"]) for r in rare) == 10
        code += '\n  def rare: ModelFrame = ModelFrame.of(\n    "y" -> ModelFrame.numeric(' + vector([str(float(r["y"])) for r in rare]) + '),\n    "exposure" -> ModelFrame.numeric(' + vector([str(float(r["exposure"])) for r in rare]) + '),\n    "group" -> ModelFrame.factor(' + vector([q(r["group"]) for r in rare]) + ')\n  ).toOption.get\n'
    (ROOT / "src/test/scala/mixeff4s/fixtures/GlmmAnalysisFrames.scala").write_text(code)

def main() -> None:
    references = json.loads((TRANCHE / "references.json").read_text())
    integrals = json.loads((TRANCHE / "integrals.json").read_text())
    assert references["schema_version"] == integrals["schema_version"] == 1
    assert set(references["fits"]) == set(integrals["integrals"]) == set(REQUIRED)
    for dataset in REQUIRED:
        fits = references["fits"][dataset]
        assert set(fits) == {"nAGQ1", "nAGQ5", "nAGQ9", "nAGQ15"}
        assert all(finite(row) for row in fits.values())
        assert all(row["n_agq"] in (1, 5, 9, 15) for row in fits.values())
        high = fits["nAGQ15"]
        integral = integrals["integrals"][dataset]
        assert finite(integral) and integral["max_relative_error_bound"] <= 1e-10
        # The expected saturated correction is exact for binomial trials and
        # Poisson counts. Residual records AGQ quadrature error, not a hidden
        # objective redefinition.
        shift_error = abs(high["full_density_normalization_shift_measured"] - high["full_density_normalization_shift_expected"])
        assert shift_error <= 2e-3, (dataset, shift_error)
    for oracle in integrals["fixed_parameter_oracles"]:
        assert finite(oracle) and oracle["relative_error_bound"] <= 1e-10
    for dataset in REQUIRED:
        assert (FRAMES / f"{dataset}.csv").is_file()

    provenance = {
        "schema_version": 1,
        "frame_sha256": {name: digest(FRAMES / f"{name}.csv") for name in REQUIRED},
        "artifact_sha256": {
            "references.json": digest(TRANCHE / "references.json"),
            "integrals.json": digest(TRANCHE / "integrals.json"),
        },
        "generator_sha256": {
            "freeze_glmm_analysis.R": digest(ROOT / "comparison" / "scripts" / "freeze_glmm_analysis.R"),
            "embed_glmm_analysis.py": digest(Path(__file__)),
        },
        "factor_coding": references["factor_coding"],
        "normalization_policy": "Retain lme4 raw -2 logLik and full conditional-density -2 log integral; compare only after measured response-only saturated correction.",
        "controls": references["controls"],
        "python": platform.python_version(),
    }
    (TRANCHE / "reference-provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    payload = json.dumps({"references": references, "integrals": integrals}, indent=2, ensure_ascii=False)
    chunks = [payload[i:i + 12000] for i in range(0, len(payload), 12000)]
    scala = """package mixeff4s.glmm

/** Generated independent R/AGQ and full-density integration witness. No R dependency in default tests. */
private[mixeff4s] object GlmmAnalysisFrozen:
  val json: String = Vector(
"""
    scala += ",\n".join(f"    \"\"\"{chunk}\"\"\"" for chunk in chunks)
    scala += "\n  ).mkString\n"
    OUT.write_text(scala)
    embed_workflow_frames()
    print(f"validated {len(REQUIRED)} AGQ datasets and {len(integrals['fixed_parameter_oracles'])} fixed-parameter integration oracles")

if __name__ == "__main__":
    main()
