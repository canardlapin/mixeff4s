#!/usr/bin/env python3
"""Copy small mixeff-rs frames and emit the portable Scala embedding."""

from __future__ import annotations

import json
import re
from pathlib import Path

RUST = Path("/Users/bbuchsbaum/code/rust/mixeff-rs")
ROOT = Path("/Users/bbuchsbaum/code/scala/mixeff4s")
REVISION = "d6b81ce8b5974b2f0e42fcf650ad6811e5202a4f"

VENDORED = [
    "dyestuff",
    "dyestuff2",
    "rail",
    "cake",
    "ergostool",
    "oats",
    "orthodont",
    "oxide",
    "machines",
    "station_season_duration",
    "singular",
]

EXTRACT_KEYS = {
    (
        "sleepstudy",
        "Reaction ~ 1 + Days + (1 + Days | Subject)",
        "Gaussian",
        "Identity",
        "ML",
    ),
    (
        "sleepstudy",
        "Reaction ~ 1 + Days + (1 + Days | Subject)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "penicillin",
        "diameter ~ 1 + (1 | plate) + (1 | sample)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "penicillin",
        "diameter ~ 1 + (1 | plate) + (1 | sample)",
        "Gaussian",
        "Identity",
        "ML",
    ),
    (
        "pastes",
        "strength ~ 1 + (1 | sample) + (1 | batch)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "dyestuff",
        "Yield ~ 1 + (1 | Batch)",
        "Gaussian",
        "Identity",
        "ML",
    ),
    (
        "dyestuff",
        "Yield ~ 1 + (1 | Batch)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "dyestuff2",
        "Yield ~ 1 + (1 | Batch)",
        "Gaussian",
        "Identity",
        "ML",
    ),
    (
        "dyestuff2",
        "Yield ~ 1 + (1 | Batch)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "rail",
        "travel ~ 1 + (1 | Rail)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "ergostool",
        "effort ~ 1 + Type + (1 | Subject)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "cake",
        "angle ~ 1 + recipe * temperature + (1 | recipe:replicate)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "oats",
        "yield ~ 1 + Variety * nitro + (1 | Block) + (1 | Block:Variety)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "orthodont",
        "distance ~ 1 + age * Sex + (1 + age | Subject)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "oxide",
        "Thickness ~ 1 + (1 | Lot) + (1 | Lot:Wafer)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "oxide",
        "Thickness ~ 1 + (1 | Lot/Wafer)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "machines",
        "score ~ 1 + Machine + (1 | Worker) + (1 | Worker:Machine)",
        "Gaussian",
        "Identity",
        "REML",
    ),
    (
        "machines",
        "score ~ 1 + Machine + (1 + Machine | Worker)",
        "Gaussian",
        "Identity",
        "REML",
    ),
}


def parse_levels(blob: str) -> list[str] | None:
    match = re.search(r"levels\s*=\s*\[(.*?)\]", blob, re.S)
    if not match:
        return None
    return re.findall(r'"([^"]*)"', match.group(1))


def parse_meta(name: str) -> dict:
    text = (RUST / "datasets" / name / "meta.toml").read_text()
    source = re.search(r'^source = "(.*)"', text, re.M).group(1)
    license_m = re.search(r'^license = "(.*)"', text, re.M)
    n_rows = int(re.search(r"^n_rows = (\d+)", text, re.M).group(1))
    columns = []
    for block in re.split(r"\n\[\[columns\]\]\n", text)[1:]:
        block = block.split("\n[[")[0]
        col_name = re.search(r'^name = "(.*)"', block, re.M).group(1)
        col_type = re.search(r'^type = "(.*)"', block, re.M).group(1)
        col = {"name": col_name, "type": col_type}
        levels = parse_levels(block)
        if levels is not None:
            col["levels"] = levels
        columns.append(col)
    return {
        "name": name,
        "n_rows": n_rows,
        "source": source,
        "license": license_m.group(1) if license_m else "",
        "mixeff_rs_revision": REVISION,
        "columns": columns,
    }


def scala_triple(text: str) -> str:
    if '"""' in text:
        raise SystemExit("CSV contains triple quotes")
    # Keep JSON/CSV bytes intact. Scala \\uXXXX in a triple-quoted string
    # would rewrite escapes; write raw Unicode instead.
    return '"""' + text + '"""'


def slim_result(record: dict, engine: str, comparable: bool) -> dict:
    out = {
        "dataset": record["dataset"],
        "formula": record["formula"],
        "family": record["family"],
        "link": record["link"],
        "estimator": record["estimator"],
        "engine": engine,
        "status": record["status"],
        "n_obs": record.get("n_obs"),
        "objective": record.get("objective"),
        "objective_definition": record.get("objective_definition"),
        "objective_comparable": comparable,
        "beta": record.get("beta"),
        "coef_names": record.get("coef_names"),
        "theta": record.get("theta"),
        "sigma": record.get("sigma"),
        "is_singular": record.get("is_singular"),
    }
    return {k: v for k, v in out.items() if v is not None}


def key_of(record: dict) -> tuple:
    return (
        record["dataset"],
        record["formula"],
        record["family"],
        record["link"],
        record["estimator"],
    )


def main() -> None:
    datasets_root = ROOT / "comparison" / "datasets"
    datasets_root.mkdir(parents=True, exist_ok=True)
    csvs = {}
    metas = {}
    for name in VENDORED:
        src = RUST / "datasets" / name
        dest = datasets_root / name
        dest.mkdir(parents=True, exist_ok=True)
        csv = (src / "data.csv").read_text()
        if not csv.endswith("\n"):
            csv += "\n"
        (dest / "data.csv").write_text(csv)
        meta = parse_meta(name)
        (dest / "meta.json").write_text(
            json.dumps(meta, indent=2, ensure_ascii=False) + "\n"
        )
        csvs[name] = csv
        metas[name] = json.dumps(meta, indent=2, ensure_ascii=False) + "\n"

    lines = [
        "package mixeff4s.comparison",
        "",
        "/** Vendored mixeff-rs frames. Regenerated by comparison/scripts/vendor_harness.py. */",
        "private[comparison] object EmbeddedFrames:",
        f'  val mixeffRsRevision: String = "{REVISION}"',
        "  val csv: Map[String, String] = Map(",
    ]
    for name in VENDORED:
        lines.append(f'    "{name}" -> {scala_triple(csvs[name])},')
    lines.append("  )")
    lines.append("  val meta: Map[String, String] = Map(")
    for name in VENDORED:
        lines.append(f'    "{name}" -> {scala_triple(metas[name])},')
    lines.append("  )")
    lines.append("")
    (ROOT / "src/test/scala/mixeff4s/comparison/EmbeddedFrames.scala").write_text(
        "\n".join(lines) + "\n"
    )

    rust = json.loads((RUST / "comparison" / "rust_results.json").read_text())
    lme4 = json.loads((RUST / "comparison" / "lme4_results.json").read_text())
    extracted = []
    for record in rust["results"]:
        if key_of(record) in EXTRACT_KEYS:
            extracted.append(slim_result(record, "mixeff-rs", True))
    for record in lme4["results"]:
        if key_of(record) in EXTRACT_KEYS:
            extracted.append(slim_result(record, "lme4", True))

    pins = [
        {
            "dataset": "sleepstudy",
            "formula": "reaction ~ 1 + days + (1 + days | subj)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "ML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 180,
            "objective": 1751.9393444636682,
            "objective_definition": "deviance",
            "objective_comparable": True,
            "beta": [251.40510484848454, 10.467285959596126],
            "coef_names": ["(Intercept)", "days"],
            "theta": [0.9292297167514472, 0.01816466496782548, 0.22264601131030412],
            "sigma": 25.591813564885108,
            "tolerances": {
                "objective": 1e-2,
                "beta": 1e-3,
                "theta": 1e-3,
                "sigma": 1e-3,
            },
        },
        {
            "dataset": "sleepstudy",
            "formula": "reaction ~ 1 + days + (1 + days | subj)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "REML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 180,
            "objective": 1743.6282719599442,
            "objective_definition": "restricted_deviance",
            "objective_comparable": True,
            "beta": [251.40510484848528, 10.467285959595493],
            "coef_names": ["(Intercept)", "days"],
            "theta": [0.9667417690560796, 0.015169059384716037, 0.2309099529619309],
            "sigma": 25.591795732317802,
            "tolerances": {
                "objective": 1e-2,
                "beta": 1e-3,
                "theta": 1e-3,
                "sigma": 1e-3,
            },
        },
        {
            "dataset": "sleepstudy",
            "formula": "reaction ~ 1 + days + (1 + days || subj)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "ML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 180,
            "objective": 1752.003255140962,
            "objective_definition": "deviance",
            "objective_comparable": True,
            "beta": [251.4051048484854, 10.467285959595674],
            "coef_names": ["(Intercept)", "days"],
            "theta": [0.9458043022417869, 0.22692740996014607],
            "sigma": 25.55613836753517,
            "is_singular": False,
            "tolerances": {
                "objective": 1e-2,
                "beta": 1e-3,
                "theta": 1e-3,
                "sigma": 1e-3,
            },
        },
        {
            "dataset": "penicillin",
            "formula": "diameter ~ 1 + (1 | plate) + (1 | sample)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "ML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 144,
            "objective": 332.1883486700085,
            "objective_definition": "deviance",
            "objective_comparable": True,
            "beta": [22.97222222222222],
            "coef_names": ["(Intercept)"],
            "theta": [1.5375939045981573, 3.219792193110907],
            "tolerances": {"objective": 1e-2, "beta": 1e-3, "theta": 1e-2},
        },
        {
            "dataset": "penicillin",
            "formula": "diameter ~ 1 + (1 | plate) + (1 | sample)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "REML",
            "engine": "mixeff-rs",
            "status": "ok",
            "n_obs": 144,
            "objective": 330.86058899126897,
            "objective_definition": "restricted_deviance",
            "objective_comparable": True,
            "beta": [22.97222222222248],
            "coef_names": ["(Intercept)"],
            "theta": [1.5396773350745998, 3.51241122181154],
            "sigma": 0.549923173720829,
            "tolerances": {
                "objective": 1e-2,
                "beta": 1e-3,
                "theta": 1e-2,
                "sigma": 1e-3,
            },
        },
        {
            "dataset": "pastes",
            "formula": "strength ~ 1 + (1 | batch / cask)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "ML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 60,
            "objective": 247.9944658624955,
            "objective_definition": "deviance",
            "objective_comparable": True,
            "beta": [60.0533333333333],
            "coef_names": ["(Intercept)"],
            "theta": [3.5269029347766856, 1.3299137410046242],
            "tolerances": {"objective": 1e-2, "beta": 1e-3, "theta": 0.09},
        },
        {
            "dataset": "contraception",
            "formula": "use_num ~ 1 + age + age2 + urban + livch + (1 | urban_dist)",
            "family": "Bernoulli",
            "link": "Logit",
            "estimator": "fast_pirls",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 1934,
            "objective": 2361.657202855648,
            "objective_definition": "profiled_fast_pirls_deviance",
            "objective_comparable": True,
            "theta": [0.5720746212924732],
            "sigma": 1.0,
            "tolerances": {"objective": 1.0, "theta": 0.01},
        },
        {
            "dataset": "oxide",
            "formula": "Thickness ~ 1 + (1 | Lot) + (1 | Lot:Wafer)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "REML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 72,
            "objective": 454.02206930988217,
            "objective_definition": "restricted_deviance",
            "objective_comparable": True,
            "beta": [2000.1527777778349],
            "coef_names": ["(Intercept)"],
            "theta": [1.6892198720377958, 3.2149146331969516],
            "sigma": 3.5453258700425447,
            "is_singular": False,
            "tolerances": {
                "objective": 1e-2,
                "beta": 1e-3,
                "theta": 1e-3,
                "sigma": 1e-3,
            },
        },
        {
            "dataset": "oxide",
            "formula": "Thickness ~ 1 + (1 | Lot/Wafer)",
            "family": "Gaussian",
            "link": "Identity",
            "estimator": "REML",
            "engine": "mixedmodels.jl",
            "status": "ok",
            "n_obs": 72,
            "objective": 454.02206930988217,
            "objective_definition": "restricted_deviance",
            "objective_comparable": True,
            "beta": [2000.1527777778349],
            "coef_names": ["(Intercept)"],
            "theta": [1.6892198720377958, 3.2149146331969516],
            "sigma": 3.5453258700425447,
            "is_singular": False,
            "tolerances": {
                "objective": 1e-2,
                "beta": 1e-3,
                "theta": 1e-3,
                "sigma": 1e-3,
            },
        },
    ]

    catalog = {
        "schema": {"name": "mixeff4s.frozen_reference", "version": 1},
        "source": {
            "mixeff_rs_revision": REVISION,
            "lme4_tool": lme4.get("tool"),
            "note": (
                "Seed pins are the numbers already checked in mixeff4s fit suites. "
                "Additional rust and lme4 rows were copied from mixeff-rs comparison "
                "artifacts at the recorded revision. No live R or Julia was run. "
                "A frozen lme4 row is not a mixeff4s claim until the scorecard says so. "
                "There are no p-values."
            ),
        },
        "results": pins + extracted,
    }
    dest = ROOT / "comparison" / "frozen" / "references.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    frozen_text = json.dumps(catalog, indent=2, ensure_ascii=False) + "\n"
    dest.write_text(frozen_text)
    (ROOT / "src/test/scala/mixeff4s/comparison/EmbeddedFrozen.scala").write_text(
        "\n".join(
            [
                "package mixeff4s.comparison",
                "",
                "/** Frozen references. Regenerated by comparison/scripts/vendor_harness.py. */",
                "private[comparison] object EmbeddedFrozen:",
                f"  val json: String = {scala_triple(frozen_text)}",
                "",
            ]
        )
        + "\n"
    )
    print(
        f"vendored {len(VENDORED)} datasets; froze {len(catalog['results'])} reference rows"
    )


if __name__ == "__main__":
    main()
