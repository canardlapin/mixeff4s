#!/usr/bin/env python3
"""Add independently frozen GLMM claims without replacing older reference records."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FORMULAS = {
    "cbpp": "proportion ~ 1 + period + (1 | herd)",
    "culcitalogreg": "predation ~ 1 + ttt + (1 | block)",
    "synthetic_poisson_offset": "count ~ 1 + x + (1 | group)",
    "synthetic_poisson_crossed": "count ~ 1 + x + (1 + x | site) + (1 | batch)",
}

def add_glmm_pins(catalog):
    rows = json.loads((ROOT / "comparison/glmm-tranche/references.json").read_text())["results"]
    for row in rows:
        joint = row["n_agq"] == 1
        pin = {
            "dataset": row["dataset"], "formula": FORMULAS[row["dataset"]],
            "family": "Binomial" if row["dataset"] == "cbpp" else "Bernoulli" if row["dataset"] == "culcitalogreg" else "Poisson",
            "link": row["link"], "estimator": "joint_laplace" if joint else "fast_pirls", "engine": "lme4",
            "status": "ok", "n_obs": row["n"],
            "objective": row["objective_minus2loglik_with_constants"],
            "objective_definition": "minus_two_loglikelihood_with_response_constants",
            "objective_comparable": joint,
            "beta": row["beta"], "coef_names": row["coef_names"], "theta": row["theta"],
            "is_singular": row["singular"],
            "tolerances": {"objective": 1e-4, "beta": 1e-3, "theta": 1e-3},
            "provenance": "comparison/glmm-tranche/provenance.json",
            "response_and_offset_contract": "cbpp: proportion=incidence/size, trials=size; Poisson: offset=log(exposure); case weights=1",
            "reference_estimator": "lme4 nAGQ=1" if joint else "lme4 nAGQ=0, documented fast divergence only",
        }
        if not joint and row["dataset"] == "synthetic_poisson_offset":
            pin["tolerances"]["theta"] = 0.02
            pin["tolerances"]["objective"] = 0.01
        key = lambda p: tuple(p[k] for k in ("dataset", "formula", "family", "link", "estimator", "engine"))
        existing = [p for p in catalog["results"] if key(p) == key(pin)]
        if existing:
            assert existing == [pin], f"conflicting GLMM pin {key(pin)}"
        else:
            catalog["results"].append(pin)
    return catalog

def main():
    path = ROOT / "comparison/frozen/references.json"
    from embed_glmm_analysis import add_glmm_analysis_pins
    catalog = add_glmm_analysis_pins(add_glmm_pins(json.loads(path.read_text())))
    text = json.dumps(catalog, indent=2, ensure_ascii=False) + "\n"
    path.write_text(text)
    chunks = [text[i:i + 12000] for i in range(0, len(text), 12000)]
    (ROOT / "src/test/scala/mixeff4s/comparison/EmbeddedFrozen.scala").write_text(
        'package mixeff4s.comparison\n\n/** Frozen references; see comparison/scripts. */\nprivate[comparison] object EmbeddedFrozen:\n  val json: String = Vector(\n' +
        ',\n'.join('    """' + chunk + '"""' for chunk in chunks) + '\n  ).mkString\n')
    scorecard = ROOT / "comparison/parity_scorecard.toml"
    content = scorecard.read_text()
    marker = "# Qualified canonical GLMM workflow claims."
    if marker not in content:
        content += "\n" + marker + "\n"
        for name, formula in FORMULAS.items():
            for joint in (True, False):
                family = "Binomial" if name == "cbpp" else "Bernoulli" if name == "culcitalogreg" else "Poisson"
                reason = ("Joint canonical Laplace; full response constants, modes, means and first-order evidence compared on JVM/JS against explicit tolPwrss=1e-10 lme4. "
                          "cbpp uses explicit trials=size; Poisson uses GlmmOptions.offset=log(exposure).") if joint else (
                          "Labelled fast-PIRLS, documented divergence from lme4 nAGQ=0. Not joint estimator parity. "
                          "Small Poisson offset fixture theta differs by about 0.014; objective and mean differences remain visible.")
                row = {"dataset": name, "formula": formula, "family": family, "link": "Log" if family == "Poisson" else "Logit",
                       "estimator": "joint_laplace" if joint else "fast_pirls", "class": "release_blocking_parity" if joint else "documented_divergence",
                       "reference": "lme4_joint_laplace" if joint else "lme4_numeric_without_objective_constants",
                       "objective_abs_tol": 0.01 if not joint and name == "synthetic_poisson_offset" else 1e-4,
                       "beta_abs_tol": 1e-3, "theta_abs_tol": 0.02 if not joint and name == "synthetic_poisson_offset" else 1e-3,
                       "reason": reason}
                content += "\n[[row]]\n" + "".join(f"{k} = {json.dumps(v)}\n" for k, v in row.items())
        scorecard.write_text(content)
    marker = "# Qualified scalar adaptive quadrature claims."
    if marker not in content:
        content += "\n" + marker + "\n"
        for name in ("cbpp", "culcitalogreg", "synthetic_poisson_offset"):
            for nodes in (5, 9, 15):
                family = "Binomial" if name == "cbpp" else "Bernoulli" if name == "culcitalogreg" else "Poisson"
                row = dict(dataset=name, formula=FORMULAS[name], family=family, link="Log" if family == "Poisson" else "Logit",
                           estimator=f"joint_agq_{nodes}", **{"class": "release_blocking_parity"}, reference="lme4_joint_agq",
                           objective_abs_tol=1e-4, beta_abs_tol=1e-3, theta_abs_tol=1e-3,
                           reason=f"Scalar AGQ{nodes}: full-density objective, coefficients and predictions compared on JVM/JS. lme4 residual-density constants corrected explicitly and independently integrated. No global-optimality or coverage claim.")
                content += "\n[[row]]\n" + "".join(f"{k} = {json.dumps(v)}\n" for k, v in row.items())
        scorecard.write_text(content)
    scala = ROOT / "src/test/scala/mixeff4s/comparison/Scorecard.scala"
    source = scala.read_text()
    start = source.index('    """|')
    end = source.index('""".stripMargin', start) + len('""".stripMargin')
    source = source[:start] + '    """|' + '\n      |'.join(content.splitlines()) + '\n      |""".stripMargin' + source[end:]
    scala.write_text(source)
    print(f"Preserved prior records; catalog now {len(catalog['results'])} records; added GLMM scorecard claims")

if __name__ == "__main__": main()
