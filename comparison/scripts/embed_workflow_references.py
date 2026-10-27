#!/usr/bin/env python3
"""Add classic workflow lme4 pins without changing existing references; regenerate portable embeddings."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def add_workflow_pins(catalog):
    path=ROOT/'comparison/workflow-tranche/lme4-workflow.json'
    if not path.exists(): return catalog
    workflow=json.loads(path.read_text())
    def key(row): return tuple(row[k] for k in ['dataset','formula','estimator','engine'])
    present={key(row) for row in catalog['results']}
    for fit in workflow['results']:
        def vector(value): return value if isinstance(value,list) else [value]
        row=dict(dataset=fit['dataset'],formula=fit['formula'],family='Gaussian',link='Identity',estimator=fit['criterion'],engine='lme4',status='ok',n_obs=fit['n'],objective=fit['objective'],objective_definition='deviance' if fit['criterion']=='ML' else 'restricted_deviance',objective_comparable=True,beta=vector(fit['beta']),coef_names=vector(fit['coef_names']),theta=vector(fit['theta']),sigma=fit['sigma'],is_singular=fit['singular'],tolerances=dict(objective=1e-3,beta=1e-3,theta=5e-3 if fit['dataset']=='pastes' else 1e-3,sigma=1e-3))
        if key(row) not in present:
            catalog['results'].append(row); present.add(key(row))
    catalog['source']['note']='Earlier pins retain their original provenance. Classic LMM workflow pins were generated with live R/lme4; scripts, package versions, frame hashes and full prediction/mode references are under comparison/workflow-tranche. Existing reference values are preserved. There are no p-values.'
    return catalog

def main():
    path=ROOT/'comparison/frozen/references.json'
    catalog=add_workflow_pins(json.loads(path.read_text()))
    text=json.dumps(catalog,indent=2,ensure_ascii=False)+'\n';path.write_text(text)
    (ROOT/'src/test/scala/mixeff4s/comparison/EmbeddedFrozen.scala').write_text('package mixeff4s.comparison\n\n/** Generated frozen reference catalog. */\nprivate[comparison] object EmbeddedFrozen:\n  val json: String = """'+text+'"""\n')
    print(f"Embedded {len(catalog['results'])} references, preserving existing numerical pins")
if __name__=='__main__': main()
