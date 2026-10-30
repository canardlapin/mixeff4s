#!/usr/bin/env python3
"""Seal or verify local inference evidence and exact source; never rerun numerical models."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'comparison/inference-tranche'


def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def check(mapping):
    for name,digest in mapping.items():assert sha(ROOT/name)==digest,name


def archive_check(path,mapping):
    with tarfile.open(path) as archive:
        for name,digest in mapping.items():
            assert hashlib.sha256(archive.extractfile(name).read()).hexdigest()==digest,name


def finalize():
    baseline=read(BASE/'starting-inputs.json')
    allow=set(read(BASE/'format-preservation.json')['existing_changed_allowlist'])
    check({p:h for p,h in baseline['inputs'].items() if p not in allow})
    assert all((ROOT/p).exists() for p in baseline['untracked'])
    check(read(BASE/'review-inputs.json'))
    assert sha(BASE/'review-inputs.json')=='93315a1768a3e4f5048bec135c6c991a68e3a004bce7cdc9b6bd15f48b2657a4'
    meta=read(BASE/'checks/full-final-03.log.meta.json');assert meta['exit_code']==0
    text=(BASE/'checks/full-final-03.log').read_text()
    assert 'Passed: Total 440, Failed 0, Errors 0, Passed 440' in text
    assert 'Passed: Total 436, Failed 0, Errors 0, Passed 436, Ignored 4' in text
    assert '[error]' not in text
    assert all('multiple main classes detected' in line for line in text.splitlines() if '[warn]' in line)
    assert 'attempts=39; unqualified=Vector(); calibration=Right' in text
    provenance=read(BASE/'reference-provenance.json')
    assert sha(BASE/'references.json')==provenance['artifacts']['references_json_sha256']
    assert sha(ROOT/provenance['generator']['path'])==provenance['generator']['sha256']
    for row in provenance['frames'].values():assert sha(ROOT/row['path'])==row['sha256']
    assert read(BASE/'R-freeze-normalization.log.meta.json')['exit_code']==0
    study=BASE/'study-final'
    check(read(study/'inputs.json'))
    archive_check(study/'source.tar.gz',read(study/'inputs.json'))
    assert read(study/'run.log.meta.json')['exit_code']==0
    assert read(BASE/'checks/study-final.log.meta.json')['exit_code']==0
    result=read(study/'results.json')
    assert len(result['summary'])==2 and all(r['requested']==16 for r in result['summary'])
    assert result['maximum_response_difference']<=1e-12
    assert len(result['records']['INFERENCE_STUDY'])==32
    assert len(result['records']['PROFILE_STUDY'])==32
    paths=[]
    for folder in ['src','.jvm/src','.js/src']:
        paths += [p for p in (ROOT/folder).rglob('*') if p.is_file()]
    paths += [ROOT/p for p in ['build.sbt','project/plugins.sbt','project/build.properties','.scalafmt.conf','AGENTS.md',
        'README.md','comparison/INFERENCE_TRANCHE.md','comparison/frozen/references.json','comparison/parity_scorecard.toml',
        'comparison/scripts/freeze_inference.R','comparison/scripts/embed_inference.py','comparison/scripts/run_inference_qualification.py']]
    paths += [ROOT/row['path'] for row in provenance['frames'].values()]
    paths += [BASE/'references.json']
    inputs={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))}
    with tarfile.open(BASE/'final-source.tar.gz','w:gz') as archive:
        for name in inputs:archive.add(ROOT/name,arcname=name)
    evidence={str(p.relative_to(ROOT)):sha(p) for p in sorted(BASE.rglob('*')) if p.is_file() and p.name!='receipt.json' and '__pycache__' not in p.parts}
    receipt=dict(schema_version=1,starting_checkpoint=baseline['checkpoint'],gale_revision='f869613cec0a89e57b6c995b0a02cf471ac7127c',
        input_sha256=inputs,evidence_sha256=evidence,
        checks=dict(jvm_passed=440,js_passed=436,js_ignored=4,glmm_workflow_null_pairs=39,
            gaussian_null_witness_pairs=7,outer_study_draws_per_runtime=16,inner_bootstrap_pairs_requested=39,
            study_summary=result['summary']),
        boundaries=['Local checks only; no hosted CI, push, release or publication',
            'Local numerical evidence does not establish finite MLE, global optimum, nominal coverage or exact-size calibration',
            'Retained coefficient/residual/scalar independent RE profiles only; GLMM and correlated component profiles refused',
            'Every null-bootstrap failure retained; any unqualified pair prevents calibrated extraction',
            'Same16 study outer draws across runtimes; current-thread JVM allocation, JS unavailable'])
    (BASE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--finalize',action='store_true');args=parser.parse_args()
    result=finalize() if args.finalize else read(BASE/'receipt.json')
    check(result['input_sha256']);check(result['evidence_sha256']);archive_check(BASE/'final-source.tar.gz',result['input_sha256'])
    print(f"Verified {len(result['input_sha256'])} inputs and {len(result['evidence_sha256'])} evidence files; exact source archive matches.")


if __name__=='__main__':main()
