#!/usr/bin/env python3
"""Finalize or verify the exact local GLMM workflow evidence, without rerunning models."""
import argparse
import collections
import hashlib
import json
import pathlib
import re
import subprocess
import tarfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
BASE = ROOT / 'comparison/glmm-analysis-tranche'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def verify_mapping(mapping):
    for name, digest in mapping.items():
        assert sha(ROOT / name) == digest, name


def verify_archive(path, hashes):
    with tarfile.open(path) as archive:
        for name, digest in hashes.items():
            assert hashlib.sha256(archive.extractfile(name).read()).hexdigest() == digest, name


def finalize():
    prior = read(BASE / 'starting-inputs.json')
    allowed = set(read(BASE / 'format-preservation.json')['changed_paths_allowlist'])
    unaffected = {p: h for p, h in prior['input_sha256'].items() if p not in allowed}
    verify_mapping(unaffected)
    assert all((ROOT / p).exists() for p in prior['remaining_untracked'])
    old = json.loads(subprocess.check_output(['git', 'show', prior['checkpoint_head'] + ':comparison/frozen/references.json'], cwd=ROOT))
    current = read(ROOT / 'comparison/frozen/references.json')
    canonical = lambda row: json.dumps(row, sort_keys=True, separators=(',', ':'))
    a = collections.Counter(map(canonical, old['results']))
    b = collections.Counter(map(canonical, current['results']))
    assert len(old['results']) == 62 and len(current['results']) == 71 and not (a - b)
    verify_mapping(read(BASE / 'review-inputs.json'))
    benchmark = BASE / 'benchmarks-final'
    verify_mapping(read(benchmark / 'input-sha256.json'))
    verify_archive(benchmark / 'source.tar.gz', read(benchmark / 'input-sha256.json'))
    verify_archive(BASE / 'bootstrap-producer-source.tar.gz', read(BASE / 'bootstrap-producer-inputs.json'))
    runs = read(benchmark / 'results.json')
    assert runs['samples'] == 960 and runs['matched_pairs'] == 480 and runs['repeats'] == 3
    assert runs['distinct_responses'] == dict(binomial=12, crossed_poisson=8, poisson_offset=12, rare_bernoulli=12)
    assert runs['maximum_response_difference'] <= 1e-12
    assert all(read(benchmark / f'run-{i}.log.meta.json')['exit_code'] == 0 for i in [1, 2, 3])
    assert read(BASE / 'checks/benchmark-final.log.meta.json')['exit_code'] == 0
    assert read(BASE / 'checks/full-final.log.meta.json')['exit_code'] == 0
    text = (BASE / 'checks/full-final.log').read_text()
    assert '[info] Passed: Total 411, Failed 0, Errors 0, Passed 411' in text
    assert '[info] Passed: Total 407, Failed 0, Errors 0, Passed 407, Ignored 3' in text
    assert 'Qualified bootstrap refits: 39/39' in text
    warnings = [line for line in text.splitlines() if '[warn]' in line]
    assert warnings == ["[warn] multiple main classes detected: run 'show discoveredMainClasses' to see the list"]
    assert '[error]' not in text
    numpy = read(BASE / 'samplers/numpy-random-source.json')
    assert sha(BASE / 'samplers/numpy-distributions.c') == numpy['distribution_source']['sha256']
    assert sha(BASE / 'samplers/NUMPY_RANDOM_LICENSE.md') == numpy['license_source']['sha256']
    assert sha(ROOT / 'src/main/resources/META-INF/licenses/numpy-random.md') == numpy['license_source']['sha256']
    provenance = read(BASE / 'reference-provenance.json')
    for name, digest in provenance['artifact_sha256'].items():
        assert sha(BASE / name) == digest
    for name, digest in provenance['generator_sha256'].items():
        assert sha(ROOT / 'comparison/scripts' / name) == digest
    for name, digest in provenance['frame_sha256'].items():
        assert sha(ROOT / 'comparison/glmm-tranche/frames' / (name + '.csv')) == digest
    bootstrap = read(BASE / 'bootstrap-r-provenance.json')['inputs_sha256']
    for path, digest in [(BASE / 'frames/binomial_workflow.csv', bootstrap['binomial_workflow_csv']),
                         (BASE / 'bootstrap-responses.csv', bootstrap['bootstrap_responses_csv']),
                         (BASE / 'bootstrap-response-provenance.json', bootstrap['bootstrap_response_provenance_json']),
                         (ROOT / 'comparison/scripts/freeze_glmm_bootstrap.R', bootstrap['generator'])]:
        assert sha(path) == digest
    source = []
    for directory in ['src', '.jvm/src', '.js/src']:
        source.extend(p for p in (ROOT / directory).rglob('*') if p.is_file())
    source += [ROOT / p for p in ['build.sbt', 'project/plugins.sbt', 'project/build.properties', '.scalafmt.conf',
                                  'AGENTS.md', 'README.md', 'comparison/GLMM_ANALYSIS_TRANCHE.md',
                                  'comparison/parity_scorecard.toml', 'comparison/frozen/references.json']]
    source += [p for p in (ROOT / 'comparison/scripts').glob('*') if p.is_file() and p.suffix in ['.py', '.R']]
    source = sorted(set(source))
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in source}
    with tarfile.open(BASE / 'final-source.tar.gz', 'w:gz') as archive:
        for name in inputs:
            archive.add(ROOT / name, arcname=name)
    evidence = {str(p.relative_to(ROOT)): sha(p) for p in sorted(BASE.rglob('*'))
                if p.is_file() and p.name != 'receipt.json'}
    result = dict(schema_version=1, starting_checkpoint=prior['checkpoint_head'],
                  gale_revision='f869613cec0a89e57b6c995b0a02cf471ac7127c',
                  input_sha256=inputs, evidence_sha256=evidence,
                  checks=dict(jvm_passed=411, js_passed=407, js_ignored=3, example_qualified_replicates=39,
                              benchmark_samples=960, benchmark_matched_pairs=480,
                              frozen_records_prior_preserved=62, frozen_records_current=71,
                              unrelated_tracked_unchanged=len(unaffected), prior_untracked_names_preserved=22),
                  boundaries=['Local checks and packages only; no hosted CI or publication',
                              'First-order qualification does not certify finite joint MLE, global optimum, Hessian or interval coverage',
                              'Benchmarks retain fit failures and strict-assessment refusals; no universal speedup claim',
                              'Bounded runtime replay, declared sampler limits and scalar AGQ support only'])
    (BASE / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--finalize', action='store_true')
    args = parser.parse_args()
    receipt = finalize() if args.finalize else read(BASE / 'receipt.json')
    verify_mapping(receipt['input_sha256'])
    verify_mapping(receipt['evidence_sha256'])
    verify_archive(BASE / 'final-source.tar.gz', receipt['input_sha256'])
    print(f"Verified {len(receipt['input_sha256'])} inputs and {len(receipt['evidence_sha256'])} evidence files; source archive exact.")


if __name__ == '__main__':
    main()
