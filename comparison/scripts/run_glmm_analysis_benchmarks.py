#!/usr/bin/env python3
"""Matched cold-start GLMM rebuild/refit runs; preserve response identity and every outcome."""
import argparse
import hashlib
import json
import os
import pathlib
import platform
import statistics
import subprocess
import tarfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=pathlib.Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--responses', type=int, default=12)
    args = parser.parse_args()
    if args.repeats < 1 or not 1 <= args.responses <= 64:
        parser.error('positive repeats and1 <= responses <=64 required')
    root = pathlib.Path(__file__).resolve().parents[2]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs = list((root / 'src').rglob('*.scala'))
    for directory in ['.jvm/src', '.js/src']:
        inputs += list((root / directory).rglob('*.scala'))
    inputs += [root / p for p in ['build.sbt', 'project/plugins.sbt', 'project/build.properties']]
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    hashes = {str(p.relative_to(root)): sha(p) for p in sorted(inputs)}
    (output / 'input-sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
    with tarfile.open(output / 'source.tar.gz', 'w:gz') as archive:
        for relative in hashes:
            archive.add(root / relative, arcname=relative)
    runs = []
    bank = {}
    maximum_response_difference = 0.0
    expected_per_platform = 12 * args.responses + 2 * min(8, args.responses)
    expected_responses = 3 * args.responses + min(8, args.responses)
    for repeat in range(args.repeats):
        log = output / f'run-{repeat+1}.log'
        command = ['python3', str(pathlib.Path.home() / '.agents/skills/lean-logs/scripts/run_logged.py'),
                   '--log', str(log), '--timeout', '900', '--head-bytes', '0', '--tail-bytes', '1600', '--',
                   'sbt', '-batch', 'mixeff4sJVM/testOnly mixeff4s.benchmark.GlmmWorkflowBenchmarkSuite',
                   'mixeff4sJS/testOnly mixeff4s.benchmark.GlmmWorkflowBenchmarkSuite']
        result = subprocess.run(command, cwd=root, env=dict(os.environ, MIXEFF4S_BENCHMARK='1',
                                MIXEFF4S_GLMM_BENCH_RESPONSES=str(args.responses)))
        if result.returncode:
            raise SystemExit(result.returncode)
        text = log.read_text()
        rows = [json.loads(line[len('GLMM_BENCH '):]) for line in text.splitlines() if line.startswith('GLMM_BENCH ')]
        responses = [json.loads(line[len('GLMM_BENCH_RESPONSE '):]) for line in text.splitlines() if line.startswith('GLMM_BENCH_RESPONSE ')]
        assert len(rows) == 2 * expected_per_platform, (log, len(rows))
        assert len(responses) == 2 * expected_responses, (log, len(responses))
        for row in responses:
            key = (row['workload'], row['index'])
            original = bank.setdefault(key, row)
            assert row['seed'] == original['seed']
            error = max(abs(a-b) for a, b in zip(row['y'], original['y']))
            maximum_response_difference = max(maximum_response_difference, error)
            assert error <= 1e-12, (key, error)
        for row in rows:
            row['response_sha256'] = hashlib.sha256(json.dumps(bank[(row['workload'], row['index'])]['y']).encode()).hexdigest()
        (output / f'run-{repeat+1}.json').write_text(json.dumps(rows, indent=2) + '\n')
        runs.append(rows)
    summary = []
    keys = sorted({(r['platform'], r['workload'], r['nAgq'], r['phase']) for r in runs[0]})
    for key in keys:
        groups = [[r for r in rows if (r['platform'], r['workload'], r['nAgq'], r['phase']) == key] for rows in runs]
        medians = [statistics.median(r['ns'] / 1e6 for r in group) for group in groups]
        alloc = [statistics.median(r['allocated_bytes'] for r in g) if all(r['allocated_bytes'] is not None for r in g) else None for g in groups]
        summary.append(dict(platform=key[0], workload=key[1], nAgq=key[2], phase=key[3],
                            max_evaluations=groups[0][0]['max_evaluations'],
                            process_medians_ms=medians, median_ms=statistics.median(medians),
                            allocated_bytes=statistics.median(alloc) if all(x is not None for x in alloc) else None,
                            fresh_assessment_ms=statistics.median(statistics.median(r['fresh_assessment_ns']/1e6 for r in g) for g in groups),
                            fit_failures=sum(r['error'] is not None for g in groups for r in g),
                            unqualified=sum(r['error'] is None and not r['qualified'] for g in groups for r in g),
                            attempts=sum(map(len, groups))))
    for relative, digest in hashes.items():
        assert sha(root / relative) == digest, f'source changed during timing: {relative}'
    distinct = {name: len({json.dumps(row['y']) for (workload, _), row in bank.items() if workload == name}) for name in sorted({k[0] for k in bank})}
    result = dict(repeats=args.repeats, responses_per_small_workload=args.responses,
                  samples=sum(map(len, runs)), matched_pairs=sum(map(len, runs))//2,
                  distinct_responses=distinct, maximum_response_difference=maximum_response_difference,
                  host=dict(platform=platform.platform(), machine=platform.machine()), summary=summary,
                  limits=['Fits include their built-in outer diagnostic; fresh inference assessment is timed separately',
                          'Crossed measured attempts have an explicit2000-evaluation budget; generation fit uses10000; small cases use10000',
                          'Each response draws new random effects and response noise; process repetitions reuse that response bank',
                          'Cold starts preserve estimator/controls; full rebuild includes frame/design compilation',
                          'Current-thread JVM allocation only; JS allocation unavailable; no peak-memory claim',
                          'Fit failures and stricter assessment failures are distinct; all observations retained'])
    (output / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f"{result['samples']} samples / {result['matched_pairs']} matched pairs; distinct responses {distinct}")
    for row in summary:
        if row['phase'] == 'refit':
            old = next(r for r in summary if (r['platform'],r['workload'],r['nAgq'],r['phase']) == (row['platform'],row['workload'],row['nAgq'],'rebuild'))
            print(row['platform'], row['workload'], row['nAgq'], round(old['median_ms'],3), '->', round(row['median_ms'],3),
                  'ms; fit failures',row['fit_failures'],'unqualified',row['unqualified'],'of',row['attempts'])


if __name__ == '__main__':
    main()
