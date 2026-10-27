#!/usr/bin/env python3
"""Summarize raw opt-in benchmark logs without dropping samples or failed fits.

Run with MIXEFF4S_BENCHMARK=1 sbt -batch
  'mixeff4sJVM/testOnly mixeff4s.benchmark.EngineBenchmarkSuite'
  'mixeff4sJS/testOnly mixeff4s.benchmark.EngineBenchmarkSuite'
Capture stdout+stderr and actual exit status with run_logged.py, then invoke:
  python3 comparison/scripts/engine_benchmark.py LOG OUTPUT_JSON
"""
import json
import pathlib
import platform
import statistics
import subprocess
import sys

source, target = map(pathlib.Path, sys.argv[1:])
samples, fits = [], []
for line in source.read_text().splitlines():
    if line.startswith('BENCH '):
        samples.append(json.loads(line[6:]))
    if line.startswith('BENCH_FIT '):
        fits.append(json.loads(line[10:]))
if not samples or not fits:
    raise SystemExit('No benchmark samples or fit witnesses')
rows = []
for key in sorted({(s['platform'], s['workload'], s['phase']) for s in samples}):
    group = [s for s in samples if (s['platform'], s['workload'], s['phase']) == key]
    rows.append(dict(platform=key[0], workload=key[1], phase=key[2], samples=len(group),
        median_ms=statistics.median(s['ns']/1e6 for s in group),
        min_ms=min(s['ns']/1e6 for s in group), max_ms=max(s['ns']/1e6 for s in group),
        median_allocated_bytes=statistics.median(s['allocated_bytes'] for s in group)
            if all(s['allocated_bytes'] is not None for s in group) else None))
cpu = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'], capture_output=True, text=True) if sys.platform == 'darwin' else None
result = dict(schema_version=1, source_log=str(source),
    host=dict(platform=platform.platform(), machine=platform.machine(),
              cpu=cpu.stdout.strip() if cpu and cpu.returncode == 0 else 'unavailable'),
    samples=samples, fits=fits, summary=rows,
    limitations=[
        'One process per runtime per log; first fit is first use of the workload, not a new process for each workload',
        'JVM allocated bytes cover the current thread; Scala.js allocation is unavailable',
        'Heap deltas include GC and do not measure peak memory',
        'JVM runs within sbt; use the same environment and independent repeat runs for strong timing claims',
        'response_batch_rebuild_5 recompiles each response; reusable refit is a distinct phase when available'
    ])
target.write_text(json.dumps(result, indent=2)+'\n')
print(f'{len(samples)} samples, {len(fits)} fit witnesses, {len(rows)} summary rows -> {target}')
for row in rows:
    if row['phase'] == 'warm_fit':
        print(f"{row['platform']} {row['workload']}: {row['median_ms']:.3f} ms, {row['median_allocated_bytes']} allocated bytes")
