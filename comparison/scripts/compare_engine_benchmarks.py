#!/usr/bin/env python3
"""Compare process-level medians; keep differing APIs/strategies explicitly named."""
import json
import re
from pathlib import Path
import statistics
import sys

before, after, output = map(Path, sys.argv[1:])
old = [json.loads(p.read_text()) for p in sorted(before.glob('run-*.json')) if re.fullmatch(r'run-\d+\.json', p.name)]
new = [json.loads(p.read_text()) for p in sorted(after.glob('run-*.json')) if re.fullmatch(r'run-\d+\.json', p.name)]
if not old or not new:
    raise SystemExit('Both sides need successful run JSON files')

def summarize(runs, platform, workload, phase):
    rows = [r for run in runs for r in run['summary'] if (r['platform'],r['workload'],r['phase']) == (platform,workload,phase)]
    if len(rows) != len(runs):
        raise ValueError(f'missing phase {platform}/{workload}/{phase}')
    medians = [r['median_ms'] for r in rows]
    allocation = [r['median_allocated_bytes'] for r in rows]
    return dict(process_medians_ms=medians, median_ms=statistics.median(medians),
                min_ms=min(medians),max_ms=max(medians),
                allocated_bytes=statistics.median(allocation) if all(a is not None for a in allocation) else None)

comparisons = []
phases = [(p,p) for p in ['first_fit','compile','workspace','objective','warm_fit','inspection','response_batch_rebuild_5']]
phases += [('warm_fit','compiled_fit'),('warm_fit','compiled_fit_gradient'),
           ('objective','objective_gradient'),
           ('response_batch_rebuild_5','response_batch_refit_independent_5'),
           ('response_batch_rebuild_5','response_batch_refit_gradient_independent_5'),
           ('response_batch_rebuild_5','response_batch_refit_gradient_warm_5')]
workloads = sorted({(r['platform'],r['workload']) for r in old[0]['summary']})
for platform,workload in workloads:
    for a,b in phases:
        first = summarize(old,platform,workload,a)
        second = summarize(new,platform,workload,b)
        comparisons.append(dict(platform=platform,workload=workload,baseline_phase=a,after_phase=b,
            before=first,after=second,time_ratio_before_over_after=first['median_ms']/second['median_ms'] if second['median_ms'] else None))
# Verify objective witnesses, independently of timing comparisons.
checks = []
for platform,workload in workloads:
    ref = [f['objective'] for run in old for f in run['fits'] if (f['platform'],f['workload'])==(platform,workload)]
    for run_index,run in enumerate(new):
        for f in run['fits']:
            if (f['platform'],f['workload'])==(platform,workload):
                error = abs(f['objective']-statistics.median(ref))
                checks.append(dict(platform=platform,workload=workload,strategy=f.get('strategy','derivative_free'),run=run_index+1,objective_difference=error,within_0_002=error<=0.002))
if not all(c['within_0_002'] for c in checks):
    raise SystemExit('Objective witness mismatch: do not qualify timing ratios')
output.write_text(json.dumps(dict(baseline_runs=len(old),after_runs=len(new),comparisons=comparisons,objective_checks=checks,
    limits=['Ratios compare process medians on this host, not universal speed guarantees','JVM current-thread allocation only; JS allocations unavailable','API and optimizer changes are separate comparison rows','WarmSequential refits are explicitly order-dependent']),indent=2)+'\n')
print(f'{len(comparisons)} comparisons and {len(checks)} objective checks -> {output}')
for phase in ['warm_fit','compiled_fit_gradient','response_batch_refit_independent_5','response_batch_refit_gradient_warm_5']:
    print('\n'+phase)
    for c in comparisons:
        if c['after_phase']==phase:
            print(f"{c['platform']:3} {c['workload']:18} {c['before']['median_ms']:9.3f} -> {c['after']['median_ms']:9.3f} ms, ratio {c['time_ratio_before_over_after']:.2f}")
