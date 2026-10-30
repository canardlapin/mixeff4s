#!/usr/bin/env python3
"""Run bounded scientific observations once per runtime; preserve failures and source identity."""
import argparse
import collections
import hashlib
import json
import os
from pathlib import Path
import statistics
import subprocess
import tarfile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[2]
    output=args.output.resolve()
    output.mkdir(parents=True,exist_ok=True)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    sources=sorted([p for folder in ['src','.jvm/src','.js/src'] for p in (root/folder).rglob('*.scala')]+[root/'build.sbt'])
    inputs={str(p.relative_to(root)):sha(p) for p in sources}
    (output/'inputs.json').write_text(json.dumps(inputs,indent=2)+'\n')
    with tarfile.open(output/'source.tar.gz','w:gz') as archive:
        for p in sources:archive.add(p,arcname=str(p.relative_to(root)))
    command=['python3',str(Path.home()/'.agents/skills/lean-logs/scripts/run_logged.py'),
        '--log',str(output/'run.log'),'--timeout','1200','--head-bytes','0','--tail-bytes','1500','--',
        'sbt','-batch','mixeff4sJVM/testOnly mixeff4s.stats.InferenceQualificationSuite',
        'mixeff4sJS/testOnly mixeff4s.stats.InferenceQualificationSuite']
    result=subprocess.run(command,cwd=root,env=dict(os.environ,MIXEFF4S_BENCHMARK='1'))
    if result.returncode:raise SystemExit(result.returncode)
    text=(output/'run.log').read_text()
    records={}
    for prefix in ['INFERENCE_STUDY','INFERENCE_STUDY_RESPONSE','PROFILE_STUDY']:
        records[prefix]=[json.loads(line[len(prefix)+1:]) for line in text.splitlines() if line.startswith(prefix+' ')]
        assert len(records[prefix])==32,(prefix,len(records[prefix]))
    banks=collections.defaultdict(list)
    for row in records['INFERENCE_STUDY_RESPONSE']:banks[row['index']].append(row)
    maximum=0
    for pair in banks.values():
        assert len(pair)==2 and pair[0]['seed']==pair[1]['seed']
        maximum=max(maximum,max(abs(a-b) for a,b in zip(pair[0]['y'],pair[1]['y'])))
    assert maximum<=1e-12
    summary=[]
    for platform in ['jvm','js']:
        rows=[r for r in records['INFERENCE_STUDY'] if r['platform']==platform]
        profiles=[r for r in records['PROFILE_STUDY'] if r['platform']==platform]
        calibrated=[r for r in rows if r['pvalue'] is not None]
        covered=[r for r in rows if r['coefficient_covered'] is not None]
        profile_covered=[r for r in profiles if r['covered'] is not None]
        count=sum(r['pvalue']<=.1 for r in calibrated)
        summary.append(dict(platform=platform,requested=16,bootstrap_inner_requested=39,
            calibration_qualified=len(calibrated),calibration_refused=16-len(calibrated),
            rejections_at_point_1=count,conditional_rejection_rate=count/len(calibrated) if calibrated else None,
            rejection_fraction_bounds_with_refusals=[count/16,(count+16-len(calibrated))/16],
            pair_attempts=sum(r['pair_attempts'] for r in rows),pair_failures=sum(r['pair_failures'] for r in rows),
            coefficient_qualified=len(covered),coefficient_refused=16-len(covered),coefficient_covered=sum(r['coefficient_covered'] for r in covered),
            coefficient_coverage_bounds_with_refusals=[sum(r['coefficient_covered'] for r in covered)/16,(sum(r['coefficient_covered'] for r in covered)+16-len(covered))/16],
            gaussian_profiles_qualified=len(profile_covered),gaussian_profiles_refused=16-len(profile_covered),gaussian_profiles_covered=sum(r['covered'] for r in profile_covered),
            gaussian_coverage_bounds_with_refusals=[sum(r['covered'] for r in profile_covered)/16,(sum(r['covered'] for r in profile_covered)+16-len(profile_covered))/16],
            median_workflow_ms=statistics.median(r['ns']/1e6 for r in rows),
            median_workflow_allocated_bytes=statistics.median(r['allocated_bytes'] for r in rows) if platform=='jvm' else None,
            median_gaussian_profile_ms=statistics.median(r['ns']/1e6 for r in profiles),
            calibration_error_codes=dict(collections.Counter(r['calibration_error'] for r in rows if r['calibration_error'])),
            coefficient_error_codes=dict(collections.Counter(r['coefficient_error'] for r in rows if r['coefficient_error'])),
            profile_error_codes=dict(collections.Counter(r['error'] for r in profiles if r['error']))))
    for name,digest in inputs.items():assert sha(root/name)==digest,f'source changed during study: {name}'
    result=dict(records=records,summary=summary,maximum_response_difference=maximum,
        interpretation=['16 predetermined outer draws per runtime; same draws across runtimes, not32 independent experiments',
            'Fitted-null plug-in calibration,39 pairs per eligible original; 90percent normal/profile intervals; true omitted slope0',
            'Refusals included in requested denominators and possible rejection bounds; conditional rates labeled separately',
            'Bounded study does not establish nominal coverage or calibration accuracy',
            'Current-thread JVM allocation; JS unavailable; oneprocess per runtime, timings include cold refits and inference'])
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
