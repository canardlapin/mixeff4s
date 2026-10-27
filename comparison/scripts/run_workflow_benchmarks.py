#!/usr/bin/env python3
"""Sequential matched workflow runs; retain raw outcomes, allocation and independent assessment cost."""
import argparse,collections,hashlib,json,os,pathlib,statistics,subprocess,platform
p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--repeats',type=int,default=3);a=p.parse_args()
if a.repeats<1:p.error('repeats must be positive')
root=pathlib.Path(__file__).resolve().parents[2];out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
inputs=list((root/'src').rglob('*.scala'))+[root/'build.sbt',root/'project/plugins.sbt',root/'project/build.properties']
for directory in ['.jvm/src','.js/src']: inputs+=list((root/directory).rglob('*.scala'))
hashes={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs}
(out/'input-sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
all_rows=[]
reference_responses={}
maximum_response_error=0.0
for index in range(a.repeats):
 log=out/f'run-{index+1}.log'
 cmd=['python3',str(pathlib.Path.home()/'.agents/skills/lean-logs/scripts/run_logged.py'),'--log',str(log),'--timeout','900','--','sbt','-batch','mixeff4sJVM/testOnly mixeff4s.benchmark.WorkflowBenchmarkSuite','mixeff4sJS/testOnly mixeff4s.benchmark.WorkflowBenchmarkSuite']
 subprocess.run(cmd,cwd=root,env=dict(os.environ,MIXEFF4S_BENCHMARK='1'),check=True)
 rows=[json.loads(line[len('FLOW_BENCH '):]) for line in log.read_text().splitlines() if line.startswith('FLOW_BENCH ')]
 assert len(rows)==480,(log,len(rows))
 responses=[json.loads(line[len('FLOW_RESPONSE '):]) for line in log.read_text().splitlines() if line.startswith('FLOW_RESPONSE ')]
 assert len(responses)==78,(log,len(responses))
 for response in responses:
  key=response['index']
  if key not in reference_responses: reference_responses[key]=response
  expected=reference_responses[key]
  assert response['seed']==expected['seed']
  assert len(response['y'])==len(expected['y'])
  maximum_response_error=max(maximum_response_error,max(abs(x-y) for x,y in zip(response['y'],expected['y'])))
  assert maximum_response_error<=1e-10,maximum_response_error
 (out/f'run-{index+1}.json').write_text(json.dumps(rows,indent=2)+'\n');all_rows.append(rows)
summaries=[]
keys=sorted({(r['platform'],r['workload'],r['strategy'],r['phase']) for r in all_rows[0]})
for key in keys:
 groups=[[r for r in rows if (r['platform'],r['workload'],r['strategy'],r['phase'])==key] for rows in all_rows]
 medians=[statistics.median(r['ns']/1e6 for r in group) for group in groups]
 alloc=[statistics.median(r['allocated_bytes'] for r in group) if all(r['allocated_bytes'] is not None for r in group) else None for group in groups]
 summaries.append(dict(platform=key[0],workload=key[1],strategy=key[2],phase=key[3],process_medians_ms=medians,median_ms=statistics.median(medians),allocated_bytes=statistics.median(alloc) if all(x is not None for x in alloc) else None,fit_failures=sum(r['error'] is not None for g in groups for r in g),unqualified=sum(r['error'] is None and not r['qualified'] for g in groups for r in g),attempts=sum(map(len,groups)),assessment_ms=statistics.median(statistics.median(r['assessment_ns']/1e6 for r in g) for g in groups)))
for name,digest in hashes.items(): assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest,f'source changed during timing: {name}'
result=dict(maximum_seeded_response_difference=maximum_response_error,runs=a.repeats,samples=sum(map(len,all_rows)),host=dict(platform=platform.platform(),machine=platform.machine()),summary=summaries,limits=['Each response independently draws new random effects and residual noise','Refit starts are cold; rebuild includes response frame construction and design compilation','Assessment time is measured separately; fit-failure and failed-assessment counts are distinct','Current-thread JVM allocation only; JS allocation unavailable; no peak-memory claim','Process medians on this host are not universal performance guarantees'])
(out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
print(f"{result['samples']} samples, {len(summaries)} summaries; raw outcomes retained at {out}")
for row in summaries:
 if row['phase']=='refit':
  before=next(r for r in summaries if (r['platform'],r['workload'],r['strategy'],r['phase'])==(row['platform'],row['workload'],row['strategy'],'rebuild'))
  print(row['platform'],row['workload'],row['strategy'],round(before['median_ms'],3),'->',round(row['median_ms'],3),'ms; failed',row['fit_failures'],'unqualified',row['unqualified'],'of',row['attempts'])
