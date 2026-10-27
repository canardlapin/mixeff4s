#!/usr/bin/env python3
"""Sequential independent sbt/Node benchmark runs, with raw output and exit receipts.

Requires the same JDK/sbt/Node and pinned Gale access as the normal build. No
performance verdict is inferred here. A baseline archive is extracted under a
new temporary directory and hash-checked before execution.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--repeats', type=int, default=3)
parser.add_argument('--baseline', type=Path)
args = parser.parse_args()
if args.repeats < 1:
    parser.error('--repeats must be positive')
root = Path(__file__).resolve().parents[2]
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
working = root
if args.baseline:
    archive_dir = args.baseline.resolve()
    working = Path(tempfile.mkdtemp(prefix='mixeff4s-baseline-'))
    with tarfile.open(archive_dir/'source.tar.gz') as archive:
        archive.extractall(working, filter='data')
    hashes = json.loads((archive_dir/'source-sha256.json').read_text())
    for name, digest in hashes.items():
        if hashlib.sha256((working/name).read_bytes()).hexdigest() != digest:
            raise SystemExit(f'baseline hash mismatch: {name}')
else:
    inputs = list((root/'src').rglob('*.scala')) + [root/'build.sbt', root/'project/plugins.sbt', root/'project/build.properties']
    for platform in ['.jvm', '.js']:
        inputs += list((root/platform/'src').rglob('*.scala'))
    hashes = {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
(output/'input-sha256.json').write_text(json.dumps(hashes, indent=2)+'\n')
(output/'run-settings.json').write_text(json.dumps(dict(working_directory=str(working),repeats=args.repeats,baseline=bool(args.baseline)),indent=2)+'\n')
env = dict(os.environ, MIXEFF4S_BENCHMARK='1')
runner = Path.home()/'.agents/skills/lean-logs/scripts/run_logged.py'
for index in range(args.repeats):
    log = output/f'run-{index+1}.log'
    command = ['python3',str(runner),'--log',str(log),'--timeout','900','--','sbt','-batch',
               'mixeff4sJVM/testOnly mixeff4s.benchmark.EngineBenchmarkSuite',
               'mixeff4sJS/testOnly mixeff4s.benchmark.EngineBenchmarkSuite']
    outcome = subprocess.run(command,cwd=working,env=env)
    if outcome.returncode:
        raise SystemExit(outcome.returncode)
    subprocess.run(['python3',str(root/'comparison/scripts/engine_benchmark.py'),str(log),str(output/f'run-{index+1}.json')],check=True)
print(f'Completed {args.repeats} independent runs; sources at {working}; receipts at {output}')
