#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
# Part of postalt_rs; see NOTICE.md for BWA/bwakit provenance.
"""Repeat sequential wall/CPU/RSS measurements; discard output, retain JSON metadata.

Uses GNU /usr/bin/time: elapsed/user/system values are seconds and maximum RSS is
KiB. Measurements include process startup and ALT-index construction. This script
neither flushes the filesystem cache nor reserves a CPU; report cache state and
host contention separately. SAM and BAM timings have different I/O workloads.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('alt')
p.add_argument('sam')
p.add_argument('--bam')
p.add_argument('--binary', default='target/release/postalt_rs')
p.add_argument('--reference', default=os.environ.get('POSTALT_REFERENCE','/opt/apps/manual/bwa-mem2/2.3.0/noarch/bin/bwa-postalt.js'))
p.add_argument('--k8', default=os.environ.get('K8','k8'))
p.add_argument('--repeats', type=int, default=3)
p.add_argument('--json', type=Path, required=True)
a = p.parse_args()
if a.repeats < 1: p.error('--repeats must be positive')
commands = {
    'k8_sam': [a.k8,a.reference,a.alt,a.sam],
    'rust_sam_compat': [a.binary,'--k8-1.2-compat',a.alt,a.sam],
    'rust_sam_default': [a.binary,a.alt,a.sam],
}
if a.bam: commands['rust_bam_2_workers'] = [a.binary,'-O','bam','-@','2',a.alt,a.bam]
result = {'platform': platform.platform(), 'reference_sha256': hashlib.sha256(Path(a.reference).read_bytes()).hexdigest(),
          'input_sam_bytes': Path(a.sam).stat().st_size, 'runs':{}, 'median':{}}
with tempfile.TemporaryDirectory() as tmp:
    timing = Path(tmp)/'time'
    # Run sequentially to avoid making competing implementations time each
    # other's CPU/I/O pressure. No internal benchmark subprocesses overlap.
    for name,cmd in commands.items():
        runs = []
        for rep in range(a.repeats):
            with open(os.devnull,'wb') as sink:
                subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o',str(timing),*cmd],stdout=sink,check=True)
            wall,user,system,rss = map(float,timing.read_text().split())
            runs.append(dict(wall_seconds=wall,user_seconds=user,system_seconds=system,max_rss_kib=int(rss)))
            print(name,rep+1,runs[-1],flush=True)
        result['runs'][name] = runs
        result['median'][name] = {k:statistics.median(r[k] for r in runs) for k in runs[0]}
# Compare identical SAM input/output modes. The BAM row additionally performs
# decompression/compression and is deliberately excluded from this speedup ratio.
result['sam_speedup'] = result['median']['k8_sam']['wall_seconds']/result['median']['rust_sam_compat']['wall_seconds']
a.json.parent.mkdir(parents=True,exist_ok=True)
a.json.write_text(json.dumps(result,indent=2)+'\n')
print('SAM speedup:',result['sam_speedup'])
