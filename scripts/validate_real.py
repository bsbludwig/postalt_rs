#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
# Part of postalt_rs; see NOTICE.md for BWA/bwakit provenance.
"""Compare SAM and BAM paths on a prepared dataset without saving large outputs.
Requires k8, the reference script, and samtools. Prints and optionally saves hashes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from compatibility import REFERENCE, repaired_reference


def digest(command, normalize=False):
    """Hash streamed output and count bytes/lines without retaining alignments.

    normalize runs samtools view -h --no-PG so BAM and SAM compare after the same
    serialization rules (notably float spelling and reference aliases). It does
    not sort or filter. Both producer and normalizer must exit successfully.
    """
    p = subprocess.Popen(command, stdout=subprocess.PIPE)
    normalizer = None
    stream = p.stdout
    if normalize:
        normalizer = subprocess.Popen(['samtools','view','--no-PG','-h','-'], stdin=p.stdout, stdout=subprocess.PIPE)
        p.stdout.close()
        stream = normalizer.stdout
    sha = hashlib.sha256()
    total = lines = 0
    while True:
        chunk = stream.read(1024*1024)
        if not chunk: break
        sha.update(chunk); total += len(chunk); lines += chunk.count(b'\n')
    stream.close()
    if normalizer and normalizer.wait(): raise RuntimeError('samtools normalization failed')
    if p.wait(): raise RuntimeError(f'command failed: {command}')
    return dict(sha256=sha.hexdigest(),bytes=total,lines=lines)

def main():
    """Check unmodified/repaired reference outputs and optional normalized BAM.

    Real-data runs omit HLA extraction, so they need no File.write oracle wrapper.
    Only default-mode sequence/quality reversal is repaired in a temporary script;
    compatibility mode executes the installed reference unchanged.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('alt')
    p.add_argument('sam')
    p.add_argument('--bam')
    p.add_argument('--binary',default='target/release/postalt_rs')
    p.add_argument('--reference',default=os.environ.get('POSTALT_REFERENCE',REFERENCE))
    p.add_argument('--k8',default=os.environ.get('K8','k8'))
    p.add_argument('--json',type=Path)
    a = p.parse_args()
    result = {}
    with tempfile.TemporaryDirectory() as temp:
        fixed = Path(temp)/'fixed.js'; fixed.write_text(repaired_reference(Path(a.reference).read_text()))
        for mode,reference,extra in [('compat',a.reference,['--k8-1.2-compat']),('default',str(fixed),[])]:
            expected = digest([a.k8,reference,a.alt,a.sam])
            actual = digest([a.binary,*extra,a.alt,a.sam])
            if actual != expected: raise AssertionError(f'{mode}: reference {expected}, Rust {actual}')
            result[mode] = actual
            print(f'PASS {mode}: {actual}',flush=True)
        if a.bam:
            expected = digest([a.binary,a.alt,a.sam],normalize=True)
            actual = digest([a.binary,'-O','bam','-@','2',a.alt,a.bam],normalize=True)
            if actual != expected: raise AssertionError(f'BAM: expected {expected}, actual {actual}')
            result['bam_normalized'] = actual
            print(f'PASS BAM: {actual}',flush=True)
    if a.json:
        a.json.parent.mkdir(parents=True,exist_ok=True)
        a.json.write_text(json.dumps(result,indent=2)+'\n')

if __name__ == "__main__": main()
