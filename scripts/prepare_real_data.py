#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
# Part of postalt_rs; see NOTICE.md for BWA/bwakit provenance.
"""Extract indexed ALT + MHC reads; restore contiguous read-end groups.

Requires samtools. Sorting cannot reconstruct the aligner's exact original order;
within each read end, original primary records precede secondary/supplementary.
No clinical reads are checked into the repository.
"""
import argparse
import itertools
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('bam')
p.add_argument('directory', type=Path)
a = p.parse_args()
a.directory.mkdir(parents=True, exist_ok=True)
stats = subprocess.check_output(['samtools', 'idxstats', a.bam], text=True)
# Select regional alignments, not all distant records with a selected QNAME.
# Include a primary MHC region so earlier primary/ALT interactions are exercised.
regions = [l.split('\t')[0] for l in stats.splitlines() if '_alt\t' in l]
regions += ['chr6:30000000-31000000']
extract = subprocess.Popen(['samtools', 'view', '--no-PG', '-u', a.bam, *regions], stdout=subprocess.PIPE)
subprocess.run(['samtools', 'sort', '--no-PG', '-n', '-@', '2', '-m', '512M', '-o', str(a.directory/'names.bam'), '-'], stdin=extract.stdout, check=True)
extract.stdout.close()
if extract.wait(): raise RuntimeError('samtools extraction failed')
view = subprocess.Popen(['samtools', 'view', '--no-PG', '-h', str(a.directory/'names.bam')], stdout=subprocess.PIPE, text=True)
with (a.directory/'input.sam').open('w') as out:
    records = (line for line in view.stdout)
    first = None
    for line in records:
        if line.startswith('@'): out.write(line)
        else: first = line; break
    # Name sorting alone can interleave mates. Split by FLAG 0xc0 and put primary
    # records before secondary/supplementary (0x900); stable ties retain samtools
    # order. This approximates aligner order, which coordinate sorting has lost.
    if first:
        for _, group in itertools.groupby(itertools.chain([first], records), lambda l: l.split('\t', 1)[0]):
            out.writelines(sorted(group, key=lambda l: (int(l.split('\t', 2)[1]) & 0xc0, bool(int(l.split('\t', 2)[1]) & 0x900))))
if view.wait(): raise RuntimeError('samtools view failed')
subprocess.run(['samtools', 'view', '--no-PG', '-b', '-o', str(a.directory/'input.bam'), str(a.directory/'input.sam')], check=True)
print(a.directory/'input.sam')
