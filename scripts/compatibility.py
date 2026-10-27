#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-only
# Part of postalt_rs; see NOTICE.md for BWA/bwakit provenance.
"""Differential tests against k8 bwa-postalt.js r985; no third-party Python packages.

--write-golden refreshes checked-in synthetic fixtures using the reference ONLY.
Normal runs compare Rust and reference byte-for-byte, including HLA FASTQ files.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = '/opt/apps/manual/bwa-mem2/2.3.0/noarch/bin/bwa-postalt.js'

def alt_line(name, pos=1001, cigar='1000M', flag=0, ctg='chr1'):
    """Construct an ALT-to-primary SAM record with one-based POS and absent bases."""
    return f'{name}\t{flag}\t{ctg}\t{pos}\t60\t{cigar}\t*\t0\t0\t*\t*\n'

ALT = ''.join([
    alt_line('altA'), alt_line('altB'), alt_line('rev', 1001, flag=16),
    alt_line('multi'), alt_line('multi', 8001), alt_line('ins', cigar='100M50I850M'),
    alt_line('del', cigar='100M50D900M'), alt_line('clip', cigar='50S900M50H'),
    alt_line('rclip', cigar='50H900M50S', flag=16),
    alt_line('skip', cigar='100M50N900M'), alt_line('unplaced', 0, '*', 4, '*'),
    alt_line('HLA-A*01', 1001), alt_line('HLA-B*01', 1501),
    alt_line('HLA-A*02', 1001, flag=16), alt_line('other', 2001, ctg='chr2'),
])
CONTIGS = ['chr1','chr2','altA','altB','rev','multi','ins','del','clip','rclip','skip','unplaced','HLA-A*01','HLA-A*02','HLA-B*01','other']
HEADER = '@HD\tVN:1.6\tSO:unsorted\n' + ''.join(f'@SQ\tSN:{c}\tLN:20000\n' for c in CONTIGS) + '@RG\tID:rg1\tSM:test\n'

def rec(name='q', ctg='altA', pos=101, cigar='100M', flag=0, mq=10, nm=0, xa=None, tags=(), seq=None, qual=None, mate='*'):
    """Build a synthetic SAM alignment; explicit sequence/quality override defaults.

    Default base/quality lengths follow query-consuming CIGAR operations. Explicit
    overrides intentionally allow missing fields and reference edge cases. This
    helper does not call Rust or derive expected results from the implementation.
    """
    n = sum(int(l) for l,op in re.findall(r'(\d+)([MIS=X])', cigar))
    if seq is None: seq = ('ACGTTGCAN' * ((n+8)//9))[:n] or '*'
    if qual is None: qual = ''.join(chr(33+i%40) for i in range(n)) or '*'
    f = [name,str(flag),ctg,str(pos),str(mq),cigar,mate,'1100' if mate != '*' else '0','-200' if mate != '*' else '0',seq,qual]
    if nm is not None: f += [f'NM:i:{nm}']
    if xa is not None: f += ['XA:Z:'+xa]
    return '\t'.join(f + list(tags)) + '\n'

def cases():
    """Return named (settings, SAM body) cases for reference-derived golden output.

    ratio selects -r; headerless omits the common header. safe_write enables only
    the test oracle's k8 File.write workaround, not a different HLA algorithm.
    Names identify preserved quirks separately from repaired reversal behavior.
    """
    base = rec
    return {
        'passthrough': ({}, base(ctg='chr1',tags=['ZZ:Z:hello world','pa:f:0.500000','RG:Z:rg1'])),
        'unmapped': ({}, base(ctg='*',pos=0,cigar='*',flag=4,seq='*',qual='*')),
        'hard_clip_reported': ({}, base(cigar='10H90M',xa='altB,+101,90M,0;')),
        'single_alt': ({}, base()),
        'single_unplaced_alt': ({}, base(ctg='unplaced')),
        'same_locus': ({}, base(xa='altB,+101,100M,1;')),
        'primary_capped': ({}, base(ctg='chr1',pos=1100,mq=20,xa='altA,+101,100M,0;')),
        'alt_floor': ({}, base(mq=40,xa='altB,+101,100M,0;chr2,+100,100M,0;')),
        'group_score_difference': ({}, base(xa='altB,+101,100M,0;chr2,+100,100M,1;')),
        'other_group_best': ({}, base(ctg='chr1',pos=1100,nm=5,mq=60,xa='altA,+101,100M,5;chr2,+100,100M,0;')),
        'tied_groups': ({}, base(ctg='chr1',pos=1100,mq=60,xa='altA,+101,100M,0;chr2,+100,100M,0;')),
        'reverse_alt': ({}, base(ctg='rev',pos=801)),
        'reverse_generated': ({}, base(xa='rev,-801,100M,0;',flag=65,mate='=',tags=['RG:Z:rg1'])),
        'reverse_reported': ({}, base(flag=16,xa='altB,+101,100M,0;')),
        'iupac_reverse': ({'safe_write':True}, base(cigar='33M',seq='WSATUGCYRKMBDHVNwsatugcyrkmbdhvnA',xa='altB,-101,33M,0;')),
        'missing_sequence': ({'safe_write':True}, base(seq='*',qual='*',xa='altB,-101,100M,0;')),
        'missing_quality': ({'safe_write':True}, base(qual='*',xa='altB,-101,100M,0;')),
        'multi_lift': ({}, base(ctg='multi',xa='altA,+101,100M,0;')),
        'insertion_projection': ({}, base(ctg='ins',pos=111,cigar='20M')),
        'deletion_projection': ({}, base(ctg='del',pos=81,cigar='50M')),
        'skip_projection': ({}, base(ctg='skip',pos=81,cigar='50M')),
        'clipped_start': ({}, base(ctg='clip',pos=31)),
        'clipped_end_quirk': ({}, base(ctg='clip',pos=901)),
        'reverse_clip': ({}, base(ctg='rclip',pos=101)),
        'hard_clip_xa': ({}, base(xa='altB,+101,5H95M,0;')),
        'nm_clamp_and_gaps': ({}, base(cigar='30M10I20M5D40M',nm=0,xa='altB,+101,30M10I20M5D40M,0;chr2,+101,100M,5;')),
        'equal_x_ignored': ({}, base(cigar='20=10X70M',xa='altB,+101,20=10X70M,0;')),
        'no_nm': ({}, base(nm=None,xa='altB,+101,100M,0;')),
        'empty_xa': ({}, base(xa=';')),
        'duplicate_tags': ({}, base(xa='altB,+101,100M,0;',tags=['om:i:7','lt:Z:old;','RG:Z:rg1','ZZ:H:00FF','ZZ:Z:duplicate'])),
        'nonoverlap_pa': ({}, base(ctg='chr2',tags=['pa:f:0.9'])+base(flag=2048)),
        'ratio_0_5': ({'ratio':0.5}, base(ctg='chr2',tags=['pa:f:0.9'])+base(flag=2048)),
        'ratio_capped': ({'ratio':2}, base(ctg='chr2',tags=['pa:f:0.9'])+base(flag=2048)),
        'ratio_negative': ({'ratio':-1}, base(ctg='chr2',tags=['pa:f:0'])+base(flag=2048)),
        'overlap_adjustment': ({}, base(ctg='chr1',pos=1100,mq=60,tags=['pa:f:0.9'])+base(flag=2048,mq=20)),
        'opposite_strand': ({}, base(ctg='chr1',pos=1100,mq=60,flag=16,tags=['pa:f:0.9'])+base(flag=2048)),
        'restore_om': ({}, base(ctg='chr1',pos=1100,mq=0,tags=['om:i:50','pa:f:0.9'])+base(flag=2048,mq=20)),
        'zero_om_not_restored': ({}, base(ctg='chr1',pos=1100,mq=0,tags=['om:i:0','pa:f:0.9'])+base(flag=2048,mq=20)),
        'current_cigar_overlap_quirk': ({}, base(ctg='chr1',pos=1000,cigar='10M',tags=['pa:f:0.9'])+base(pos=51,cigar='100M',flag=2048)),
        'read_end_boundary': ({}, base(ctg='chr2',flag=65,tags=['pa:f:0.9'])+base(flag=129)),
        'name_boundary': ({}, base(name='q1',ctg='chr2',tags=['pa:f:0.9'])+base(name='q2')),
        'hla_primary_hard_clip': ({}, base(ctg='chr1',pos=1100,cigar='10H90M')),
        'hla_union_first_record': ({}, base(ctg='chr2',flag=81)+base(pos=501,flag=2129)),
        'multiple_lifts_no_pa_reduction': ({}, base(ctg='chr2',tags=['pa:f:0.1'])+base(ctg='multi',flag=2048)),
        'transitive_overlap': ({}, base(ctg='chr1',pos=1100,xa='altA,+151,100M,0;altB,+201,100M,0;')),
        'touching_not_overlap': ({}, base(xa='altB,+201,100M,0;')),
        'header_only': ({}, ''),
        'headerless': ({'headerless':True}, base()),
        'final_line_no_newline': ({}, base().rstrip('\n')),
    }

def randomized(seed, count):
    """Generate reproducible template groups exercising interacting r985 branches.

    Each template has contiguous read-end groups, possibly with several records
    and competing XA hits. CIGAR/score combinations probe the postalt algorithm;
    generated XA supplementary records are not guaranteed to satisfy every BAM
    length invariant, so this corpus is compared as SAM.
    """
    rng = random.Random(seed)
    cigars = ['100M','10S90M','90M10S','30M5I65M','30M5D70M','30M5N70M','10H90M','30=10X60M','10M2P90M']
    ctgs = CONTIGS
    lines = []
    for i in range(count):
        for end in rng.sample([0,64,128], rng.randint(1,3)):
            for j in range(rng.randint(1,3)):
                ctg = rng.choice(ctgs)
                pos = rng.randint(1,1200) if ctg not in ['chr1','chr2'] else rng.randint(900,2200)
                xa = ''.join(f'{rng.choice(ctgs)},{rng.choice(["+","-"])}{rng.randint(1,2200)},{rng.choice(cigars)},{rng.randint(0,15)};' for _ in range(rng.randint(0,7)))
                tags = [f'pa:f:{rng.choice([0,0.4,0.99,1,1.1])}','RG:Z:rg1']
                if rng.random() < .1: tags.append(f'om:i:{rng.randint(0,60)}')
                lines.append(rec(name=f'random{i}',ctg=ctg,pos=pos,cigar=rng.choice(cigars),flag=end|rng.choice([0,16])|(2048 if j else 0),mq=rng.randint(0,60),nm=rng.randint(0,15),xa=xa,tags=tags))
    return ''.join(lines)

def run(cmd, alt, sam, settings, directory, prefix):
    """Run an oracle or replacement and return SAM bytes plus HLA file bytes.

    Clear this prefix's previous FASTQs before each run so absent/empty results
    cannot accidentally compare against stale files. Raise on any subprocess
    failure; never accept partial stdout from a crashed k8 run as golden output.
    """
    outprefix = directory/prefix
    for path in directory.glob(prefix+'.*.fq'): path.unlink()
    args = [*cmd, *([] if settings.get('no_hla') else ['-p', str(outprefix)]), '-r', str(settings.get('ratio',1)), str(alt), str(sam)]
    proc = subprocess.run(args, capture_output=True)
    if proc.returncode: raise RuntimeError(f'{args}: {proc.stderr.decode(errors="replace")}')
    fq = {p.name.removeprefix(prefix+'.'): p.read_bytes() for p in directory.glob(prefix+'.*.fq')}
    return proc.stdout, fq

def safe_file_write(text):
    """Wrap only k8's file-writing transport to avoid its string-buffer overflow.

    A Uint8Array copies identical ASCII FASTQ bytes to the ArrayBuffer overload.
    Do not use k8_encode here: its string conversion has the same missing-NUL-room
    issue. The original script's HLA selection/formatting remains unchanged.
    This wrapper is test-only and never modifies the installed reference file.
    """
    # k8 1.2 File.write(string) allocates len bytes but V8 also writes a NUL.
    # The ArrayBuffer overload avoids that heap overflow, with identical bytes.
    wrapper = "var original_write = File.prototype.write; File.prototype.write = function(s) { var b = new Uint8Array(s.length); for (var i=0; i<s.length; ++i) b[i]=s.charCodeAt(i); return original_write.call(this, b.buffer); };\n"
    return wrapper + text.removeprefix('#!/bin/env k8\n')

def repaired_reference(text):
    """Repair the oracle's two reversal expressions, leaving its algorithm intact.

    Match the chosen default behavior: use the reference IUPAC table, preserve '*'
    and unknown characters, and reverse quality independently. Assertions make
    an unexpected reference version fail rather than silently skip the repair.
    """
    # Only replace Bytes-based reversal; all lift/group/MAPQ/HLA code is untouched.
    old = "aux.set(t[9], 0); aux.revcomp(); rs = aux.toString();"
    new = "rs = t[9] == '*' ? '*' : t[9].split('').reverse().map(function(c) { var a='WSATUGCYRKMBDHVNwsatugcyrkmbdhvn', b='WSTAACGRYMKVHDBNwstaacgrymkvhdbn'; var i=a.indexOf(c); return i < 0 ? c : b.charAt(i); }).join('');"
    assert old in text
    text = text.replace(old,new)
    old = "aux.set(t[10],0); aux.reverse(); rq = aux.toString();"
    assert old in text
    return text.replace(old, "rq = t[10].split('').reverse().join('');")

def main():
    """Compare both modes, or explicitly regenerate independent golden fixtures.

    Failures retain input/actual/expected files under target for inspection. The
    manifest hashes the original reference and records case-specific workarounds;
    temporary repaired/wrapped scripts are never substituted for that provenance.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', default=str(ROOT/'target/release/postalt_rs'))
    p.add_argument('--reference', default=os.environ.get('POSTALT_REFERENCE',REFERENCE))
    p.add_argument('--k8', default=os.environ.get('K8','k8'))
    p.add_argument('--write-golden',action='store_true')
    p.add_argument('--random',type=int,default=5000)
    p.add_argument('--seed',type=int,default=985)
    a = p.parse_args()
    manifest = []
    suite = cases()
    if not a.write_golden: suite['randomized'] = ({'safe_write':True},randomized(a.seed,a.random))
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        alt = temp/'alt.sam'; alt.write_text(ALT)
        fixed = temp/'fixed.js'; fixed.write_text(safe_file_write(repaired_reference(Path(a.reference).read_text())))
        safe = temp/'safe.js'; safe.write_text(safe_file_write(Path(a.reference).read_text()))
        for name,(settings,body) in suite.items():
            sam = temp/'input.sam'; sam.write_text(('' if settings.get('headerless') else HEADER)+body)
            oracle,fastq = run([a.k8,str(safe) if settings.get('safe_write') else a.reference],alt,sam,settings,temp,'js')
            corrected,corrected_fq = run([a.k8,str(fixed)],alt,sam,settings,temp,'fixed')
            # Expected files always originate from k8. Rust is invoked only in
            # comparison mode, so a regression cannot rewrite its own oracle.
            if a.write_golden:
                d = ROOT/'tests/fixtures'/name; d.mkdir(parents=True,exist_ok=True)
                for path in d.glob('*.fq'): path.unlink()
                (d/'input.sam').write_bytes(sam.read_bytes()); (d/'expected.sam').write_bytes(corrected); (d/'k8.sam').write_bytes(oracle)
                for file,data in fastq.items(): (d/file).write_bytes(data)
                manifest.append({'name':name, **settings})
            else:
                actual,fq = run([a.binary,'--k8-1.2-compat'],alt,sam,settings,temp,'rs')
                if actual != oracle or fq != fastq:
                    failure = ROOT/'target/compatibility-failure'/name; failure.mkdir(parents=True,exist_ok=True)
                    (failure/'input.sam').write_bytes(sam.read_bytes()); (failure/'expected.sam').write_bytes(oracle); (failure/'actual.sam').write_bytes(actual)
                    for file,data in fastq.items(): (failure/('expected.'+file)).write_bytes(data)
                    for file,data in fq.items(): (failure/('actual.'+file)).write_bytes(data)
                    raise AssertionError(f'{name}: differs; inspect {failure}')
                actual,fq = run([a.binary],alt,sam,settings,temp,'rsfixed')
                if actual != corrected or fq != corrected_fq:
                    failure = ROOT/'target/compatibility-failure'/name; failure.mkdir(parents=True,exist_ok=True)
                    (failure/'input.sam').write_bytes(sam.read_bytes()); (failure/'expected.sam').write_bytes(corrected); (failure/'actual.sam').write_bytes(actual)
                    raise AssertionError(f'{name}: corrected mode differs; inspect {failure}')
                print(f'PASS {name} (both modes)')
    if a.write_golden:
        (ROOT/'tests/fixtures/alt.sam').write_text(ALT)
        (ROOT/'tests/fixtures/manifest.json').write_text(json.dumps({'reference':'bwa-postalt.js r985','sha256':hashlib.sha256(Path(a.reference).read_bytes()).hexdigest(),'cases':manifest},indent=2)+'\n')
        print(f'Wrote {len(manifest)} independent k8 golden cases')
    else: print(f'All {len(suite)} differential cases passed (random seed {a.seed}, {a.random} templates)')

if __name__ == '__main__': main()
