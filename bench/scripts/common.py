"""Shared helpers for the benchmark workflow (imported by the scripts, not by the Snakefile)."""
import gzip, re

def read_simmap(path, chrs):
    """Per-chromosome (pos, sex-averaged cM) arrays from a Bherer-style simmap."""
    import numpy as np
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            c, pos, m, f = line.split()[:4]
            if c in chrs: out.setdefault(c, []).append((int(pos), (float(m) + float(f)) / 2))
    return {c: np.array(v) for c, v in out.items()}

def mapped_ranges(simmap_path, chrs):
    """[first, last] map position per chromosome, clipped to the simulated length.  ped-sim transmits
    nothing outside the map, so truth and genotypes are reference-homozygous there."""
    rng = {}
    with open(simmap_path) as fh:
        next(fh)
        for line in fh:
            c, pos = line.split()[:2]; pos = int(pos)
            if c not in chrs: continue
            lo, hi = rng.get(c, (pos, pos)); rng[c] = (min(lo, pos), max(hi, pos))
    return {c: (lo, min(hi, chrs[c])) for c, (lo, hi) in rng.items()}

def read_roh_bed(fn, samples, min_len=0):
    iv = {s: [] for s in samples}; cur = None
    for line in open(fn):
        if line.startswith('track'): cur = re.search(r'Ind: (\S+)', line).group(1); iv.setdefault(cur, [])
        elif line.strip():
            f = line.split(); a, b = int(f[1]), int(f[2])
            if b - a >= min_len: iv[cur].append((f[0].replace('chr', ''), a, b))
    return iv

def read_rohan(fn, min_len=0):
    """ROHan .hmmrohl.gz: '#ROH_ID CHROM BEGIN END ROH_LENGTH VALIDATED_SITES'."""
    iv = []
    with gzip.open(fn, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            f = line.split()
            if len(f) < 4: continue
            c, a, b = f[1].replace('chr', ''), int(f[2]), int(f[3])
            if b - a >= min_len: iv.append((c, a, b))
    return iv

def clip(iv, rng):
    out = []
    for c, a, b in iv:
        if c not in rng: continue
        lo, hi = rng[c]; a, b = max(a, lo), min(b, hi)
        if b > a: out.append((c, a, b))
    return out

def merge(iv):
    out = []
    for c, a, b in sorted(iv):
        if out and out[-1][0] == c and a <= out[-1][2]: out[-1] = (c, out[-1][1], max(b, out[-1][2]))
        else: out.append((c, a, b))
    return out

def overlap(A, B):
    tot = 0; by = {}
    for c, a, b in B: by.setdefault(c, []).append((a, b))
    for c in by: by[c].sort()
    for c, a, b in A:
        for a2, b2 in by.get(c, []):
            if b2 < a: continue
            if a2 > b: break
            tot += max(0, min(b, b2) - max(a, a2))
    return tot

def length(iv): return sum(b - a for _, a, b in iv)
