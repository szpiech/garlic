"""Base-level scoring of every ROH call set against ped-sim HBD truth and against the 'attainable' set
(pedigree HBD plus the segments GARLIC calls on the true genotypes, which are background founder
autozygosity).  One row per (method, depth, sample); pooled summary per (method, depth)."""
import sys, os
sys.path.insert(0, snakemake.params.scriptdir)
import pandas as pd, numpy as np
from common import mapped_ranges, read_roh_bed, read_rohan, clip, merge, overlap, length

cfg = snakemake.config; P = snakemake.params
CHRS = {k: int(v) for k, v in P.chrs.items()}
MIN = int(cfg['score']['min_len'])
rng = mapped_ranges(cfg['tools']['simmap'], CHRS)
genome = sum(hi - lo + 1 for lo, hi in rng.values())
samples = [l.strip() for l in open(snakemake.input.samples) if l.strip()]
truth = {s: [] for s in samples}
for l in open(snakemake.input.truth):
    c, a, b, s = l.split()
    if s in truth: truth[s].append((c, int(a), int(b)))
truth = {s: clip(v, rng) for s, v in truth.items()}
ceiling = {s: clip(v, rng) for s, v in read_roh_bed(snakemake.input.ceiling, samples, MIN).items()}
attain = {s: merge(truth[s] + ceiling[s]) for s in samples}

rows = []
def add(method, depth, calls):
    for s in samples:
        if s not in calls: continue
        cl = clip(calls[s], rng)
        rows.append(dict(method=method, depth=depth, sample=s, pedigree=s[:2], truth_bp=length(truth[s]), called_bp=length(cl),
                         ov_truth=overlap(cl, truth[s]), ov_attain=overlap(cl, attain[s]), ov_ceiling=overlap(cl, ceiling[s]),
                         ceiling_bp=length(ceiling[s]), froh_truth=length(truth[s]) / genome, froh_called=length(cl) / genome,
                         froh_ceiling=length(ceiling[s]) / genome))
add('garlic_true_genotypes', 'inf', ceiling)
for fn in snakemake.input.garlic:                      # roh/<arm>_d<depth>.roh.bed
    stem = os.path.basename(fn).replace('.roh.bed', ''); arm, d = stem.rsplit('_d', 1)
    add('garlic_' + arm, int(d), read_roh_bed(fn, samples, MIN))
for est in ('mid', 'min', 'max'):
    by = {}
    for fn in snakemake.input.rohan:                   # roh/rohan/<sample>_d<depth>/out.mid.hmmrohl.gz
        sd = os.path.basename(os.path.dirname(fn)); s, d = sd.rsplit('_d', 1)
        f2 = fn.replace('out.mid.', f'out.{est}.')
        if os.path.exists(f2): by.setdefault(int(d), {})[s] = read_rohan(f2, MIN)
    for d, calls in by.items(): add('rohan' if est == 'mid' else f'rohan_{est}', d, calls)
df = pd.DataFrame(rows)
df.to_csv(snakemake.output.per_individual, index=False)
g = df.groupby(['method', 'depth'], sort=False)
summ = pd.DataFrame({'n_individuals': g['sample'].size(),
                     'precision_vs_pedigree': g.ov_truth.sum() / g.called_bp.sum(),
                     'precision_vs_attainable': g.ov_attain.sum() / g.called_bp.sum(),
                     'recall_pedigree': g.ov_truth.sum() / g.truth_bp.sum(),
                     'recall_ceiling_calls': g.ov_ceiling.sum() / g.ceiling_bp.sum(),
                     'froh_called_mean': g.froh_called.mean(), 'froh_truth_mean': g.froh_truth.mean(),
                     'froh_mae_vs_ceiling': g.apply(lambda x: (x.froh_called - x.froh_ceiling).abs().mean(), include_groups=False)}).reset_index()
summ.round(4).to_csv(snakemake.output.summary, index=False)
print(summ.round(3).to_string(index=False), file=sys.stderr)
