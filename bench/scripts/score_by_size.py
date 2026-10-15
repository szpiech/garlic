"""Recall and precision by ROH size class for every caller (GARLIC arms, PLINK, ROHan) at every depth.

Truth is truth/ibd_all.bed.gz: every within-individual IBD tract >= 10 kb with its TMRCA (chr start end sample
tmrca span, 1-based inclusive), clipped to the genetic-map span.  Size classes follow GARLIC's --size-bounds:
0.1-0.5, 0.5-1 and >1 Mb (calls under 0.1 Mb are reported for precision only).
  recall (class k)    bases of truth tracts of length k (>= min_truth, TMRCA <= max_tmrca) covered by calls of any length
  recall_any_age      the same with no TMRCA limit
  precision (class k) bases of calls of length k inside truth tracts >= min_truth of any age
  precision_10kb      the same against tracts >= 10 kb
  frac_<TMRCA class>  composition of the called bases of class k by the TMRCA of the tract they lie in
An arm is scored over the individuals it has output for (all individuals for GARLIC and PLINK; ROHan per BAM).
An empty GARLIC .roh.bed (no automatic cutoff) is a result: zero calls."""
import sys, os, re, gzip
import numpy as np, pandas as pd

S = snakemake
MIN_T, MAX_T = float(S.params.min_truth), float(S.params.max_tmrca)
CLASSES = [(0, 1e5, '<0.1 Mb'), (1e5, 5e5, '0.1-0.5 Mb'), (5e5, 1e6, '0.5-1 Mb'), (1e6, np.inf, '>1 Mb')]
TBINS = [-np.inf, 20, 100, 1000, np.inf]; TLAB = ['tmrca_le20', 'tmrca_20_100', 'tmrca_100_1000', 'tmrca_gt1000']
samples = [l.strip() for l in open(S.input.samples) if l.strip()]
mapped = {}
for l in open(S.input.mapped):
    c, a, b = l.split()[:3]; mapped[c.replace('chr', '')] = (int(a) + 1, int(b))

def clip(df):
    df = df[df.chr.isin(mapped)].copy()
    lo = df.chr.map(lambda c: mapped[c][0]); hi = df.chr.map(lambda c: mapped[c][1])
    df['start'] = np.maximum(df.start, lo); df['end'] = np.minimum(df.end, hi)
    df = df[df.end >= df.start].copy(); df['len'] = df.end - df.start + 1
    return df

def size_class(L):
    return pd.cut(L, [c[0] for c in CLASSES] + [np.inf], labels=[c[2] for c in CLASSES], right=False).astype(str)

def merge_iv(a):
    if len(a) == 0: return np.zeros((0, 2), dtype=np.int64)
    a = a[np.argsort(a[:, 0])]; out = [a[0].copy()]
    for s, e in a[1:]:
        if s <= out[-1][1] + 1: out[-1][1] = max(out[-1][1], e)
        else: out.append(np.array([s, e]))
    return np.array(out, dtype=np.int64)

def covered(segs, ref):
    """bases of each segment in segs (n,2) covered by the merged intervals ref (m,2)"""
    out = np.zeros(len(segs), dtype=np.int64)
    if len(ref) == 0 or len(segs) == 0: return out
    rs, re_ = ref[:, 0], ref[:, 1]
    for i, (s, e) in enumerate(segs):
        lo = np.searchsorted(re_, s); hi = np.searchsorted(rs, e, side='right')
        if hi > lo: out[i] = np.sum(np.minimum(e, re_[lo:hi]) - np.maximum(s, rs[lo:hi]) + 1)
    return out

# ------------------------------------------------------------------ truth
truth = pd.read_csv(S.input.truth, sep='\t', header=None, names=['chr', 'start', 'end', 'ind', 'tmrca', 'span'], dtype={'chr': str})
truth = clip(truth[truth.ind.isin(samples)])
truth['size_class'] = size_class(truth.len); truth['tclass'] = pd.cut(truth.tmrca, TBINS, labels=TLAB).astype(str)
REF = {}                                     # (ind, chr) -> merged reference intervals
for (s, c), g in truth.groupby(['ind', 'chr']):
    REF[(s, c)] = {'ge_min': merge_iv(g[g.len >= MIN_T][['start', 'end']].to_numpy()),
                   'ge_10kb': merge_iv(g[['start', 'end']].to_numpy()),
                   **{t: merge_iv(g[g.tclass == t][['start', 'end']].to_numpy()) for t in TLAB}}
EMPTY = np.zeros((0, 2), dtype=np.int64)

# ------------------------------------------------------------------ call sets
def read_garlic(fn):
    rows = []; ind = None
    for l in open(fn):
        if l.startswith('track'): ind = re.search(r'Ind: (\S+)', l).group(1); continue
        p = l.split()
        if len(p) >= 3: rows.append((ind, p[0].replace('chr', ''), int(p[1]) + 1, int(p[2])))   # BED half-open
    return rows

def read_plink(fn):
    rows = []
    with open(fn) as fh:
        hdr = fh.readline().split(); ix = {k: hdr.index(k) for k in ('IID', 'CHR', 'POS1', 'POS2')}
        for l in fh:
            p = l.split()
            if p: rows.append((p[ix['IID']], p[ix['CHR']].replace('chr', ''), int(p[ix['POS1']]), int(p[ix['POS2']])))
    return rows

def read_rohan(fn, ind):
    rows = []
    with gzip.open(fn, 'rt') as fh:
        for l in fh:
            if l.startswith('#'): continue
            p = l.split()
            if len(p) >= 4: rows.append((ind, p[1].replace('chr', ''), int(p[2]), int(p[3])))
    return rows

GARLIC_INPUT = {'calls': 'genotype calls', 'pl': 'genotype calls, --gl-type PL', 'beagle': 'ANGSD likelihoods',
                'beaglee4': 'ANGSD likelihoods, --error 1e-4', 'beaglee2': 'ANGSD likelihoods, --error 1e-2'}
arms = []                                    # (method, input, depth, rows, individuals)
for fn in S.input.garlic:
    stem = os.path.basename(fn).replace('.roh.bed', '')
    if stem == 'truth': arms.append(('GARLIC', 'true genotypes', 'true', read_garlic(fn), samples))
    else:
        a, d = stem.rsplit('_d', 1); arms.append(('GARLIC', GARLIC_INPUT.get(a, a), d, read_garlic(fn), samples))
for fn in S.input.plink:
    h, src = re.match(r'h(\d+)_(\S+)\.hom$', os.path.basename(fn)).groups()
    arms.append((f'PLINK (window-het {h})', 'true genotypes' if src == 'truth' else 'genotype calls',
                 'true' if src == 'truth' else src[1:], read_plink(fn), samples))
by = {}
for fn in S.input.rohan:                    # roh/rohan/<sample>_d<depth>/out.mid.hmmrohl.gz
    s, d = os.path.basename(os.path.dirname(fn)).rsplit('_d', 1)
    for est in ('mid', 'min', 'max'):
        f2 = fn.replace('out.mid.', f'out.{est}.')
        if os.path.exists(f2):
            r = by.setdefault((est, d), ([], [])); r[0].extend(read_rohan(f2, s)); r[1].append(s)
for (est, d), (rows, inds) in sorted(by.items()):
    arms.append(('ROHan' if est == 'mid' else f'ROHan ({est})', 'BAM', d, rows, inds))

# ------------------------------------------------------------------ scoring
summary, per_ind, by_tmrca = [], [], []
for method, inp, depth, rows, inds in arms:
    calls = clip(pd.DataFrame(rows, columns=['ind', 'chr', 'start', 'end']).astype({'chr': str}) if rows else
                 pd.DataFrame({'ind': pd.Series(dtype=str), 'chr': pd.Series(dtype=str), 'start': pd.Series(dtype=int), 'end': pd.Series(dtype=int)}))
    calls = calls[calls.ind.isin(inds)].copy(); calls['size_class'] = size_class(calls.len)
    for col in ['ov_min', 'ov_10kb'] + TLAB: calls[col] = 0
    for (s, c), g in calls.groupby(['ind', 'chr']):
        ref = REF.get((s, c)); segs = g[['start', 'end']].to_numpy()
        calls.loc[g.index, 'ov_min'] = covered(segs, ref['ge_min'] if ref else EMPTY)
        calls.loc[g.index, 'ov_10kb'] = covered(segs, ref['ge_10kb'] if ref else EMPTY)
        for t in TLAB: calls.loc[g.index, t] = covered(segs, ref[t] if ref else EMPTY)
    tr = truth[truth.ind.isin(inds) & (truth.len >= MIN_T)].copy(); tr['ov'] = 0
    merged_calls = {k: merge_iv(g[['start', 'end']].to_numpy()) for k, g in calls.groupby(['ind', 'chr'])}
    for (s, c), g in tr.groupby(['ind', 'chr']):
        tr.loc[g.index, 'ov'] = covered(g[['start', 'end']].to_numpy(), merged_calls.get((s, c), EMPTY))
    eff = tr[tr.tmrca <= MAX_T]; n = len(inds); key = dict(cohort=S.params.cohort, method=method, input=inp, depth=depth)
    for _, _, lab in CLASSES + [(0, 0, 'all >=0.1 Mb')]:
        if lab == 'all >=0.1 Mb':
            te, ta, cc = eff, tr, calls[calls.len >= 1e5]
        else:
            te, ta, cc = eff[eff.size_class == lab], tr[tr.size_class == lab], calls[calls.size_class == lab]
        cl = cc.len.sum()
        summary.append(dict(**key, size_class=lab, n_individuals=n,
                            truth_mb_per_ind=te.len.sum() / n / 1e6, recall=te.ov.sum() / te.len.sum() if te.len.sum() else np.nan,
                            truth_any_age_mb_per_ind=ta.len.sum() / n / 1e6, recall_any_age=ta.ov.sum() / ta.len.sum() if ta.len.sum() else np.nan,
                            n_calls=len(cc), called_mb_per_ind=cl / n / 1e6,
                            precision=cc.ov_min.sum() / cl if cl else np.nan, precision_10kb=cc.ov_10kb.sum() / cl if cl else np.nan,
                            **{f'frac_{t}': cc[t].sum() / cl if cl else np.nan for t in TLAB},
                            frac_no_tract=(cl - cc[TLAB].sum().sum()) / cl if cl else np.nan))
    for s in inds:
        for _, _, lab in CLASSES:
            te = eff[(eff.ind == s) & (eff.size_class == lab)]; cc = calls[(calls.ind == s) & (calls.size_class == lab)]
            per_ind.append(dict(**key, sample=s, size_class=lab, truth_bp=te.len.sum(), truth_covered_bp=te.ov.sum(),
                                called_bp=cc.len.sum(), called_in_truth_bp=cc.ov_min.sum()))
    g = tr.groupby(['size_class', 'tclass']).agg(truth_bp=('len', 'sum'), covered_bp=('ov', 'sum')).reset_index()
    for _, r in g.iterrows():
        by_tmrca.append(dict(**key, size_class=r.size_class, tmrca_class=r.tclass, truth_mb_per_ind=r.truth_bp / n / 1e6,
                             recall=r.covered_bp / r.truth_bp if r.truth_bp else np.nan))
pd.DataFrame(summary).round(4).to_csv(S.output.by_size, index=False)
pd.DataFrame(by_tmrca).round(4).to_csv(S.output.by_tmrca, index=False)
pd.DataFrame(per_ind).to_csv(S.output.per_individual, index=False)
print(f"{len(arms)} call sets scored; {len(truth)} truth tracts", file=sys.stderr)
