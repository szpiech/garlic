"""Any-age truth for scoring by size class: every within-individual IBD tract >= min_span with its TMRCA.
Output truth/ibd_all.bed.gz: chr start end sample tmrca span (1-based inclusive), one row per tract
(a tract is a stretch with a single most recent common ancestor of the individual's two haplotypes).

population cohort: tskit ibd_segments within the individual's two sample nodes of the founder trees.
pedigree cohort: ped-sim's --bp file gives, for each haplotype, the founder haplotype it copies along the
  chromosome.  Where both haplotypes copy the same founder haplotype the tract is pedigree HBD (TMRCA = the
  number of generations back to that founder); elsewhere the two founder haplotypes are IBD wherever the
  founder trees put them under one MRCA (TMRCA = node time + generations back to the founders).
  Founder haplotype k is haplotype k % 2 of the founder on line k // 2 of the --founder_ids file; this is
  checked against the genotypes (each simulated haplotype must carry its founder's alleles), and the HBD
  part is checked against ped-sim's own HBD segments."""
import sys, gzip, re, collections
import numpy as np, tskit

S = snakemake
CHRS = {k: int(v) for k, v in S.params.chrs.items()}
MIN = float(S.params.min_span)
samples = [l.strip() for l in open(S.input.samples) if l.strip()]
trees = {c: tskit.load(f) for c, f in zip(CHRS, S.input.trees)}
log = open(S.log[0], 'w')
out = gzip.open(S.output[0], 'wt'); n_out = 0

def emit(c, a, b, s, t):
    global n_out
    if b - a + 1 >= MIN: out.write(f"{c}\t{a}\t{b}\t{s}\t{t:.1f}\t{b - a + 1}\n"); n_out += 1

def gen(name):
    m = re.search(r'_g(\d+)-', name); return int(m.group(1)) if m else 1

if S.params.cohort != 'pedigree':
    for c, ts in trees.items():
        for s in samples:
            a, b = ts.individual(int(s[1:])).nodes
            res = ts.ibd_segments(within=[a, b], min_span=MIN, store_segments=True)
            for lst in res.values():
                for seg in lst:
                    emit(c, int(seg.left) + 1, min(int(seg.right), CHRS[c]), s, ts.node(seg.node).time)
else:
    ids = [l.split() for l in open(S.input.ids) if l.strip()]            # ped-sim founder id, input VCF sample
    fhap = lambda k: (ids[k // 2][1], k % 2, gen(ids[k // 2][0]))       # (input sample F###, haplotype, generation)
    # --bp: <sample> s<sex> h<hap> <chr>|<start> <founder_hap>:<end> ... per chromosome
    bp = collections.defaultdict(dict)
    for line in open(S.input.bp):
        f = line.split(); s, h = f[0], int(f[2][1]); cur = None; segs = {}
        for tok in f[3:]:
            if '|' in tok: cur, st = tok.split('|'); st = int(st); segs[cur] = []
            else:
                k, e = map(int, tok.split(':'))
                if cur in CHRS and st <= CHRS[cur]: segs[cur].append((st, min(e, CHRS[cur]), k))   # clip to the simulated length
                st = e + 1
        bp[s][h] = segs
    missing = [s for s in samples if s not in bp]
    assert not missing, f"no --bp rows for {missing[:3]}"

    # check the founder-haplotype mapping against genotypes at every 25th founder site
    fvcf = {}; isamp = {}
    with gzip.open(S.input.fvcf, 'rt') as fh:
        for line in fh:
            if line.startswith('##'): continue
            if line.startswith('#'): cols = line.rstrip().split('\t'); fidx = {n: i for i, n in enumerate(cols)}; continue
            break
    with gzip.open(S.input.fvcf, 'rt') as fh:
        i = 0
        for line in fh:
            if line[0] == '#': continue
            i += 1
            if i % 25: continue
            p = line.rstrip().split('\t'); fvcf[(p[0], int(p[1]))] = p
    agree = tot = 0
    with gzip.open(S.input.ivcf, 'rt') as fh:
        for line in fh:
            if line.startswith('##'): continue
            if line.startswith('#'): icols = line.rstrip().split('\t'); continue
            p = line.rstrip().split('\t'); key = (p[0], int(p[1]))
            if key not in fvcf: continue
            fp = fvcf[key]
            for j, s in enumerate(icols[9:]):
                if s not in samples: continue
                gt = p[9 + j].split('|')
                for h in (0, 1):
                    for st, e, k in bp[s][h].get(p[0], []):
                        if st <= key[1] <= e:
                            fs, fh_, _ = fhap(k); fgt = fp[fidx[fs]].split('|')
                            agree += fgt[fh_] == gt[h]; tot += 1; break
    conc = agree / max(tot, 1)
    print(f"founder-haplotype mapping: {agree}/{tot} alleles agree ({conc:.5f})", file=log, flush=True)
    assert tot > 1000 and conc > 0.999, f"founder-haplotype mapping check failed: {conc:.4f} over {tot} alleles"

    hbd_bp = 0; pair_cache = {}
    for c, ts in trees.items():
        nodes = {}
        for s in samples:
            g_ind = gen(s); h0 = bp[s][0].get(c, []); h1 = bp[s][1].get(c, [])
            # intersect the two haplotypes' founder segments
            i = j = 0
            while i < len(h0) and j < len(h1):
                a = max(h0[i][0], h1[j][0]); b = min(h0[i][1], h1[j][1]); ka, kb = h0[i][2], h1[j][2]
                if b >= a:
                    fa, fb = fhap(ka), fhap(kb)
                    if ka == kb:
                        emit(c, a, b, s, g_ind - fa[2]); hbd_bp += b - a + 1
                    else:
                        na = ts.individual(int(fa[0][1:])).nodes[fa[1]]; nb = ts.individual(int(fb[0][1:])).nodes[fb[1]]
                        key = (c, min(na, nb), max(na, nb))
                        if key not in pair_cache:
                            res = ts.ibd_segments(between=[[na], [nb]], store_segments=True)
                            pair_cache[key] = sorted((int(sg.left) + 1, int(sg.right), ts.node(sg.node).time) for lst in res.values() for sg in lst)
                        back = g_ind - (fa[2] + fb[2]) / 2
                        for l, r, t in pair_cache[key]:
                            l2, r2 = max(l, a), min(r, b)
                            if r2 >= l2: emit(c, l2, min(r2, CHRS[c]), s, t + back)
                if h0[i][1] < h1[j][1]: i += 1
                else: j += 1
    # truth_hbd.bed rows are clipped at the chromosome end, so a segment that starts beyond it has end < start
    ped_hbd = sum(max(0, int(l.split()[2]) - int(l.split()[1]) + 1) for l in open(S.input.hbd) if l.strip())
    print(f"HBD from --bp: {hbd_bp} bp; ped-sim HBD segments: {ped_hbd} bp", file=log, flush=True)
    assert abs(hbd_bp - ped_hbd) <= 0.01 * max(ped_hbd, 1), "HBD from --bp disagrees with ped-sim's HBD segments"
out.close()
print(f"{n_out} within-individual IBD tracts (span >= {MIN:g}, any age)", file=log)
