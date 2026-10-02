"""Truth for a population cohort: within-individual IBD segments of the simulated tree sequence whose
coalescence time is at most `ibd_max_time` generations ago -- runs of homozygosity by descent from the
population's own demography, with no pedigree.  Output: truth_hbd.bed (chr start end sample), 1-based
inclusive like ped-sim's .seg, clipped to the chromosome."""
import sys, tskit
cfg = snakemake.config
CHRS = {k: int(v) for k, v in snakemake.params.chrs.items()}
T = float(cfg['population_cohort']['ibd_max_time'])
min_span = float(cfg['population_cohort'].get('ibd_min_span', 100000))
samples = [l.strip() for l in open(snakemake.input.samples) if l.strip()]
# founder sample F%03d is individual index int(name[1:]) -> nodes 2i, 2i+1
n = 0
with open(snakemake.output[0], 'w') as fo:
    for c in CHRS:
        ts = tskit.load(f"{snakemake.params.treedir}/chr{c}.trees")
        for s in samples:
            i = int(s[1:]); a, b = ts.individual(i).nodes
            res = ts.ibd_segments(within=[a, b], max_time=T, min_span=min_span, store_segments=True)
            segs = []
            for pair, lst in res.items():
                for seg in lst: segs.append((int(seg.left) + 1, min(int(seg.right), CHRS[c])))
            for l, r in sorted(segs):
                if r > l: fo.write(f"{c}\t{l}\t{r}\t{s}\n"); n += 1
print(f"{n} within-individual IBD segments (coalescence <= {T:g} generations, span >= {min_span:g})", file=sys.stderr)
