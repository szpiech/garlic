"""HBD rows of ped-sim's --seg_haps .seg file -> truth_hbd.bed (chr start end sample), clipped to the
simulated chromosome lengths."""
import collections
CHRS = {k: int(v) for k, v in snakemake.params.chrs.items()}
hbd = collections.defaultdict(list)
with open(snakemake.input.seg) as fh:
    for line in fh:
        f = line.split()
        # s1 h1 s2 h2 chr start end type gstart gend glen
        if len(f) >= 8 and f[7] == 'HBD' and f[4] in CHRS:
            hbd[f[0]].append((f[4], int(f[5]), min(int(f[6]), CHRS[f[4]])))
with open(snakemake.output[0], 'w') as fo:
    for s, segs in hbd.items():
        for c, a, b in sorted(segs): fo.write(f"{c}\t{a}\t{b}\t{s}\n")
