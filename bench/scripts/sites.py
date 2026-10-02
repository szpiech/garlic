"""Shared ascertained site list (founder MAF >= threshold) in the forms ANGSD, bcftools and GARLIC need,
plus the mapped-region BED used to restrict every input and call."""
import sys
sys.path.insert(0, snakemake.params.scriptdir)
from common import mapped_ranges
cfg = snakemake.config
CHRS = {k: int(v) for k, v in snakemake.params.chrs.items()}
maf = float(cfg['ascertainment_maf'])
rng = mapped_ranges(cfg['tools']['simmap'], CHRS)
o = snakemake.output
n = nm = 0
with open(snakemake.input.sites) as fi, open(o.sites, 'w') as fo, open(o.bed, 'w') as fb, open(o.alleles, 'w') as fa, open(o.mapped_bed, 'w') as fmb:
    for line in fi:
        c, p, r, a, f = line.split(); f = float(f)
        if min(f, 1 - f) < maf: continue
        fo.write(f"{c}\t{p}\t{r}\t{a}\n"); fb.write(f"{c}\t{int(p)-1}\t{p}\n"); fa.write(f"{c}\t{p}\t{r},{a}\n"); n += 1
        lo, hi = rng[c]
        if lo <= int(p) <= hi: fmb.write(f"{c}\t{int(p)-1}\t{p}\n"); nm += 1
with open(o.regions, 'w') as fh:
    for c in CHRS: fh.write(c + "\n")
with open(o.mapped, 'w') as fh:
    for c, (lo, hi) in rng.items(): fh.write(f"{c}\t{lo-1}\t{hi}\n")
print(f"{n} ascertained sites, {nm} inside the mapped region", file=sys.stderr)
