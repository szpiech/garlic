"""Restrict an ANGSD beagle file to the mapped region (marker = chr_pos, split at the last underscore)."""
import gzip
rng = {}
for line in open(snakemake.input.mapped):
    c, a, b = line.split(); rng[c] = (int(a) + 1, int(b))
with gzip.open(snakemake.input.beagle, 'rt') as fi, gzip.open(snakemake.output[0], 'wt', compresslevel=4) as fo:
    fo.write(next(fi))
    for line in fi:
        c, pos = line.split('\t', 1)[0].rsplit('_', 1)
        lo, hi = rng.get(c, (1, 0))
        if lo <= int(pos) <= hi: fo.write(line)
