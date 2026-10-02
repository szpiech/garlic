"""ped-sim def file, and the map / interference files restricted to the simulated chromosomes
(ped-sim requires the interference file to list exactly the map's chromosomes)."""
cfg = snakemake.config
CHRS = set(snakemake.params.chrs)
n = int(cfg['n_per_pedigree'])
with open(snakemake.output.deffile, 'w') as fh:
    for name, lines in cfg['pedigrees'].items():
        fh.write(f"def {name} {n} {len(lines)}\n")
        for l in lines: fh.write(l + "\n")
with open(cfg['tools']['simmap']) as fi, open(snakemake.output.simmap, 'w') as fo:
    fo.write(next(fi))
    for line in fi:
        if line.split()[0] in CHRS: fo.write(line)
with open(cfg['tools']['intf']) as fi, open(snakemake.output.intf, 'w') as fo:
    for line in fi:
        if line.split()[0] in CHRS: fo.write(line)
