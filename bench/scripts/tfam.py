"""TFAM (for --gl-beagle) and --pop file (for --vcf) with every sample in one population."""
names = [l.strip() for l in open(snakemake.input[0]) if l.strip()]
with open(snakemake.output.tfam, 'w') as fh:
    for s in names: fh.write(f"bench\t{s}\t0\t0\t0\t-9\n")
with open(snakemake.output.popfile, 'w') as fh:
    for s in names: fh.write(f"{s}\tbench\n")
