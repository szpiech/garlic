"""msprime founder haplotypes -> founders.vcf.gz (phased, VCF POS = 1-based), ref.fa (random sequence
carrying the ancestral allele at every variant), sites.tsv (chr pos ref alt founder_alt_frequency)."""
import sys, gzip
sys.path.insert(0, snakemake.params.scriptdir)
import numpy as np, msprime
from common import read_simmap

cfg = snakemake.config
CHRS = {k: int(v) for k, v in snakemake.params.chrs.items()}

def demography(F):
    """One population; config founders.demography selects its history (times in generations ago).
    constant:   Ne throughout.
    bottleneck: Ne_now until t_start, Ne_bottleneck from t_start to t_end, Ne before t_end.
    small:      a small constant population (Ne_small)."""
    model = F.get('demography', 'constant')
    dem = msprime.Demography()
    if model == 'constant':
        dem.add_population(name='P', initial_size=int(F['Ne']))
    elif model == 'bottleneck':
        dem.add_population(name='P', initial_size=int(F['Ne']))
        dem.add_population_parameters_change(time=float(F['bottleneck_start']), initial_size=int(F['bottleneck_Ne']), population='P')
        dem.add_population_parameters_change(time=float(F['bottleneck_end']), initial_size=int(F['Ne']), population='P')
    elif model == 'small':
        dem.add_population(name='P', initial_size=int(F['small_Ne']))
    else:
        raise ValueError(f"unknown demography {model}")
    return dem
F = cfg['founders']; SEED = int(F['seed']); NHAP = int(F['n_haplotypes'])
rng = np.random.default_rng(SEED)
maps = read_simmap(cfg['tools']['simmap'], CHRS)
names = [f"F{i:03d}" for i in range(NHAP // 2)]
bases = np.array(list("ACGT"))
with open(snakemake.output.ref, 'w') as ref_fa, gzip.open(snakemake.output.vcf, 'wt') as vcf, open(snakemake.output.sites, 'w') as sites:
    vcf.write("##fileformat=VCFv4.2\n")
    for c, L in CHRS.items(): vcf.write(f"##contig=<ID={c},length={L}>\n")
    vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
    for c, L in CHRS.items():
        m = maps[c]; m = m[m[:, 0] < L]
        assert len(m) >= 2, f"chr{c}: genetic map has no positions below {L}"
        pos = np.concatenate([[0], m[:, 0], [L]]); cm = np.concatenate([[0], m[:, 1], [m[-1, 1]]])
        pos, idx = np.unique(pos, return_index=True); cm = cm[idx]
        rate = np.diff(cm) / 100.0 / np.diff(pos); rate = np.where(np.isfinite(rate) & (rate >= 0), rate, 0.0)
        ts = msprime.sim_ancestry(samples=NHAP // 2, demography=demography(F), recombination_rate=msprime.RateMap(position=pos, rate=rate),
                                  sequence_length=L, random_seed=SEED + int(c), ploidy=2)
        ts = msprime.sim_mutations(ts, rate=float(F['mu']), random_seed=SEED + 100 + int(c), model=msprime.BinaryMutationModel())
        ts.dump(f"{snakemake.params.treedir}/chr{c}.trees")
        seq = rng.choice(bases, size=L)
        ref_fa.write(f">{c}\n")
        for k in range(0, L, 60): ref_fa.write("".join(seq[k:k + 60]) + "\n")
        G = ts.genotype_matrix(); positions = ts.sites_position.astype(int)      # 0-based
        keep = np.ones(len(positions), bool); keep[1:] &= np.diff(positions) > 0; keep &= positions < L - 1
        G = G[keep]; positions = positions[keep]
        ref_alleles = seq[positions]
        alt_alleles = np.array([rng.choice([b for b in "ACGT" if b != r]) for r in ref_alleles])
        n = 0
        for i, p0 in enumerate(positions):
            g = G[i]; p = p0 + 1
            if g.sum() == 0 or g.sum() == len(g): continue
            gt = "\t".join(f"{g[2*j]}|{g[2*j+1]}" for j in range(len(g) // 2))
            vcf.write(f"{c}\t{p}\t{c}_{p}\t{ref_alleles[i]}\t{alt_alleles[i]}\t.\tPASS\t.\tGT\t{gt}\n")
            sites.write(f"{c}\t{p}\t{ref_alleles[i]}\t{alt_alleles[i]}\t{g.mean():.5f}\n"); n += 1
        print(f"chr{c}: {n} segregating sites", file=sys.stderr)
