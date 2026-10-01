#!/usr/bin/env python3
"""GARLIC low-coverage benchmark: simulated ROH truth, Illumina reads, three
GARLIC inputs and ROHan, scored against ped-sim HBD segments.

Stages (each idempotent, re-runnable):
  founders   msprime founder haplotypes on chr20-22 (b37 lengths, Bherer sex-averaged map)
             -> founders.vcf.gz, ref.fa (ancestral allele at every site), sites.tsv
  pedsim     ped-sim pedigrees with planted inbreeding -> truth HBD segments, true genotypes VCF
  reads      per individual: two haplotype FASTAs, ART paired-end reads at each depth
  map        bwa-mem2 -> sorted BAM
  gl         ANGSD GL (beagle) + MAF at the shared site list; bcftools calls at the same sites
  roh        GARLIC x {true genotypes, bcftools calls, calls+GQ, beagle}; ROHan on BAMs
  score      base-level precision/recall vs truth, FROH error -> bench_scores.csv

All paths under /data/workspace/garlic-lowcov/bench by default.
"""
import argparse, gzip, os, subprocess, sys, json, random, shutil, time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = Path(os.environ.get('BENCH_ROOT', '/data/workspace/garlic-lowcov/bench'))
P    = Path('/data/workspace/garlic-lowcov')
GARLIC = P / 'bin/garlic'
ROHAN  = P / 'tools/rohan/bin/rohan'
PEDSIM = Path('/data/workspace/rohml/bin/ped-sim')
SIMMAP = Path('/data/workspace/rohml/maps/refined_mf_X.simmap')
INTF   = Path('/data/workspace/rohml/tools/ped-sim/interfere/nu_p_campbell_X.tsv')

CHRS = {'20': 63025520, '21': 48129895, '22': 51304566}   # GRCh37
DEPTHS = [1, 2, 4, 8]
ROHAN_DEPTHS = [1, 4]
SMALL = os.environ.get('BENCH_SMALL') == '1'        # smoke test: 24 Mb of chr21/22, one pedigree copy, two depths
if SMALL:
    CHRS = {'21': 24000000, '22': 24000000}; DEPTHS = [1, 4]   # the b37 map starts at 9.4 / 16.1 Mb on these
NE = 20000; MU = 1.25e-8; SEED = 2026
N_FOUNDERS_HAP = 400          # founder haplotypes to simulate (200 diploids; ped-sim draws from these)
READ_LEN = 150; FRAG = 350; FRAG_SD = 40
ERR = 0.001                    # garlic --error

def sh(cmd, log=None, check=True, env=None):
    t0 = time.time()
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)
    if log:
        with open(log, 'a') as fh:
            fh.write(f"$ {cmd}\n{r.stdout}{r.stderr}[exit {r.returncode}, {time.time()-t0:.0f}s]\n")
    if check and r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}\n{r.stderr[-2000:]}")
    return r

def done(flag): return (ROOT / 'flags' / flag).exists()
def mark(flag): (ROOT / 'flags').mkdir(parents=True, exist_ok=True); (ROOT / 'flags' / flag).touch()

# ----------------------------------------------------------------------------
def read_simmap():
    """Sex-averaged cM per chromosome from the Bherer simmap (chr, pos, male_cM, female_cM)."""
    import numpy as np
    out = {}
    with open(SIMMAP) as fh:
        next(fh)
        for line in fh:
            c, pos, m, f = line.split()[:4]
            if c not in CHRS: continue
            out.setdefault(c, []).append((int(pos), (float(m) + float(f)) / 2))
    return {c: np.array(v) for c, v in out.items()}

def stage_founders():
    if done('founders'): return
    import msprime, numpy as np
    rng = np.random.default_rng(SEED)
    maps = read_simmap()
    (ROOT / 'founders').mkdir(parents=True, exist_ok=True)
    ref_fa = open(ROOT / 'founders/ref.fa', 'w')
    vcf = gzip.open(ROOT / 'founders/founders.vcf.gz', 'wt')
    sites = open(ROOT / 'founders/sites.tsv', 'w')
    sample_names = [f"F{i:03d}" for i in range(N_FOUNDERS_HAP // 2)]
    vcf.write("##fileformat=VCFv4.2\n")
    for c, L in CHRS.items(): vcf.write(f"##contig=<ID={c},length={L}>\n")
    vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")
    bases = np.array(list("ACGT"))
    nsites_total = 0
    for c, L in CHRS.items():
        m = maps[c]; m = m[m[:, 0] < L]
        assert len(m) >= 2, f'chr{c}: genetic map has no positions below {L}'
        pos = np.concatenate([[0], m[:, 0], [L]]); cm = np.concatenate([[0], m[:, 1], [m[-1, 1]]])
        pos, idx = np.unique(pos, return_index=True); cm = cm[idx]
        rate = np.diff(cm) / 100.0 / np.diff(pos)
        rate = np.where(np.isfinite(rate) & (rate >= 0), rate, 0.0)
        ratemap = msprime.RateMap(position=pos, rate=rate)
        ts = msprime.sim_ancestry(samples=N_FOUNDERS_HAP // 2, population_size=NE, recombination_rate=ratemap,
                                  sequence_length=L, random_seed=SEED + int(c), ploidy=2)
        ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED + 100 + int(c), model=msprime.BinaryMutationModel())
        # Reference: random sequence, ancestral allele at every variant.  Written in 60-col lines.
        seq = rng.choice(bases, size=L)
        ref_fa.write(f">{c}\n")
        G = ts.genotype_matrix()                      # sites x haplotypes (0/1)
        positions = ts.sites_position.astype(int)          # 0-based; VCF POS = positions + 1
        # one variant per position: keep the first at duplicate positions
        keep = np.ones(len(positions), bool)
        keep[1:] &= np.diff(positions) > 0
        keep &= positions < L - 1
        G = G[keep]; positions = positions[keep]
        ref_alleles = seq[positions]
        alt_alleles = np.array([rng.choice([b for b in "ACGT" if b != r]) for r in ref_alleles])
        for k in range(0, L, 60): ref_fa.write("".join(seq[k:k + 60]) + "\n")
        for i, p0 in enumerate(positions):
            g = G[i]; p = p0 + 1
            if g.sum() == 0 or g.sum() == len(g): continue
            gt = "\t".join(f"{g[2*j]}|{g[2*j+1]}" for j in range(len(g) // 2))
            vcf.write(f"{c}\t{p}\t{c}_{p}\t{ref_alleles[i]}\t{alt_alleles[i]}\t.\tPASS\t.\tGT\t{gt}\n")
            sites.write(f"{c}\t{p}\t{ref_alleles[i]}\t{alt_alleles[i]}\t{g.mean():.5f}\n")
            nsites_total += 1
        print(f"[founders] chr{c}: {len(positions)} sites", flush=True)
    ref_fa.close(); vcf.close(); sites.close()
    # ANGSD refuses a .fai that is not strictly newer than the FASTA
    sh(f"cd {ROOT}/founders && sleep 2 && samtools faidx ref.fa && touch ref.fa.fai && bwa-mem2 index ref.fa > bwa_index.log 2>&1", log=ROOT / 'founders/cmd.log')
    print(f"[founders] {nsites_total} segregating sites written", flush=True)
    mark('founders')

# ----------------------------------------------------------------------------
PEDS = {  # name: (generations, def lines) -- offspring of full sibs, first cousins, second cousins, outbred
    # generation lines: <gen> <#print> <#branches> [branch specs]; "1:1_2" = branch 1's parents are branches 1 and 2 of the previous generation
    'fs': (3, ["1 0 1", "2 0 2", "3 1 1 1:1_2"]),
    'c1': (4, ["1 0 1", "2 0 2", "3 0 2", "4 1 1 1:1_2"]),
    'c2': (5, ["1 0 1", "2 0 2", "3 0 2", "4 0 2", "5 1 1 1:1_2"]),
    'ob': (2, ["1 0 1", "2 1 1"]),
}
N_PER_PED = 1 if SMALL else 5

def stage_pedsim():
    if done('pedsim'): return
    d = ROOT / 'pedsim'; d.mkdir(parents=True, exist_ok=True)
    with open(d / 'bench.def', 'w') as fh:
        for name, (gens, lines) in PEDS.items():
            fh.write(f"def {name} {N_PER_PED} {gens}\n")
            for l in lines: fh.write(l + "\n")
    # map restricted to the simulated chromosomes; male cM on X irrelevant here
    with open(SIMMAP) as fi, open(d / 'map.simmap', 'w') as fo:
        fo.write(next(fi))
        for line in fi:
            if line.split()[0] in CHRS: fo.write(line)
    # interference file must list exactly the chromosomes in the map
    with open(INTF) as fi, open(d / 'intf.tsv', 'w') as fo:
        for line in fi:
            if line.split()[0] in CHRS: fo.write(line)
    sh(f"cd {d} && {PEDSIM} -d bench.def -m map.simmap -i {ROOT}/founders/founders.vcf.gz -o bench --intf intf.tsv --seg_haps --keep_phase --founder_ids --err_rate 0 --miss_rate 0 > pedsim.log 2>&1",
       log=d / 'cmd.log')
    # ped-sim writes plain gzip (or plain text); re-compress with bgzip so it can be indexed
    if (d / 'bench.vcf').exists():
        sh(f"cd {d} && bgzip -f bench.vcf", log=d / 'cmd.log')
    else:
        sh(f"cd {d} && zcat bench.vcf.gz | bgzip -c > bench.bgz.vcf.gz && mv bench.bgz.vcf.gz bench.vcf.gz", log=d / 'cmd.log')
    sh(f"cd {d} && bcftools index -f -t bench.vcf.gz", log=d / 'cmd.log')
    # sample list of the printed (final-generation) individuals
    sh(f"cd {d} && bcftools query -l bench.vcf.gz > samples.txt", log=d / 'cmd.log')
    names = [l.strip() for l in open(d / 'samples.txt')]
    print(f"[pedsim] {len(names)} samples: {names[:6]} ...", flush=True)
    # truth: HBD rows of the .seg file, per sample
    hbd = {}
    nrow = 0
    with open(d / 'bench.seg') as fh:
        for line in fh:
            f = line.split(); nrow += 1
            if nrow <= 3: print("[pedsim] seg row:", line.rstrip(), flush=True)
            # --seg_haps format: s1 h1 s2 h2 chr start end type gstart gend glen
            if len(f) >= 8 and f[7] == 'HBD':
                hbd.setdefault(f[0], []).append((f[4], int(f[5]), min(int(f[6]), CHRS[f[4]])))
    if not hbd:  # fall back to "same sample, both haplotypes"
        with open(d / 'bench.seg') as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 8 and f[0] == f[2]:
                    hbd.setdefault(f[0], []).append((f[4], int(f[5]), int(f[6])))
    with open(d / 'truth_hbd.bed', 'w') as fo:
        for s, segs in hbd.items():
            for c, a, b in sorted(segs): fo.write(f"{c}\t{a}\t{b}\t{s}\n")
    print(f"[pedsim] truth HBD segments: {sum(len(v) for v in hbd.values())} over {len(hbd)} samples", flush=True)
    mark('pedsim')

# ----------------------------------------------------------------------------
def samples():
    return [l.strip() for l in open(ROOT / 'pedsim/samples.txt')]

def _hap_fastas(s):
    """Two haplotype FASTAs for sample s from the ped-sim VCF (phased) and the reference."""
    out = ROOT / f'reads/{s}'; out.mkdir(parents=True, exist_ok=True)
    if (out / 'hap1.fa').exists() and (out / 'hap2.fa').exists(): return
    sh(f"cd {out} && bcftools view -s {s} -Oz -o {s}.vcf.gz {ROOT}/pedsim/bench.vcf.gz && bcftools index -f {s}.vcf.gz && "
       f"bcftools consensus -f {ROOT}/founders/ref.fa -H 1 {s}.vcf.gz > hap1.fa 2>cons1.log && "
       f"bcftools consensus -f {ROOT}/founders/ref.fa -H 2 {s}.vcf.gz > hap2.fa 2>cons2.log", log=out / 'cmd.log')

def _reads(args):
    s, depth = args
    out = ROOT / f'reads/{s}/d{depth}'; out.mkdir(parents=True, exist_ok=True)
    if (out / 'done').exists(): return s, depth
    # half the depth from each haplotype; ART writes paired FASTQ
    half = depth / 2.0
    for h in (1, 2):
        sh(f"cd {out} && art_illumina -ss HS25 -i ../hap{h}.fa -p -l {READ_LEN} -f {half} -m {FRAG} -s {FRAG_SD} "
           f"-na -q -rs {SEED + depth * 7 + h} -o h{h}_ > art{h}.log 2>&1", log=out / 'cmd.log')
    sh(f"cd {out} && cat h1_1.fq h2_1.fq | gzip -1 > r1.fq.gz && cat h1_2.fq h2_2.fq | gzip -1 > r2.fq.gz && rm -f h?_?.fq", log=out / 'cmd.log')
    (out / 'done').touch()
    return s, depth

def stage_reads(threads):
    if done('reads'): return
    ss = samples()
    with ProcessPoolExecutor(min(threads, len(ss))) as ex:
        list(ex.map(_hap_fastas, ss))
    jobs = [(s, d) for s in ss for d in DEPTHS]
    with ProcessPoolExecutor(threads) as ex:
        for s, d in ex.map(_reads, jobs): print(f"[reads] {s} {d}x", flush=True)
    mark('reads')

def _map(args):
    s, depth, t = args
    out = ROOT / f'reads/{s}/d{depth}'
    if (out / 'aln.bam.bai').exists(): return s, depth
    sh(f"cd {out} && bwa-mem2 mem -t {t} -R '@RG\\tID:{s}_d{depth}\\tSM:{s}' {ROOT}/founders/ref.fa r1.fq.gz r2.fq.gz 2>bwa.log "
       f"| samtools sort -@2 -o aln.bam - 2>sort.log && samtools index aln.bam", log=out / 'cmd.log')
    return s, depth

def stage_map(threads):
    if done('map'): return
    ss = samples()
    jobs = [(s, d, 4) for s in ss for d in DEPTHS]
    with ProcessPoolExecutor(max(1, threads // 4)) as ex:
        for s, d in ex.map(_map, jobs): print(f"[map] {s} {d}x", flush=True)
    mark('map')

# ----------------------------------------------------------------------------
def stage_gl(threads):
    """ANGSD genotype likelihoods at the shared site list, one run per depth over all samples;
    bcftools calls at the same sites."""
    if done('gl'): return
    g = ROOT / 'gl'; g.mkdir(parents=True, exist_ok=True)
    ss = samples()
    # Shared ascertainment: founder sites with MAF >= 0.05 in the founder sample.
    with open(ROOT / 'founders/sites.tsv') as fi, open(g / 'sites.txt', 'w') as fo, open(g / 'sites.bed', 'w') as fb, open(g / 'alleles.tsv', 'w') as fa:
        n = 0
        for line in fi:
            c, p, r, a, f = line.split(); f = float(f)
            if min(f, 1 - f) >= 0.05:
                fo.write(f"{c}\t{p}\t{r}\t{a}\n"); fb.write(f"{c}\t{int(p)-1}\t{p}\n"); fa.write(f"{c}\t{p}\t{r},{a}\n"); n += 1
    print(f"[gl] {n} ascertained sites", flush=True)
    with open(g / 'regions.txt', 'w') as fh:
        for c in CHRS: fh.write(c + "\n")
    # forced-allele file for bcftools call -C alleles: chr pos ref,alt
    sh(f"cd {g} && bgzip -f alleles.tsv && tabix -f -s1 -b2 -e2 alleles.tsv.gz", log=g / 'cmd.log')
    sh(f"cd {g} && angsd sites index sites.txt", log=g / 'cmd.log')
    for d in DEPTHS:
        with open(g / f'bams_d{d}.txt', 'w') as fh:
            for s in ss: fh.write(f"{ROOT}/reads/{s}/d{d}/aln.bam\n")
    # ANGSD (8 threads each) and bcftools (1 thread each) for every depth run concurrently
    jobs = [(d, 'angsd') for d in DEPTHS] + [(d, 'bcftools') for d in DEPTHS]
    with ProcessPoolExecutor(len(jobs)) as ex:
        for d, tool in ex.map(_gl_one, jobs): print(f"[gl] {tool} depth {d}x done", flush=True)
    mark('gl')

def _gl_one(args):
    d, tool = args
    g = ROOT / 'gl'
    if tool == 'angsd':
        if not (g / f'd{d}.beagle.gz').exists():
            sh(f"cd {g} && angsd -b bams_d{d}.txt -ref {ROOT}/founders/ref.fa -out d{d} -sites sites.txt -rf regions.txt "
               f"-GL 1 -doGlf 2 -doMajorMinor 3 -doMaf 1 -doCounts 1 -dumpCounts 4 -minMapQ 20 -minQ 20 -baq 1 -C 50 -P 8 > angsd_d{d}.log 2>&1",
               log=g / f'cmd_d{d}.log')
    else:
        if not (g / f'd{d}.calls.vcf.gz').exists():
            sh(f"cd {g} && bcftools mpileup -f {ROOT}/founders/ref.fa -T sites.bed -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 -b bams_d{d}.txt -Ou 2>mp_d{d}.log "
               f"| bcftools call -m -C alleles -T alleles.tsv.gz -Oz -o d{d}.calls.vcf.gz 2>call_d{d}.log && bcftools index -f d{d}.calls.vcf.gz", log=g / f'cmd_d{d}.log')
    return d, tool

# ----------------------------------------------------------------------------
def mapped_ranges():
    """Per-chromosome [first, last] position of the genetic map used by ped-sim.  Outside it
    ped-sim transmits nothing (truth and genotypes are reference-homozygous), so every input
    and every call is restricted to it; in the human assemblies these flanks are the
    acrocentric arms and centromeres, which real analyses mask as well."""
    rng = {}
    with open(ROOT / 'pedsim/map.simmap') as fh:
        next(fh)
        for line in fh:
            c, pos = line.split()[:2]; pos = int(pos)
            lo, hi = rng.get(c, (pos, pos)); rng[c] = (min(lo, pos), max(hi, pos))
    return {c: (lo, min(hi, CHRS[c])) for c, (lo, hi) in rng.items() if c in CHRS}

def write_mapped_inputs():
    g = ROOT / 'gl'; rng = mapped_ranges()
    with open(g / 'mapped.bed', 'w') as fh:
        for c, (lo, hi) in rng.items(): fh.write(f"{c}\t{lo-1}\t{hi}\n")
    for d in DEPTHS:
        out = g / f'd{d}.map.beagle.gz'
        if not out.exists():
            with gzip.open(g / f'd{d}.beagle.gz', 'rt') as fi, gzip.open(out, 'wt', compresslevel=4) as fo:
                fo.write(next(fi))
                for line in fi:
                    c, pos = line.split('\t', 1)[0].rsplit('_', 1); pos = int(pos)
                    lo, hi = rng[c]
                    if lo <= pos <= hi: fo.write(line)
        if not (g / f'd{d}.map.calls.vcf.gz').exists():
            sh(f"cd {g} && bcftools view -T mapped.bed -Oz -o d{d}.map.calls.vcf.gz d{d}.calls.vcf.gz && bcftools index -f d{d}.map.calls.vcf.gz", log=g / 'cmd.log')
    if not (g / 'sites.mapped.bed').exists():
        with open(g / 'sites.bed') as fi, open(g / 'sites.mapped.bed', 'w') as fo:
            for line in fi:
                c, a, b = line.split(); lo, hi = rng[c]
                if lo <= int(b) <= hi: fo.write(line)
    return rng

def write_tfam(path, names):
    with open(path, 'w') as fh:
        for s in names: fh.write(f"bench\t{s}\t0\t0\t0\t-9\n")
    with open(str(path).replace('.tfam', '.pop'), 'w') as fh:
        for s in names: fh.write(f"{s}\tbench\n")

def stage_roh(threads):
    if done('roh'): return
    r = ROOT / 'roh'; r.mkdir(parents=True, exist_ok=True)
    ss = samples()
    write_tfam(r / 'bench.tfam', ss)
    write_mapped_inputs()
    # --no-centromere: a synthetic reference has no assembled gap and no built-in table.
    common = f"--no-centromere --winsize 60 --error {ERR} --size-bounds 500000 1000000 --force --quiet --threads {min(threads, 8)}"
    vcfc = common + " --pop bench.pop"
    # 1. True genotypes (ceiling).  Restrict to the ascertained sites.
    if not (r / 'truth.roh.bed').exists():
        sh(f"cd {r} && bcftools view -T {ROOT}/gl/sites.mapped.bed -s {','.join(ss)} -Oz -o truth.vcf.gz {ROOT}/pedsim/bench.vcf.gz && "
           f"{GARLIC} --vcf truth.vcf.gz {vcfc} --out truth > truth.stdout 2> truth.stderr", log=r / 'cmd.log')
    for d in DEPTHS:
        # 2. bcftools hard calls, run-wide --error
        if not (r / f'calls_d{d}.roh.bed').exists():
            sh(f"cd {r} && {GARLIC} --vcf {ROOT}/gl/d{d}.map.calls.vcf.gz {vcfc} --out calls_d{d} > calls_d{d}.stdout 2> calls_d{d}.stderr", log=r / 'cmd.log', check=False)
        # 3. bcftools calls with per-genotype PL (bcftools writes PL, not GQ; the per-genotype-error path the manual warns about)
        if not (r / f'pl_d{d}.roh.bed').exists():
            sh(f"cd {r} && {GARLIC} --vcf {ROOT}/gl/d{d}.map.calls.vcf.gz --gl-type PL {vcfc.replace('--error '+str(ERR), '')} --out pl_d{d} > pl_d{d}.stdout 2> pl_d{d}.stderr", log=r / 'cmd.log', check=False)
        # 4. beagle likelihoods, at the shared --error and at two others (sensitivity)
        if not (r / f'beagle_d{d}.roh.bed').exists():
            sh(f"cd {r} && {GARLIC} --gl-beagle {ROOT}/gl/d{d}.map.beagle.gz --tfam bench.tfam {common} --out beagle_d{d} > beagle_d{d}.stdout 2> beagle_d{d}.stderr", log=r / 'cmd.log', check=False)
        for e, tag in ((1e-4, 'e4'), (1e-2, 'e2')):
            if not (r / f'beagle{tag}_d{d}.roh.bed').exists():
                sh(f"cd {r} && {GARLIC} --gl-beagle {ROOT}/gl/d{d}.map.beagle.gz --tfam bench.tfam {common.replace('--error '+str(ERR), '--error '+str(e))} --out beagle{tag}_d{d} > beagle{tag}_d{d}.stdout 2> beagle{tag}_d{d}.stderr", log=r / 'cmd.log', check=False)
        print(f"[roh] garlic depth {d}x done", flush=True)
    mark('roh')

def stage_rohan(threads):
    """ROHan per BAM (slow: parallel over samples). Separate stage so it can be added once the build works."""
    if done('rohan'): return
    if not ROHAN.exists():
        print("[rohan] binary missing, skipping", flush=True); return
    ss = samples()
    # ROHan's local-heterozygosity step takes ~2-3 h per 162 Mb BAM with two threads here, so it is
    # run at the two depths where the GARLIC inputs disagree most (1x) and agree (4x): 40 BAMs.
    rohan_depths = [d for d in ROHAN_DEPTHS if d in DEPTHS]
    jobs = [(s, d) for d in rohan_depths for s in ss]
    with ProcessPoolExecutor(max(1, threads // 2)) as ex:
        for s, d in ex.map(_rohan, jobs): print(f"[rohan] {s} {d}x", flush=True)
    mark('rohan')

def _rohan(args):
    s, d = args
    out = ROOT / f'roh/rohan/{s}_d{d}'; out.mkdir(parents=True, exist_ok=True)
    if (out / 'done').exists(): return s, d
    sh(f"cd {out} && {ROHAN} -t 2 --size 1000000 --seed {SEED} -o out {ROOT}/founders/ref.fa {ROOT}/reads/{s}/d{d}/aln.bam > rohan.log 2>&1", log=out / 'cmd.log', check=False)
    (out / 'done').touch()
    return s, d

# ----------------------------------------------------------------------------
def stage_score():
    import pandas as pd, numpy as np
    r = ROOT / 'roh'
    truth = pd.read_csv(ROOT / 'pedsim/truth_hbd.bed', sep='\t', names=['chr', 'start', 'end', 'sample'], dtype={'chr': str})
    ss = samples()
    rng = mapped_ranges()
    genome = sum(hi - lo + 1 for lo, hi in rng.values())      # F_ROH denominator = mapped length

    def clip(iv):
        out = []
        for c, a, b in iv:
            if c not in rng: continue
            lo, hi = rng[c]; a, b = max(a, lo), min(b, hi)
            if b > a: out.append((c, a, b))
        return out

    def intervals_from_roh_bed(fn):
        iv = {}
        cur = None
        for line in open(fn):
            if line.startswith('track'):
                import re; cur = re.search(r'Ind: (\S+)', line).group(1); iv.setdefault(cur, [])
            elif line.strip():
                f = line.split(); iv.setdefault(cur, []).append((f[0].replace('chr', ''), int(f[1]), int(f[2])))
        return iv

    def intervals_from_rohan(s, d):
        fn = ROOT / f'roh/rohan/{s}_d{d}/out.mid.hmmrohl.gz'
        if not fn.exists(): return None
        iv = []
        with gzip.open(fn, 'rt') as fh:
            for line in fh:
                f = line.split()
                if not f or f[0].startswith('#') or f[0] == 'chr' or f[0] == 'ROH_ID': continue
                try: iv.append((f[1].replace('chr', ''), int(f[2]), int(f[3])))
                except Exception: continue
        return iv

    def overlap(a, b):
        """Total overlap in bp between two interval lists (chr,start,end)."""
        tot = 0
        by = {}
        for c, s_, e in b: by.setdefault(c, []).append((s_, e))
        for c in by: by[c].sort()
        for c, s_, e in a:
            for s2, e2 in by.get(c, []):
                if e2 < s_: continue
                if s2 > e: break
                tot += max(0, min(e, e2) - max(s_, s2))
        return tot

    rows = []
    def score(method, depth, calls):
        for s in ss:
            t_all = clip([(c, a, b) for c, a, b, smp in truth[truth['sample'] == s][['chr', 'start', 'end', 'sample']].itertuples(index=False)])
            cl_all = calls.get(s, []) if calls is not None else None
            if cl_all is None: continue
            cl_all = clip(cl_all)
            for min_len in (0, 1000000):
                t = [x for x in t_all if x[2] - x[1] >= min_len]; cl = [x for x in cl_all if x[2] - x[1] >= min_len]
                tlen = sum(b - a for _, a, b in t); clen = sum(b - a for _, a, b in cl)
                ov = overlap(cl, t)
                rows.append(dict(method=method, depth=depth, sample=s, pedigree=s[:2], min_len=min_len,
                                 n_truth=len(t), n_called=len(cl), truth_bp=tlen, called_bp=clen, overlap_bp=ov,
                                 precision=ov / clen if clen else np.nan, recall=ov / tlen if tlen else np.nan,
                                 froh_truth=tlen / genome, froh_called=clen / genome))
    if (r / 'truth.roh.bed').exists(): score('garlic_true_genotypes', 'inf', intervals_from_roh_bed(r / 'truth.roh.bed'))
    for d in DEPTHS:
        for m in ['calls', 'pl', 'beagle', 'beaglee4', 'beaglee2']:
            fn = r / f'{m}_d{d}.roh.bed'
            if fn.exists(): score(f'garlic_{m}', d, intervals_from_roh_bed(fn))
        if ROHAN.exists():
            calls = {}
            for s in ss:
                iv = intervals_from_rohan(s, d)
                if iv is not None: calls[s] = iv
            if calls: score('rohan', d, calls)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / 'bench_scores.csv', index=False)
    d1 = df[df.min_len == 1000000]
    g = d1.groupby(['method', 'depth'])
    summ = pd.DataFrame({'precision_pooled': g.overlap_bp.sum() / g.called_bp.sum(), 'recall_pooled': g.overlap_bp.sum() / g.truth_bp.sum(),
                         'froh_called': g.froh_called.mean(), 'froh_truth': g.froh_truth.mean()})
    print(summ.round(3).to_string())

# ----------------------------------------------------------------------------
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stages', nargs='*', default=['founders', 'pedsim', 'reads', 'map', 'gl', 'roh', 'rohan', 'score'])
    ap.add_argument('--threads', type=int, default=32)
    a = ap.parse_args()
    ROOT.mkdir(parents=True, exist_ok=True)
    for st in a.stages:
        t0 = time.time()
        {'founders': stage_founders, 'pedsim': stage_pedsim, 'reads': lambda: stage_reads(a.threads),
         'map': lambda: stage_map(a.threads), 'gl': lambda: stage_gl(a.threads), 'roh': lambda: stage_roh(a.threads),
         'rohan': lambda: stage_rohan(a.threads), 'score': stage_score}[st]()
        print(f"== stage {st} done in {time.time()-t0:.0f}s", flush=True)
