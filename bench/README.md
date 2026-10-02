# Low-coverage ROH benchmark

Snakemake workflow used to benchmark the `low-cov` branch on simulated low-coverage
sequence data. It is written for the lab office box (tool paths in `config.yaml`) and is
kept here for provenance; it is not part of the build.

```
conda activate /data/workspace/garlic-lowcov/env-bench
snakemake -s bench/Snakefile --cores 32                                      # full benchmark
snakemake -s bench/Snakefile --cores 16 --config small=true root=<dir>       # 24 Mb smoke test (~4 min)
snakemake -s bench/Snakefile --cores 32 --config run_rohan=false             # GARLIC arms only
snakemake -s bench/Snakefile -n --rerun-triggers mtime                       # what would run
```

Pipeline (one rule per box; `pedsim` is a checkpoint because the sample names come from it):

```
founders ─┬─> index_reference ──────────────────────────┐
          └─> pedsim_inputs ─> pedsim ─┬─> truth         │
                                       ├─> tfam          │
                                       └─> hap_fasta ─> reads ─> map ─┬─> bamlist ─┬─> angsd ─> mapped_beagle ─> garlic_beagle (x3 --error)
founders ─> sites ─> index_sites ─────────────────────────────────────┘            └─> bcftools_call ─> mapped_calls ─┬─> garlic_calls
                                                                       └─> rohan                                      └─> garlic_pl
pedsim + sites ─> truth_vcf ─> garlic_truth (ceiling)
everything ─> score ─> bench_scores.csv, bench_summary.csv
```

1. `founders` – msprime (Ne 20 000, mu 1.25e-8, Bherer sex-averaged map) founder
   haplotypes on chr20-22 (GRCh37 lengths); a random reference with the ancestral allele
   at every site; `founders.vcf.gz`, `sites.tsv`.
2. `pedsim` – ped-sim pedigrees with planted consanguinity (config `pedigrees`: offspring
   of full sibs, first cousins, second cousins, outbred; five each); `--seg_haps` HBD
   segments are the truth; `--keep_phase --err_rate 0 --miss_rate 0`; output re-bgzipped.
3. `hap_fasta`, `reads`, `map` – per-haplotype FASTA (`bcftools consensus`), ART HiSeq
   2500 150 bp paired reads at 1/2/4/8x total depth, bwa-mem2 + samtools sort.
4. `sites` – one shared ascertained site list (founder MAF >= 0.05) in the forms ANGSD,
   bcftools and GARLIC need, plus `mapped.bed`: ped-sim transmits nothing outside the
   genetic map, so every input and call is restricted to the map's span (118 of 162 Mb).
5. `angsd` / `bcftools_call` – `-GL 1 -doGlf 2 -doMaf 1` beagle likelihoods per depth;
   `mpileup | call -m -C alleles` hard calls with PL at the same sites.
6. `garlic_*` – GARLIC on true genotypes (ceiling), on hard calls with `--error`, on hard
   calls with `--gl-type PL` (bcftools emits PL, not GQ; file prefix `pl_`), and on the
   beagle file with `--gl-beagle` at the shared `--error` and at the `error_sweep` values;
   all with `--no-centromere --winsize 60 --size-bounds 500000 1000000`. A run that
   cannot place the automatic cutoff leaves an empty `.roh.bed`; that is a result.
7. `rohan` – ROHan per BAM at `rohan_depths` (4-6 h per 162 Mb BAM with two threads;
   skipped when the binary is absent or `run_rohan=false`).
8. `score` – base-level precision/recall of runs >= 1 Mb against pedigree HBD and against
   the attainable set (pedigree HBD plus the segments called on the true genotypes, which
   are background founder autozygosity), F_ROH error; per individual and pooled.

`bench.py` is the earlier stage-by-stage driver with the same output layout; the numbers in
the manual's benchmark section were produced with it, and a directory it produced can be
resumed by the Snakefile (the dry run against such a directory reuses founders, pedigrees,
reads, alignments and ROHan, and regenerates the site list and everything downstream of it
because two of the `sites` outputs did not exist under the old layout).
