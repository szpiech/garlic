# Low-coverage ROH benchmark

Snakemake workflow used to benchmark the `low-cov` branch of GARLIC on simulated
low-coverage sequence data. It simulates individuals with known runs of homozygosity (ROH),
simulates reads at 1x, 2x, 4x and 8x, calls ROH with GARLIC, PLINK 1.9 and ROHan, and scores
every call set against the simulated truth. Tool paths in the configuration files are those of
the lab office server; edit the `tools` block for another machine.

The full description (software, versions, install sources, every command and argument, and a
flowchart) is in `garlic_lowcov_benchmark_pipeline.docx`.

## Install

```
conda env create -f envs/env-bench.yaml      # Snakemake, msprime, tskit, samtools, bcftools, bwa-mem2, ART, ANGSD, ...
conda env create -f envs/env-plink.yaml      # PLINK v1.90b6.21
```

ped-sim (github.com/williamslab/ped-sim), ROHan (github.com/grenaud/ROHan) and GARLIC
(`low-cov` branch, `cd src && make`) are built from source. The Bherer et al. (2017) map is
converted to `refined_mf_X.simmap` with the commands in the ped-sim README, and the
interference file is `interfere/nu_p_campbell_X.tsv` from the ped-sim repository. Put all of
these paths in the `tools` block of `config.yaml` (the population configs inherit it).

## Run

```
conda activate env-bench
snakemake -s bench/Snakefile --cores 32                                              # pedigree cohort
snakemake -s bench/Snakefile --cores 32 --configfile bench/config_pop_bottleneck.yaml   # bottleneck population
snakemake -s bench/Snakefile --cores 32 --configfile bench/config_pop_small.yaml        # Ne 2,000 population
snakemake -s bench/Snakefile --cores 16 --config small=true root=<dir>                  # 24-Mb smoke test
snakemake -s bench/Snakefile --cores 32 --config run_rohan=true                         # also ROHan (off by default)
snakemake -s bench/Snakefile -n                                                         # dry run
```

`root` in each configuration file is the output directory. ROHan takes 4.5 to 6 h per BAM with
two threads, so it runs only with `run_rohan=true`, at the depths in `rohan_depths`.

To reuse a directory made by the earlier stage driver `bench.py` (kept in the repository), add
`--rerun-triggers mtime`. After changing a rule, rerun it and everything downstream with
`-R <rule>`; deleting a rule's outputs does not trigger a rerun when the score tables are newer.

## Rules, in order

| rule | software | function |
|---|---|---|
| `founders` | msprime, tskit | 200 diploid founders on chr20-22 (GRCh37 lengths), tree sequences, random reference |
| `pedsim_inputs`, `pedsim` | ped-sim or bcftools | pedigree cohort (offspring of full sibs, first and second cousins, outbred; five each) or 20 unrelated individuals |
| `truth` | Python, tskit | pedigree HBD, or IBD with TMRCA <= 100 generations (used by `score`) |
| `pedsim_bp` | ped-sim | rerun with the logged seed and `--bp`; must reproduce `bench.seg` exactly |
| `truth_all` | Python, tskit | every within-individual IBD tract >= 10 kb with its TMRCA (used by `score_by_size`) |
| `index_reference`, `hap_fasta`, `reads`, `map` | samtools, bcftools, ART, bwa-mem2 | haplotype FASTA, 150-bp paired reads, sorted BAM |
| `sites`, `index_sites`, `bamlist` | Python, ANGSD, htslib | one shared site list (founder MAF >= 0.05) inside the genetic-map span |
| `angsd` | ANGSD | beagle genotype likelihoods (`-GL 1 -doGlf 2`) |
| `bcftools_call`, `gt_check` | bcftools | `call -m -G - -P 0`: maximum-likelihood genotypes with PL; check of GT against argmin(PL) |
| `mapped_beagle`, `mapped_calls`, `tfam`, `truth_vcf` | Python, bcftools | inputs restricted to the map span |
| `garlic_truth`, `garlic_calls`, `garlic_pl`, `garlic_beagle` | GARLIC | true genotypes, genotype calls, `--gl-type PL`, `--gl-beagle` at three `--error` values |
| `plink_roh` | PLINK 1.9 | `--homozyg` with the Ceballos et al. (2018) settings, `--homozyg-window-het` 3 and 1 |
| `rohan` | ROHan | BAM at `rohan_depths` (only with `run_rohan=true`) |
| `score` | Python | runs >= 1 Mb: recall vs `truth`, precision vs the attainable set, FROH error |
| `score_by_size` | Python | recall and precision by size class (0.1-0.5, 0.5-1, >1 Mb) against `truth_all` |

## Genotype calls

`bcftools call -m` by default multiplies each sample's likelihoods by a Hardy-Weinberg prior over
all samples, and at one or two reads that prior decides the genotype (70% of heterozygous calls
at 1x in the pedigree cohort had PL favouring a homozygote). `-G -` removes that prior but leaves
the substitution-rate prior (`-P`), which turned 10% of two-or-more-read heterozygotes into 0/0.
`-G - -P 0` gives the maximum-likelihood genotype; `gl/d<depth>.gtcheck.tsv` tabulates GT against
argmin(PL) by read count.

## Outputs (under `root`)

- `bench_summary.csv`, `bench_scores.csv`: runs >= 1 Mb, pooled and per individual.
- `bench_scores_by_size.csv`: per caller, input, depth and size class: `recall` (truth tracts
  >= 100 kb with TMRCA <= 1000 generations), `recall_any_age`, `precision` (called bases inside
  tracts >= 100 kb of any age), `precision_10kb`, and the TMRCA composition of the called bases.
- `bench_recall_by_size_tmrca.csv`: recall by size and TMRCA class.
- `bench_scores_by_size_ind.csv`: the same bases per individual.
- `roh/`: GARLIC `.roh.bed`, PLINK `plink/*.hom`, ROHan `rohan/<sample>_d<depth>/`.

A GARLIC run that cannot place its automatic LOD cutoff leaves an empty `.roh.bed`; that is a
result (no calls), not a pipeline error. PLINK reads the shared site list with no further MAF or
HWE filter.
