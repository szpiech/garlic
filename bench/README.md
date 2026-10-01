# Low-coverage ROH benchmark

`bench.py` is the driver used to benchmark the `low-cov` branch on simulated
low-coverage sequence data. It is written for the lab office box (paths at the
top of the file) and is kept here for provenance; it is not part of the build.

Pipeline (one stage per function, each idempotent and resumable):

1. `founders` – msprime (Ne 20 000, mu 1.25e-8, Bherer sex-averaged map) founder
   haplotypes on chr20-22 (GRCh37 lengths); a random reference sequence with the
   ancestral allele at every site; `founders.vcf.gz`, `sites.tsv`.
2. `pedsim` – ped-sim pedigrees with planted consanguinity: offspring of full sibs
   (`fs`), first cousins (`c1`), second cousins (`c2`) and outbred (`ob`), five
   each; `--seg_haps` HBD segments are the truth; `--err_rate 0 --miss_rate 0`.
3. `reads` – per-haplotype FASTA (`bcftools consensus`) and ART HiSeq 2500
   150 bp paired reads at 1, 2, 4 and 8x total depth.
4. `map` – bwa-mem2 + samtools sort.
5. `gl` – one shared ascertained site list (founder MAF >= 0.05); ANGSD
   `-GL 1 -doGlf 2 -doMaf 1` beagle likelihoods per depth; `bcftools mpileup |
   call -m -C alleles` hard calls with PL at the same sites.
6. `roh` – GARLIC on true genotypes (ceiling), on hard calls with `--error`, on
   hard calls with `--gl-type PL`, and on the beagle file with `--gl-beagle`;
   all with `--no-centromere --winsize 60 --size-bounds 500000 1000000`.
7. `rohan` – ROHan per BAM (optional; skipped if the binary is absent).
8. `score` – base-level precision/recall of called ROH against HBD truth and
   FROH error, per individual and method -> `bench_scores.csv`.

`BENCH_SMALL=1` runs a 24 Mb smoke test (chr21/22, one pedigree copy, 1x and 4x).
