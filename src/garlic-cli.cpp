#include "garlic-cli.h"
#include "garlic-data.h"   //MISSING: the sentinel a --lod-cutoff must stay above
#include <iostream>
#include <fstream>
#include <utility>
#include "garlic-platform.h"   //garlicMkdir, garlicIsPathSep
#include <cerrno>

const string VERSION = "2.1.0";

const string PREAMBLE = "\ngarlic v" + VERSION + " -- a program to call runs of homozygosity in genetic data.\n\
Source code and binaries can be found at <https://www.github.com/szpiech/garlic>.\n\
\n\
Citations:\n\
\n\
A Blant, et al. (2017) bioRxiv, doi: 10.1101/177352\n\
ZA Szpiech, et al. (2017) Bioinformatics, doi: 10.1093/bioinformatics/btx102\n\
TJ Pemberton, et al. (2012) AJHG, 91: 275–292\n";

const string ARG_OVERLAP_FRAC = "--overlap-frac";
const double DEFAULT_OVERLAP_FRAC = 0.25;
const string HELP_OVERLAP_FRAC = "The fraction of the windows covering a SNP that must be\n\
\tabove the LOD cutoff for the SNP to be in a run. Runs are therefore at\n\
\tleast WINSIZE*OVERLAP_FRAC SNPs long. If 0, garlic uses the smallest\n\
\tsensible value, 1/winsize.";

const string ARG_AUTO_OVERLAP_FRAC = "--auto-overlap-frac";
const bool DEFAULT_AUTO_OVERLAP_FRAC = false;
const string HELP_AUTO_OVERLAP_FRAC = "Choose --overlap-frac from marker density (see\n\
\t--auto-overlap-coef).";

const string ARG_OUTFILE = "--out";
const string DEFAULT_OUTFILE = "outfile";
const string HELP_OUTFILE = "The base name for all output files.";

const string ARG_THREADS = "--threads";
const int DEFAULT_THREADS = 1;
const string HELP_THREADS = "The number of threads for the LOD scores, the LD calculations\n\
\tunder --weighted, the LOD score KDE and ROH assembly.";

const string ARG_ERROR = "--error";
const double DEFAULT_ERROR = -1;
const string HELP_ERROR = "The assumed genotyping error rate. Required. With --gl-beagle\n\
\tit is the chance that a site inside a run looks non-autozygous for any\n\
\treason (see --gl-beagle).";

const string ARG_WINSIZE = "--winsize";
const int DEFAULT_WINSIZE = 0;
const string HELP_WINSIZE = "The window size, in SNPs, over which LOD scores are calculated.\n\
\tWith --auto-winsize it is where the search starts.";

const string ARG_WINSIZE_MULTI = "--winsize-multi";
const int DEFAULT_WINSIZE_MULTI = -1;
const string HELP_WINSIZE_MULTI = "Several window sizes, in SNPs. garlic writes\n\
\tthe LOD score KDE for each and calls no ROH. With --auto-winsize it\n\
\tchooses the best size from the list.";

const string ARG_AUTO_WINSIZE = "--auto-winsize";
const bool DEFAULT_AUTO_WINSIZE = false;
const string HELP_AUTO_WINSIZE = "Choose the window size automatically. With --weighted it\n\
\tis set from SNP density (see --auto-winsize-coef). Otherwise the search\n\
\tstarts at --winsize, which must be given, and grows by\n\
\t--auto-winsize-step SNPs until the LOD score density meets\n\
\t--auto-winsize-threshold, up to --max-winsize.";

const string ARG_AUTO_WINSIZE_STEP = "--auto-winsize-step";
const int DEFAULT_AUTO_WINSIZE_STEP = 10;
const string HELP_AUTO_WINSIZE_STEP = "The step size, in SNPs, of the --auto-winsize search.";

const string ARG_MAX_GAP = "--max-gap";
const int DEFAULT_MAX_GAP = 200000;
const string HELP_MAX_GAP = "No LOD score window is calculated across a gap between two loci\n\
\tlarger than this many bp.";

const string ARG_RESAMPLE = "--resample";
const int DEFAULT_RESAMPLE = 0;
const string HELP_RESAMPLE = "The number of resamples for estimating allele frequencies. If\n\
\t0, garlic uses the frequencies calculated from the data.";

const string ARG_TPED = "--tped";
const string DEFAULT_TPED = "none";
const string HELP_TPED = "A TPED file holding the map and the genotypes.";

const string ARG_TFAM = "--tfam";
const string DEFAULT_TFAM = "none";
const string HELP_TFAM = "A TFAM file holding population and individual IDs.";

const string ARG_VCF = "--vcf";
const string DEFAULT_VCF = "none";
const string HELP_VCF = "A VCF file (plain, gzipped or BGZF) to read genotypes from\n\
\tinstead of --tped. Biallelic SNVs with a GT field are used. Multiallelic\n\
\tsites, indels and symbolic ALTs are skipped, and the counts are\n\
\treported. Sample names come from the #CHROM line, so --tfam is not\n\
\taccepted. Use --pop for population labels.\n\
\t\n\
\tSites must be grouped by chromosome. Ploidy must be 2, except on a\n\
\tchromosome declared with --sex-system or --sex-chr, where a haploid call\n\
\tis read as hemizygous. With --phased every call must use '|'.";

const string ARG_VCF_PASS_ONLY = "--vcf-pass-only";
const bool DEFAULT_VCF_PASS_ONLY = false;
const string HELP_VCF_PASS_ONLY = "Skip sites whose FILTER is neither PASS nor '.'. Without\n\
\tthis flag such sites are kept and their number is reported, because\n\
\tFILTER records the variant caller's thresholds rather than garlic's.";

const string ARG_POP = "--pop";
const string DEFAULT_POP = "none";
const string HELP_POP = "A file mapping sample ID to population, in place of the labels\n\
\tin the TFAM. Whitespace separated, '#' comments and blank lines ignored,\n\
\t.gz accepted:\n\
\t\n\
\t  <sample_id>  <population>  [sex]\n\
\t\n\
\tThe column order is the reverse of a TFAM's first two columns. The\n\
\toptional third column is sex in PLINK coding (1 male, 2 female, 0 or -9\n\
\tunknown). Unlike in a TFAM, any other value is an error.\n\
\t\n\
\tRows are matched by sample ID, not by order. Every sample in the data\n\
\tneeds a row. Extra rows are ignored and counted, so one file can be used\n\
\tfor a subset of the cohort.";

const string ARG_BEAGLE = "--gl-beagle";
const string DEFAULT_BEAGLE = "none";
const string HELP_BEAGLE = "Genotype likelihoods in beagle format, as written by ANGSD\n\
\t-doGlf 2: a header of marker, allele1, allele2 and three columns per\n\
\tindividual, then one row per site with P(reads | allele1/allele1),\n\
\tP(reads | allele1/allele2) and P(reads | allele2/allele2). This is the\n\
\tinput for low-coverage sequence data. The per-site term of the LOD score\n\
\tis averaged over the three genotypes, weighted by the likelihoods, so an\n\
\tuncertain site contributes little evidence and a site with no reads (a\n\
\tflat triplet) contributes exactly 0.\n\
\t\n\
\tReplaces --tped and --vcf. --tfam is still required and is matched to\n\
\tthe beagle columns by position (the order of ANGSD's -b bam list).\n\
\t--error is still required. It is the chance that a site inside a run\n\
\tlooks non-autozygous for any reason, while the likelihoods carry the\n\
\tread-level uncertainty, so --gl-type is refused. Allele frequencies are\n\
\testimated by EM over the likelihoods, and --freq-file overrides them.\n\
\t\n\
\tRun ANGSD on an ascertained site list (-sites, or at least -SNP_pval).\n\
\tGive likelihoods, not -doPost or -doGeno posteriors, which already\n\
\tinclude a Hardy-Weinberg prior. --weighted and --cm still need --map.\n\
\tSee the manual for its limits at low depth.";

const string ARG_TGLS = "--tgls";
const string DEFAULT_TGLS = "none";
const string HELP_TGLS = "A tgls file with one per-genotype value per individual, in the\n\
\trow order of the TPED. See --gl-type for what the values mean. For PL\n\
\tand GL this is not the VCF convention.";

const string ARG_GL_TYPE = "--gl-type";
const string DEFAULT_GL_TYPE = "none";
const string HELP_GL_TYPE = "Which genotype-quality field to use. With --vcf it is read\n\
\tfrom the FORMAT column. With --tgls it names the form of the single\n\
\tvalue per genotype in the .tgls file. The two are not equivalent.\n\
\t\n\
\tWith --vcf (recommended):\n\
\tGQ\terror = 10^(-GQ/10), as the VCF specification defines it.\n\
\tPL\tThe whole PL array is converted to a posterior:\n\
\t  \t  P(g) proportional to 10^(-PL_g/10),  error = 1 - P(called)/sum.\n\
\t  \tA VCF normalises PL so that the called genotype's value is 0 for\n\
\t  \tevery call, and the information is in how much worse the\n\
\t  \talternatives are.\n\
\tGL\tThe same, via PL = -10*GL.\n\
\t  \tA missing field is an error that names whichever of GQ, PL and GL\n\
\t  \tthe file does carry, and '.' for a called genotype is an error.\n\
\t\n\
\tWith --tgls (legacy, one value per genotype):\n\
\tGQ\tAs above. The only form a VCF's values can be used in directly.\n\
\tPL\terror = 1 - 10^(-PL/10), where 10^(-PL/10) is P(genotype correct).\n\
\t  \tThis is not the VCF PL field. Integer PLs give errors of 0, 0.21,\n\
\t  \t0.37, 0.50, ..., and an error of 0.001 is PL = 0.00435. VCF PLs\n\
\t  \tgiven here are rejected. Use --vcf --gl-type PL instead.\n\
\tGL\terror = 1 - 10^GL, with the same caveat.\n\
\t\n\
\tBelow roughly 4x, use --gl-beagle. Each form above gives one error rate\n\
\tper genotype, evaluated at the called genotype, so an uncertain\n\
\thomozygous call still counts as evidence for autozygosity. See\n\
\t\"Low-coverage sequence data\" in the manual.";

const string ARG_MAP = "--map";
const string DEFAULT_MAP = "none";
const string HELP_MAP = "A scaffold genetic map. Sites not in the file are interpolated,\n\
\tand sites outside its bounds are dropped. Required by --weighted and\n\
\t--cm. Format: <chr> <snpid> <genetic pos> <physical pos>. Negative\n\
\tpositions are an error.";

const string ARG_WEIGHTED = "--weighted";
const bool DEFAULT_WEIGHTED = false;
const string HELP_WEIGHTED = "Weight the LOD scores by LD and by the probabilities of no\n\
\trecombination and no mutation (wLOD). Requires --map.";

const string ARG_RAW_LOD = "--raw-lod";
const bool DEFAULT_RAW_LOD = false;
const string HELP_RAW_LOD = "Also write the LOD scores to gzip-compressed files.";

const string ARG_LOD_CUTOFF = "--lod-cutoff";
const double DEFAULT_LOD_CUTOFF = -999999;
const string HELP_LOD_CUTOFF = "A LOD score cutoff above which windows are called\n\
\tautozygous, applied to every population. By default it is chosen\n\
\tautomatically from the KDE of the LOD scores.";

const string ARG_BOUND_SIZE = "--size-bounds";
const double DEFAULT_BOUND_SIZE = -1;
const string HELP_BOUND_SIZE = "The boundaries between ROH size classes,\n\
\tin increasing order. By default they are chosen automatically with a\n\
\tGaussian mixture model (see --nclust).";

const string ARG_TPED_MISSING = "--tped-missing";
const char DEFAULT_TPED_MISSING = '0';
const string HELP_TPED_MISSING = "The one-character missing-data code in a TPED.";

const string ARG_FREQ_FILE = "--freq-file";
const string DEFAULT_FREQ_FILE = "none";
const string HELP_FREQ_FILE = "A file of allele frequencies for every variant, with one\n\
\tfrequency column per population headed by its label, or a single column\n\
\tthat applies to every population:\n\
\tCHR SNP POS ALLELE FREQ\n\
\t<chr> <locus ID> <pos> <allele> <freq>\n\
\tBy default the frequencies are calculated from the data. The .freq.gz\n\
\tthat a run writes can be given back here.";

const string ARG_FREQ_ONLY = "--freq-only";
const bool DEFAULT_FREQ_ONLY = false;
const string HELP_FREQ_ONLY = "Write the frequency file from the data and exit. Uses little\n\
\tmemory.";

const string ARG_STRATIFY_CUTOFF = "--stratify-cutoff";
const bool DEFAULT_STRATIFY_CUTOFF = false;
const string HELP_STRATIFY_CUTOFF = "Choose the automatic LOD cutoff within strata of\n\
\tindividuals that carry the same evidence per window, so that the shallow\n\
\tmembers of a cohort of mixed depth are not over-called. Individuals are\n\
\tsorted by the location of their non-autozygous LOD mode, neighbours\n\
\twhose modes differ by less than a quarter of the smaller spread share a\n\
\tstratum, and the pooled rule is applied within each stratum. A cohort of\n\
\tuniform depth is one stratum and gets exactly the pooled cutoff. The\n\
\tper-individual mode, spread, stratum and cutoff are written to\n\
\t<out>.ind.tsv. Needs the automatic cutoff.";

const string ARG_KDE_SUBSAMPLE = "--kde-subsample";
const int DEFAULT_KDE_SUBSAMPLE = 0;
const string HELP_KDE_SUBSAMPLE = "The number of individuals sampled at random for the LOD\n\
\tscore KDE. If the population has fewer, all are used. The default, 0,\n\
\tuses every individual, so the cutoff depends on the data alone. A\n\
\tpositive value saves memory but makes the cutoff depend on --seed.";

extern const string ARG_LD_SUBSAMPLE = "--ld-subsample";
extern const int DEFAULT_LD_SUBSAMPLE = 0;
extern const string HELP_LD_SUBSAMPLE = "The number of individuals sampled at random for the LD\n\
\tcalculations of --weighted. If the population has fewer, all are used. 0\n\
\tor less uses every individual, which takes longer.";

const string ARG_BUILD = "--build";
const string DEFAULT_BUILD = "none";
const string HELP_BUILD = "The genome build: hg18, hg19, hg38 or t2t-chm13\n\
\t(T2T-CHM13v2.0, UCSC hs1). It supplies the centromere locations, the\n\
\tchromosome lengths and the human pseudoautosomal regions.";

const string ARG_CENTROMERE_FILE = "--centromere";
const string DEFAULT_CENTROMERE_FILE = "none";
const string HELP_CENTROMERE_FILE = "A file of custom centromere boundaries, <chr> <start>\n\
\t<end>. One of --build, --centromere or --no-centromere is required.";

const string ARG_M = "--M";
const int DEFAULT_M = 7;
const string HELP_M = "The expected number of meioses since a recent common ancestor, for\n\
\t--weighted.";

const string ARG_MU = "--mu";
const double DEFAULT_MU = 1e-9;
const string HELP_MU = "The mutation rate per bp per generation, for --weighted.";

const string ARG_PHASED = "--phased";
const bool DEFAULT_PHASED = false;
const string HELP_PHASED = "The data are phased. Under --weighted, LD is then computed as\n\
\tr2 rather than hr2, which uses more memory. With --vcf every call must\n\
\tuse '|'. Has no effect on the LOD scores without --weighted.";

const string ARG_NCLUST = "--nclust";
const int DEFAULT_NCLUST = 3;
const string HELP_NCLUST = "The number of Gaussian mixture components used to classify ROH\n\
\tlengths.";

const string ARG_CM = "--cm";
const bool DEFAULT_CM = false;
const string HELP_CM = "Report ROH lengths and size classes in cM. ROH boundaries are still\n\
\treported in bp. Requires --map.";

const string ARG_KDE_THINNING = "--no-kde-thinning";
const bool DEFAULT_KDE_THINNING = false;
const string HELP_KDE_THINNING = "Deprecated. The same as --kde-thin-step 1, which sends\n\
\tevery LOD score window to the KDE and can increase runtime considerably.";

const string ARG_MAX_WINSIZE = "--max-winsize";
const int DEFAULT_MAX_WINSIZE = 1000;
const string HELP_MAX_WINSIZE = "The largest window size that --auto-winsize will try.";

const string ARG_DUMP_DOCS = "--dump-docs";
const string DEFAULT_DUMP_DOCS = "__none";
const string HELP_DUMP_DOCS = "Write the command line reference to stdout and exit: txt\n\
\tfor the README block, tex for the manual's appendix. Used by 'make\n\
\tdocs'.";

const string ARG_LOAD_PARAMS = "--load-params";
const string DEFAULT_LOAD_PARAMS = "__none";
const string HELP_LOAD_PARAMS = "Read flag values from a <out>.params.json written by an\n\
\tearlier run. Flags on the command line take precedence, so a run can be\n\
\trepeated with one parameter changed.";

const string ARG_QUIET = "--quiet";
const bool DEFAULT_QUIET = false;
const string HELP_QUIET = "Suppress the progress bar and informational messages. Errors and\n\
\twarnings still go to stderr and to <out>.error.";

const string ARG_VERBOSE = "--verbose";
const bool DEFAULT_VERBOSE = false;
const string HELP_VERBOSE = "Show the progress bar even when stderr is not a terminal. By\n\
\tdefault it is shown only on a terminal, because it is drawn with\n\
\tbackspaces.";

const string ARG_NO_CENTROMERE = "--no-centromere";
const bool DEFAULT_NO_CENTROMERE = false;
const string HELP_NO_CENTROMERE = "Treat every chromosome as having no assembly gap,\n\
\tinstead of requiring --build or --centromere. For assemblies with no\n\
\tcentromere gap, and for references with no built-in table. Cannot be\n\
\tcombined with --build or --centromere.";

const string ARG_AUTOSOMES_ONLY = "--autosomes-only";
const bool DEFAULT_AUTOSOMES_ONLY = false;
const string HELP_AUTOSOMES_ONLY = "Drop every chromosome that is not an autosome before\n\
\tcalling: X, Y, Z, W and the mitochondrion by name, and 23, 24, 25 and 26\n\
\twhen --build names a human assembly. A hemizygous chromosome analysed as\n\
\tan autosome is called as one run along its whole length. Without a human\n\
\t--build, bare 23-26 are refused rather than dropped, because they are\n\
\tautosomes in a species with that many chromosomes. See --sex-system.";

const string ARG_SEX_SYSTEM = "--sex-system";
const string DEFAULT_SEX_SYSTEM = "unset";
const string HELP_SEX_SYSTEM = "Which sex is heterogametic, and so which individuals are\n\
\themizygous on the shared sex chromosome:\n\
\t  xy    individuals coded 1 (male) in the TFAM or --pop are heterogametic\n\
\t  zw    individuals coded 2 (female) are heterogametic\n\
\t  none  every chromosome is diploid in every individual, and sex\n\
\t        chromosomes are not reported\n\
\tPLINK codes 1 male and 2 female in every species, so in a ZW species the\n\
\thomogametic sex is coded 1. Leaving the flag unset is not the same as\n\
\t'none'. If a sex chromosome is present and the flag is unset, garlic\n\
\tstops rather than analysing it as an autosome.";

const string ARG_SEX_CHR = "--sex-chr";
const string HELP_SEX_CHR = "The chromosome that is diploid in the\n\
\thomogametic sex and hemizygous in the heterogametic one (X under xy, Z\n\
\tunder zw). Needed only when it has another name, e.g. --sex-system zw\n\
\t--sex-chr LG12. Requires --sex-system.\n\
\tDefault: the conventional name for --sex-system";

const string ARG_SEX_CHR_DEGENERATE = "--sex-chr-degenerate";
const string HELP_SEX_CHR_DEGENERATE = "The chromosome carried only by the\n\
\theterogametic sex (Y under xy, W under zw). It is dropped before\n\
\tcalling, since no individual can carry a run of homozygosity on it.\n\
\tDefault: the conventional name for --sex-system";

const string ARG_FROH_DENOM = "--froh-denominator";
const string DEFAULT_FROH_DENOM = "analyzed";
const string HELP_FROH_DENOM = "What <out>.froh.tsv divides the autozygous total\n\
\tby.\n\
\tanalyzed (default): the span the markers cover, from the first to the\n\
\tlast locus of each analysed chromosome, less the centromere and any\n\
\tregion excluded with --par. It moves with the marker set, so FROH from\n\
\ttwo different arrays are not directly comparable.\n\
\tassembly: the full length of each analysed chromosome, as published FROH\n\
\tdenominators use. Lengths come from --chr-lengths, --build or a VCF's\n\
\t##contig header, in that order, and a chromosome with no length stops\n\
\tthe run. Not available with --cm.\n\
\tDefault: " + DEFAULT_FROH_DENOM;

const string ARG_CHR_LENGTHS = "--chr-lengths";
const string DEFAULT_CHR_LENGTHS = "_NONE";
const string HELP_CHR_LENGTHS = "A file of chromosome lengths in two whitespace-separated\n\
\tcolumns, chromosome and length. A samtools faidx .fai is accepted as it\n\
\tis. Needed for --froh-denominator assembly on a non-human assembly or on\n\
\tTPED input, where nothing in the data states a length. Takes precedence\n\
\tover --build and a VCF's ##contig header.\n\
\tDefault: none";

const string ARG_SEXCHR_LOD_CUTOFF = "--sexchr-lod-cutoff";
const double DEFAULT_SEXCHR_LOD_CUTOFF = -999999;
const string HELP_SEXCHR_LOD_CUTOFF = "Call the shared sex chromosome at this LOD score\n\
\tcutoff instead of the autosomal one. Every run logs the cutoff that the\n\
\tsex chromosome's own windows would give, which is where a value for this\n\
\tcomes from. The autosomal cutoff is the default because one chromosome's\n\
\twindows from part of a cohort often have no second minimum in the\n\
\tdensity.\n\
\tDefault: the autosomal cutoff";

const string ARG_SEXCHR_WINSIZE = "--sexchr-winsize";
const int DEFAULT_SEXCHR_WINSIZE = 0;
const string HELP_SEXCHR_WINSIZE = "Use windows of this many SNPs on the shared sex\n\
\tchromosome instead of the autosomal size. The overlap fraction is\n\
\tunchanged. A LOD score cutoff is specific to the window size it was\n\
\testimated at, so --sexchr-lod-cutoff is required as well.\n\
\tDefault: the autosomal window size";

const string ARG_PAR = "--par";
const string HELP_PAR = "Pseudoautosomal regions on the shared sex\n\
\tchromosome, as <chr>:<start>-<end>, comma or space separated. The loci\n\
\tinside are dropped, since the regions are diploid in both sexes and too\n\
\tshort to call. Any number of regions per chromosome is allowed.\n\
\tCoordinates are in the input's assembly, and a region cannot be named by\n\
\tchromosome code, because PLINK's 25 is an ordinary autosome in most\n\
\tspecies. --build supplies the human X regions from the Genome Reference\n\
\tConsortium's definition (see centromeres/par_regions.txt). Naming any\n\
\tregion here replaces them, and --par none turns them off.\n\
\tDefault: the regions for --build, or none without one";

const string ARG_EXCLUDE_REGIONS = "--exclude-regions";
const string DEFAULT_EXCLUDE_REGIONS = "none";
const string HELP_EXCLUDE_REGIONS = "A BED file of regions whose sites are dropped before\n\
\tanything is computed, such as a mappability or repeat mask or segmental\n\
\tduplications. Mapping bias in such regions turns heterozygotes into\n\
\thomozygotes, which looks like autozygosity. Standard BED (0-based\n\
\thalf-open <chrom> <chromStart> <chromEnd>, further columns ignored,\n\
\ttrack/browser/# lines skipped), gzipped or not. Chromosomes not in the\n\
\tdata are ignored. The dropped sites are removed from the frequency file\n\
\tand from every window, but the FROH denominator is not reduced (unlike\n\
\t--par). Applied on every chromosome. Use --par for pseudoautosomal\n\
\tregions.";

const string ARG_PAR_FILE = "--par-file";
const string DEFAULT_PAR_FILE = "none";
const string HELP_PAR_FILE = "The same regions read from a file, one per line as <chr>\n\
\t<start> <end>, with '#' comments and blank lines ignored and .gz\n\
\taccepted. Unlike --centromere, repeated rows for a chromosome add a\n\
\tregion. May be combined with --par.";

const string ARG_HET_RATE_BOUNDS = "--het-rate-bounds";
const string HELP_HET_RATE_BOUNDS = "<lo> <hi>: heterozygosity on the\n\
\tshared sex chromosome below which an individual reads as hemizygous and\n\
\tabove which it reads as diploid. In between it is ambiguous, and garlic\n\
\tneither calls ROH there for that individual nor counts its alleles. Used\n\
\tto check each recorded sex and to infer the missing ones. The bounds are\n\
\tfixed rather than fitted, because a cohort of one sex has no bimodality\n\
\tto fit.";

const string ARG_HAPLOID_CHR = "--haploid-chr";
const string HELP_HAPLOID_CHR = "Chromosomes haploid in every individual\n\
\t(mitochondrion, chloroplast). Dropped before calling.\n\
\tDefault: names matching M or MT, plus 26 under a human --build";

const string ARG_CHR = "--chr";
const string HELP_CHR = "Analyse only these chromosomes, e.g. --chr chr1\n\
\tchr2 chrX. Names are matched after normalisation, so '1' and 'chr1' are\n\
\tthe same. Naming a chromosome that is not in the data is an error.\n\
\tDefault: all chromosomes in the input";

const string ARG_OUTDIR = "--outdir";
const string DEFAULT_OUTDIR = "";
const string HELP_OUTDIR = "Write all output files to this directory, which is created if\n\
\tneeded. The same as prefixing --out with the path.";

const string ARG_POOL = "--pool-populations";
const bool DEFAULT_POOL = false;
const string HELP_POOL = "Analyse every individual as one population, ignoring\n\
\tthe labels, as garlic did before per-population analysis. Output names\n\
\tare then <out>.roh.bed rather than <out>.<POP>.roh.bed.\n\
\t\n\
\tAllele frequencies are pooled over everyone, which inflates\n\
\theterozygosity relative to any one population and biases the LOD scores.\n\
\tUse it when the labels are not populations (column 1 of a TFAM is the\n\
\tPLINK family ID) or to reproduce an earlier run.";

const string ARG_FROH = "--froh";
const bool DEFAULT_FROH = false;
const string HELP_FROH = "Also write <out>.froh.tsv, with each individual's autozygous\n\
\ttotal and fraction by ROH size class. The denominator is stated in the\n\
\theader.";

const string ARG_AUTO_WINSIZE_THRESHOLD = "--auto-winsize-threshold";
const double DEFAULT_AUTO_WINSIZE_THRESHOLD = 0.50;
const string HELP_AUTO_WINSIZE_THRESHOLD = "The smoothness criterion at which the\n\
\t--auto-winsize search stops. Lower values demand a smoother LOD score\n\
\tdensity, and so a larger window.";

const string ARG_KDE_POINTS = "--kde-points";
const int DEFAULT_KDE_POINTS = 512;
const string HELP_KDE_POINTS = "The number of grid points at which the LOD score density is\n\
\tevaluated. This also sets the resolution of the LOD cutoff.";

const string ARG_KDE_CUT = "--kde-cut";
const double DEFAULT_KDE_CUT = 3;
const string HELP_KDE_CUT = "Extend the density grid this many bandwidths beyond the\n\
\trange of the data.";

const string ARG_MODE_SPAN = "--mode-smooth-span";
const int DEFAULT_MODE_SPAN = 20;
const string HELP_MODE_SPAN = "The number of adjacent grid points used to smooth the\n\
\tdensity when finding its modes. Must be smaller than --kde-points.";

const string ARG_AUTO_WINSIZE_COEF = "--auto-winsize-coef";
const string HELP_AUTO_WINSIZE_COEF = "<slope> <intercept>: coefficients\n\
\tof winsize = slope*log(density) + intercept, used by --auto-winsize with\n\
\t--weighted. The defaults are an empirical fit to human SNP-array data.\n\
\tChange them for other ascertainment schemes or species.\n\
\tDefault: 8.3235 138.0521";

const string ARG_AUTO_OVERLAP_COEF = "--auto-overlap-coef";
const string HELP_AUTO_OVERLAP_COEF = "<slope> <intercept>: coefficients\n\
\tof overlap% = slope*log(density) + intercept, used by\n\
\t--auto-overlap-frac. The defaults are an empirical fit to human\n\
\tSNP-array data.\n\
\tDefault: 6.375 63.888";

const string ARG_GMM_MAX_ITER = "--gmm-max-iter";
const int DEFAULT_GMM_MAX_ITER = 1000;
const string HELP_GMM_MAX_ITER = "The maximum number of EM iterations when fitting ROH size\n\
\tclasses.";

const string ARG_GMM_TOL = "--gmm-tol";
const double DEFAULT_GMM_TOL = 1e-5;
const string HELP_GMM_TOL = "The convergence tolerance of the ROH size-class EM fit.";

const string ARG_KDE_THIN_STEP = "--kde-thin-step";
const int DEFAULT_KDE_THIN_STEP = 0;
const string HELP_KDE_THIN_STEP = "Send every Nth LOD score window to the KDE. 0, the\n\
\tdefault, uses the window size, so the windows used do not overlap. 1\n\
\tuses every window. --no-kde-thinning is a deprecated alias for 1.";

const string ARG_VERSION = "--version";
const bool DEFAULT_VERSION = false;
const string HELP_VERSION = "Print the version and exit.";

const string ARG_FORCE = "--force";
const bool DEFAULT_FORCE = false;
const string HELP_FORCE = "Overwrite existing output files. Without it garlic will not\n\
\toverwrite an existing <out>.roh.bed.";

const string ARG_SEED = "--seed";
const int DEFAULT_SEED = 0;
const string HELP_SEED = "The random seed. Affects --kde-subsample, --ld-subsample and\n\
\t--resample. If 0 (the default), a seed is drawn and written to the log\n\
\tso that the run can be repeated exactly.";

const string ARG_FEATURES = "--features";
const string DEFAULT_FEATURES = "none";
const string HELP_FEATURES = "A file classifying sites, in four whitespace-separated\n\
\tcolumns:\n\
\t  <chr> <pos> <allele> <class>\n\
\tLines beginning '#' are comments. Class labels are arbitrary text, and\n\
\tthere may be any number of them. A site may appear more than once, so\n\
\ttwo annotation schemes can classify the same variant, but the same class\n\
\ttwice at one site is an error. Chromosome names are matched ignoring\n\
\tcase and a 'chr' prefix.\n\
\t\n\
\tFor each individual garlic counts the classified genotypes inside runs\n\
\tof each size class and outside them, and writes <out>.counts.tsv.\n\
\tGenotypes are read from the input the run was called from, or from\n\
\t--tped-counting or --vcf-counting if given.";

const string ARG_FEATURE_TPED = "--tped-counting";
const string DEFAULT_FEATURE_TPED = "none";
const string HELP_FEATURE_TPED = "A TPED to count classified genotypes from, instead of\n\
\tthe file the run was called from, for example when the calls come from\n\
\tarray data and the classified variants are sequenced. The two need not\n\
\tshare a single site. Requires --tfam-counting and --features.\n\
\t\n\
\tThese genotypes take no part in calling and are not filtered, since\n\
\tmonomorphic and rare sites are what is being counted.";

const string ARG_ROH_FILE = "--roh-file";
const string DEFAULT_ROH_FILE = "none";
const string HELP_ROH_FILE = "Count classified genotypes against the calls in an existing\n\
\tgarlic .roh.bed instead of calling ROH. --tped and --vcf are not used.\n\
\tGive the genotypes with --tped-counting and --tfam-counting or with\n\
\t--vcf-counting, and the classification with --features.\n\
\t\n\
\tA .roh.bed does not record which chromosomes were analysed or what the\n\
\tsize class letters mean, so no site can be marked UNASSESSED and the\n\
\tsize class boundaries are reported as unknown. Counting inside the run\n\
\tthat makes the calls knows both.\n\
\t\n\
\tgarlic 1.1.6a was released with two chromStart conventions under one\n\
\tversion string, so the convention is inferred from the file and\n\
\treported.";

const string ARG_FEATURE_VCF = "--vcf-counting";
const string DEFAULT_FEATURE_VCF = "none";
const string HELP_FEATURE_VCF = "A VCF to count classified genotypes from, instead of\n\
\tthe file the run was called from. Sample names come from its #CHROM\n\
\tline. Requires --features.\n\
\t\n\
\tUnlike the VCF a run is called from, indels and multiallelic sites are\n\
\tread here. The classified allele is matched against REF and every ALT,\n\
\tso a classification that is not a biallelic SNV still counts.";

const string ARG_FEATURE_TFAM = "--tfam-counting";
const string DEFAULT_FEATURE_TFAM = "none";
const string HELP_FEATURE_TFAM = "The TFAM naming the samples in --tped-counting, in\n\
\tcolumn order. Individuals are matched to the ROH calls by ID, and one\n\
\tpresent in only one of the two files is reported and left out of the\n\
\ttable.";


//Written by the Makefile so the stamp cannot go stale.  Absent when building
//outside the Makefile or from a tarball with no .git, hence the guard.
#if defined(__has_include)
#  if __has_include("garlic-version.h")
#    include "garlic-version.h"
#  endif
#endif
#ifndef GARLIC_GIT_SHA
#define GARLIC_GIT_SHA "unknown"
#endif

param_t *getCLI(int argc, char *argv[], int &status)
{
	param_t *params = new param_t;
	params->setPreamble(PREAMBLE);
	params->addFlag(ARG_OVERLAP_FRAC, DEFAULT_OVERLAP_FRAC, "", HELP_OVERLAP_FRAC);
	params->addFlag(ARG_AUTO_OVERLAP_FRAC, DEFAULT_AUTO_OVERLAP_FRAC, "", HELP_AUTO_OVERLAP_FRAC);
	params->addFlag(ARG_OUTFILE, DEFAULT_OUTFILE, "", HELP_OUTFILE);
	params->addFlag(ARG_THREADS, DEFAULT_THREADS, "", HELP_THREADS);
	params->addFlag(ARG_ERROR, DEFAULT_ERROR, "", HELP_ERROR);
	params->addFlag(ARG_WINSIZE, DEFAULT_WINSIZE, "", HELP_WINSIZE);
	params->addFlag(ARG_MAX_GAP, DEFAULT_MAX_GAP, "", HELP_MAX_GAP);
	params->addFlag(ARG_RESAMPLE, DEFAULT_RESAMPLE, "", HELP_RESAMPLE);
	params->addFlag(ARG_TPED, DEFAULT_TPED, "", HELP_TPED);
	params->addFlag(ARG_TFAM, DEFAULT_TFAM, "", HELP_TFAM);
	params->addFlag(ARG_VCF, DEFAULT_VCF, "", HELP_VCF);
	params->addFlag(ARG_VCF_PASS_ONLY, DEFAULT_VCF_PASS_ONLY, "", HELP_VCF_PASS_ONLY);
	params->addFlag(ARG_POP, DEFAULT_POP, "", HELP_POP);
	params->addFlag(ARG_TGLS, DEFAULT_TGLS, "", HELP_TGLS);
	params->addFlag(ARG_BEAGLE, DEFAULT_BEAGLE, "", HELP_BEAGLE);
	params->addFlag(ARG_GL_TYPE, DEFAULT_GL_TYPE, "", HELP_GL_TYPE);
	params->addFlag(ARG_MAP, DEFAULT_MAP, "", HELP_MAP);
	params->addFlag(ARG_WEIGHTED, DEFAULT_WEIGHTED, "", HELP_WEIGHTED);
	params->addFlag(ARG_RAW_LOD, DEFAULT_RAW_LOD, "", HELP_RAW_LOD);
	params->addListFlag(ARG_BOUND_SIZE, DEFAULT_BOUND_SIZE, "", HELP_BOUND_SIZE);
	params->addFlag(ARG_LOD_CUTOFF, DEFAULT_LOD_CUTOFF, "", HELP_LOD_CUTOFF);
	params->addFlag(ARG_TPED_MISSING, DEFAULT_TPED_MISSING, "", HELP_TPED_MISSING);
	params->addFlag(ARG_FREQ_FILE, DEFAULT_FREQ_FILE, "", HELP_FREQ_FILE);
	params->addFlag(ARG_FREQ_ONLY, DEFAULT_FREQ_ONLY, "", HELP_FREQ_ONLY);
	params->addListFlag(ARG_WINSIZE_MULTI, DEFAULT_WINSIZE_MULTI, "", HELP_WINSIZE_MULTI);
	params->addFlag(ARG_KDE_SUBSAMPLE, DEFAULT_KDE_SUBSAMPLE , "", HELP_KDE_SUBSAMPLE);
	params->addFlag(ARG_STRATIFY_CUTOFF, DEFAULT_STRATIFY_CUTOFF, "", HELP_STRATIFY_CUTOFF);
	params->addFlag(ARG_LD_SUBSAMPLE, DEFAULT_LD_SUBSAMPLE , "", HELP_LD_SUBSAMPLE);
	params->addFlag(ARG_AUTO_WINSIZE, DEFAULT_AUTO_WINSIZE, "", HELP_AUTO_WINSIZE);
	params->addFlag(ARG_AUTO_WINSIZE_STEP, DEFAULT_AUTO_WINSIZE_STEP, "", HELP_AUTO_WINSIZE_STEP);
	params->addFlag(ARG_BUILD, DEFAULT_BUILD, "", HELP_BUILD);
	params->addFlag(ARG_CENTROMERE_FILE, DEFAULT_CENTROMERE_FILE, "", HELP_CENTROMERE_FILE);
	params->addFlag(ARG_M, DEFAULT_M, "", HELP_M);
	params->addFlag(ARG_MU, DEFAULT_MU, "", HELP_MU);
	params->addFlag(ARG_PHASED, DEFAULT_PHASED, "", HELP_PHASED);
	params->addFlag(ARG_NCLUST, DEFAULT_NCLUST, "", HELP_NCLUST);
	params->addFlag(ARG_CM, DEFAULT_CM, "", HELP_CM);
	params->addFlag(ARG_KDE_THINNING, DEFAULT_KDE_THINNING, "", HELP_KDE_THINNING);
	params->addFlag(ARG_SEED, DEFAULT_SEED, "", HELP_SEED);
	params->addFlag(ARG_MAX_WINSIZE, DEFAULT_MAX_WINSIZE, "", HELP_MAX_WINSIZE);
	params->addFlag(ARG_VERSION, DEFAULT_VERSION, "", HELP_VERSION);
	params->addFlag(ARG_FORCE, DEFAULT_FORCE, "", HELP_FORCE);
	params->addFlag(ARG_KDE_THIN_STEP, DEFAULT_KDE_THIN_STEP, "", HELP_KDE_THIN_STEP);
	params->addFlag(ARG_DUMP_DOCS, DEFAULT_DUMP_DOCS, "", HELP_DUMP_DOCS);
	params->addFlag(ARG_LOAD_PARAMS, DEFAULT_LOAD_PARAMS, "", HELP_LOAD_PARAMS);
	params->addFlag(ARG_QUIET, DEFAULT_QUIET, "", HELP_QUIET);
	params->addFlag(ARG_VERBOSE, DEFAULT_VERBOSE, "", HELP_VERBOSE);
	params->addFlag(ARG_NO_CENTROMERE, DEFAULT_NO_CENTROMERE, "", HELP_NO_CENTROMERE);
	params->addFlag(ARG_AUTOSOMES_ONLY, DEFAULT_AUTOSOMES_ONLY, "", HELP_AUTOSOMES_ONLY);
	params->addFlag(ARG_SEX_SYSTEM, DEFAULT_SEX_SYSTEM, "", HELP_SEX_SYSTEM);
	params->addListFlag(ARG_SEX_CHR, "_NONE", "", HELP_SEX_CHR);
	params->addListFlag(ARG_SEX_CHR_DEGENERATE, "_NONE", "", HELP_SEX_CHR_DEGENERATE);
	params->addListFlag(ARG_HAPLOID_CHR, "_NONE", "", HELP_HAPLOID_CHR);
	params->addListFlag(ARG_HET_RATE_BOUNDS, 0.02, "", HELP_HET_RATE_BOUNDS);
	params->addFlag(ARG_FROH_DENOM, DEFAULT_FROH_DENOM, "", HELP_FROH_DENOM);
	params->addFlag(ARG_CHR_LENGTHS, DEFAULT_CHR_LENGTHS, "", HELP_CHR_LENGTHS);
	params->addFlag(ARG_SEXCHR_LOD_CUTOFF, DEFAULT_SEXCHR_LOD_CUTOFF, "", HELP_SEXCHR_LOD_CUTOFF);
	params->addFlag(ARG_SEXCHR_WINSIZE, DEFAULT_SEXCHR_WINSIZE, "", HELP_SEXCHR_WINSIZE);
	params->addListFlag(ARG_PAR, "_NONE", "", HELP_PAR);
	params->addFlag(ARG_PAR_FILE, DEFAULT_PAR_FILE, "", HELP_PAR_FILE);
	params->addFlag(ARG_EXCLUDE_REGIONS, DEFAULT_EXCLUDE_REGIONS, "", HELP_EXCLUDE_REGIONS);
	params->addListFlag(ARG_CHR, "_ALL", "", HELP_CHR);
	params->addFlag(ARG_OUTDIR, DEFAULT_OUTDIR, "", HELP_OUTDIR);
	params->addFlag(ARG_POOL, DEFAULT_POOL, "", HELP_POOL);
	params->addFlag(ARG_FROH, DEFAULT_FROH, "", HELP_FROH);
	params->addFlag(ARG_AUTO_WINSIZE_THRESHOLD, DEFAULT_AUTO_WINSIZE_THRESHOLD, "", HELP_AUTO_WINSIZE_THRESHOLD);
	params->addFlag(ARG_KDE_POINTS, DEFAULT_KDE_POINTS, "", HELP_KDE_POINTS);
	params->addFlag(ARG_KDE_CUT, DEFAULT_KDE_CUT, "", HELP_KDE_CUT);
	params->addFlag(ARG_MODE_SPAN, DEFAULT_MODE_SPAN, "", HELP_MODE_SPAN);
	params->addListFlag(ARG_AUTO_WINSIZE_COEF, 0.0, "", HELP_AUTO_WINSIZE_COEF);
	params->addListFlag(ARG_AUTO_OVERLAP_COEF, 0.0, "", HELP_AUTO_OVERLAP_COEF);
	params->addFlag(ARG_GMM_MAX_ITER, DEFAULT_GMM_MAX_ITER, "", HELP_GMM_MAX_ITER);
	params->addFlag(ARG_GMM_TOL, DEFAULT_GMM_TOL, "", HELP_GMM_TOL);
	params->addFlag(ARG_FEATURES, DEFAULT_FEATURES, "", HELP_FEATURES);
	params->addFlag(ARG_FEATURE_TPED, DEFAULT_FEATURE_TPED, "", HELP_FEATURE_TPED);
	params->addFlag(ARG_FEATURE_TFAM, DEFAULT_FEATURE_TFAM, "", HELP_FEATURE_TFAM);
	params->addFlag(ARG_FEATURE_VCF, DEFAULT_FEATURE_VCF, "", HELP_FEATURE_VCF);
	params->addFlag(ARG_ROH_FILE, DEFAULT_ROH_FILE, "", HELP_ROH_FILE);

	//---- how --help and the generated reference are grouped ----------------
	//
	//Ordered roughly as a run meets them -- what you read in, what genome it
	//is, how windows and LOD scores are formed, how runs are called and
	//classified, what comes out -- rather than alphabetically, which put
	//--M next to --auto-overlap-coef and told a newcomer nothing.  Flags are
	//sorted within each group, so the editorial judgement is the grouping and
	//nothing else.
	//
	//Adding a flag means naming it here too.  That is enforced, not trusted:
	//a flag in no group makes --help say so and makes 'make docs' fail, and
	//the docs job in CI runs 'make docs', so it cannot reach a release.
	{
		vector< pair< string, vector<string> > > cats;
		vector<string> f;

		f.clear();
		f.push_back(ARG_TPED);        f.push_back(ARG_TFAM);
		f.push_back(ARG_VCF);         f.push_back(ARG_VCF_PASS_ONLY);
		f.push_back(ARG_BEAGLE);
		f.push_back(ARG_TGLS);        f.push_back(ARG_GL_TYPE);
		f.push_back(ARG_MAP);         f.push_back(ARG_FREQ_FILE);
		f.push_back(ARG_TPED_MISSING);
		cats.push_back(make_pair(string("Input files"), f));

		f.clear();
		f.push_back(ARG_POP);         f.push_back(ARG_POOL);
		cats.push_back(make_pair(string("Populations"), f));

		f.clear();
		f.push_back(ARG_BUILD);       f.push_back(ARG_CENTROMERE_FILE);
		f.push_back(ARG_NO_CENTROMERE); f.push_back(ARG_CHR);
		f.push_back(ARG_CHR_LENGTHS); f.push_back(ARG_AUTOSOMES_ONLY);
		f.push_back(ARG_SEX_SYSTEM);  f.push_back(ARG_SEX_CHR);
		f.push_back(ARG_SEX_CHR_DEGENERATE); f.push_back(ARG_HAPLOID_CHR);
		f.push_back(ARG_PAR);         f.push_back(ARG_PAR_FILE);
		f.push_back(ARG_EXCLUDE_REGIONS);
		cats.push_back(make_pair(string("Genome and chromosomes"), f));

		f.clear();
		f.push_back(ARG_WINSIZE);     f.push_back(ARG_AUTO_WINSIZE);
		f.push_back(ARG_AUTO_WINSIZE_STEP); f.push_back(ARG_AUTO_WINSIZE_THRESHOLD);
		f.push_back(ARG_AUTO_WINSIZE_COEF); f.push_back(ARG_MAX_WINSIZE);
		f.push_back(ARG_WINSIZE_MULTI); f.push_back(ARG_WEIGHTED);
		f.push_back(ARG_M);           f.push_back(ARG_MU);
		f.push_back(ARG_ERROR);       f.push_back(ARG_MAX_GAP);
		f.push_back(ARG_PHASED);
		cats.push_back(make_pair(string("LOD scores and windows"), f));

		f.clear();
		f.push_back(ARG_LOD_CUTOFF);  f.push_back(ARG_OVERLAP_FRAC);
		f.push_back(ARG_AUTO_OVERLAP_FRAC); f.push_back(ARG_AUTO_OVERLAP_COEF);
		f.push_back(ARG_SEXCHR_LOD_CUTOFF); f.push_back(ARG_SEXCHR_WINSIZE);
		f.push_back(ARG_HET_RATE_BOUNDS); f.push_back(ARG_LD_SUBSAMPLE);
		cats.push_back(make_pair(string("Calling runs of homozygosity"), f));

		f.clear();
		f.push_back(ARG_STRATIFY_CUTOFF);
		f.push_back(ARG_KDE_SUBSAMPLE); f.push_back(ARG_KDE_POINTS);
		f.push_back(ARG_KDE_CUT);     f.push_back(ARG_KDE_THIN_STEP);
		f.push_back(ARG_KDE_THINNING); f.push_back(ARG_MODE_SPAN);
		cats.push_back(make_pair(string("Choosing the LOD cutoff automatically"), f));

		f.clear();
		f.push_back(ARG_BOUND_SIZE); f.push_back(ARG_NCLUST);
		f.push_back(ARG_CM);          f.push_back(ARG_GMM_MAX_ITER);
		f.push_back(ARG_GMM_TOL);
		cats.push_back(make_pair(string("Size classification"), f));

		f.clear();
		f.push_back(ARG_FEATURES);    f.push_back(ARG_FEATURE_TPED);
		f.push_back(ARG_FEATURE_TFAM); f.push_back(ARG_FEATURE_VCF);
		f.push_back(ARG_ROH_FILE);
		cats.push_back(make_pair(string("Counting classified genotypes"), f));

		f.clear();
		f.push_back(ARG_OUTFILE);     f.push_back(ARG_OUTDIR);
		f.push_back(ARG_FROH);        f.push_back(ARG_FROH_DENOM);
		f.push_back(ARG_FREQ_ONLY);   f.push_back(ARG_RAW_LOD);
		f.push_back(ARG_FORCE);
		cats.push_back(make_pair(string("Output"), f));

		f.clear();
		f.push_back(ARG_THREADS);     f.push_back(ARG_SEED);
		f.push_back(ARG_RESAMPLE);    f.push_back(ARG_LOAD_PARAMS);
		f.push_back(ARG_QUIET);       f.push_back(ARG_VERBOSE);
		f.push_back(ARG_VERSION);     f.push_back(ARG_HELP);
		f.push_back(ARG_HELP_SHORT);  f.push_back(ARG_DUMP_DOCS);
		cats.push_back(make_pair(string("Run control"), f));

		params->setHelpCategories(cats);
	}

	//A bare invocation is a usage error, not a successful no-op run.
	if (argc < 2)
	{
		params->printHelp();
		delete params;
		status = PARAM_ERROR;
		return NULL;
	}

	status = params->parseCommandLine(argc, argv);

	if (status == PARAM_OK && params->isFlagSet(ARG_DUMP_DOCS))
	{
		bool ok = params->writeHelpDoc(cout, params->getStringFlag(ARG_DUMP_DOCS));
		delete params;
		status = ok ? PARAM_HELP : PARAM_ERROR;
		return NULL;
	}

	//Applied after parsing so the command line wins, and before main reads
	//any value.
	if (status == PARAM_OK && params->isFlagSet(ARG_LOAD_PARAMS))
	{
		if (!params->loadFlagsJSON(params->getStringFlag(ARG_LOAD_PARAMS)))
		{
			delete params;
			status = PARAM_ERROR;
			return NULL;
		}
	}

	if (status == PARAM_OK && params->getBoolFlag(ARG_VERSION))
	{
		cout << "garlic v" << VERSION << " (" << GARLIC_GIT_SHA << ")\n";
		delete params;
		status = PARAM_HELP;
		return NULL;
	}

	if (status != PARAM_OK)
	{
		delete params;
		return NULL;
	}
	return params;
}

void writeParamsJSON(string file, param_t *params, vector< pair<string,string> > &resolved)
{
	//ios::binary: see writeROHData in garlic-roh.cpp.
	ofstream out(file.c_str(), ios::binary);
	if (out.fail())
	{
		LOG.err("ERROR: Failed to open", file);
		return;
	}
	out << "{\n";
	out << "  \"garlic_version\": \"" << VERSION << "\",\n";
	out << "  \"git_sha\": \"" << GARLIC_GIT_SHA << "\",\n";
	out << "  \"resolved\": {\n";
	for (unsigned int i = 0; i < resolved.size(); i++)
	{
		out << "    \"" << resolved[i].first << "\": " << resolved[i].second;
		if (i + 1 < resolved.size()) out << ",";
		out << "\n";
	}
	out << "  },\n";
	vector<string> given = params->setFlags();
	out << "  \"set\": [";
	for (unsigned int i = 0; i < given.size(); i++)
	{
		if (i) out << ", ";
		out << "\"" << given[i] << "\"";
	}
	out << "],\n";
	out << "  \"flags\": {\n";
	params->writeFlagsJSON(out, "    ");
	out << "  }\n";
	out << "}\n";
	out.close();
	LOG.log("Effective parameters:", file);
	return;
}

//Refuse to clobber a previous run's calls unless asked to.
bool checkOutfileClobber(string outfile, bool force, string suffix)
{
	if (force) return false;
	string roh = outfile + suffix;
	ifstream probe(roh.c_str());
	if (probe.good())
	{
		probe.close();
		LOG.err("ERROR: Output file already exists:", roh);
		LOG.err("\tPass --force to overwrite, or choose another --out.");
		return true;
	}
	return false;
}

bool checkSeed(int seed){
	if(seed < 0){
		LOG.err("ERROR: Random seed must be >= 0 (0 means draw one automatically).");
		return true;
	}
	return false;
}

bool checkPopFile(string popfile, string tpedfile, string vcffile, string beaglefile){
	if (popfile.compare(DEFAULT_POP) == 0) return false;
	if (tpedfile.compare(DEFAULT_TPED) == 0 && vcffile.compare(DEFAULT_VCF) == 0
	        && beaglefile.compare(DEFAULT_BEAGLE) == 0){
		LOG.err("ERROR: --pop needs input data; give --tped, --vcf or --gl-beagle.");
		return true;
	}
	return false;
}

bool checkMaxWinsize(int maxWinsize, int winsize){
	if(maxWinsize < winsize){
		LOG.err("ERROR: --max-winsize must be >= the starting --winsize.");
		return true;
	}
	return false;
}

bool checkCM(string mapfile, bool CM){
	if(CM && mapfile.compare(DEFAULT_MAP) == 0){
		LOG.err("ERROR: Must provide mapfile if you wish to construct ROH in genetic map units.");
		return true;
	}
	else return false;
}
bool checkNCLUST(int nclust){
	if(nclust <= 0){
		LOG.err("ERROR: Must choose positive number for number of GMM clusters.");
		return true;
	}
	else return false;
}

bool checkM(int M){
	if(M <= 0){
		LOG.err("ERROR: M must be an integer > 0.");
		return true;
	}
	else return false;
}

bool checkMU(double mu){
	if(mu <= 0 || mu >= 1){
		LOG.err("ERROR: mu must be between 0 and 1.");
		return true;
	}
	else return false;
}


bool checkBuild(string BUILD)
{
	if (BUILD.compare("hg18") != 0 &&
	        BUILD.compare("hg19") != 0 &&
	        BUILD.compare("hg38") != 0 &&
	        BUILD.compare("t2t-chm13") != 0 &&
	        BUILD.compare(DEFAULT_BUILD) != 0)
	{
		LOG.err("ERROR: Must choose hg18, hg19, hg38 or t2t-chm13 for build version.");
		return true;
	}
	return false;
}


bool checkBuildAndCentromereFile(string BUILD, string centromereFile, bool NO_CENTROMERE) {
	bool haveBuild = (BUILD.compare(DEFAULT_BUILD) != 0);
	bool haveFile = (centromereFile.compare(DEFAULT_CENTROMERE_FILE) != 0);

	if (NO_CENTROMERE)
	{
		//Saying both "there is no gap" and "here is the gap" is a contradiction,
		//so reject it rather than silently preferring one.
		if (haveBuild || haveFile) {
			LOG.err("ERROR: --no-centromere cannot be combined with --build or --centromere.");
			return true;
		}
		return false;
	}

	if (!haveBuild && !haveFile) {
		LOG.err("ERROR: Must choose hg18/hg19/hg38/t2t-chm13 for build version, provide a custom");
		LOG.err("ERROR: centromere file, or pass --no-centromere if the assembly has no gap.");
		return true;
	}
	return false;
}

bool checkMultiWinsizes(vector<int> &multiWinsizes, bool &WINSIZE_EXPLORE, bool wasSet)
{
	if (wasSet)
	{
		for (unsigned int i = 0; i < multiWinsizes.size(); i++)
		{
			if (multiWinsizes[i] <= 0)
			{
				//cerr << "ERROR: SNP window sizes must be > 1.\n";
				LOG.err("ERROR: SNP window sizes must be > 1.");
				return true;
			}
		}
		WINSIZE_EXPLORE = true;
	}
	return false;
}

bool checkAutoFreq(string freqfile, bool FREQ_ONLY, bool &AUTO_FREQ)
{
	if (freqfile.compare(DEFAULT_FREQ_FILE) != 0)
	{
		AUTO_FREQ = false;
		if (FREQ_ONLY)
		{
			//cerr << "ERROR: Specifying both " << ARG_FREQ_ONLY << " and " << ARG_FREQ_FILE << " accomplishes nothing useful.\n";
			LOG.err("ERROR: Specifying both", ARG_FREQ_ONLY, false);
			LOG.err(" and", ARG_FREQ_FILE, false);
			LOG.err(" accomplishes nothing useful.");
			return true;
		}
	}
	return false;
}

bool checkAutoWinsizeStep(int auto_winsize_step) {
	if (auto_winsize_step <= 0) {
		LOG.err("ERROR: Step size for automatic window selection must be positive.");
		return true;
	}
	return false;
}

bool checkAutoWinsize(bool WINSIZE_EXPLORE, bool AUTO_WINSIZE)
{
	//Check if both AUTO_WINSIZE and WINSIZE_EXPLORE are set
	//If so, exit with error.
	if (WINSIZE_EXPLORE && AUTO_WINSIZE)
	{
		//cerr << "ERROR: Must set only one of " << ARG_WINSIZE_MULTI << " and " << ARG_AUTO_WINSIZE << ".\n";
		LOG.err("ERROR: Must set only one of", ARG_WINSIZE_MULTI, false);
		LOG.err(" and", ARG_AUTO_WINSIZE);
		return true;
	}
	return false;
}

bool checkAutoCutoff(double LOD_CUTOFF, bool &AUTO_CUTOFF, bool wasSet)
{
	//Was keyed on LOD_CUTOFF != -999999, so -999999 was an unusable value.
	if (wasSet) {
		AUTO_CUTOFF = false;
		//A window that must not be called carries MISSING (-9999), and the
		//assembler decides by `window >= cutoff`.  A cutoff at or below that
		//sentinel therefore calls the windows whose whole purpose is to be
		//uncallable: the tails of every chromosome, the windows spanning an
		//assembly gap, and -- since the sex-chromosome work -- every window
		//of an individual that is not diploid there.  No real LOD score is
		//anywhere near it.
		if (LOD_CUTOFF <= double(MISSING)) {
			LOG.err("ERROR:", ARG_LOD_CUTOFF, false);
			LOG.err(" must be greater than", MISSING, false);
			LOG.err(", which is the value that marks a window as uncallable.");
			return true;
		}
	}
	return false;
}

bool checkFROHDenominator(const string &mode, int &denom, bool CM, bool FROH,
                          bool chrLengthsSet)
{
	if (mode.compare("analyzed") == 0)      denom = FROH_ANALYZED;
	else if (mode.compare("assembly") == 0) denom = FROH_ASSEMBLY;
	else
	{
		LOG.err("ERROR:", ARG_FROH_DENOM, false);
		LOG.err(" must be analyzed or assembly, not", mode, false);
		LOG.err(".");
		return true;
	}

	//Nothing to divide: --froh-denominator describes a file that is not being
	//written, and silently accepting it would leave the user believing they
	//had asked for something.
	if (denom == FROH_ASSEMBLY && !FROH)
	{
		LOG.err("ERROR:", ARG_FROH_DENOM, false);
		LOG.err(string(" assembly describes ") + ARG_FROH + ", which was not given.");
		return true;
	}

	//A chromosome length is a number of base pairs and the denominator under
	//--cm is a number of centimorgans.  The equivalent would be a total map
	//length per chromosome, which is a different table and a different flag.
	if (denom == FROH_ASSEMBLY && CM)
	{
		LOG.err("ERROR:", ARG_FROH_DENOM, false);
		LOG.err(string(" assembly cannot be used with ") + ARG_CM + ": a chromosome's length in");
		LOG.err("\tbase pairs is not its length in centimorgans.");
		return true;
	}

	//The file is a source of lengths, not a mode.  Reading it and then not
	//using it would be the quieter failure.
	if (chrLengthsSet && denom != FROH_ASSEMBLY)
	{
		LOG.err("ERROR:", ARG_CHR_LENGTHS, false);
		LOG.err(string(" supplies the lengths that ") + ARG_FROH_DENOM + " assembly divides by,");
		LOG.err("\tand nothing else reads them. Give that too, or drop this.");
		return true;
	}
	return false;
}

bool checkSexChrEstimates(double cutoff, bool cutoffSet, int winsize, bool winsizeSet,
                          bool WINSIZE_EXPLORE, bool FREQ_ONLY)
{
	if (!cutoffSet && !winsizeSet) return false;

	const string which = (cutoffSet && winsizeSet)
	                     ? (ARG_SEXCHR_LOD_CUTOFF + " and " + ARG_SEXCHR_WINSIZE)
	                     : (cutoffSet ? ARG_SEXCHR_LOD_CUTOFF : ARG_SEXCHR_WINSIZE);

	if (FREQ_ONLY) {
		LOG.err("ERROR:", which, false);
		LOG.err(string(" concern ROH calling and ") + ARG_FREQ_ONLY + " calls none.");
		return true;
	}

	//The exploration reports a table over window sizes, and every row of it
	//would need its own sex-chromosome cutoff to mean anything.
	if (WINSIZE_EXPLORE) {
		LOG.err("ERROR:", which, false);
		LOG.err(string(" cannot be used with ") + ARG_WINSIZE_MULTI + ".");
		return true;
	}

	if (cutoffSet && cutoff <= double(MISSING)) {
		LOG.err("ERROR:", ARG_SEXCHR_LOD_CUTOFF, false);
		LOG.err(" must be greater than", MISSING, false);
		LOG.err(", which is the value that marks a window as uncallable.");
		return true;
	}

	if (winsizeSet) {
		if (winsize <= 0) {
			LOG.err("ERROR:", ARG_SEXCHR_WINSIZE, false);
			LOG.err(" must be a positive number of SNPs.");
			return true;
		}
		//A LOD score cutoff belongs to the window size it was estimated at:
		//more SNPs in a window means more terms in the sum, and the scores
		//move with it.  Carrying the autosomal cutoff onto windows of another
		//size would give a plausible-looking but meaningless set of calls, so
		//it is refused rather than warned about.
		if (!cutoffSet) {
			LOG.err("ERROR:", ARG_SEXCHR_WINSIZE, false);
			LOG.err(string(" changes what a LOD score on that chromosome means, so the "
			               "autosomal cutoff no longer applies to it. Give ") +
			        ARG_SEXCHR_LOD_CUTOFF + " as well.");
			return true;
		}
	}
	return false;
}

bool checkBoundSizes(vector<double> &boundSizes, bool &AUTO_BOUNDS, bool wasSet){

	if(!wasSet){
		return false;
	}
	else {
		AUTO_BOUNDS = false;

		for(unsigned int i = 0; i < boundSizes.size(); i++){		
			if (boundSizes[i] <= 0){
				LOG.err("ERROR: User provided size boundaries must be positive.");
				return true;
			}
			if(i > 0){
				if (boundSizes[i] <= boundSizes[i-1]){
					LOG.err("ERROR: User provided size boundaries must be in strictly increasing order.");
					return true;
				}
			}
		}
	}
	return false;
}

bool checkRequiredFiles(string tpedfile, string tfamfile, string vcffile, string tglsfile, string beaglefile,
                        string rohfile)
{
	bool haveTped = (tpedfile.compare(DEFAULT_TPED) != 0);
	bool haveVcf  = (vcffile.compare(DEFAULT_VCF)   != 0);
	bool haveTfam = (tfamfile.compare(DEFAULT_TFAM) != 0);
	bool haveRoh  = (rohfile.compare(DEFAULT_ROH_FILE) != 0);
	bool haveBeagle = (beaglefile.compare(DEFAULT_BEAGLE) != 0);

	if ((haveTped && haveVcf) || (haveBeagle && (haveTped || haveVcf)))
	{
		LOG.err("ERROR: --tped, --vcf and --gl-beagle are alternative sources of the same data; give one.");
		return true;
	}
	if (haveBeagle)
	{
		//Genotype likelihoods carry no sample metadata, and a .tgls file would
		//be a second, incompatible statement of per-genotype uncertainty.
		if (haveRoh)
		{
			LOG.err("ERROR:", ARG_ROH_FILE, false);
			LOG.err(" reads calls from a file, so no ROH are called and", ARG_BEAGLE, false);
			LOG.err(" is not used.");
			return true;
		}
		if (!haveTfam)
		{
			LOG.err("ERROR:", ARG_BEAGLE, false);
			LOG.err(" needs --tfam for individual IDs, population labels and sex, one row per");
			LOG.err("\tbeagle individual in the same order.");
			return true;
		}
		if (tglsfile.compare(DEFAULT_TGLS) != 0)
		{
			LOG.err("ERROR: --tgls and --gl-beagle are alternative statements of genotype uncertainty; give one.");
			return true;
		}
		return false;
	}
	//--roh-file calls nothing, so the calling input is not merely unnecessary
	//there: giving it would say the run should call ROH and read them from a
	//file at the same time.
	if (haveRoh && (haveTped || haveVcf))
	{
		LOG.err("ERROR:", ARG_ROH_FILE, false);
		LOG.err(" reads calls from a file, so no ROH are called and", haveTped ? ARG_TPED : ARG_VCF, false);
		LOG.err(" is not used.");
		LOG.err("\tGive the genotypes to count with --tped-counting/--tfam-counting or --vcf-counting.");
		return true;
	}
	if (haveRoh) return false;
	if (!haveTped && !haveVcf)
	{
		LOG.err("ERROR: Must provide genotypes: --tped with --tfam, --vcf, or --gl-beagle with --tfam.");
		return true;
	}
	if (haveVcf && haveTfam)
	{
		//A VCF names its own samples, so a TFAM could only disagree with it.
		LOG.err("ERROR: --vcf takes sample names from its #CHROM line, so --tfam is not used.");
		LOG.err("\tUse --pop to supply population labels.");
		return true;
	}
	if (haveTped && !haveTfam)
	{
		LOG.err("ERROR: Must provide both a tped and a tfam file.");
		return true;
	}
	if (haveVcf && tglsfile.compare(DEFAULT_TGLS) != 0)
	{
		//Not merely redundant: reading PL or GL from a .tgls file cannot be
		//correct, because that format holds ONE value per genotype and a VCF
		//normalises the called genotype's value to exactly 0.  With --vcf the
		//FORMAT column has the whole array, which is the point of --gl-type
		//there.  Silently ignoring --tgls would leave the user believing the
		//file was used.
		LOG.err("ERROR: --tgls and --vcf are alternative sources of genotype qualities.");
		LOG.err("\tWith --vcf, --gl-type reads GQ/PL/GL from the FORMAT column instead.");
		return true;
	}
	return false;
}

bool checkMapFile(string mapfile, bool WEIGHTED){
	if(mapfile.compare(DEFAULT_MAP) == 0 && WEIGHTED){
		LOG.err("ERROR: Weighted LOD score method requires a map file.");
		return true;
	}
	return false;
}

bool checkThreads(int numThreads)
{
	if (numThreads <= 0)
	{
		//cerr << "ERROR: Number of threads must be > 0.\n";
		LOG.err("ERROR: Number of threads must be > 0.");
		return true;
	}
	return false;
}

bool checkError(double error, string tglsfile, bool wasSet, bool haveVCFLikelihoods)
{
	if (!wasSet)
	{
		//--vcf with --gl-type reads per-genotype error rates from the FORMAT
		//column, so it is a third source of them alongside --error and --tgls.
		if (tglsfile.compare(DEFAULT_TGLS) == 0 && !haveVCFLikelihoods) {
			LOG.err("ERROR: --error must be given, or per-genotype error rates must come");
			LOG.err("ERROR: from a TGLS file (--tgls) or a VCF FORMAT field (--vcf --gl-type).");
			return true;
		}
		return false;
	}
	if (error <= 0 || error >= 1)
	{
		LOG.err("ERROR: Genotype error rate must be > 0 and < 1.");
		return true;
	}
	return false;
}

//mkdir -p, so --outdir can name a nested path.
//Is this a path component there is nothing to create for -- the root, . or
//.., or a bare drive letter such as C:?  mkdir fails on each with a misleading
//message.  Applied to already-trimmed components.
static bool isPathStub(const string &p, bool (*isSep)(char)){
	if(p.empty() || p == "." || p == "..") return true;
	if(p.size() == 1 && isSep(p[0])) return true;
	if(p.size() == 2 && p[1] == ':') return true;
	return false;
}

bool makeOutdir(string dir){
	if(dir.empty()) return false;
	string partial;
	for(unsigned int i = 0; i < dir.size(); i++){
		partial += dir[i];
		if(garlicIsPathSep(dir[i]) || i + 1 == dir.size()){
			string p = partial;
			if(p.size() > 1 && garlicIsPathSep(p[p.size()-1])) p.erase(p.size()-1);
			if(isPathStub(p, garlicIsPathSep)) continue;
			if(garlicMkdir(p.c_str()) != 0 && errno != EEXIST){
				LOG.err("ERROR: Could not create output directory:", p);
				return true;
			}
			//EEXIST alone does not mean the path is usable: a regular file of
			//that name gives EEXIST as well, and garlic then aborted when the
			//writers could not open their outputs (exit 134, no diagnosis).
			if(!garlicIsDir(p.c_str())){
				LOG.err("ERROR: Output path exists but is not a directory:", p);
				return true;
			}
		}
	}

	//Postcondition, checked once on the whole path rather than inferred from the
	//loop above.  The loop's per-component checks depend on mkdir's and stat's
	//errno conventions agreeing across platforms, and on MinGW a path whose
	//parent is a regular file still came back as success.  What actually has to
	//be true when this returns is simply that `dir` is now a usable directory,
	//so test exactly that.
	string full = dir;
	while(full.size() > 1 && garlicIsPathSep(full[full.size()-1])) full.erase(full.size()-1);
	if(!isPathStub(full, garlicIsPathSep) && !garlicIsDir(full.c_str())){
		LOG.err("ERROR: Output directory could not be created:", full);
		return true;
	}
	return false;
}

bool checkAutoWinsizeThreshold(double t){
	if(t <= 0){ LOG.err("ERROR: --auto-winsize-threshold must be > 0."); return true; }
	return false;
}

bool checkKDEPoints(int m){
	if(m < 16){ LOG.err("ERROR: --kde-points must be >= 16."); return true; }
	return false;
}

bool checkKDECut(double c){
	if(c <= 0){ LOG.err("ERROR: --kde-cut must be > 0."); return true; }
	return false;
}

bool checkModeSpan(int s){
	if(s < 2){ LOG.err("ERROR: --mode-smooth-span must be >= 2."); return true; }
	return false;
}

bool checkCoefPair(vector<double> &coef, string flag){
	if(coef.size() != 2){
		LOG.err("ERROR: " + flag + " takes exactly two values: <slope> <intercept>.");
		return true;
	}
	return false;
}

bool checkGMMParams(int maxIter, double tol){
	if(maxIter < 1){ LOG.err("ERROR: --gmm-max-iter must be >= 1."); return true; }
	if(tol <= 0){ LOG.err("ERROR: --gmm-tol must be > 0."); return true; }
	return false;
}

bool checkKDEThinStep(int step){
	if(step < 0){
		LOG.err("ERROR: --kde-thin-step must be >= 0 (0 means use the window size).");
		return true;
	}
	return false;
}

//--size-bounds and --nclust are alternatives: supplying bounds skips the GMM
//that --nclust configures.  This used to be silent.
void warnBoundsOverridesNclust(bool boundsSet, bool nclustSet){
	if(boundsSet && nclustSet){
		LOG.err("WARNING: --size-bounds was given, so --nclust is ignored.");
		LOG.err("\tThe GMM that --nclust configures only runs when size boundaries are");
		LOG.err("\tchosen automatically.");
	}
}

bool checkGLType(string TYPE, string tglsfile, string beaglefile)
{
	if ( TYPE.compare("GQ") != 0 && TYPE.compare("GL") != 0 && TYPE.compare("PL") != 0 && tglsfile.compare(DEFAULT_TGLS) != 0 ) {
		LOG.err("ERROR: Must choose GQ/GL/PL for genotype likelihood format or provide a single error rate with --error.");
		return true;
	}
	//A beagle file IS the likelihood; --gl-type names a per-genotype quality
	//form that would be read from somewhere else.  Refuse rather than ignore.
	if ( beaglefile.compare(DEFAULT_BEAGLE) != 0 && TYPE.compare(DEFAULT_GL_TYPE) != 0 ) {
		LOG.err("ERROR: --gl-type has no meaning with --gl-beagle; the triplet is the likelihood.");
		return true;
	}
	return false;
}

bool checkWinsize(int winsize, bool WINSIZE_EXPLORE, bool AUTO_WINSIZE, bool WEIGHTED, bool FREQ_ONLY)
{
	//--freq-only never forms a window, so do not make the user invent a
	//--winsize just to get past this check.
	if (FREQ_ONLY) return false;
	if (winsize <= 1){
		if(!WINSIZE_EXPLORE && !(AUTO_WINSIZE && WEIGHTED)){
			LOG.err("ERROR: SNP window size must be > 1. If using --auto-winsize, this is the starting value.");
			return true;
		}
	}
	return false;
}

bool checkMaxGap(int MAX_GAP)
{
	if (MAX_GAP < 0){
		LOG.err("ERROR: Max gap must be > 0.");
		return true;
	}
	else if (MAX_GAP < 1000){
		LOG.err("WARNING: max gap set very low:", MAX_GAP);
	}
	return false;
}

bool checkOverlapFrac(double OVERLAP_FRAC){
	if(OVERLAP_FRAC < 0 || OVERLAP_FRAC > 1){
		LOG.err("ERROR: Overlap fraction must be >= 0 and <= 1.");
		return true;
	}
	return false;
}