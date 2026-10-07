"""Replace bcftools' prior-based GT by the maximum-likelihood genotype from PL, and report how often they differ.

bcftools call -m multiplies each sample's genotype likelihoods by a population prior; at 1-2 reads the
likelihoods barely separate a homozygote from a heterozygote, so the prior decides and intermediate-frequency
sites are called heterozygous even where the reads favour a homozygote.  A het called that way inside a true
ROH costs log10(eps)-scale LOD.  Here GT := argmin(PL) (ties -> missing), DP by AD sum if present.
usage: python ml_genotypes.py in.vcf.gz out.vcf.gz stats.tsv"""
import gzip, sys
from collections import Counter
inp, out, stats = sys.argv[1:4]
cnt = Counter()
GT = {0: '0/0', 1: '0/1', 2: '1/1'}
with gzip.open(inp, 'rt') as f, gzip.open(out, 'wt') as o:
    for line in f:
        if line.startswith('#'):
            o.write(line); continue
        p = line.rstrip('\n').split('\t'); fmt = p[8].split(':')
        if 'GT' not in fmt or 'PL' not in fmt or ',' in p[4]:
            o.write(line); continue
        gi, pi = fmt.index('GT'), fmt.index('PL'); ai = fmt.index('AD') if 'AD' in fmt else -1
        for j in range(9, len(p)):
            s = p[j].split(':')
            gt = s[gi].replace('|', '/')
            try: pl = [int(v) for v in s[pi].split(',')]
            except (ValueError, IndexError): cnt[('none', gt, 'missingPL')] += 1; continue
            try: dp = sum(int(v) for v in s[ai].split(',')) if ai >= 0 else -1
            except ValueError: dp = -1
            m = min(pl); arg = [k for k, v in enumerate(pl) if v == m]
            ml = GT[arg[0]] if len(arg) == 1 else './.'
            dpc = str(dp) if 0 <= dp <= 3 else ('4+' if dp > 3 else 'NA')
            cnt[(dpc, gt, ml)] += 1
            s[gi] = ml; p[j] = ':'.join(s)
        o.write('\t'.join(p) + '\n')
with open(stats, 'w') as fs:
    fs.write('reads\tbcftools_gt\tml_gt\tn\n')
    for (a, b, c), n in sorted(cnt.items()): fs.write(f'{a}\t{b}\t{c}\t{n}\n')
