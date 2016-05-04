This README gives an overview of basic commands for running GARLIC with the example data provided.

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --auto-winsize --out example

The example data is derived from human genotype chip data with hg18 coordinates (--build hg18), and a priori we do not know how large to make our window so we use the built-in window size selection algorithm (--auto-winsize).  All output files will be named starting with "example" (--out example).  This produces the following files:

example.error
example.log
example.60SNPs.kde
example.freq.gz
example.roh.bed

If the automatic window size selection algorithm fails, you can output the KDEs of the LOD score distribution for multiplw window sizes (without calling ROH) by using the --winsize-multi argument, i.e.

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize-multi 30 40 50 60 70 80 90 --out example

This will generate KDEs for your inspection.  Once you've chosen a window size, you should rerun garlic specifying that size.  For example, if you choose a window size of 60 SNPs, then you would run

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize 60 --out example

If you already know what LOD score cutoff to use (say you are analyzing more individuals from a previously studies population), you can use the --lod-cutoff argument, i.e. if your LOD score cutoff is known to be 2.5 then you would run

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize 60 --out example --lod-cutoff 2.5

If you already know what size thresholds to use for size classificaation (say you are analyzing more individuals from a previously studies population), you can use the --size-bounds argument, i.e. if your size thresholds are known to be 500000 and 1000000 for the boundaries between short/med and med/long, respectively, then you would run

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize 60 --out example --lod-cutoff 2.5 --size-bounds 500000 1000000

Other command line arguments are described in the manual.
