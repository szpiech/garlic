**NOTE TO WINDOWS USERS: In order to successfully run GARLIC you must have ann_figtree_version.dll and figtree.dll in the same folder as garlic.exe 

This README gives an overview of basic commands for running GARLIC with the example data provided.

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --auto-winsize --winsize 30 --out example --error 0.001

The example data is derived from human genotype chip data with hg18 coordinates (--build hg18), and a priori we do not know how large to make our window so we use the built-in window size selection algorithm (--auto-winsize).  --auto-winsize searches upward from the value given by --winsize, so a starting size is required; earlier versions of this file omitted it, and the command failed.  All output files will be named starting with "example" (--out example).  On this data the search selects a window of 50 SNPs, so this produces:

example.log
example.50SNPs.kde
example.freq.gz
example.roh.bed
example.params.json

example.error is created only if something is actually written to it, so a clean run does not produce one.

example.params.json records the effective value of every command line flag together with the values the run settled on -- the selected window size, LOD score cutoff, size class boundaries and random seed.  Passing it back with --load-params repeats the run; flags given on the command line override the file, so you can repeat a run with one parameter changed.

The output files shipped alongside this README come from the --winsize 60 command further down, not from the --auto-winsize command above, which is why the KDE file here is example.60SNPs.kde rather than example.50SNPs.kde.  Automatic LOD cutoff selection is deterministic, so these files reproduce exactly: cutoff 1.88468, size class boundaries 542417 and 1.77283e+06, 15601 ROH.

If the automatic window size selection algorithm fails, you can output the KDEs of the LOD score distribution for multiple window sizes (without calling ROH) by using the --winsize-multi argument, i.e.

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize-multi 30 40 50 60 70 80 90 --out example --error 0.001

This will generate KDEs for your inspection.  Once you've chosen a window size, you should rerun garlic specifying that size.  For example, if you choose a window size of 60 SNPs, then you would run

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize 60 --out example --error 0.001

If you already know what LOD score cutoff to use (say you are analyzing more individuals from a previously studies population), you can use the --lod-cutoff argument, i.e. if your LOD score cutoff is known to be 2.5 then you would run

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize 60 --out example --lod-cutoff 2.5 --error 0.001

If you already know what size thresholds to use for size classificaation (say you are analyzing more individuals from a previously studies population), you can use the --size-bounds argument, i.e. if your size thresholds are known to be 500000 and 1000000 for the boundaries between short/med and med/long, respectively, then you would run

garlic --tped example.tped.gz --tfam example.tfam --build hg18 --winsize 60 --out example --lod-cutoff 2.5 --size-bounds 500000 1000000 --error 0.001

ROH calls land in example.roh.bed.  Following the BED specification the start coordinate is 0-based and the end coordinate is exclusive, so column 3 minus column 2 equals the ROH length in column 5.  Versions before this one wrote a 1-based start, so coordinates from older output are shifted by one base relative to these.

Adding --froh writes example.froh.tsv as well: one row per individual per size class, giving the ROH count, the total autozygous length and FROH.  The denominator is stated in the file's header, because there is no single conventional choice and values are not comparable across denominators.

Other command line arguments are listed by --help, in the README, and in the manual.

To run the wLOD, you must provide a map file and give the --weighted flag:

garlic --tped example.tped.gz --tfam example.tfam --map example.map.gz --weighted --build hg18 --winsize 60 --out example --error 0.001 
