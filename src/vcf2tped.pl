#!/usr/bin/perl -w

# DEPRECATED.  garlic reads a VCF directly:
#
#     garlic --vcf in.vcf.gz --pop pops.txt ...
#
# and that path is not merely more convenient, it is more correct.  This script
# writes a tped plus a tgls file, and a tgls file holds ONE genotype-quality
# value per genotype.  A VCF normalises PL and GL so the CALLED genotype's
# value is exactly 0, so that single value carries no information -- it is 0
# for every call, confident or not.  garlic now rejects such a file rather than
# computing with it (it used to produce a LOD of -16 per heterozygote, or, for
# an all-integer PL file, call nothing at all).
#
# --vcf --gl-type PL/GL reads the whole FORMAT array and computes the
# posterior, which is the only correct use of a normalised PL.  --gl-type GQ is
# the one form this script's output can still be used with, because GQ needs no
# normalisation.
#
# It also writes the population column as a copy of the sample ID, which is
# what --pop exists to fix.
#
# Kept so existing workflows keep running.  It will print this notice on every
# invocation.

use strict;

print STDERR "\n";
print STDERR "WARNING: vcf2tped.pl is DEPRECATED.  garlic reads a VCF directly:\n";
print STDERR "WARNING:   garlic --vcf $ARGV[0] --pop <popfile> ...\n" if defined $ARGV[0];
print STDERR "WARNING: The tgls file this writes holds one value per genotype, which cannot\n";
print STDERR "WARNING: represent a VCF's PL or GL (those are normalised to 0 for the called\n";
print STDERR "WARNING: genotype).  garlic rejects such a file.  Use --vcf --gl-type PL instead;\n";
print STDERR "WARNING: --gl-type GQ is the one form this output still works with.\n";
print STDERR "WARNING: See the README for details.\n";
print STDERR "\n";


if ($#ARGV < 0){
	print STDERR "./vcf2tped.pl <vcf file>\n" and die;
}

my $vcffile = $ARGV[0];
my $tpedfile = $vcffile;
my $tfamfile = $vcffile;

unless($tpedfile =~ s/\.vcf/\.tped/g){
	print STDERR "$vcffile not recognized as a vcf file based on file name extension.\n" and die;
}

unless($tfamfile =~ s/\.vcf/\.tfam/g){
	print STDERR "$vcffile not recognized as a vcf file based on file name extension.\n" and die;
}


open(FIN,"<",$vcffile) or die "$vcffile $!";
open(TPED, ">", $tpedfile) or die "$tpedfile $!";
while(defined(my $line = <FIN>)){
	chomp $line;
	if($line =~ m/^#CHROM/){
		open(TFAM, ">",$tfamfile) or die "$tfamfile $!";
		my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split(/\s+/,$line);
		for my $id (@data){
			print TFAM "0\t$id\t0\t0\t0\t0\n";
		}
		close(TFAM);
		next;		
	}
	elsif($line =~ m/^#/){
		next;
	}

	my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split(/\s+/,$line);

	print TPED "$CHROM\t$ID\t0\t$POS";
	for my $dat (@data){
		if($dat =~ m/((\d|\.)(\/|\|)(\d|\.))(:.+)?/){
			if($2 eq '0'){
				print TPED "\t$REF";
			}
			elsif($2 eq '1'){
				print TPED "\t$ALT";
			}
			elsif($2 eq '.'){
				print TPED "\t0";	
			}

			if($4 eq '0'){
				print TPED "\t$REF";
			}
			elsif($4 eq '1'){
				print TPED "\t$ALT";
			}
			elsif($4 eq '.'){
				print TPED "\t0";	
			}

		}
	}
	print TPED "\n";
}
close(FIN);
close(TPED);




