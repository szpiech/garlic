#!/usr/bin/perl -w 

use strict;

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




