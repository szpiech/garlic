#!/usr/bin/perl -w 

use strict;

my $file = "centromere_hg38.txt";

my %start;
my %end;

open(FIN,"<",$file) or die $!;
for my $line (<FIN>){
	chomp $line;
	if (! ($line =~ m/^#/) ){
		my ($bin,$chrom,$start,$end,$name) = split(/\s+/,$line);
		if(exists $start{$chrom}){
			if($start{$chrom} > $start){
				$start{$chrom} = $start;
			}
		}
		else{
			$start{$chrom} = $start;
		}

		if(exists $end{$chrom}){
			if($end{$chrom} < $end){
				$end{$chrom} = $end;
			}
		}
		else{
			$end{$chrom} = $end;
		}
	}
}
close(FIN);

for my $chrom (keys %start){
	print "gapStart[\"${chrom}\"] = ",$start{$chrom},";\n";
	print "gapEnd[\"${chrom}\"] = ",$end{$chrom},";\n";
}