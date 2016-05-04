#!/usr/bin/perl -w

use strict;
use POSIX qw/floor/;
use POSIX qw/ceil/;

my $startchr = 22;

if ( $#ARGV < 4 ) {
    print STDERR
        "./count_features_in_roh.pl <feature file> <roh file> <tped/vcf file> <num chr> <out file>\n";
    print STDERR
        "If a TPED file is given, a TFAM file is also expected with the same basename.\n";
    print STDERR
        "Support for gzipped files depends on having gunzip installed.\n";
    die;
}

my $FEATUREFILE = $ARGV[0];
my $ROHFILE     = $ARGV[1];
my $GENFILE     = $ARGV[2];
my $numchr      = $ARGV[3];
my $OUTFILE     = $ARGV[4];
my $FILETYPE;
my $TFAMFILE;
my $ISZIP = 0;
my @genfilelist;

if ( $GENFILE =~ m/.*\.vcf(\.gz)?$/ ) {
    $FILETYPE = "VCF";
    if ( defined $1 ) {
        $ISZIP = 1;
    }
}
elsif ( $GENFILE =~ m/.*\.tped(\.gz)?$/ ) {
    $FILETYPE = "TPED";
    if ( defined $1 ) {
        $ISZIP = 1;
    }
    $TFAMFILE = $GENFILE;
    $TFAMFILE =~ s/\.tped/\.tfam/g;
}
else {
    print STDERR
        "ERROR: $GENFILE not recognzied as vcf or tped (based on file name extension).\n";
    die;
}

if ( $GENFILE =~ m/(.+)?chr\d+(.+)?/ ) {
    my $front = $1;
    my $back  = $2;

    for ( my $i = $startchr; $i <= $numchr; $i++ ) {
        my $filename = "${front}chr${i}${back}";
        $filename =~ s/^\s+//g;
        $filename =~ s/\s+$//g;
        push( @genfilelist, $filename );
    }
}
else {
    print STDERR 'ERROR: Could not match pattern m/(.+)?chr\d+(.+)?/ to ',
        $GENFILE, ' to generate file list.', "\n";
    print STDERR "\t",
        'Please name your files *chr{num}*, i.e. data.chr1.tped, data.chr2.tped, etc.',
        "\n"
        and die;
}

my %effect;
my %effectTypes;
my @effectList;
open( FIN, "<", $FEATUREFILE ) or die $!;
while ( defined( my $line = <FIN> ) ) {
    chomp $line;
    my ( $chrpos, $ref, $alt, $effect ) = split( /\s+/, $line );
    my ( $chr, $pos ) = split( /:/, $chrpos );

    #print $chr, " ", $pos, "\n";
    $effect{$chr}{$pos}{$alt} = $effect;
    $effectTypes{$effect} = 1;
}
close(FIN);

@effectList = keys %effectTypes;

my %ROH;    #popind -> chr -> array of arrays

print STDERR "Reading $ROHFILE...\n";
open( ROH, "<$ROHFILE" ) or die $!;

my $pop;
my $ind;

#my %ind2pop;
while ( defined( my $rohline = <ROH> ) ) {
    chomp $rohline;

    if ( $rohline =~ m/^track .+Ind: (.+) Pop:(.+) ROH.+/ ) {
        $ind = $1;
        $pop = $2;

        #       $ind2pop{$ind} = $pop;
    }
    else {
        my ( $chr, $start, $end, $class, $size, $junk )
            = split( /\s+/, $rohline, 6 );
        my @tmp = ( $start, $end - 1, $class );
        push( @{ $ROH{$ind}{$chr} }, \@tmp );
    }
}
close(ROH);

my %counts;    #indID -> ROH class -> feature = count

#my %pr_damaging;
#my %po_damaging;
#my %benign;

my @indlist;
##If VCF file
#load individuals, all [should be] the same so just take first vcf file in the list
if ( $FILETYPE eq "VCF" ) {
    my $tempfile = $genfilelist[0];
    if ($ISZIP) {
        open( FIN, "gunzip -c $tempfile |" ) or die "$tempfile $!";
    }
    else {
        open( FIN, "<", $tempfile ) or die "$tempfile $!";
    }
    print STDERR "Loading individual list.\n";

    my ($junk1, $junk2, $junk3, $junk4, $junk5,
        $junk6, $junk7, $junk8, $junk9
    );
    while ( defined( my $line = <FIN> ) ) {
        chomp $line;
        if ( $line =~ m/^##/ ) {
            next;
        }
        elsif ( $line =~ m/^#CHROM/ ) {
            (   $junk1, $junk2, $junk3, $junk4, $junk5,
                $junk6, $junk7, $junk8, $junk9, @indlist
            ) = split( /\s+/, $line );

            for my $ind (@indlist) {
                for my $f (@effectList) {
                    $counts{$ind}{'A'}{$f}    = 0;
                    $counts{$ind}{'B'}{$f}    = 0;
                    $counts{$ind}{'C'}{$f}    = 0;
                    $counts{$ind}{'NONE'}{$f} = 0;
                }
            }
        }
        else {
            last;
        }
    }
    close(FIN);
}
##Otherwise TPED/TFAM
elsif ( $FILETYPE eq "TPED" ) {
    my $tempfile = $TFAMFILE;
    if ($ISZIP) {
        open( FIN, "gunzip -c $tempfile |" ) or die "$tempfile $!";
    }
    else {
        open( FIN, "<", $tempfile ) or die "$tempfile $!";
    }

    while ( defined( my $line = <FIN> ) ) {
        chomp $line;
        my ( $fid, $iid, $junk ) = split( /\s+/, $line, 3 );
        push( @indlist, $iid );
    }

    for my $ind (@indlist) {
        for my $f (@effectList) {
            $counts{$ind}{'A'}{$f}    = 0;
            $counts{$ind}{'B'}{$f}    = 0;
            $counts{$ind}{'C'}{$f}    = 0;
            $counts{$ind}{'NONE'}{$f} = 0;
        }
    }

    close(FIN);
}

for ( my $chr = $startchr; $chr <= $numchr; $chr++ ) {
    my $chrstr = "chr$chr";
    my $file   = $genfilelist[ $chr - $startchr ];
    open( FIN, "<", $file ) or die $!;
    print STDERR "${chrstr}\n";
    my $count = 0;
    if ( $FILETYPE eq "VCF" ) {
        while ( my $line = <FIN> ) {
            chomp $line;
            if ( $line =~ m/^#/ ) {
                next;
            }

            my ($c, $pos,  $rsid, $ref,    $alt,
                $q, $pass, $info, $format, @genotypes
            ) = split( /\s+/, $line );

            if ( $count % 100000 == 0 ) {
                print STDERR "$chrstr $pos\n";
            }
            $count++;

            if ( exists $effect{$chrstr}{$pos} ) {

                my $functionalAllele;
                my $functionalAlleleATCG;
                if ( exists $effect{$chrstr}{$pos}{$ref} ) {
                    $functionalAllele     = 0;
                    $functionalAlleleATCG = $ref;
                }
                elsif ( exists $effect{$chrstr}{$pos}{$alt} ) {
                    $functionalAllele     = 1;
                    $functionalAlleleATCG = $alt;
                }
                else {
                    print STDERR
                        "Neither $ref nor $alt are in the feature file, but $chrstr:${pos} is.";
                    next;
                }

                #print STDERR "$chr $pos ", join(' ',@genotypes), "\n";
                for ( my $i = 0; $i <= $#indlist; $i++ ) {
                    my ( $a1, $a2 ) = split( /\//, $genotypes[$i] );
                    if ( $a1 eq '.' ) {
                        next;
                    }
                    if ( $a1 eq $functionalAllele and $a1 eq $a2 ) {
                        my $class
                            = hitsInterval(
                            \@{ $ROH{ $indlist[$i] }{$chrstr} }, $pos );

                        #print STDERR "$class\n";
                        my $effect
                            = $effect{$chrstr}{$pos}{$functionalAlleleATCG};

                        if ( $class eq '0' ) {
                            $counts{ $indlist[$i] }{'NONE'}{$effect}++;
                        }
                        else {
                            $counts{ $indlist[$i] }{$class}{$effect}++;
                        }

                    }
                }
            }
        }
    }
    elsif ( $FILETYPE eq "TPED" ) {
        while ( my $line = <FIN> ) {
            chomp $line;
            my ( $c, $rsid, $gpos, $pos, @genotypes ) = split( /\s+/, $line );
            if ( $count % 100000 == 0 ) {
                print STDERR "$chrstr $pos\n";
            }
            $count++;

            if ( exists $effect{$chrstr}{$pos} ) {

                for ( my $i = 0; $i <= $#indlist; $i++ ) {
                    my ( $a1, $a2 )
                        = ( $genotypes[ 2 * $i ], $genotypes[ 2 * $i + 1 ] );
                    if ( $a1 eq '0' ) {
                        next;
                    }
                    if ( exists $effect{$chrstr}{$pos}{$a1} and $a1 eq $a2 ) {
                        my $class
                            = hitsInterval(
                            \@{ $ROH{ $indlist[$i] }{$chrstr} }, $pos );

                        #print STDERR "$class\n";
                        my $effect = $effect{$chrstr}{$pos}{$a1};

                        if ( $class eq '0' ) {
                            $counts{ $indlist[$i] }{'NONE'}{$effect}++;
                        }
                        else {
                            $counts{ $indlist[$i] }{$class}{$effect}++;
                        }

                    }
                }
            }
        }
    }
    close(FIN);
}

open( FOUT, ">", $OUTFILE ) or die $!;
my @ROHSizeClasses = ( "A", "B", "C", "NONE" );

for my $f (sort @effectList) {
    for my $class (@ROHSizeClasses) {
        print FOUT "${f}${class} ";
    }
}
print FOUT "\n";

for my $ind (@indlist) {

    #   if(exists $ind2pop{$ind}){
    #       my $pop = $ind2pop{$ind};
    print FOUT "$ind";
    for my $f (sort @effectList) {
        for my $class (@ROHSizeClasses) {
            print FOUT " ", $counts{$ind}{$class}{$f};
        }
    }

    print FOUT "\n";
}

close(FOUT);

#takes an array of intervals and a query
#returns the 3rd element of the interval if query is inside one of the intervals
#returns 0 otherwise
sub hitsInterval {
    my @intervals = @{ $_[0] };
    my $query     = $_[1];

    #print $padding, "\n";
    my $numIntervals = @intervals;

    if ( $numIntervals == 0 ) {
        return 0;
    }

    #print STDERR "number of intervals: $numIntervals\n";

    #print "Starting search at ", floor($numIntervals/2), "\n";

    my @range = ( 0, $numIntervals - 1 );
    my $index = floor( ( $range[1] - $range[0] ) / 2 );
    my $result;
    my $inside = 0;

    do {
        #print "$query\n";
        #print "Index: $index\n";
        #print "Range: ", $range[0], " ", $range[1], "\n";

        if ($query <= $intervals[$index][1] &&   #$query is inside an interval
            $query >= $intervals[$index][0]
            )
        {
#print STDERR "$query is inside [",$intervals[$index][0], ",",$intervals[$index][1] ,"]\n";
            $result = 0;
            $inside = 1;
        }
        elsif ( $query > $intervals[$index][1] ) {

            #print STDERR "$query > ", $intervals[$index][1], "\n";
            if ( $index == $numIntervals - 1 ) {
                $result = 0;
            }
            elsif ( $query < $intervals[ $index + 1 ][0] ) {
                $result = 0;

#print STDERR "Off target query: ";
#print STDERR "[",$intervals[$index][0], ",",$intervals[$index][1] ,"]";
#print STDERR " < $query < ";
#print STDERR "[",$intervals[$index+1][0], ",",$intervals[$index+1][1] ,"]\n";
            }
            else {
                $result   = 1;
                $range[0] = $index;
                $index    = ceil( ( $range[1] - $range[0] ) / 2 ) + $range[0];

 #print "New index = ", $index, "\n";
 #print "[",$intervals[$index][0], ",",$intervals[$index][1] ,"]\n";# and die;
            }
        }
        elsif ( $query < $intervals[$index][0] ) {

            #print STDERR "$query < ", $intervals[$index][0], "\n";
            if ( $index == 0 ) {
                $result = 0;
            }
            elsif ( $query > $intervals[ $index - 1 ][1] ) {
                $result = 0;

  #print STDERR "Off target query: ";
  #print STDERR "[",$intervals[$index-1][0], ",",$intervals[$index-1][1] ,"]";
  #print STDERR " < $query < ";
  #print STDERR "[",$intervals[$index][0], ",",$intervals[$index][1] ,"]\n";
            }
            else {
                $result = -1;
                $range[1] = $index;
                $index = floor( ( $range[1] - $range[0] ) / 2 ) + $range[0];
            }
        }
    } while ( $result != 0 );

    #print "$result\n";
    #die;
    #print STDERR $result, "\n";

    if ( $inside == 0 ) {
        return $inside;
    }

    return $intervals[$index][2];

}
