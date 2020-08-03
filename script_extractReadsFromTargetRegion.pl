#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options + usage statement
my %opts;
getopts('b:r:f:h', \%opts); # Values in %opts

if ($opts{'h'} || !$opts{'b'} || !$opts{'r'} || !$opts{'f'}) {
print "
	Required:
		-b      Input BAM file
		-r      Region of interest: chr1:1000000-1100000
		-f	fastq file (gzipped)

	Optional:
		-h      This helpful help

	\n";
	exit 0;
}

my $inBam = $opts{'b'};
my $targetRegion = $opts{'r'};
my $fastqFile = $opts{'f'};

# Generate sam file
print "Making target.sam file for reads in $targetRegion\n";
`samtools view $inBam $targetRegion > target.sam`;

# Open SAM and identify read IDs that have been split
print "Getting readIDs from target.sam\n";
my(%reads);
open INF, "target.sam";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;

	$reads{$F[0]} = 1;
}
close INF;

## Get sequence of read
print "Getting original sequence data from $fastqFile\n";
my $output;
open INF, "gunzip -c $fastqFile |" or die "Can't open fastq file $fastqFile: $!";
while (<INF>) {
	my($readID) = $_ =~ /^\@(\S+)/;
	my $seq = <INF>;
	my $foo = <INF>;
	my $qual = <INF>;

	next unless $reads{$readID} == 1;

	chomp($seq);
	$output .= ">$readID\n$seq\n";
}
close INF;

#$targetRegion =~ s/(\:,\-)/_/g;

open OUTF,">isolatedReads.$targetRegion.fa";
print OUTF $output;
close OUTF;

`rm -f target.sam`;
