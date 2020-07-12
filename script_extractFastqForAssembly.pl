#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options + usage statement
my %opts;
getopts('b:r:h', \%opts); # Values in %opts

if ($opts{'h'} || !$opts{'b'} || !$opts{'r'}) {
print "
	Required:
		-b      Input BAM file
		-r      Region of interest: chr1:1000000-1100000

	Optional:
		-h      This helpful help

		-f	fastq file (gzipped)

	\n";
	exit 0;
}

my $inBam = $opts{'b'};
my $targetRegion = $opts{'r'};

# Generate sam file
`samtools view $inBam $targetRegion > target.sam`;

# Open SAM and identify read IDs that have been split
my(%output,%count);
open INF, "target.sam";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;

	my $readID = $F[0];
	my $mapkey = $F[1];
	my $chr = $F[2];
	my $pos = $F[3];
	my $qual = $F[4];
	my $CIGAR = $F[5];

	my($startSplit) = $CIGAR =~ /^(\d+)(S|H)\S+/;
	my($endSplit)   = $CIGAR =~ /(\d+)(S|H)$/;

	print "$readID\t$startSplit\t$endSplit\n$CIGAR\n\n";

	$output{$readID} .= "$chr,$pos,$qual,$mapkey,$startSplit,$endSplit|";
	$count{$readID}++ if $startSplit > 100 || $endSplit > 100;
}
close INF;

if ($opts{'f'}) {
	open INF, "gunzip -c $opts{'f'} |" or die "Can't open fastq file $opts{'f'}: $!";
	while (<INF>) {
		my($readID) = $_ =~ /^\@(\S+)/;
		my $seq = <INF>;
		my $foo = <INF>;
		my $qual = <INF>;
	}
	close INF;
}

foreach my $readID (keys %count) {
	#next if $count{$readID} == 1;

	print "$readID\n";
	my $data = $output{$readID};
	foreach my $foo (split /\|/, $data) {
		$foo =~ s/\,/\t/g;
		print " - $foo\n";
	}
	print "\n";
}
