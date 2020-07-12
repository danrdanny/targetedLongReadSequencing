#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options + usage statement
my %opts;
getopts('f:h', \%opts); # Values in %opts

if ($opts{'h'} || !$opts{'f'}) {
print "
	Required:
		-f	fastq file

	Optional:
		-h      This helpful help

	\n";
	exit 0;
}

my $file = $opts{'f'}; 
my $catCmd = "cat";
   $catCmd = "zcat" if $opts{'f'} =~ /\.gz$/;

sub executeCommand {
	my $r = `$_[0]\n`;
	chomp($r);
	return($r);
}

my %data;

print "Length\treadCount\tbases\taveReadLen\n";

foreach my $readLen ("0", "1000","3000","5000") {
	$data{"reads_$readLen"}	= executeCommand("$catCmd $file | awk \'{if(NR%4==2) print length(\$1)}\' | awk '{ if (\$1 >= $readLen) print \$1 }' | wc -l");
	$data{"bases_$readLen"}	= executeCommand("$catCmd $file | awk \'{if(NR%4==2) print length(\$1)}\' | awk '{ if (\$1 >= $readLen) print \$1 }' | awk '{sum+=\$1} END {printf \"%.0f\\n\", sum}'");

	$data{"aveReadLen_$readLen"} = sprintf("%0.0f", $data{"bases_$readLen"} / $data{"reads_$readLen"});

	print "$readLen\t$data{\"reads_$readLen\"}\t$data{\"bases_$readLen\"}\t$data{\"aveReadLen_$readLen\"}\n";
}

print "\n";
