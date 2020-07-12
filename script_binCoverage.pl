#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options + usage statement
my %opts;
getopts('b:c:hp:', \%opts); # Values in %opts

if ( $opts{'h'} || !$opts{'b'} || !$opts{'c'} || !$opts{'p'} ) {
print "
        Required:
                -b	bam file you want to calculate coverage for
		-c	chromosome to calc coverage for, format is chr<num>
			pass 'all' for all chromosomes
		-p	prefix for output coverage file

        Optional:
                -h      This helpful help

        \n";
        exit 0;
}

my $sample = $opts{'b'};
my $samplePrefix = $opts{'p'};

my @chr;
if ($opts{'c'} eq "all") {
	foreach my $chrNum (1..22) {
		push(@chr,"chr$chrNum");
	}

	push(@chr,"chrX");
} else {
	push(@chr,"$opts{'c'}");
}

my $tmpScriptDir = "tmpScriptDirectory";
`mkdir -p $tmpScriptDir`;

my(@pids);
print "[".localtime(time)."] Running samtools depth on each chr in \@chr.\n";
foreach my $chr (@chr) {
	next unless $chr =~ /^chr/;

	die "could not fork" unless defined(my $pid = fork);
	unless ($pid) { #child execs
        	my $cmd = "samtools depth -a -Q 0 -r $chr $sample > $tmpScriptDir/$chr.depthOfCoverage.tsv";
		`$cmd`;
		exit 1;
	}
    	push @pids, $pid; #parent stores children's pids
}

#wait for all children to finish
print "[".localtime(time)."] Waiting for samtools depth programs to complete.\n";
for my $pid (@pids) {
	waitpid $pid, 0;
}
print "[".localtime(time)."] All samtools depth programs done.\n";

## ------------ ##
## Bin coverage ##
## ------------ ##

my(@pids);
print "[".localtime(time)."] Binning coverage for each chr in \@chr.\n";
foreach my $chr (@chr) {
	next unless $chr =~ /^chr/;

	my $inf = "$tmpScriptDir/$chr.depthOfCoverage.tsv";

	if ($chr =~ /(chr1|chr2|chr3|chr4)/) {
		binCoverage($inf,$chr);
	} else {
		die "could not fork" unless defined(my $pid = fork);
		unless ($pid) { #child execs
			binCoverage($inf,$chr);

			exit 1;
		}
    		push @pids, $pid; #parent stores children's pids
		sleep 90;
	}
}

#wait for all children to finish
print "[".localtime(time)."] Waiting for binning programs to complete.\n";
for my $pid (@pids) {
	waitpid $pid, 0;
}
print "[".localtime(time)."] All binning programs done.\n";

## ------------------------------------------ ##
## Consolidate data, delete intermeidate data ##
## ------------------------------------------ ##

my $finalOutputHeader = "chr\tpos\tcount";
my $finalCoverageFile = "$samplePrefix.coverage.10kbwindow.tsv";
`echo \"$finalOutputHeader\" > $finalCoverageFile`;

print "[".localtime(time)."] \n";
print "[".localtime(time)."] Consolidating data.\n";
foreach my $chr (@chr) {
	next unless $chr =~ /^chr/;
	my $windowFile = "$tmpScriptDir/$chr.depthOfCoverage.10kbWindows.tsv";

	if (!-e $windowFile) {
		print "[".localtime(time)."] WARNING: $chr depth of coverage file missing\n";
	}

	`cat $tmpScriptDir/$chr.depthOfCoverage.10kbWindows.tsv >> $finalCoverageFile`;
}

# Delete data


# Done

sub binCoverage {
	my($inf,$chr) = ($_[0],$_[1]);

	my(%depthOfCov,%totalDepth,%totalDepthCount,$lastChr,%chrMin,%chrMax);
	open INF,"$inf" or die "Can't open $inf: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;
		next unless $F[2] > 0;
		next if $_ =~ /^#/;
		my $chr = $F[0];
		my $pos = $F[1];

		$depthOfCov{$chr}{$pos} = $F[2];
		$totalDepth{$chr}+= $F[2];
		$totalDepthCount{$chr}++;
		$chrMin{$chr} = $F[1] if $F[1] < $chrMax{$chr} || !$chrMin{$chr};
		$chrMax{$chr} = $F[1] if $F[1] > $chrMax{$chr};
	}
	close INF;
	#print "[".localtime(time)."] Data collected, analyzing.\n";

	my $output;
	foreach my $chr (sort keys %totalDepth) {
		my $aveDepth = sprintf("%0.1f", $totalDepth{$chr} / $totalDepthCount{$chr});
		my $min         = 0; #$chrMin{$chr};
		my $max         = $min + 10000;
		my $step        = 10000;

		until ($min >= $chrMax{$chr}) {
			my($totalDepth,$count);
			foreach my $id ($min..$max) {
				next unless $depthOfCov{$chr}{$id}; #> 0; #id < $min || $id >= $max;
				$totalDepth += $depthOfCov{$chr}{$id};
				$count++;
			}

			my($currAveDepth,$ratio,$logDepth);
			if ($count == 0) {
				$currAveDepth = 0;
			} else {
				$currAveDepth = $totalDepth / $count;
				$currAveDepth = sprintf("%0.0f", $currAveDepth);
			}

			$output .= "$chr\t$min\t$currAveDepth\n";
	
			$min += $step;
			$max += $step;
		}
	}

	open OUTF,">$tmpScriptDir/$chr.depthOfCoverage.10kbWindows.tsv";
	print OUTF $output;
	close OUTF;
}
