#!/usr/bin/perl

use strict;

open INF,"$ARGV[0]";
while (<INF>) {
	chomp($_);
	next if $_ =~ /^\#/;

	my($gnomadAF) = $_ =~ /gnomad_AF\=(\S+)\;/;

	my(@F) = split /\t/, $_;
	if (!$gnomadAF) {
		print "$F[0]\t$F[1]\t.\t$F[3]\t$F[4]\tna\t.\t.\n";
	} else {
		next if $gnomadAF > 0.02;

		print "$F[0]\t$F[1]\t.\t$F[3]\t$F[4]\t$gnomadAF\t.\t.\n";
	}
}
