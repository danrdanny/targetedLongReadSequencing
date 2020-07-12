#!/usr/bin/perl

use strict;
use Getopt::Std;
use String::ShellQuote qw( shell_quote );

## Command-line options + usage statement
my %opts;
getopts('t:s:adhmnqv', \%opts); # Values in %opts

if ($opts{'h'} || !$opts{'t'}) {
print "
	Required:
		-t      Threads
		-s      location of sample sheet, toml format:
				# General toml format: key = \"value\" pairs, [section names], and # comments

				sampleName 	= \"sampleName\" 	# will become <newName>.<readVer>.bam
				fastqFile 	= [\"fastq 1\", etc]  	# can be list or single file
				readVer		= \"version\" 		# example: guppy361
				genesTargeted 	= [\"gene1\", etc] 	# pulls coordinate data from hg38.genePos.tsv

	Optional:
		-h	This helpful help

		Coverage/Read Statistics:
		-q	Calculate fastq file stats for individual files and combined (forked)
		-d	Calculate depth of coverage from bam file (requires bam file, forked)
		#-c	Calculate coverage of target regions (requires bam file and coverage file from option d)

		Align:
		-a	Align with minimap2 (fast, does not fork)
		-n	Align with ngmlr (slow, does not fork)

		Variants:
		-v	Call variants with longshot, annotate with VEP, and put low freq variants through spliceai
		-m	Call structural variants with Sniffles and SVIM

		Future/Other:
		#-?	Assemble target region from sample sheet
		#-?	Make methylation-specific bam files
		#-?	Phase target regions (either methylation or based on longshot/whatshap)
	\n";
	exit 0;
}

## ------------------------------------------------------- ##
## Define and collect variables used throughout the script ##
## ------------------------------------------------------- ##

my(%fastq,$sampleName,$readVer,%geneList,$geneCount);

my $sampleSheet;
open INF,"$opts{'s'}" or die "can't open $opts{'s'}: $!\n";
while (<INF>) {
	next if $_ =~ /^#/;
	chomp($_);
	$sampleSheet .= "$_\n";
}
close INF;

my($sampleName) = $sampleSheet =~ /sampleName\s\=\s\"(\S+)\"/;
my($readVer)    = $sampleSheet =~ /readVer\s\=\s\"(\S+)\"/;
my($fastqData)  = $sampleSheet =~ /fastqFile\s\=\s\[(.+)\]/;

$sampleSheet =~ s/\,//g;
$sampleSheet =~ s/\n/,/g;
my($geneData)   = $sampleSheet =~ /genesTargeted\s\=\s\[(.+)\]/;

foreach my $foo (split /\,/, $fastqData) {
	$foo =~ s/\"|\s+//g;
	next unless $foo =~ /[a-z]/i;

	$fastq{"count"}++;
	$fastq{"fileList"} .= "$foo, ";
	$fastq{"file"}{$foo} = 1;
}
$fastq{"fileList"} =~ s/\,\s$//;

print "Please confirm\n
Sample name: 	$sampleName
Read version:	$readVer
Fastq list: 	$fastq{'fileList'}\n";

foreach my $foo (split /\,/, $geneData) {
	next unless $foo =~ /[a-z]/i;
	my(@D) = split /\t/, $foo;
	
	$D[1] =~ s/\s+//g;
	$D[2] =~ s/\s+//g;
	$D[3] =~ s/\s+//g;
	$D[4] =~ s/\s+//g;
	
	$geneList{$D[0]}{"gene"} = "";
	$geneList{$D[0]}{"targetRegion"} = "";
	$geneList{$D[0]}{"regionToAnalyze"} = "$D[1]:$D[2]-$D[3]";
	$geneList{$D[0]}{"chr"} = "$D[1]";
	++$geneCount;

	print "\t$D[0]\t$geneList{$D[0]}{'regionToAnalyze'}\n";
}

print "Total genes targeted:	$geneCount\n";
print "\n";


## ----------------------------------------------------------------------------------- ##
## Use varaibles from above to define additional variables used throughout the program ##
## ----------------------------------------------------------------------------------- ##

my $out_fastqFileStats	= "$sampleName.$readVer.fastq.stats";
my $out_finalBam	= "$sampleName.$readVer.bam";
my $out_finalBam_ngmlr	= "$sampleName.$readVer.ngmlr.bam";
#my $out_longshotBam	= "$sampleName.$readVer.longshot.$geneListForFileName.bam";
my $out_longshotBamH1	= "$sampleName.$readVer.longshot.H1.bam";
my $out_longshotBamH2	= "$sampleName.$readVer.longshot.H2.bam";
my $out_depthOfCoverage	= "$sampleName.$readVer.coverage.10kb.tsv";

## Global Data
my $refGenome	= "/n/dat/hg38/hg38.no_alt.fa";
my $ontmmifile	= "/n/dat/hg38/hg38.no_alt.ont.mmi";
my $threads	= $opts{'t'};
#my $logDir	= "/n/projects/targetedLongReadSeq/";
#my $logfile 	= getLogFileName($sampleName);

my $vcfAnnoTemplateLoc = "/n/projects/targetedLongReadSeq/template.vcfAnno.conf.toml";
my $spliceaiTemplateLoc = "/n/projects/targetedLongReadSeq/template.spliceai.in";

## ------------------------ ##
## Collect fastq file stats ##
## ------------------------ ##

if ($opts{'q'}) {
	print "[".localtime(time)."] RUN: Calculating fastq stats for $sampleName (forked)\n";

	my $pid = fork();
	die "unable to fork: $!" unless defined($pid);
	if (!$pid) {
		my %fastqStats;
		foreach my $fastqFile (keys %{$fastq{'file'}}) {
			my $output = "File\tLength\treadCount\tbases\taveReadLen\n";

			foreach my $readLen ("0","1000","3000","5000") {
				my $reads = executeCommand("zcat $fastqFile | awk \'{if(NR%4==2) print length(\$1)}\' | awk '{ if (\$1 >= $readLen) print \$1 }' | wc -l");
				chomp($reads);
				my $bases = executeCommand("zcat $fastqFile | awk \'{if(NR%4==2) print length(\$1)}\' | awk '{ if (\$1 >= $readLen) print \$1 }' | awk '{sum+=\$1} END {printf \"%.0f\\n\", sum}'");
				chomp($bases);

				$fastqStats{"$readLen-reads"} += $reads;
				$fastqStats{"$readLen-bases"} += $bases;

				my $aveReadLen = sprintf("%0.0f", $bases / $reads);
				$output .= "$fastqFile\t$readLen\t$reads\t$bases\t$aveReadLen\n";
			}

			open OUTF,">fastqStats.$fastqFile";
			print OUTF $output;
			close OUTF;
		}

		my $output = "File\tLength\treadCount\tbases\taveReadLen\n";
		foreach my $readLen ("0","1000","3000","5000") {
			my $aveReadLen = sprintf("%0.0f", $fastqStats{"$readLen-bases"} / $fastqStats{"$readLen-reads"});
			$output .= "all\t$readLen\t".$fastqStats{"$readLen-reads"}."\t".$fastqStats{"$readLen-bases"}."\t$aveReadLen\n";

			open OUTF,">fastqStats.all";
			print OUTF $output;
			close OUTF;
		}

		print "[".localtime(time)."] DONE: Calculating fastq stats for $sampleName (forked)\n";
		exit 1;
	}
}

## --------------- ##
## Align to genome ##
## --------------- ##

if ($opts{'a'} || $opts{'n'}) {

	my $tempFastq = "tmp.fastq.gz";
	foreach my $fastqFile (keys %{$fastq{'file'}}) {
		executeCommand("cat $fastqFile >> $tempFastq");
	}

	# align with minimap2
	if ($opts{'a'} && !-e $out_finalBam) {
		print "[".localtime(time)."] RUN: Aligning $sampleName with minimap2.\n";
	
		executeCommand("minimap2 -t $threads -a $ontmmifile $tempFastq | samtools sort -\@$threads -o $out_finalBam -");
		executeCommand("samtools index $out_finalBam");

		print "[".localtime(time)."] DONE: Aligning $sampleName with minimap2.\n";
	} elsif ($opts{'a'}) {
		print "[".localtime(time)."] WARN: Did not align $sampleName with minimap2 because $out_finalBam exists.\n";
	}

	# align with ngmlr
	if ($opts{'n'} && !-e $out_finalBam_ngmlr) {
		print "[".localtime(time)."] RUN: Aligning $sampleName with ngmlr.\n";

		executeCommand("ngmlr -x ont -t $threads -r $refGenome -q $tempFastq | samtools sort -\@$threads -o $out_finalBam_ngmlr -");
		executeCommand("samtools index $out_finalBam_ngmlr");

		print "[".localtime(time)."] DONE: Aligning $sampleName with ngmlr.\n";
	} elsif ($opts{'n'}) {
		print "[".localtime(time)."] WARN: Did not align $sampleName with ngmlr because $out_finalBam_ngmlr exists.\n";
	}

	executeCommand("rm -f $tempFastq");
}

## --------------------------- ##
## Calculate depth of coverage ##
## --------------------------- ##

if ($opts{'d'} && -e $out_finalBam) {
	print "[".localtime(time)."] RUN: Calculating depth of coverage for $sampleName based on $out_finalBam (forked)\n";

	my $pid = fork();
	die "unable to fork: $!" unless defined($pid);
	if (!$pid) {
		my $cmd = "perl /n/projects/targetedLongReadSeq/script_calculateCoverage.pl -b $out_finalBam -c all -p $sampleName";
		executeCommand($cmd);

		print "[".localtime(time)."] DONE: Calculating depth of coverage for $sampleName (forked)\n";
		exit 1;
	}
} elsif ($opts{'d'}) {
	print "[".localtime(time)."] WARN: Unable to calculate depth of coverage for $sampleName as bam file does not exist ($out_finalBam)\n";
}

## -------------------------------- ##
## Generate per-gene coverage plots ##
## -------------------------------- ##


## ------------------------------------- ##
## Convert bam files to view methylation ##
## ------------------------------------- ##

#~/bin/nanopore-methylation-utilities/mtsv2bedGraph.py -i ./S004-VARS2.guppy361.chr611.methylationCalls.tsv -g /n/dat/hg38/hg38.no_alt.fa | sort -k1,1 -k2,2n | bgzip > methylation.bed.gz. (quick)
#~/bin/nanopore-methylation-utilities/convert_bam_for_methylation.py -b ./S004-VARS2.guppy361.allReads.bam -c ./methylation.bed.gz -f /n/dat/hg38/hg38.no_alt.fa | samtools sort -o S004.meth.bam (slow)

## ------------------------ ##
## Call Structural Variants ##
## ------------------------ ##

if ($opts{'m'}) { 
	my $out_sniffles_minimap2 = "$sampleName-out_sv_sniffles_minimap2.vcf";
	my $out_sniffles_ngmlr    = "$sampleName-out_sv_sniffles_ngmlr.vcf";

	my $cmd = "mkdir -p sniffles";
	executeCommand($cmd);

	if (-e $out_finalBam) {

		# assume it does not have the MD tag
		my $bamForSniffles = "$sampleName.forSniffles.MDtag.minimap2.bam";
		if (!-e $bamForSniffles) {
			print "[".localtime(time)."] RUN: adding MD tag to $out_finalBam.\n";
			executeCommand("samtools calmd -\@$threads -b $out_finalBam $refGenome > $bamForSniffles");
			executeCommand("samtools index $bamForSniffles");
		}

		print "[".localtime(time)."] RUN: running sniffles on $bamForSniffles.\n";
		executeCommand("sniffles -m $bamForSniffles -v sniffles/$out_sniffles_minimap2");

		print "[".localtime(time)."] RUN: running svim on $bamForSniffles.\n";
		executeCommand("svim alignment --read_names --sequence_alleles --insertion_sequences --sample $sampleName svim_minimap2 $out_finalBam $refGenome");

		executeCommand("rm -f $bamForSniffles $bamForSniffles.bai");

	} else {
		print "[".localtime(time)."] WARN: Cannot run sniffles or svim on $out_finalBam as it does not exist.\n";
	}

	if (-e $out_finalBam_ngmlr) {

		# assume it does not have the MD tag
		my $bamForSniffles = "$sampleName.forSniffles.MDtag.ngmlr.bam";
		if (!-e $bamForSniffles) {
			print "[".localtime(time)."] RUN: adding MD tag to $out_finalBam.\n";
			executeCommand("samtools calmd -\@$threads -b $out_finalBam_ngmlr $refGenome > $bamForSniffles");
			executeCommand("samtools index $bamForSniffles");
		}

		print "[".localtime(time)."] RUN: running sniffles on $bamForSniffles.\n";
		executeCommand("sniffles -m $bamForSniffles -v sniffles/$out_sniffles_ngmlr");

		print "[".localtime(time)."] RUN: running svim on $out_finalBam_ngmlr.\n";
		executeCommand("svim alignment --read_names --sequence_alleles --insertion_sequences --sample $sampleName svim_ngmlr $out_finalBam_ngmlr $refGenome");

		executeCommand("rm -f $bamForSniffles $bamForSniffles.bai");

	} else {
		print "[".localtime(time)."] WARN: Cannot run sniffles or svim on $out_finalBam_ngmlr as it does not exist.\n";
	}
		
}


## ------------- ##
## Call Variants ##
## ------------- ##

if ($opts{'v'}) { 
	my $cmd = "mkdir -p variants";
	executeCommand($cmd);

	## first grab the vcfAnno toml template
	my $vcfAnnoTemplate;
	open INF,"$vcfAnnoTemplateLoc";
	while (<INF>) {
		$vcfAnnoTemplate .= $_;
	}
	close INF;

	## make a new out_forSpliceAI file with the correct header
	my $out_vcfForSpliceAI = "variants/spliceai.in";
	my $out_spliceAIFinalOut = "variants/spliceai.out";
	my $out_variantListForSpliceAI;
	my $cmd = "cp $spliceaiTemplateLoc ./variants/spliceai.in";
	executeCommand($cmd);

	## for every gene of interest run longshot and annotate
	foreach my $gene (keys %geneList) {
		my $out_longshotBam   = "variants/$sampleName.$readVer.longshot.$gene.bam";
		my $out_longshotVCF   = "variants/$sampleName.$readVer.longshot.$gene.vcf";
		my $out_longshotVCFAF = "variants/$sampleName.$readVer.longshot.$gene.AF.VEP.vcf";

		my $longshotCmd = "longshot -A --strand_bias_pvalue_cutoff -y 30 -r $geneList{$gene}{'regionToAnalyze'} -O $out_longshotBam --bam $out_finalBam --ref $refGenome --out $out_longshotVCF"; 
		executeCommand($longshotCmd);

		# need to make new vcfAnno file
		my $chr = $geneList{$gene}{"chr"};
		   $chr =~ s/chr//;
		my $tmp_vcfAnnoTemplate = $vcfAnnoTemplate;
		$tmp_vcfAnnoTemplate =~ s/\<chr\>/$chr/;

		open OUTF,">variants/vcfAnno.conf.toml";
		print OUTF $tmp_vcfAnnoTemplate;
		close OUTF;

		my $vcfannoCmd = "vcfanno variants/vcfAnno.conf.toml $out_longshotVCF > $out_longshotVCFAF"; 
		executeCommand($vcfannoCmd);
		
		open INF,"$out_longshotVCFAF";
		while (<INF>) {
			chomp($_);
			next if $_ =~ /^\#/;

			my($gnomadAF) = $_ =~ /\;gnomad_AF\=(\S+)\;/;
			next if $gnomadAF > 0.02;

			$gnomadAF = "-" if !$gnomadAF;
			my(@F) = split /\t/, $_;
			$out_variantListForSpliceAI .= "$F[0]\t$F[1]\t.\t$F[3]\t$F[4]\t$gnomadAF\t.\t.\n";
		}
	}
	#$out_variantListForSpliceAI =~ s/\n$//;
	$out_variantListForSpliceAI =~ s/chr//g;

	open OUTF,">>$out_vcfForSpliceAI";
	print OUTF $out_variantListForSpliceAI;
	close OUTF;

	# once done, run spliceai
	my $spliceaiCmd = "spliceai -I $out_vcfForSpliceAI -O $out_spliceAIFinalOut -R $refGenome -A grch38";
	executeCommand($spliceaiCmd);

}


## -------------------------------------- ##
## Summarize variants from target regions ##
## -------------------------------------- ##

if ($opts{'z'}) {
# per gene
# SNVs
# spliceAI
# cat spliceai.out | grep -v "#" | awk -F"|" '{ if ($3 > .2 || $4 > .2 || $5 > .2 || $6 > .2) print}'


# SVIM
#   | awk '$1 ~ /chr6/' | awk '{if($2>30900000) print}' | awk '{if($2<30930000) print}' | awk '{if($6>6) print}' |


	#foreach my $gene (keys %geneList) {
		
		#open INF,"$out_longshotVCFAF";
		#while (<INF>) {
			#chomp($_);
			#next if $_ =~ /^\#/;
#
			#my($gnomadAF) = $_ =~ /\;gnomad_AF\=(\S+)\;/;
			#next if $gnomadAF > 0.02;
#
			#$gnomadAF = "-" if !$gnomadAF;
			#my(@F) = split /\t/, $_;
			#$out_variantListForSpliceAI .= "$F[0]\t$F[1]\t.\t$F[3]\t$F[4]\t$gnomadAF\t.\t.\n";
		#}
	#}
	##$out_variantListForSpliceAI =~ s/\n$//;
	#$out_variantListForSpliceAI =~ s/chr//g;
#
	#open OUTF,">>$out_vcfForSpliceAI";
	#print OUTF $out_variantListForSpliceAI;
	#close OUTF;
}

## ----------- ##
## Subroutines ##

## ----------- ##
## Subroutines ##
## ----------- ##

sub getLogFileName {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$year += 1900;
	$mon += 1;
	$mon = "0$mon" if $mon =~ /^\d$/; # this is because I'm OCD and can't stand when the log file names don't look nice
	$mday = "0$mday" if $mday =~ /^\d$/;
	$hour = "0$hour" if $hour =~ /^\d$/;
	$min = "0$min" if $min =~ /^\d$/;
	$sec = "0$sec" if $sec =~ /^\d$/;
	my $logfile = "log.$_[0].$year-$mon-$mday\_$hour\:$min\:$sec";
	return($logfile);
}

sub executeCommand {
	#open LOGF, ">>$logDir/$logfile";
	print "[".localtime()."] CMD: $_[0]\n";
	#print LOGF "[".localtime()."] CMD: $_[0]\n";
	#close LOGF;
	my $output = `$_[0]`;
	return($output);
}

#sub logData {
#	open LOGF, ">>$logDir/$logfile";
#	print LOGF "[".localtime()."] LOG: $_[0]\n";
#	close LOGF;
#	return 1;
#}
