#!/usr/bin/perl

use strict;
use Getopt::Std;

## Command-line options + usage statement
my %opts;
getopts('at:f:hds', \%opts); # Values in %opts

if ($opts{'h'} || !$opts{'t'}) {
print "
	Required:
		-t      Threads
		-f      location of sample sheet, toml format:
				# General toml format: key = \"value\" pairs, [section names], and # comments

				fastqLoc = [\"fastq 1\", \"fastq 2\", etc]  #
				readVer	= \"version\" # example: guppy361
				newName = \"newName\" # will become <newName>.<readVer>.bam

				genesToPhase = [\"gene1\", \"gene2\", etc] # pulls coordinate data from hg38.genePos.tsv
								        # and calls SNPs/phases +/- 1kb of gene

	Optional:
		-h      This helpful help
		-d	Calculate depth of coverage (forked)
		-s	Calculate fastq file stats (forked)

		#-a	Assemble target region from sample sheet

	\n";
	exit 0;
}

## Define and collect variables used throughout the script
my($fastqLoc,$fastqCount,$sampleName,$readVer,%geneList);

#open INF,"$opts{'f'}" or die "can't open $opts{'f'}: $!\n";
#while (<INF>) {
	#chomp($_);

	#($sampleName) = $_ =~ /^newName\s\=\s\"(\S+)\"$/ if $_ =~ /^newName/;
	#($fastqFile)  = $_ =~ /^fastqLoc\s\=\s\[\(.+)\]$/ if $_ =~ /^fastqLoc/;
	#($readVer)    = $_ =~ /^readVer\s\=\s\"(\S+)\"$/ if $_ =~ /^readVer/;

	#if ($_ =~ /(fastqLoc|genesToPhase)/) {
		#my($chunk) = $_ =~ /\[\"(.+)\]/;
		#my(@D) = split /\"/, $chunk;

		#foreach my $foo (@D) {
			#$foo =~ s/(\"|\,|\s+)//g;
			#next unless $foo =~ /[a-z]/i;

			#++$fastqCount        if $_ =~ /^fastqLoc/;
			#$fastqList{$foo} = 1 if $_ =~ /^fastqLoc/;
			#$geneList{$foo}  = 1 if $_ =~ /^genesToPhase/;
		#}
	#}
#}
#close INF;

my $sampleName = "S007";
my $readVer = "guppy361";

## Use varaibles from above to define additional variables used throughout the program
my $out_fastqFileStats	= "$sampleName.$readVer.fastq.stats";
my $out_finalBam	= "$sampleName.$readVer.bam";
#my $out_longshotBam	= "$sampleName.$readVer.longshot.$geneListForFileName.bam";
my $out_longshotBamH1	= "$sampleName.$readVer.longshot.H1.bam";
my $out_longshotBamH2	= "$sampleName.$readVer.longshot.H2.bam";
my $out_depthOfCoverage	= "$sampleName.$readVer.coverage.10kb.tsv";

my %regionsOfInterest;
#   $regionsOfInterest{"ALMS1"} = "";
#   $regionsOfInterest{"NPHP4"} = "";
#   $regionsOfInterest{"VARS2"} = "chr6,30900000,30940000";

my %geneList;
# For S005
   #$geneList{"HSPG2"}{"fullGenePos"} = "chr1:21822244-21937310";
   #$geneList{"HERC2"}{"fullGenePos"} = "chr15:28111037-28322173";
   #$geneList{"SCN4A"}{"fullGenePos"} = "chr17:63938554-63972918";
   #$geneList{"LIFR"}{"fullGenePos"}  = "chr5:38474668-38606290";
   #$geneList{"CLCN1"}{"fullGenePos"} = "chr7:143316111-143352083";
   #$geneList{"HSPG2"}{"chr"} = "chr1";
   #$geneList{"HERC2"}{"chr"} = "chr15";
   #$geneList{"SCN4A"}{"chr"} = "chr17";
   #$geneList{"LIFR"}{"chr"}  = "chr5";
   #$geneList{"CLCN1"}{"chr"} = "chr7";

# For S007
  $geneList{"chr7target"}{"fullGenePos"} = "chr7:104000000-129000000";
  $geneList{"chr7target"}{"chr"} = "chr7";

## Global Data
my $refGenome	= "/n/dat/hg38/hg38.no_alt.fa";
my $ontmmifile	= "/n/dat/hg38/hg38.no_alt.ont.mmi";
my $sampleSheet	= $opts{'s'};
my $threads	= $opts{'t'};
#my $logDir	= "/n/projects/targetedLongReadSeq/";
#my $logfile 	= getLogFileName($sampleName);


## Collect fastq file stats
#if ($opts{'d'}) {

## Align to genome

#executeCommand("minimap2 -t $threads -a $ontmmifile $fastqFile | samtools sort -@$threads -o $out_finalBam -");
#executeCommand("samtools index $out_finalBam");

## Generate coverage plots

# Convert bam for methylation
#~/bin/nanopore-methylation-utilities/mtsv2bedGraph.py -i ./S004-VARS2.guppy361.chr611.methylationCalls.tsv -g /n/dat/hg38/hg38.no_alt.fa | sort -k1,1 -k2,2n | bgzip > methylation.bed.gz. (quick)
#~/bin/nanopore-methylation-utilities/convert_bam_for_methylation.py -b ./S004-VARS2.guppy361.allReads.bam -c ./methylation.bed.gz -f /n/dat/hg38/hg38.no_alt.fa | samtools sort -o S004.meth.bam (slow)

my $vcfAnnoTemplateLoc = "/n/projects/targetedLongReadSeq/template.vcfAnno.conf.toml";

if ($opts{'a'}) { # a for annotate
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
	my $cmd = "cp /n/projects/targetedLongReadSeq/template.spliceai.in ./variants/spliceai.in";
	executeCommand($cmd);

	## for every gene of interest run longshot and annotate
	foreach my $gene (keys %geneList) {
		my $out_longshotBam   = "variants/$sampleName.$readVer.longshot.$gene.bam";
		my $out_longshotVCF   = "variants/$sampleName.$readVer.longshot.$gene.vcf";
		my $out_longshotVCFAF = "variants/$sampleName.$readVer.longshot.$gene.AF.VEP.vcf";

		my $longshotCmd = "longshot -A --strand_bias_pvalue_cutoff -y 30 -r $geneList{$gene}{'fullGenePos'} -O $out_longshotBam --bam $out_finalBam --ref $refGenome --out $out_longshotVCF"; 
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

			if (!$gnomadAF) {
				#print "No match:\n$_\n\n";
			}

			next if $gnomadAF > 0.02;

			#print "$gnomadAF\n";

			$gnomadAF = "-" if !$gnomadAF;
			my(@F) = split /\t/, $_;
			$out_variantListForSpliceAI .= "$F[0]\t$F[1]\t.\t$F[3]\t$F[4]\t$gnomadAF\t.\t.\n";
		}
		#$out_variantListForSpliceAI =~ s/\n$//;
		$out_variantListForSpliceAI =~ s/chr//g;

		open OUTF,">>$out_vcfForSpliceAI";
		print OUTF $out_variantListForSpliceAI;
		close OUTF;
	}

	# once done, run spliceai
	my $spliceaiCmd = "spliceai -I $out_vcfForSpliceAI -O $out_spliceAIFinalOut -R $refGenome -A grch38";
	executeCommand($spliceaiCmd);

}


#

## subroutines used
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
