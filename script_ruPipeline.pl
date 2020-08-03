#!/usr/bin/perl

use strict;
use Getopt::Std;
#use String::ShellQuote qw( shell_quote );

## Command-line options + usage statement
my %opts;
getopts('t:s:adhmnrqv', \%opts); # Values in %opts

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
		-r	Report variants (list options ... )

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
my($fastqData,$geneData,$analysisRegion);

foreach my $foo (split /\]/, $sampleSheet) {
	if ($foo =~ /fastqFile/) {
		($fastqData) = $foo =~ /fastqFile\s\=\s\[(.+)/;
	} elsif ($foo =~ /genesTargeted/) {
		$foo =~ s/(\,|\])//g;
		$foo =~ s/\n/,/g;
		($geneData) = $foo =~ /genesTargeted\s\=\s\[(.+)/;
	} elsif ($foo =~ /analysisRegion/) {
		$foo =~ s/(\,|\])//g;
		$foo =~ s/\n/,/g;
		($analysisRegion) = $foo =~ /analysisRegion\s\=\s\[(.+)/;
	}
}

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
Fastq list: 	$fastq{'fileList'}
\n";

## Genes targeted
print "Genes targeted:\n";
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
	$geneList{$D[0]}{"regionStart"} = "$D[2]";
	$geneList{$D[0]}{"regionEnd"}   = "$D[3]";
	$geneList{$D[0]}{"chr"} = "$D[1]";
	++$geneCount;

	print "\t$D[0]\t$geneList{$D[0]}{'regionToAnalyze'}\n";
}
print "Total: $geneCount\n";
print "\n";

## Analysis region
print "Analysis regions:\n";
foreach my $foo (split /\,/, $analysisRegion) {
	next unless $foo =~ /[a-z]/i;
	my(@D) = split /\t/, $foo;
	
	$D[0] =~ s/\s+//g;
	$D[1] =~ s/\s+//g;
	$D[2] =~ s/\s+//g;
	$D[3] =~ s/\s+//g;
	
	push @{$geneList{'analysisRegion'}}, "$D[0]|$D[1]:$D[2]-$D[3]";
	print "\t$D[0]\t$D[1]:$D[2]-$D[3]\n";
}
print "\n";
print "\n";


## ----------------------------------------------------------------------------------- ##
## Use varaibles from above to define additional variables used throughout the program ##
## ----------------------------------------------------------------------------------- ##

my $out_fastqFileStats	= "$sampleName.$readVer.fastq.stats";
my $out_finalBam	= "$sampleName.$readVer.bam";
my $out_finalBam_ngmlr	= "$sampleName.$readVer.ngmlr.bam";
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
	my $snifflesPath = "/home/dem/bin/Sniffles-1.0.12/bin/sniffles-core-1.0.11/sniffles";

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
		executeCommand("$snifflesPath -m $bamForSniffles -v sniffles/$out_sniffles_minimap2");

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
		executeCommand("$snifflesPath -m $bamForSniffles -v sniffles/$out_sniffles_ngmlr");

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

	## Make each output directory
	my $longshotDir = "variant_longshot";
	my $clairDir    = "variant_clair";
	my $medakaDir   = "variant_medaka";

	foreach my $blah ($longshotDir, $clairDir, $medakaDir) {
		executeCommand("mkdir -p $blah");
	}

	## Grab the vcfAnno toml template
	my $vcfAnnoTemplate;
	open INF,"$vcfAnnoTemplateLoc";
	while (<INF>) {
		$vcfAnnoTemplate .= $_;
	}
	close INF;

	## subroutine to grab low freq alleles in each round of variant calling
	my(%variantDetail);
	sub filterVariantVCFFile {
		my($vcf,$caller) = ($_[0],$_[1]);

		open INF,"$vcf";
		while (<INF>) {
			chomp($_);
			next if $_ =~ /^\#/;
			next if $_ =~ /LowQual/ && $caller eq "clair"; # clair calls everything, so the low quality calls are not valuable

			my($gnomadAF) = $_ =~ /gnomad_AF\=(\S+)\;/;
			next if $gnomadAF > 0.02;

			$gnomadAF = 0 if !$gnomadAF;
			my(@F) = split /\t/, $_;
			#$out_variantListForSpliceAI .= "$F[0]\t$F[1]\t.\t$F[3]\t$F[4]\t$gnomadAF\t$caller\t.\n";

			my $key = "$F[0]|$F[1]|$F[3]|$F[4]";
#print "$F[0]\t$F[1]\t$F[3]\t$F[4]\t$gnomadAF\t$caller\n";
			$variantDetail{$key}{"caller"} .= "$caller,";
			$variantDetail{$key}{"AF"} = $gnomadAF;

#print "\thash_caller: $variantDetail{$key}{'caller'}\n";
#print "\thash_af:     $variantDetail{$key}{'AF'}\n";
		}
		return 1;
	}

	## for every gene of interest run longshot, medaka, and clair then annotate their vcf files
	foreach my $aligner ("minimap2","ngmlr") {
		my $variantBam = $out_finalBam; # this is the minimap2 bam
		   $variantBam = $out_finalBam_ngmlr if $aligner eq "ngmlr";

		foreach my $gene (keys %geneList) {
			# need to make new vcfAnno file
			my $chr = $geneList{$gene}{"chr"};
		   	$chr =~ s/chr//;
			my $tmp_vcfAnnoTemplate = $vcfAnnoTemplate;
			$tmp_vcfAnnoTemplate =~ s/\<chr\>/$chr/;

			open OUTF,">vcfAnno.conf.toml";
			print OUTF $tmp_vcfAnnoTemplate;
			close OUTF;

			## -------- ##
			## LONGSHOT ##
			## -------- ##

			my $out_longshotBam   = "$longshotDir/$sampleName.$readVer.longshot.$aligner.$gene.bam";
			my $out_longshotVCF   = "$longshotDir/$sampleName.$readVer.longshot.$aligner.$gene.vcf";
			my $out_longshotVCFAF = "$longshotDir/$sampleName.$readVer.longshot.$aligner.$gene.AF.VEP.vcf";

			my $longshotCmd = "longshot -A --strand_bias_pvalue_cutoff -y 30 -r $geneList{$gene}{'regionToAnalyze'} -O $out_longshotBam --bam $variantBam --ref $refGenome --out $out_longshotVCF"; 
			executeCommand($longshotCmd);
			executeCommand("vcfanno vcfAnno.conf.toml $out_longshotVCF > $out_longshotVCFAF");
			filterVariantVCFFile($out_longshotVCFAF,"longshot_$aligner");
		
			## ------ ##
			## MEDAKA ##
			## ------ ##

			my $out_medakaVCFAF = "$medakaDir/$sampleName.$readVer.medaka.$aligner.$gene.AF.VEP.vcf";
	
			my $cmd  = "\#!/bin/bash\nsource /home/dem/miniconda3/etc/profile.d/conda.sh\n";
		   	$cmd .= "source activate medaka\n";
		   	$cmd .= "medaka_variant -i $variantBam -f $refGenome -o $medakaDir/$aligner/$gene -t $threads -r $geneList{$gene}{'regionToAnalyze'}\n";

			open OUTF,">medaka.sh";
			print OUTF $cmd;
			close OUTF;

			executeCommand("bash medaka.sh");
			executeCommand("vcfanno vcfAnno.conf.toml $medakaDir/$aligner/$gene/round_1.vcf > $out_medakaVCFAF");
			filterVariantVCFFile($out_medakaVCFAF,"medaka_$aligner");

			## ----- ##
			## CLAIR ##
			## ----- ##

			my $out_clairVCFAF = "$clairDir/$sampleName.$readVer.clair.$aligner.$gene.AF.VEP.vcf";
	
			my $cmd  = "\#!/bin/bash\nsource /home/dem/miniconda3/etc/profile.d/conda.sh\n";
		   	$cmd .= "source activate clair-env\n";
		   	$cmd .= "clair.py callVarBam --chkpnt_fn /home/dem/ont/model --bam_fn $variantBam --ref_fn $refGenome --call_fn $clairDir/$gene.$aligner.vcf --sampleName $gene --pysam_for_all_indel_bases --threads $threads --qual 100 --ctgName chr$chr --ctgStart $geneList{$gene}{'regionStart'} --ctgEnd $geneList{$gene}{'regionEnd'}";
	
			open OUTF,">clair.sh";
			print OUTF $cmd;
			close OUTF;

			executeCommand("bash clair.sh");
			executeCommand("vcfanno vcfAnno.conf.toml $clairDir/$gene.$aligner.vcf > $out_clairVCFAF");
			filterVariantVCFFile($out_clairVCFAF,"clair_$aligner");
		}
	}

	## Remove temporary files
	executeCommand("rm -f medaka.sh clair.sh vcfAnno.conf.toml");

	## Run spliceai
	my $out_variantListForSpliceAI;
	foreach my $foo (keys %variantDetail) {
		my($chr,$pos,$ref,$alt) = $foo =~ /(\S+)\|(\d+)\|(\S+)\|(\S+)/;
		my $af      = $variantDetail{$foo}{"AF"};
		my $callers = $variantDetail{$foo}{"caller"};

		$callers =~ s/\,$//g;
		my $callerCount;
		foreach (split /\,/, $callers) { ++$callerCount; }

		##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
		##INFO=<ID=CC,Number=1,Type=Integer,Description="Caller Count">
		##INFO=<ID=CL,Number=1,Type=String,Description="Caller List">
		#CHROM	POS	ID	REF	ALT	INFO

		$chr =~ s/chr//;
		$out_variantListForSpliceAI .= "$chr\t$pos\t.\t$ref\t$alt\tAF=$af;CC=$callerCount;CL=$callers\n";
	}

	## make a new out_forSpliceAI file with the correct header
	my $out_vcfForSpliceAI = "spliceai.in";
	my $out_spliceAIFinalOut = "spliceai.out";
	executeCommand("cp $spliceaiTemplateLoc ./spliceai.in");

	open OUTF,">>$out_vcfForSpliceAI";
	print OUTF $out_variantListForSpliceAI;
	close OUTF;

	executeCommand("spliceai -I $out_vcfForSpliceAI -O $out_spliceAIFinalOut -R $refGenome -A grch38");
}


## -------------------------------------- ##
## Summarize variants from target regions ##
## -------------------------------------- ##

if ($opts{'r'}) {
	my $spliceaiCutoff	= 0.5;
	my(%geneTargets,%coverage);

	my $in_spliceai		= "spliceai.in";
	my $out_spliceai	= "spliceai.out";
	my $out_sniffles_minimap2 = "sniffles/$sampleName\-out_sv_sniffles_minimap2.vcf";
	my $out_sniffles_ngmlr 	  = "sniffles/$sampleName\-out_sv_sniffles_ngmlr.vcf";
	my $out_svim_minimap2 	= "svim_minimap2/variants.vcf";
	my $out_svim_ngmlr    	= "svim_ngmlr/variants.vcf";

	open INF,"/n/projects/targetedLongReadSeq/sequencingTargets/hg38.genePos.tsv" or die "can't open file: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;

		next unless $geneList{$F[1]};
		$geneTargets{$F[1]} = "$F[0]:$F[3]-$F[4]";
	}
	close INF;

	# targetRegion and analysisRegion are different. analysisRegion may be deletions within a target region
	foreach my $i ( 0 .. $#{ $geneList{'analysisRegion'} } ) {
		my($region,$chr,$start,$end) = $geneList{'analysisRegion'}[$i] =~ /(\S+)\|(\S+)\:(\d+)\-(\d+)/;
		$geneTargets{$region} = "$chr:$start-$end";
	}

	## Coverage ##
	if (-e "$sampleName.coverage.1kbwindow.tsv") {
		foreach my $gene (keys %geneTargets) {
			my($chr,$start,$end) = $geneTargets{$gene} =~ /(\S+)\:(\d+)\-(\d+)/;
			`cat $sampleName.coverage.1kbwindow.tsv | grep $chr > tmp.out`;

			my($totalCov,$covCount);
			open INF,"tmp.out" or die "can't open file: $!";
			while (<INF>) {
				chomp($_);
				my(@D) = split /\t/, $_;
				next unless $D[1] >= $start && $D[1] <= $end;
				++$covCount;
				$totalCov += $D[2];
			}
			close INF;
			my $aveCov = sprintf("%0.0f", $totalCov/$covCount);
			$coverage{$gene} = $aveCov;
			#print "$gene\t$aveCov\n";
			`rm -f tmp.out`;
		}
	} else {
		print "Can't calc coverage as $sampleName.Coverage.1kbwindow.tsv is missing\n";
	}

	print "##_SUMMARY_DATA\n";
	print "Gene\tCoverage\tgeneRegion\tchr\tgeneStart\tgeneEnd\ttargetRegion\n";
	foreach my $gene (keys %geneTargets) {
		my($chr,$start,$end) = $geneTargets{$gene} =~ /(\S+)\:(\d+)\-(\d+)/;
		print "$gene\t$coverage{$gene}\t$geneTargets{$gene}\t$chr\t$start\t$end\t$geneList{$gene}{'regionToAnalyze'}\n";
	}
	print "##\n";


	## spliceAI output ##

		# 2. SNVs -> spliceAI result + AF + caller + coverage
 		# indels -> spliceAI result + AF + caller + coverage
 		# cat spliceai.out | grep -v "#" | awk -F"|" '{ if ($3 > .2 || $4 > .2 || $5 > .2 || $6 > .2) print}'

	print "##_SpliceAI\n";

	my %spliceaiin;
	open INF,"$in_spliceai";
	while (<INF>) {
		next if $_ =~ /^#/;
		chomp($_);
		my(@F) = split /\t/, $_;
		next unless $F[5] =~ /AF/;
		my($callers) = $F[5] =~ /CL\=(\S+)/;
		my($af)      = $F[5] =~ /AF\=(\S+)\;CC/;

		$spliceaiin{$F[0]}{$F[1]}{$F[3]}{$F[4]} = "$af\t$callers";
	}
	close INF;

	open INF,"$out_spliceai";
	while (<INF>) {
		next if $_ =~ /^#/;
		chomp($_);
		my(@F) = split /\t/, $_;
		next unless $F[7] =~ /SpliceAI/;
		my($spliceAIout) = $F[7] =~ /SpliceAI\=(.+)/;

		# only report variants within the region or gene of interest
		my $skip = 1;
		foreach my $gene (keys %geneTargets) {
			my($chr,$start,$end) = $geneTargets{$gene} =~ /(\S+)\:(\d+)\-(\d+)/;
			next unless $chr eq "chr$F[0]";
			$skip = 0 if $F[1] >= $start && $F[1] <= $end;
		}
		next if $skip == 1;
	
		foreach my $foo (split /\,/, $spliceAIout) {
			my(@D) = split /\|/, $foo;
			my $gene = $D[1];
			my $skip = 1;

			#next unless $geneList{$D[1]};

			$skip = 0 if $D[2] >= $spliceaiCutoff;
			$skip = 0 if $D[3] >= $spliceaiCutoff;
			$skip = 0 if $D[4] >= $spliceaiCutoff;
			$skip = 0 if $D[5] >= $spliceaiCutoff;

			next if $skip == 1;

			print "$gene\tchr$F[0]:$F[1]\t$_\t$spliceaiin{$F[0]}{$F[1]}{$F[3]}{$F[4]}\n";
		}
	}
	close INF;
	print "##\n"; 
	print "##_SV_Data\n"; 

	## SVs ##

	foreach my $gene (keys %geneList) {
		print "##_SV_data_for_$gene\n";
		#my $chr = $geneList{$gene}{'chr'};
		#my $targetStart = $geneList{$gene}{'regionStart'};
		#my $targetEnd   = $geneList{$gene}{'regionEnd'};
		my($chr,$geneStart,$geneEnd) = $geneTargets{$gene} =~ /(\S+)\:(\d+)\-(\d+)/;

		foreach my $svFile ($out_sniffles_minimap2,$out_sniffles_ngmlr,$out_svim_minimap2,$out_svim_ngmlr) {
			
			open INF,"$svFile" or die "can't open $svFile: $!";
			while (<INF>) {
				next if $_ =~ /^#/;
				chomp($_);
				my(@F) = split /\t/, $_;
				next unless $chr eq $F[0];
				next if $F[5] <= 2 && $svFile =~ /svim/;
	
				foreach my $pos ($geneStart..$geneEnd) {
					print "$gene\t$chr:$pos\t$svFile\t$_\n" if $pos == $F[1];
				}
			}
			close INF;
		}
	}
 
	print "##\n";
	print "##_Done_with_SV_data\n";
	print "##\n";
}


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
