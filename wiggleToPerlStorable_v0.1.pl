#!/usr/bin/perl -w
$|++;
use strict;
use File::Path;
use Time::HiRes qw( time );
use Storable;
use File::Basename;
use File::Spec::Functions qw(rel2abs);
use File::Copy;

######################################################################################################################################################
#
#	Description
#		This is a perl script to convert a wiggle pileup to perl storables folder
#
#	Input
#		--wigglePath=				path; path of the wiggle file;
#		--refFastaPath=				reference genome sequence, used for getting the name of the contig, as the reference to look for contig data in the pileup files;
#		--outDir=					output directory; default = ./wiggleMerger/
#
#	Output
#
#	Usage
#		./wiggleToPerlStorable_v0.1.pl --wigglePath=/Volumes/B_MPro2TB/NGS/results/plasmodium/full_Spliced/E99_ctpsm_s7_20hr/IGVTools/finalSAM.combined.sorted.wig.gz --refFastaPath=/Volumes/A_MPro2TB/softwareForNGS/myPerlScripts/pipeLines/marathon/resources/genome/Pf3D7_01_v3_JAN2012_withMitoPlstd.fa
#
#	History:
#
#		V0.1
#			debut, but built based on wiggleMerger_v0.1.pl;
#
######################################################################################################################################################
#
#==========================================================Main body starts==========================================================================#

#----------Read parameters ----------#
my ($wigglePath, $refFastaPath, $outDir) = &readParameters();
printCMDLogOrFinishMessage("CMDLog");

#----------Read contig names
my $cntgLenHsh_ref = readContigLengthFromFasta();

#----------Read wiggle
storeWigToPerlStorable($wigglePath, $cntgLenHsh_ref, $outDir);

printCMDLogOrFinishMessage("finishMessage");

exit;

#========================================================= Main body ends ===========================================================================#

########################################################################## readParameters
sub readParameters {
	
	my $dirPath = dirname(rel2abs($0));
	my $outDir = "$dirPath/wiggleToPerlStorable/";

	foreach my $param (@ARGV) {
		if ($param =~ m/--wigglePath=/) {$wigglePath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--refFastaPath=/) {$refFastaPath = substr ($param, index ($param, "=")+1);}
		elsif ($param =~ m/--outDir=/) {$outDir = substr ($param, index ($param, "=")+1);} 
	}
	
	#---check the files
	open (TEST, "$refFastaPath") || die "Can't open refFastaPath\n"; close TEST;
	open (TEST, "$wigglePath") || die "Can't open wigglePath\n"; close TEST;

	chop $outDir if ($outDir =~ m/\/$/); #---remove the last slash
	system "mkdir -pm 777 $outDir";
	
	return ($wigglePath, $refFastaPath, $outDir);
}
########################################################################## printCMDLogOrFinishMessage
sub printCMDLogOrFinishMessage {

	my $CMDLogOrFinishMessage = $_[0];
	
	if ($CMDLogOrFinishMessage eq "CMDLog") {
		#---open a log file if it doesnt exists
		my $scriptNameXext = $0;
		$scriptNameXext =~ s/\.\w+$//;
		open (CMDLOG, ">>$scriptNameXext.cmd.log.txt"); #---append the CMD log file
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		my $dirPath = dirname(rel2abs($0));
		my ($scriptName, $callScriptPath, $scriptSuffix) = fileparse($0, qr/\.[^.]*/);

		print CMDLOG "[".$runTime."]\t"."$dirPath/$scriptName$scriptSuffix ".(join " ", @ARGV)."\n";
		close CMDLOG;
		print "\n=========================================================================\n";
		print "$0 starts running at [$runTime]\n";
		print "=========================================================================\n\n";

	} elsif ($CMDLogOrFinishMessage eq "finishMessage") {
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
		my $runTime = sprintf "%04d-%02d-%02d %02d:%02d", $year+1900, $mon+1,$mday,$hour,$min;	
		print "\n=========================================================================\n";
		print "$0 finished running at [$runTime]\n";
		print "=========================================================================\n\n";
	}
	
}
########################################################################## readContigLengthFromFasta
sub readContigLengthFromFasta {

	#print "Reading $refFastaPath for contig names.\n";
	open (INFILE, $refFastaPath);
	my (%cntgLenHsh, $seqName, $length, $seq);

	chomp (my $curntLine = <INFILE>); #get the first line
	while (my $nextLine = <INFILE>) {
		chomp $nextLine;
			
		#---Only two types of line in current line, the header or seq
		if ($curntLine =~ m/^>/) {#-- header line
			$seqName = $curntLine;
			$seqName =~ s/>//g; #---remove space
			$cntgLenHsh{$seqName} = 0;

		} else {#--seq line
			$cntgLenHsh{$seqName} = $cntgLenHsh{$seqName} + length ($curntLine);
		}
			
		#---check if next line has a > or that's the end of file
		if ($nextLine =~ m/^>/) {
			$seqName = $nextLine;
			$seqName =~ s/>//g; #---remove space
			$cntgLenHsh{$seqName} = 0;

		} elsif (eof(INFILE)) {#---this is the last line
			$cntgLenHsh{$seqName} = $cntgLenHsh{$seqName} + length ($nextLine);
		}
		
		#---next line becomes current line
		$curntLine = $nextLine;
	}
	close INFILE;

	my $contigNum = keys %cntgLenHsh;

	#print "Totally $contigNum contig names stored.\n";
	
	#---print to check
	#foreach my $cntg (sort {$a cmp $b} keys %cntgLenHsh) {print $cntg."\t".$cntgLenHsh{$cntg}."\n";}

	return (\%cntgLenHsh);
}
########################################################################## readPileupFromIndividualSamples
sub storeWigToPerlStorable {
	
	my ($wigglePath, $cntgLenHsh_ref, $outDir) = @_;

	my %cntgLenHsh = %{$cntgLenHsh_ref};

	#---check for pigz
	my $compresser = "gzip";
	my $stdout = `pigz --version 2>&1`;
	if ($stdout !~ m/command not found/) {
		$compresser = 'pigz';
	}
	
	#----create empty ary
	my ($wiggleName, $wiggleDir, $wiggleSuffix) = fileparse($wigglePath, qr/\.[^.]*/);
	my $cntgCovStroableDir = "$outDir/$wiggleName";
	system "mkdir -pm 777 $cntgCovStroableDir";
	
	my %indivSampleCntgCovIdxHsh = ();
	foreach my $cntg (keys %cntgLenHsh) {
		print "Creating empty storabe for $cntg         \r";
		my $cntgLen = $cntgLenHsh{$cntg};
		my @cntgCovAry = ();
		foreach (1..$cntgLen) {
			push @cntgCovAry, undef;
		}
		my $cntgCovStroableName = "$cntg.ary.pls";
		my $cntgCovStroablePath = "$cntgCovStroableDir/$cntgCovStroableName";
		$indivSampleCntgCovIdxHsh{$cntg} = $cntgCovStroableName;
		store(\@cntgCovAry, "$cntgCovStroablePath");
	}
	my $indivSampleCntgCovIdxHshPath = "$cntgCovStroableDir/index.hsh.pls";
	store(\%indivSampleCntgCovIdxHsh, "$indivSampleCntgCovIdxHshPath");

	#------read the wig	
	if ($wigglePath =~ m/\.gz$/) {
		open (WIGFILE, "$compresser -dc $wigglePath |");
	} elsif ($wigglePath =~ m/\.wig$/) {
		open (WIGFILE, "$wigglePath");
	}
	
	my $cntg = "initialize";
	my @cntgCovAry = ();
	my @mergeCntgCovAry = ();
	my $cntgCovStroablePath = "";
	my $mergeCntgCovStroablePath = "";
	my $span = 1;#---default value = 1
	my $procPos = 0;
	my $cntgProc = 0;
	while (my $theLine = <WIGFILE>) {
		chomp $theLine;
		next if ($theLine =~ m/^[\#|track type]/);
		
		if (($theLine =~ m/^variableStep/) or (eof(WIGFILE))) {
			if ($cntg ne "initialize") {#---not the first contig, report the previous contig, include the eof(WIGFILE) 
				print "Writing $cntg storable                                \r";
				store(\@cntgCovAry, "$cntgCovStroablePath");
				system ("$compresser -f $cntgCovStroablePath"); #---compress right after install
				$cntgProc++;
			}
			
			if ($theLine =~ m/^variableStep/) {#----new contig, include the first contig but not eof(WIGFILE)
				if ($theLine =~ m/chrom=(\S+) span=(\d+)/) {
					$cntg = $1;
					$span = $2;
					
					$cntgCovStroablePath = "$cntgCovStroableDir/$indivSampleCntgCovIdxHsh{$cntg}";

					print "Retrieving $cntg storable                          \r";
					@cntgCovAry = @{retrieve($cntgCovStroablePath)};

				} else {
					die "wiggle file variableStep line doesn't have the chrom=(\\S+) span=(\\d+) format\n";
				}
			}
			
		} else {
			my ($pos, $samplePlusCov, $sampleMinusCov) = split /\s+/, $theLine;
			$procPos++;
			die "the wiggle file is not strand specific\n" if not defined $sampleMinusCov;
			my $index = $pos - 1;
			$cntgCovAry[$index] = join ",", (int $samplePlusCov, int $sampleMinusCov);
			print "Processed $cntgProc contig and $procPos positions        \r" if (not($procPos % 1000));
		}
	}
	close WIGFILE;
	print "Read wiggle.......................done\n";
	
}
