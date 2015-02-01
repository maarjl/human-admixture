#!/usr/bin/perl
#
use Getopt::Long;
use strict;

#################
## Read in command line arguments

sub usageHelp{
  print "Usage ./neaverage.pl -o outputfile [-l inputrecombfilelist] [-v] [<inputfiles>]\n";
  print "	<inputfiles> are the the .EMprobs.out files produced by ChromoPainter required if -l omitted).\n";
  print "	<inputrecombfilelist> is a two column file, for weighing estimates by recombination distance.\n	The first containing recombination maps,\n";
  print "	and the second containing .EMprobs.out files,\n";
  print "	Example:\nchromosome1.recombfile chromosome1.chromopainter.EMprobs.out\nchromosome2.recombfile chromosome2.chromopainter.EMprobs.out\n";
  print "	outputfile will contain correct ChromoPainter switches for \n	Ne and mutation rates, i.e.\n";
  print "		\"-n <Ne> -M <mutationrate>\"\n";
  print "	which you should specify verbatim in future runs of ChromoPainter.\n";
  print "	(don't forget to OMIT the EM estimation arguments\n";
  print "	     such as \"-i 10 -in\")\n\n";
  print "	-v: Enable verbose mode\n\n";
  print "	EXAMPLE USAGE: (when all EM files are given equal weight)\n";
  print "		./neaverage.pl -o neest.txt <chromopainterroot>*EMprobs.out\n";
  print "	where <chromopainterroot>*EMprobs.out are all files produced by\n	ChromoPainter run in E-M mode, e.g.\n";
  print "		chromopainter -i 10 -in -iM -g <phasefile>\n			-r <recombfile> -o <chromopainterroot><index>\n";
  print "	and there are potentially multiple runs, indexed by <index>.\n";
  print "	EXAMPLE USAGE: (when different EM files have different weight)\n";
  print "		for i in `seq 1 22`; do echo \"my.chromosome1.recombfile my.chromosome1.\"\n";
  print "		./neaverage.pl -o neest.txt -l emfilelist.txt\n";
  die("\n");
}
sub usageHelpMessageDie{
  print "$_[0]\n";
  usageHelp();
  exit;
}

sub trim($)
{
        my $string = shift;
        $string =~ s/^\s+//;
        $string =~ s/\s+$//;
        return $string;
}


my @emfiles;  ## List of EM files to use
my @recfiles; ## List of recombination maps
my @fileweights;  ## List of EM files to use
my $verbose = 0; # quiet = 1 - verbose
my $help = 0; # detailed help
my $outputfile = "";# outputfile location
my $reclistfile = "";# recombination file location
GetOptions ('help' => \$help, 
	    'verbose' => \$verbose,
         'o|out=s' => \$outputfile,
		 'l|list=s' => \$reclistfile);


if($outputfile eq "" || ($reclistfile eq "" && scalar(@ARGV)==0)) {
  print "Outputfile muist be specified, and either a recombinationlist file or other files on the command line\n";
  $help=1;
}
if($help==1) {
	usageHelpMessageDie();
}

## Set up the outputfile (in case this crashes)
open OUTPUTFILE, ">", $outputfile or die "Cannot write to outputfile $outputfile!\n$!";

## construct the list files to process
if($reclistfile eq ""){
		for(my $i=0 ; $i< scalar(@ARGV) ; $i++){
			push @emfiles, $ARGV[$i];
			push @fileweights, 1.0/scalar(@ARGV);
		}
}else{
  open RECLIST, $reclistfile or die "Error: list file $reclistfile doesn't exist!\n$!";
  while(my $line=trim(<RECLIST>)) {
	my @tarr=split(/\s+/, $line);
	if(scalar(@tarr)>=2){
		push @recfiles, @tarr[0];
		push @emfiles, @tarr[1];
	}
  }
  close RECLIST;
}

## Construct recombination weights, if applicable
if(scalar(@recfiles)>0){
	foreach(@recfiles){
		open RECFILE, $_ or die "Error: recombination file $_ doesn't exist!\n$!";
		my $lastbploc=-1;
		my $lastrecrate=-1;
		my $recdistance=0;
		while(my $line=trim(<RECFILE>)) {
			my @tarr=split(/\s+/, $line);
			my $bpdist = $tarr[0] - $lastbploc;
			if($bpdist<0 || $lastbploc<0) { $bpdist = 0; }
			$recdistance += $lastrecrate * $bpdist;
			$lastrecrate=$tarr[1];
			$lastbploc=$tarr[0];
		}
		close RECFILE;
		push @fileweights, $recdistance;
	}
	my $totaldistance = 0;
	foreach(@fileweights){ $totaldistance+= $_; }
	for(my $i=0; $i<scalar(@fileweights); $i++) { $fileweights[$i]/=$totaldistance; }
	if($verbose) {
		for(my $i=0; $i<scalar(@fileweights); $i++) { print ("File $i ($recfiles[$i]) has weight $fileweights[$i]\n"); }
	}
}

## Now we search each EMfile for the final iteration FOR EACH INDIVIDUAL. So we get many observations per file (potentially)
my @allne;
my @allmu;
my @allweight;
for (my $file=0;$file < scalar(@emfiles); $file++){
	if($verbose){
		print "Processing file $file ($emfiles[$file]) with weight $fileweights[$file]\n";
	}
	open TMPFILE, $emfiles[$file] or die "Error: EMfile $emfiles[$file] doesn't exist!\n$!";
	my @tmpfile=<TMPFILE>;
	my $tmpfilesize =scalar(@tmpfile);
	my $lastline=@tmpfile[0];
	my @lastlineS = split(/\s+/, trim($lastline));
	my $numcols=0; # number of columns detected in the EM file
	for(my $line=1;$line<$tmpfilesize;++$line){
		my $nextline=@tmpfile[$line];
		my @nextlineS = split(/\s+/, trim($nextline));
		if(@nextlineS[0] eq "IND"){
			$numcols=scalar(@lastlineS);
			push(@allne,$lastlineS[$numcols - 2]);
			push(@allmu,$lastlineS[$numcols - 1]);
			push(@allweight,$fileweights[$file]);
			if($verbose){
				print "Found Ne=$lastlineS[$numcols - 2], Mu=$lastlineS[$numcols - 1]\n";
			}
		}
		@lastlineS=@nextlineS;
		$lastline=$nextline;
	}
	close TMPFILE;
	push(@allne,$lastlineS[$numcols - 2]);
	push(@allmu,$lastlineS[$numcols - 1]);
	push(@allweight,$fileweights[$file]);
	if($verbose){
		print "Found Ne=$lastlineS[$numcols - 2], Mu=$lastlineS[$numcols - 1]\n";
	}
}

###############
## Now we just compute the averages (which are now straight weighted averages)
my $avene=0;
my $avemu=0;
my $weightsum=0;
my $numvalues=scalar(@allne);
for(my $i=0; $i<scalar(@allne); $i++){
	$avene+=$allne[$i]*$allweight[$i];
	$avemu+=$allmu[$i]*$allweight[$i];
	$weightsum+=$allweight[$i];
}
$avene=$avene/$weightsum;
$avemu=$avemu/$weightsum;
#print("Considering $numvalues values and obtaining Ne = $avene, and global mutation rate mu = $avemu \n");


######### PRINT THE MEAN
print OUTPUTFILE "-n $avene -M $avemu\n";
close OUTPUTFILE;

