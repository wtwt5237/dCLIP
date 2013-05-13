#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use File::Temp qw/tempfile/;

##############  read options  #########################

my ($file1,$file2,$pair,$min1,$min2,$temp_dir,$dir,$filter,$step,$mut_type,$max,$precision,$help);

$min1=5; # default values
$min2=5;
$file1="";
$file2="";
$pair="";
$temp_dir=".";
$dir=".";
$step=5;
$filter=10;
$mut_type="T2C";
$max=10;
$precision=0.001;
$help=0;

GetOptions 
(
  'f1=s'=>\$file1,
  'f2=s'=>\$file2,
  'pair=s'=>\$pair,
  'm1=s'=>\$min1,
  'm2=s'=>\$min2,
  'temp=s'=>\$temp_dir,
  'dir=s'=>\$dir,
  'step=s'=>\$step,
  'filter=s'=>\$filter,
  'mut=s'=>\$mut_type,
  'max=s'=>\$max,
  'pre=s'=>\$precision,
  'h'=>\$help
);

if ($help==1) {help();}
if ($file1 eq "" || $file2 eq "") {die "Missing one or both file name(s)!\n"} # check parameters
if ($filter<=1) {die "The filter parameter should be set to be larger than 1!\n";} # process filter parameter
$filter*=$step;

################  check directory  #####################

unless (-d $temp_dir) {mkdir $temp_dir or die "Can't create temp directory!\n";}
unless (-d $dir) {mkdir $dir or die "Can't create output directory!\n";}

################  preprocess data  #####################

print "Preprocessing data\n";

my ($fh,$tempfile1,$tempfile2,$tempfile_p1,$tempfile_p2);

my $perl_path=abs_path($0);
$perl_path=~s/dCLIP\.pl//;
$perl_path="\"".$perl_path."\"";

$tempfile_p1="0";
$tempfile_p2="0";

($fh,$tempfile1)=tempfile(DIR=>$temp_dir,SUFFIX=>".tmp");
($fh,$tempfile2)=tempfile(DIR=>$temp_dir,SUFFIX=>".tmp");

if ($pair eq "")
{
  system("perl ".$perl_path."/cluster.pl \"".$file1."\" \"".$file2."\" ".$min1." ".$min2." \"".$tempfile1."\" ".$step)==0 || die "Clustering tags failed!\n";
  system("perl ".$perl_path."/preprocess.pl ".$step." \"".$tempfile1."\" ".$mut_type." \"".$tempfile2."\"")==0 || die "Writing raw basecoverage file failed!\n";
}else
{
  ($fh,$tempfile_p1)=tempfile(DIR=>$temp_dir,SUFFIX=>".tmp");
  ($fh,$tempfile_p2)=tempfile(DIR=>$temp_dir,SUFFIX=>".tmp");
 
  print "  Merge mates for file1\n"; 
  system("perl ".$perl_path."/merge_pair.pl ".$pair." \"".$file1."\" \"".$tempfile_p1."\"")==0 || die "Merging on file1 failed!\n";
  print "  Merge mates for file2\n";
  system("perl ".$perl_path."/merge_pair.pl ".$pair." \"".$file2."\" \"".$tempfile_p2."\"")==0 || die "Merging on file2 failed!\n";

  system("perl ".$perl_path."/cluster_p.pl \"".$tempfile_p1."\" \"".$tempfile_p2."\" ".$min1." ".$min2." \"".$tempfile1."\" ".$step)==0 || die "Clustering tags failed!\n";   
  system("perl ".$perl_path."/preprocess_p.pl ".$step." \"".$tempfile1."\" ".$mut_type." \"".$tempfile2."\"")==0 || die "Writing raw basecoverage file failed!\n"; 
}

################  initialization  ######################

print "Initialize Hidden Markov Model\n";

my ($M_file,$emission_file,$region_file,$gamma_file,$xis_file,$a_file,$fb_file,$data_file,$files);

($fh,$M_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".M"); # define the file names for diskcache to write to
($fh,$emission_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".emission");
($fh,$region_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".region");
($fh,$gamma_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".gamma");
($fh,$xis_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".xis");
($fh,$a_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".a");
($fh,$fb_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".fb");
($fh,$data_file)=tempfile(DIR=>$temp_dir,SUFFIX=>".data");

$files=join("\" \"",($M_file,$emission_file,$region_file,$gamma_file,$xis_file,$a_file,$fb_file,$data_file));
$files="\"".$files."\"";

system("perl ".$perl_path."/initialize.pl \"".$tempfile2."\" ".$step." ".$filter." ".$files)==0 || die "Initialize Hidden Markov Model failed!\n";

################  HMM  #################################

print "Running Hidden Markov Model\n";

system("perl ".$perl_path."/hmm.pl ".$max." ".$precision." ".$files)==0 || die "Running Hidden Markov Model failed!\n"; 

################  Viterbi  #############################

print "Running Viterbi algorithm\n";

# part of the xis table is rewritten to store result data

system("perl ".$perl_path."/viterbi.pl ".$files)==0 || die "Running Viterbi algorithm failed!\n";

################  write output  ########################

print "Writing output\n";

system("perl ".$perl_path."/output.pl \"".$dir."\" ".$step." \"".$tempfile2."\" ".$files)==0 || die "Writing output failed!\n";

###############  help documents  #######################

sub help
{
  print <<'HELP';

dCLIP: Differential CLIP-Seq datasets analysis tool
version: 1.0
Usage: perl dCLIP.pl [options]
 
Input:
  -f1	The SAM format file of the first condition.
  -f2	The SAM format file of the second condition.
  -pair	If the aligned SAM format files are from single-end experiments, leave this option unset. For paired-end files, set this option to the suffix of the names of forward reads and backward reads. For example, "F3,F5-RNA". 
  -m1	The minimum number of tags for the first condition. All tags from both conditions are pooled, collapsed and overlapped to form clusters. Only clusters with at least m1 tags of the first condition or m2 tags of the second condition will be considered. Default: 5.
  -m2   The minimum number of tags for the second condition. Default: 5.
Directory:
  -temp	The temporary directory to store intermediate files. Default: ".".
  -dir	The folder to store final output files. Default: ".".
Parameters:
  -step	The step size of profiling tag intensities. This controls the resolution of the Hidden Markov Model. Default: 5.
  -filter A filter value used for defining regions with significant binding in both conditions. A higher value will be more conservative in calling differential regions. Should be set >1.  Default: 10.
  -mut	The mutant type(s) of the marker mutations. Can be any one or combination (separated by comma) of "T2C","T2A",...,"A2G","Del","Ins". For example, "T2C,A2G" will include T-to-C and A-to-G mutations as marker mutations. "all" will include all types of mutations. Default: "T2C".
  -max	The maximum number of iterations allowed for the Hidden Markov Model. Default: 10.
  -pre	The precision of the criterion for convergence. Default: 0.001.
Help:
  -h	Print this usage message.

Output: 
  The program will produce 10 output files in the output folder.
  8 of these files named like "File1_Mutant_neg.bedgraph" are bedGraph format files storing the total or mutant tag counts on the + or - strand in the first or second condition. For example, "File1_Mutant_neg.bedgraph" stores the mutant tag count on the - strand in the first condition as a bedGraph format file. The resolution of these files is 1bp.
  "dCLIP_output.txt" stores the detailed information of the raw data and Hidden Markov Model inference results at the resolution of the step bp. 
  "dCLIP_summary.bed" stores the summary of the inference results. This file can be uploaded to UCSC Genome Browser. Red bars are regions where condition 2 has stronger binding than condition 1, while green bars are regions with equal binding strength for both conditions and blue bars are regions where condition 1 has stronger binding than condition 2.

For more information, please refer to the readme file.

HELP
  exit;
}

################  clean up  ############################

END
{
  if ($help==0)
  {
    if ($tempfile_p1 ne "0") 
    {
      unlink($tempfile_p1);
      unlink($tempfile_p2);
    }
    unlink($tempfile1);
    unlink($tempfile2);
    unlink($M_file);
    unlink($emission_file);
    unlink($region_file);
    unlink($gamma_file);
    unlink($xis_file);
    unlink($a_file);
    unlink($fb_file);
    unlink($data_file);
  }
}

exit;
