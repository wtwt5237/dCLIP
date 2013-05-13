#$!lusr/bin/perl
use strict;
use warnings;
use PDL;
use PDL::DiskCache;

my (@files)=@ARGV;
my $cache;
my ($ave_l,$ave_c,$ave_r,$sd_l,$sd_c,$sd_r,$n);
my ($M,$emission,$region,$gamma,$xis,$a,$fb,$delta,$phi,$states,$prob);
my ($starts,$ends,$starts_wo,$ends_wo);
my ($temp,$i);

##########  read data from disk  ####################

print "  Read data from disk\n";

open(FILE,$files[7]); # read scalars
($ave_l,$ave_c,$ave_r,$sd_l,$sd_c,$sd_r,$n)=split(" ",<FILE>);
close(FILE);

@files=@files[0..6];
$cache=diskcache(\@files,{rw=>1}); # cache matrices
($M,$emission,$region,$gamma,$xis,$a,$fb)=@{$cache};

#########  check inf values  ########################

if (sum($gamma) eq "inf" || sum($xis) eq "inf" || sum($fb) eq "inf" || sum($a) eq "inf") 
{
  die "  Bad values generated!";
}

#########  initialize delta and phi  ################

print "  Initialize Viterbi algorithm\n";

$starts=$region->slice("0,:"); # coordinates
$starts=where $starts,$starts>=0;
$starts_wo=$region->slice("1,:");
$starts_wo=where $starts_wo,$starts_wo>=0;
$ends=$region->slice("2,:");
$ends=where $ends,$ends>=0;
$ends_wo=$region->slice("3,:");
$ends_wo=where $ends_wo,$ends_wo>=0;
$ends_wo=$ends_wo->slice("-1:0:-1"); # $ends_wo is in reverse order

$delta=zeroes(3,$n);
$delta->dice_axis(1,$starts).=$gamma->dice_axis(1,$starts)*$emission->dice_axis(1,$starts);
$phi=zeroes(3,$n);
$states=zeroes($n);
$prob=zeroes($n);

########  dynamic programming  ######################

print "  Viterbi algorithm\n";

$a=transpose($a);

foreach $i ($starts_wo->list)
{
  $temp=$a*$delta->slice(":,".($i-1));
  $phi->slice(":,$i").=$temp->maximum_ind;
  $temp=$temp->maximum*$emission->slice(":,$i");
  $temp/=sum($temp);
  $delta->slice(":,$i").=$temp;
}

$states->index($ends).=$delta->dice_axis(1,$ends)->maximum_ind;

foreach $i ($ends_wo->list)
{
  $states->slice("$i").=at($phi,at($states,$i+1),$i+1);
}

$prob->slice(":").=$gamma->index2d($states->slice(":"),pdl((0..$n-1)));

#########  store results  ##############################

# part of the xis table is rewritten to store results

$xis->slice("(0),:").=$states;
$xis->slice("(1),:").=$prob;





