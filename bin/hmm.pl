#!/usr/bin/perl
use strict;
use warnings;
use PDL;
use PDL::DiskCache;

my ($max,$err,@files)=@ARGV;
my $cache;
my ($ave_l,$ave_c,$ave_r,$sd_l,$sd_c,$sd_r,$n);
my ($M,$emission,$region,$gamma,$xis,$a,$fb,$for,$back);
my ($scales,$iterat,$a_pre,$i,$j,$temp1,$temp2);
my ($starts,$ends,$starts_wo,$ends_wo);

##########  read data from disk  ####################

print "  Read data from disk\n";

open(FILE,$files[7]); # read scalars
($ave_l,$ave_c,$ave_r,$sd_l,$sd_c,$sd_r,$n)=split(" ",<FILE>);
close(FILE);

@files=@files[0..6];
$cache=diskcache(\@files,{rw=>1}); # cache matrices
($M,$emission,$region,$gamma,$xis,$a,$fb)=@{$cache};

#########  set initial values  #######################

print "  Set initial values\n";

$iterat=0;
$scales=zeroes($n);
$a_pre=zeroes(3,3);

$starts=$region->slice("0,:"); # coordinates
$starts=where $starts,$starts>=0;
$starts_wo=$region->slice("1,:");
$starts_wo=where $starts_wo,$starts_wo>=0;
$ends=$region->slice("2,:");
$ends=where $ends,$ends>=0;
$ends_wo=$region->slice("3,:");
$ends_wo=where $ends_wo,$ends_wo>=0;
$ends_wo=$ends_wo->slice("-1:0:-1"); # $ends_wo is in reverse order

#########  start HMM  ################################

$|=1;
print "  Start iterations\n    ";

while (sum(abs($a-$a_pre)>$err)>0 && $iterat<$max)
{
  ### update $a ###

  $a_pre.=$a;

  for $i (0..2)
  {
    for $j (0..2)
    {
      $temp1=$gamma->index2d($i,$ends_wo); # the sum of probability from state i
      $temp2=$xis->index2d($i*3+$j,$ends_wo);  # the sum of probability from state i to state j
      $a->slice("($j),($i)").=sum($temp2)/sum($temp1); # from the second index of $a to the first index of $a
    }
  }

  ### update forward ###

  $for=$fb->slice("0:2,:");

  $for->dice_axis(1,$starts).=$gamma->dice_axis(0,pdl(0,1,2))->dice_axis(1,$starts)*$emission->dice_axis(0,pdl(0,1,2))->dice_axis(1,$starts); 
  # calculate forward prob at start positions
  $scales->index($starts).=1/$for->dice_axis(1,$starts)->sumover; # calculate scaling parameters at start positions
  $for->dice_axis(1,$starts)->xchg(0,1)*=$scales->index($starts); # scale

  foreach $i ($starts_wo->list)
  {
    $temp1=$emission->slice(":,$i")*($for->slice(":,".($i-1)) x $a);
    $temp2=1/sum($temp1);
    $temp1*=$temp2;
    $for->slice(":,$i").=$temp1; # update value
    $scales->slice("$i").=$temp2; # scale
  } 

  ### update backward ###

  $back=$fb->slice("3:5,:");
  $a=transpose($a);  

  $back->dice_axis(1,$ends)->xchg(0,1).=$scales->index($ends); # calculate backward prob at end sites

  foreach $i ($ends_wo->list)
  {
    $temp1=($emission->slice(":,".($i+1))*$back->slice(":,".($i+1))) x $a;
    $temp1*=$scales->slice("$i");
    $back->slice(":,$i").=$temp1;
  }

  $a=transpose($a);

  ### update gamma ###

  $gamma.=$for*$back;
  $gamma->xchg(0,1)/=$scales;

  ### update xis ###

  for $i (0..2)
  {
    for $j (0..2)
    {
      $xis->slice("(".($i*3+$j)."),0:-2").=at($a,$j,$i)*$emission->slice("($j),1:-1")*$for->slice("($i),0:-2")*$back->slice("($j),1:-1");
    }
  }
   
  ### advance iteration ###

  $iterat++;
  print ">";
}

print "\n";
























