#!/usr/bin/perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::LiteF;
use PDL::Stats::GLM;
use PDL::Stats::Basic;
use PDL::DiskCache;

my ($file,$step,$filter,@files)=@ARGV;
my ($region_id,$tag1,$tag2,$mut1,$mut2);
my ($tag1_subset,$tag2_subset,$A,$M,$M_subset,%fit,$fit,$cons);
my ($sd_l,$sd_r,$sd_c,$ave_l,$ave_r,$ave_c,$prod,$prob,$p,$temp,$min);
my ($m1,$m2,$mu);
my ($emission,$n,$line,$pi);
my ($i,$j,@starts,@ends,$region);
my ($gamma,$xis,$a,$fb,$cache);
 
#############  read and collapse data  ###########################

print "  Read basecoverage file\n";

($region_id,$tag1,$mut1,$tag2,$mut2)=rcols $file,0,[4..$step+3],[$step+4..3+2*$step],[4+2*$step..3+3*$step],[4+3*$step..3+4*$step],{PERLCOLS=>[0],COLSEP=>"\t"} or die "Reading basecoverage file failed!\n";

$tag1=$tag1->xchg(0,1)->sumover;
$mut1=$mut1->xchg(0,1)->sumover;
$tag2=$tag2->xchg(0,1)->sumover;
$mut2=$mut2->xchg(0,1)->sumover;

#############  normalize tag density  ############################

print "  Normalize tag density data\n";

$tag1+=$filter; # log transform
$tag2+=$filter;
$tag1=log(inplace $tag1);
$tag2=log(inplace $tag2);

my $subset=($tag1>log(2*$filter))+($tag2>log(2*$filter)); # take out common peaks
$subset=($subset==2);
if (sum($subset)==0 || sum(1-$subset)==0) {die "Sequencing is too deep or not deep enough!\n";}
$tag1_subset=where $tag1,$subset;
$tag2_subset=where $tag2,$subset;

$A=$tag1_subset+$tag2_subset; # fit regression line
$M=$tag1_subset-$tag2_subset;
%fit=$M->ols($A,{const=>1,plot=>0});
$fit=$fit{"b"};

$A=$tag1+$tag2; # normalize M
$M=$tag1-$tag2;
$A=$A*at($fit,0)+at($fit,1);
$M-=$A; # M now contains the normalized M values

###############  fit model  ########################################

print "  Fit mixture model\n";

$M_subset=where $M,(1-$subset); # fit sd
$sd_c=($M_subset->stdv)*0.9;
$sd_l=$sd_c; 
$sd_r=$sd_c;

$n=$M->nelem; # assign some values
$pi=3.1415926;
$cons=-1/(2*$sd_c**2);
$prod=(at($M x transpose($M),0,0)/$n-$sd_c**2)/2;
if ($prod<0) {die "Model fitting failed! Use a larger filter parameter!\n";}
$prob=zeroes(100);
$j=0; # this is a flag

foreach $i (1..100) # grid approximation
{
  $mu=($i-50)/100+3*$sd_c;
  $p=$prod/$mu**2;
  if ($mu<=0 || $p>=0.5) 
  {
    $prob->slice($i-1).=-1e300;
    next;
  }
  $j=1;
  $temp=$p*exp($cons*($M-$mu)**2)+$p*exp($cons*($M+$mu)**2)+(1-2*$p)*exp($cons*$M**2);
  $temp=log($temp*0.9999999+0.00000005);
  $prob->slice($i-1).=sum($temp);
}

if ($j==0)
{
  print "Fitting failed! Will just use p=0.05 and mu=0.5!\n";
  $mu=0.5;
}else
{
  $i=maximum_ind($prob);
  $mu=(at($i,0,0)-49)/100+3*$sd_c;
}

$ave_l=-$mu; # assign mu values
$ave_r=$mu;
$ave_c=0;

# left component, state 1, 1<2, M<0
# center component, state 2, 1=2, M=0
# right component, state 3, 1>2, M>0

###############  initialize all matrices  ##########################

print "  Initialize matrices\n";

### initialize emission table ###

$emission=zeroes($n,3); # define emission table, the columns corresponds to left, middle and right component in the mixture model
$emission+=$M; # add x value
$temp=pdl[$ave_l,$ave_c,$ave_r]; # minus mean
$emission=$emission->xchg(0,1);
$emission-=$temp;
$emission->inplace->power(2,0); # power 2 
$temp=pdl[-2*$sd_l**2,-2*$sd_c**2,-2*$sd_r**2]; # devide by -2*sigma^2
$emission/=$temp;
$emission->inplace->exp; # exponential
$temp=pdl[sqrt(2*$pi)*$sd_l,sqrt(2*$pi)*$sd_c,sqrt(2*$pi)*$sd_r];
$emission/=$temp;

if (min($emission)<1e-300) {print "    Warning: very small probability values exist! They may cause the computation to go out of the dynamic range!\n";}

### initialize region table ###

@starts=();
$starts[0]=0;
@ends=();

# region table has four columns
# the first column stores start sites
# the second column stores non-start sites
# the third column stores end sites
# the fourth column stores non-end sites

foreach $i (0..$n-2)
{
  if ($$region_id[$i]!=$$region_id[$i+1])
  {
    push @starts,$i+1;
    push @ends,$i;
  }
}

$ends[$#ends+1]=$n-1;

$region=zeroes(4,$n);
$region->inplace->yvals;
$region->slice("0,:").=-1;
$region->index2d(0,pdl[@starts]).=pdl @starts;
$region->index2d(1,pdl[@starts]).=-1;
$region->slice("2,:").=-1;
$region->index2d(2,pdl[@ends]).=pdl @ends;
$region->index2d(3,pdl[@ends]).=-1;

### initialize gamma table  ###

$gamma=zeroes(3,$n);
$gamma.=1/3;

### initialize xis table  ###

$xis=zeroes(9,$n);
# 0->0,0->1,0->2,1->0,1->1,1->2,2->0,2->1,2->2

for $i (0..2)
{
  for $j (0..2)
  {
    $xis->slice("(".($i*3+$j)."),0:".($n-2)).=$gamma->slice("(".$i."),0:".($n-2))*$gamma->slice("(".$j."),1:".($n-1));
  }
}

###  initialize transition kernel  ###

$a=zeroes(3,3);
$a->diagonal(0,1)++;

###  initialize forward-backward table  ###

$fb=zeroes(6,$n);

###############  tie data to file  #################################

# @files contains the file names for $M_file,$emission_file,$region_file,$gamma_file,$xis_file,$a_file,$data_file

print "  Write data to disk\n";

open(FILE,">".$files[7]); # cache scalars
print FILE join(" ",($ave_l,$ave_c,$ave_r,$sd_l,$sd_c,$sd_r,$n));
close(FILE);

@files=@files[0..6];
$cache=diskcache(\@files,{rw=>1}); # cache matrices
$cache->[0]=$M;
$cache->[1]=$emission;
$cache->[2]=$region;
$cache->[3]=$gamma;
$cache->[4]=$xis;
$cache->[5]=$a;
$cache->[6]=$fb;



