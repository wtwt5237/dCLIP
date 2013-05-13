#!/usr/bin/perl
use strict;
use warnings;
use POSIX qw/floor/;
use List::Util qw/min max/;
$SIG{__WARN__}=sub{die "SAM file format error!\n"};

# $step is the step resolution desired
# $file_cluster is the combined cluster file
# $mut_type is the type of mutation wanted, bracketed by "" and separated by ",", e.g. "T2C,A2G"
# $file_count is the output file

my ($step,$file_cluster,$mut_type,$file_count)=@ARGV;
my ($chr,$strand,$tag_number);
my (@items,@starts,@ends,@muts1,@muts2,@bins1,@bins2,@bins_mut1,@bins_mut2,@file);
my ($start,$end,$i,$j);

@bins1=();
@bins2=();
@bins_mut1=();
@bins_mut2=();

print("  Calculate basecoverage\n");

open(FILE_IN,$file_cluster) or die "read cluster file ".$file_cluster." failed!\n";
open(FILE_OUT,"> ".$file_count) or die "can't write to file ".$file_count."!\n";

my $region_id=1;

while (<FILE_IN>)
{
  ($chr,$strand,$tag_number)=split("\t",$_);
  
  @muts1=();
  @muts2=();
  @starts=();
  @ends=();
  @file=();

  foreach $i (1..$tag_number) # read tag and mutant information
  {
    @items=split("\t",<FILE_IN>);
    if ($items[1]>0)
    {
      push @file,1;
      read_tag($mut_type,$strand,$items[1],$items[2],$items[3],$items[4],\@starts,\@ends,\@muts1);
    }else
    {
      push @file,-1;
      read_tag($mut_type,$strand,-$items[1],$items[2],$items[3],$items[4],\@starts,\@ends,\@muts2);
    }
  }
    
  $start=$step*floor(min(@starts)/$step);  # define cluster region
  $end=$step*(floor(max(@ends)/$step)+1)-1;

  foreach ($start..$end) {$bins1[$_-$start]=0;} # calculate base coverage on each base
  foreach ($start..$end) {$bins2[$_-$start]=0;}
  foreach ($start..$end) {$bins_mut1[$_-$start]=0;} # calcualte mutant base coverage on each base
  foreach ($start..$end) {$bins_mut2[$_-$start]=0;} 

  foreach $i (0..($tag_number-1)) 
  {
    if ($file[$i]==1)
    {
      foreach (($starts[$i]-$start)..($ends[$i]-$start)) {$bins1[$_]++;}
    }else
    {
      foreach (($starts[$i]-$start)..($ends[$i]-$start)) {$bins2[$_]++;}
    }
  }

  map {$bins_mut1[$_-$start]++;} @muts1;
  map {$bins_mut2[$_-$start]++;} @muts2;

  foreach $i (0..(($end-$start+1)/$step-1)) # print output
  {
    print FILE_OUT $region_id."\t".$chr."\t".($i*$step+$start)."\t".$strand."\t";
      
    # print information for file 1
    print FILE_OUT join("\t",@bins1[($i*$step)..($i*$step+$step-1)])."\t";
    print FILE_OUT join("\t",@bins_mut1[($i*$step)..($i*$step+$step-1)])."\t";      

    # print information for file 2
    print FILE_OUT join("\t",@bins2[($i*$step)..($i*$step+$step-1)])."\t";
    print FILE_OUT join("\t",@bins_mut2[($i*$step)..($i*$step+$step-1)])."\n";     
  }

  $region_id++;
}

close(FILE_IN);
close(FILE_OUT);

sub read_tag # read tag position and mutant position
{
  my ($mut_type,$strand,$pos,$CIGAR,$seq,$MD,$starts_ref,$ends_ref,$muts_ref)=@_;  
  my (@items,$len,@mut_pos);

  $len=0;
  while ($CIGAR=~/([0-9]+)([MD])/g) {$len+=$1;}
  push @$starts_ref,$pos;
  push @$ends_ref,($pos+$len-1); # record position information       
 
  @mut_pos=();
  read_mut($mut_type,$CIGAR,$seq,$MD,\@mut_pos,$strand);
  map {$_+=($pos-1);} @mut_pos;
  push @$muts_ref,@mut_pos; # record mutation information
}

sub read_mut # read mutatant position
{
  my ($mut_type,$CIGAR,$seq,$MD,$mut_pos_ref,$strand)=@_;
  my @mut_pos=();
  my $ref_pos=0;
  my $tag_pos=0;
  my ($regex,$match,$temp);

  if ($CIGAR=~/([0-9]+)S.*[0-9]+M/) {$seq=substr($seq,$1);} # offset soft-clipping
  $CIGAR=~s/[0-9]+H//g; # offset hard-clipping

  while ($CIGAR=~/([0-9]+)([MDI])/g)
  {
    if ($2 eq "M") 
    {
      $ref_pos+=$1;
      $tag_pos+=$1;
    }elsif ($2 eq "I")
    {
      {if ($mut_type=~/Ins|all/) {push @$mut_pos_ref,($ref_pos+1);}}
      substr($seq,$tag_pos,$1)="";
      $tag_pos+=$1;
    }else
    { 
      {if ($mut_type=~/Del|all/) {push @$mut_pos_ref,($ref_pos+1);}}
      $ref_pos+=$1;
    }
  }  

  $ref_pos=0;
  $tag_pos=0;

  while ($MD=~/([0-9]+|[ACGTN]|\^[ACGTN]+)/g) 
  {
    $match=$1;
    if ($match=~/[0-9]+/) 
    {
      $ref_pos+=$match;
      $tag_pos+=$match;      
    }elsif ($match=~/^[ACGTN]$/)
    {
      $ref_pos+=1;
      $temp=substr($seq,$tag_pos,1);

      if ($strand eq "-")
      {
        $match=transform($match);
        $temp=transform($temp);
      }

      $temp=$match."2".$temp."|all";
      $regex=qr/$temp/;
      $tag_pos+=1;
      if ($mut_type=~$regex) {push @$mut_pos_ref,$ref_pos;}
    }else
    {
      $ref_pos+=length($match)-1;
    }
  }
}

sub transform # negative strand to positive strand
{
  my $base=$_[0];
  
  if ($base eq "A")
  {
    $base="T";
  }elsif ($base eq "T")
  {
    $base="A";
  }elsif ($base eq "C")
  {
    $base="G";
  }elsif ($base eq "G")
  {
    $base="C";
  }

  return($base);
}
