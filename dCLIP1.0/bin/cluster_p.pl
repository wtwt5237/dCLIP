#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max min/;
$SIG{__WARN__}=sub{die "SAM file format error!\n"};

my ($file_m1,$file_m2,$min1,$min2,$file_cluster,$step)=@ARGV;
my %tags;
my ($temp,$chr,@pos_start,@cluster,$i,$end,$strand,$m1,$m2);

##################  cluster tags  #########################

print("  Cluster tags\n");

read_tags(\%tags,$file_m1,"");
read_tags(\%tags,$file_m2,"-");

sub read_tags
{
  my ($tags_ref,$file,$prefix)=@_;
  my ($len_f,$len_r,$start,$end,$count,@items,$temp);

  open(FILE_IN,$file) or die "read file ".$file." failed!\n";

  while (<FILE_IN>) # read and compute length of each tag
  {
    @items=split("\t",$_);

    $len_f=0;
    while ($items[4]=~/([0-9]+)([MD])/g) {$len_f+=$1;}    
    $len_r=0;
    while ($items[8]=~/([0-9]+)([MD])/g) {$len_r+=$1;}
  
    if ($items[6]=~/-1/) {$items[6]=$len_f;}
    if ($items[10]=~/-1/) {$items[10]=$len_r."\n";}

    $start=min($items[3],$items[7]);
    $end=max($items[3]+$len_f-1,$items[7]+$len_r-1);  

    $count=()=$items[4]=~/(I)/g;
    $count+=()=$items[6]=~/(\^[A-Z]+|[ATCG])/g;
    $count+=()=$items[8]=~/(I)/g;
    $count+=()=$items[10]=~/(\^[A-Z]+|[ATCG])/g;

    $temp=$tags{$items[1]."_".$items[2]}->{$start}->{$prefix.($end-$start+1)};
 
    # The first hash key is chr and strand
    # The second hash key is start
    # The third hash key is file_length
    # The hash value stores mutant count, tag name, file_pos_f, CIGAR_f, seq_f, MD_f, pos_r, CIGAR_r, seq_r, MD_r 

    if (exists ${$temp}[0])
    {
      if (abs(${$temp}[0])<$count) {@$temp=($count,$items[0],$prefix.$items[3],@items[4..10]);} 
    }else
    {
      $tags{$items[1]."_".$items[2]}->{$start}->{$prefix.($end-$start+1)}=[$count,$items[0],$prefix.$items[3],@items[4..10]];
    }
  }

  close(FILE_IN);
}

##################  write clusters  #########################

print("  Write clusters\n");

open(FILE_OUT,">".$file_cluster) or die "write to temporary cluster file failed!\n";

foreach $temp (keys %tags) # form clusters
{
  ($chr,$strand)=split("_",$temp);

  @pos_start=keys %{$tags{$temp}};
  @pos_start=sort {$a <=> $b} @pos_start; # @pos_start contains sorted start sites
  map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[0]}}; # @cluster contains the information of the tags that form one single cluster
  $end=$pos_start[0]-1+max(map {abs($_)} keys %{$tags{$temp}->{$pos_start[0]}});

  for $i (1..$#pos_start)
  {
    if ($end+$step>=$pos_start[$i])
    {
      $end=max($end,$pos_start[$i]-1+max(map {abs($_)} keys %{$tags{$temp}->{$pos_start[$i]}}));
      map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[$i]}};
    }else
    {
      $m1=0;
      $m2=0;
      map {if (${$_}[2]>0) {$m1++;} else {$m2++;}} @cluster;

      if ($m1>=$min1 || $m2>=$min2) 
      {
        print FILE_OUT $chr."\t".$strand."\t".($#cluster+1)."\n";
        map {print FILE_OUT join("\t",@{$_}[1..9]);} @cluster;
      }
      @cluster=();
      map {push @cluster,$_} values %{$tags{$temp}->{$pos_start[$i]}};
      $end=$pos_start[$i]-1+max(map {abs($_)} keys %{$tags{$temp}->{$pos_start[$i]}});
    }
  }

  $m1=0;
  $m2=0;
  map {if (${$_}[2]>0) {$m1++;} else {$m2++;}} @cluster;

  if ($m1>=$min1 || $m2>=$min2) 
  {
    print FILE_OUT $chr."\t".$strand."\t".($#cluster+1)."\n";
    map {print FILE_OUT join("\t",@{$_}[1..9]);} @cluster;
  }
  @cluster=();
}

close(FILE_OUT);

exit;




