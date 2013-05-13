#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/max/;
$SIG{__WARN__}=sub{die "SAM file format error!\n"};

my ($file1,$file2,$min1,$min2,$file_cluster,$step)=@ARGV;
my ($temp,$count,@items,%tags,$len,$chr,@pos_start,$end,@cluster,$i,$strand,$m1,$m2);
my $md_field=-1;

##################  search for MD field  #########################

print("  Locate MD field\n");

open(FILE_IN,$file1) or die "read alignment file ".$file1." failed!\n";

while (<FILE_IN>) 
{
  if ($_=~/(.*)MD\:Z\:/)
  {
    $i=()=$1=~/(\t)/g;
    if ($i>10) 
    {
      $md_field=$i;
      last;
    }
  }
}

if ($md_field<0) {die "no md field found in sam file!\n";}
close(FILE_IN);

##################  read tags  ######################################

print("  Read raw alignment file 1\n");
read_tags(\%tags,$file1,"");
print("  Read raw alignment file 2\n");
read_tags(\%tags,$file2,"-");

sub read_tags # read and compute length of each tag in each file
{
  my ($tags_ref,$file,$prefix)=@_;
  $i=0;

  open(FILE_IN,$file) or die "read alignment file ".$file." failed!\n";

  while (<FILE_IN>)
  {
    if ($_=~/^@/) {next;}
    @items=split("\t",$_);
    if ($items[5]=~/N/ || !($items[1]==0 || $items[1]==16)) {next;} # delete gapped mapping and unmapped reads
    
    $i++;
    if ($i % 10000000==0) {print "    Read in ".$i." reads\n";}

    $len=0;
    while ($items[5]=~/([0-9]+)([MD])/g) {$len+=$1;} # the mapped region length on the reference genome   

    $count=()=$items[5]=~/(I)/g;
    $count+=()=$items[$md_field]=~/(\^[A-Z]+|[ATCG])/g; # count mutant numbers on each tag
   
    $temp=$tags_ref->{$items[1]."_".$items[2]}->{$items[3]}->{$prefix.$len}; 

    if ($items[$md_field]=~/MD\:Z\:(.*)$/) # fix MD field if it does not exist
    {
      $items[$md_field]=$1;
    }else
    {
      $items[$md_field]=$len;
    }

    if (exists ${$temp}[0])
    {
      if (${$temp}[0]<$count) {@$temp=($count,$items[0],$prefix.$items[3],$items[5],$items[9],$items[$md_field]);} 
      # strand and chr are stored as the first hash key
      # pos is stored as the second hash key
      # file and length are stored as the third hash key
      # mutant count, tag name, file_pos, CIGAR, sequence, MD are stored as the hash value
    }else
    {
      $tags_ref->{$items[1]."_".$items[2]}->{$items[3]}->{$prefix.$len}=[$count,$items[0],$prefix.$items[3],$items[5],$items[9],$items[$md_field]];
    }
  }

  close(FILE_IN);
}

###################  write to temp file  #############################

print("  Write cluster file\n");

open(FILE_OUT,">".$file_cluster) or die "write to temporary cluster file failed!\n";

foreach $temp (keys %tags) # form clusters
{
  ($strand,$chr)=split("_",$temp);
  if ($strand==0) 
  {
    $strand="+";
  }else
  {
    $strand="-";
  }

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
      $m1=0; # count number of tags from each file in one cluster
      $m2=0; 
      map {if (${$_}[2]>0) {$m1++;} else {$m2++;}} @cluster;

      if ($m1>=$min1 || $m2>=$min2) 
      {
        print FILE_OUT $chr."\t".$strand."\t".($#cluster+1)."\n";
        map {print FILE_OUT ${$_}[1]."\t".${$_}[2]."\t".${$_}[3]."\t".${$_}[4]."\t".${$_}[5]."\n";} @cluster;
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
    map {print FILE_OUT ${$_}[1]."\t".${$_}[2]."\t".${$_}[3]."\t".${$_}[4]."\t".${$_}[5]."\n";} @cluster;
  }
  @cluster=();
}

close(FILE_OUT);
