#!/usr/bin/perl
use strict;
use warnings;
use PDL;
use PDL::DiskCache;

my ($dir,$step,$file,@files)=@ARGV;
my $cache;
my ($region_id,$chr,$strand,$pos,$tag1,$mut1,$tag2,$mut2);
my ($M,$emission,$region,$gamma,$xis,$a,$fb);
my ($pre_id,$i,$j,$pre_state,@state,$start,$end,@color,$strength);

##########  read data from disk  ####################

print "  Read data from disk\n";

@files=@files[0..6];
$cache=diskcache(\@files,{ro=>1}); # cache matrices
($M,$emission,$region,$gamma,$xis,$a,$fb)=@{$cache};

($region_id,$chr,$pos,$strand,$tag1,$mut1,$tag2,$mut2)=rcols $file,0,1,2,3,[4..$step+3],[$step+4..2*$step+3],[4+2*$step..3+3*$step],[3*$step+4..4*$step+3],{PERLCOLS=>[0,1,2,3],COLSEP=>"\t"};

#########  write bedGragh files  ######################

print "  Write bedGraph files\n";

write_bg($chr,$pos,$dir,$tag1,$step,$strand,"+","File1_Total_pos");
write_bg($chr,$pos,$dir,$mut1,$step,$strand,"+","File1_Mutant_pos");
write_bg($chr,$pos,$dir,$tag2,$step,$strand,"+","File2_Total_pos");
write_bg($chr,$pos,$dir,$mut2,$step,$strand,"+","File2_Mutant_pos");
write_bg($chr,$pos,$dir,$tag1,$step,$strand,"-","File1_Total_neg");
write_bg($chr,$pos,$dir,$mut1,$step,$strand,"-","File1_Mutant_neg");
write_bg($chr,$pos,$dir,$tag2,$step,$strand,"-","File2_Total_neg");
write_bg($chr,$pos,$dir,$mut2,$step,$strand,"-","File2_Mutant_neg");

sub write_bg
{
  my ($chr,$pos,$dir,$wig,$step,$strand,$pn,$file)=@_;
  my ($i,$j,$tmp_chr,$tmp_pos,$tmp_count);
  my ($pre_chr,$pre_start,$pre_end,$pre_count);

  $pre_chr="";  
  $pre_start=-1;
  $pre_end=-1;
  $pre_count=-1;
  $i=-1;
  $j=-1;

  open(FILE,">".$dir."/".$file.".bedgraph");
  
  foreach $tmp_count ($wig->xchg(0,1)->list) 
  {
    $j++;
    if ($j % $step==0)
    {
      $i++;
      $tmp_chr=$chr->[$i];
      $tmp_pos=$pos->[$i];
      $j=0;
    }else
    {
      $tmp_pos++;
    } # get new coordinates
    
    if ($tmp_count==0 || $strand->[$i] ne $pn) {next;} # skip zeroes and the wrong strands

    if ($pre_chr ne $tmp_chr || $pre_end!=$tmp_pos-1 || $pre_count!=$tmp_count) # concatenate neighboring positions with the same count
    {
      if ($pre_chr eq "")
      {
        print FILE "track type=bedGraph name=\"".$file."\" description=\"".$file."\"\n"; # write header
      }else
      {
        print FILE $pre_chr." ".$pre_start." ".($pre_end+1)." ".$pre_count."\n";
      }
       
      $pre_chr=$tmp_chr; 
      $pre_start=$tmp_pos; 
      $pre_end=$tmp_pos; 
      $pre_count=$tmp_count; 
    }else
    {   
      $pre_end=$tmp_pos;
    }
  } 

  print FILE $pre_chr." ".$pre_start." ".($pre_end+1)." ".$pre_count."\n";

  close(FILE);
}

#########  write full HMM output  ########################

print "  Write full HMM output file\n";

$tag1=$tag1->xchg(0,1)->sumover;
$tag2=$tag2->xchg(0,1)->sumover;
$mut1=$mut1->xchg(0,1)->sumover;
$mut2=$mut2->xchg(0,1)->sumover;

# the other parts of xis table will also be rewritten to contain output information for the sake of convenience

$xis->slice("(2),:").=$M;
$xis->slice("(3),:").=$tag1;
$xis->slice("(4),:").=$mut1;
$xis->slice("(5),:").=$tag2;
$xis->slice("(6),:").=$mut2;

open(FILE,">".$dir."/dCLIP_output.txt");
print FILE "id\tchrom\tstrand\tposition\tstate\tprobability\tdifferential\ttag1\tmut1\ttag2\tmut2\n";

foreach $i (0..$M->nelem-1)
{
  print FILE $region_id->[$i]."\t".$chr->[$i]."\t".$strand->[$i]."\t".$pos->[$i]."\t";
  print FILE join("\t",$xis->slice("0:6,($i)")->list)."\n";
}

close(FILE);

############  write summary HMM output  ##################

print "  Write summary HMM output file\n";

@state=$xis->slice("(0),:")->list;
$pre_state=$state[0];
$pre_id=$region_id->[0];
$i=0;
$j=1;
$strength=0;
@color=("255,0,0","0,255,0","0,0,255");

open(FILE,">".$dir."/dCLIP_summary.bed");
print FILE "track name=\"dCLIP_summary\" description=\"dCLIP_summary\" itemRgb=\"On\"\n";

while ($j<=$#state)
{
  if ($region_id->[$j]==$pre_id && $state[$j]==$pre_state)
  {
    $strength+=at($xis,5-$pre_state,$j);
    $j++;
  }else
  {
    $start=$pos->[$i];
    $end=$pos->[$j-1]+$step-1;
    if ($pre_state==1) {$strength=0;}
    print FILE join("\t",$chr->[$i],$start,$end,$pre_id,int($strength/($j-$i)/$step*10)/10,$strand->[$i],$start,$end,$color[$pre_state])."\n";
    $i=$j;
    $strength=at($xis,5-$state[$j],$j);
    $j++;
    $pre_id=$region_id->[$i];
    $pre_state=$state[$i];
  }
}

$start=$pos->[$i]; 
$end=$pos->[$j-1]+$step-1;
print FILE join("\t",$chr->[$i],$start,$end,$pre_id,int($strength/($j-$i)/$step*10)/10,$strand->[$i],$start,$end,$color[$pre_state])."\n";

close(FILE);








