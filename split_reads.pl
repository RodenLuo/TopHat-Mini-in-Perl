#!/usr/bin/perl
die "usage: perl split_reads.pl <input fastq> <output splitted fastq>\n" unless @ARGV == 2;

$input = $ARGV[0];
$output = $ARGV[1];

open Unaligned,$input;
open Output,">$output";
#@unaligned = <Unaligned>;

#$line is a line in sam file
while($line = <Unaligned>){
  chomp $line;
  $read[0] = $line;
  $read[9] = <Unaligned>;
  $_ = <Unaligned>;
  $read[10] = <Unaligned>;
  #$read[0]	read name
  #$read[9]	sequence
  #$read[10]	quality
  #if(length($read[9]) < 75){next;}

  @sub_read = unpack("(A25)*",$read[9]);
  @sub_quality = unpack("(A25)*",$read[10]);

  $c = 0;

  while($c<@sub_read-1){
    $start = 1+25*$c;
    print Output "$read[0]_tophat_${start}_tophat_\n$sub_read[$c]\n+\n$sub_quality[$c]\n";
    $c += 1;
  }
  $start = 1+25*$c;
  if (length($sub_read[$c]) >= 15){print Output "\@$read[0]_tophat_${start}_tophat_\n$sub_read[$c]\n+\n$sub_quality[$c]\n";}
}
