#!/usr/bin/env perl
use List::Util qw(max min sum);
die "usage: perl splice.pl <splitted_srt_name.sam> <output.sam> <ref.tdx> <ref.fa> [max gap]\n" unless (@ARGV == 4 or @ARGV == 5);

$input = $ARGV[0];
$output = $ARGV[1];
$ref_tdx = $ARGV[2];
$ref_fa = $ARGV[3];
# $ARGV[5] is used

open Reads,$input;
open Out,">$output";
open Out_Splice, ">$output.splice_junction.txt";

#hg19.fa part
open Tdx,$ref_tdx;
open Fa,$ref_fa;
@fa = <Fa>;

#store all chr start line in %tdx
%tdx = ();
while($line = <Tdx>){
  chomp $line;
  ($tdx_coo,$tdx_chr) = split(/:>/,$line);
  $tdx{$tdx_chr} = $tdx_coo;
}


$max_gap = 20000;
if($ARGV[5]){$max_gap = $ARGV[5];}
$min_gap = 70;

$line2 = <Reads>;
while($line1 = $line2){
  %reads_tophat=();#存下四个切割了的reads，键是切割点
  ($read_no,$_,$tophat_no,$_) = split(/_/,$line1);
  $reads_tophat{$tophat_no} = $line1;

  %reads_hash=();#存下每个sam记录，键是align到的染色体和坐标
  @r = split(/\t/,$line1);
  $reads_hash{"$r[2]\t$r[3]"} = $line1;

  #继续往下，只要是同一个100bp read里面来的，就继续放到哈希里面去
  while($line2 = <Reads> and $line2 =~ /^$read_no/){
    ($_,$_,$tophat_no,$_) = split(/_/,$line2);
    $reads_tophat{$tophat_no} = $line2;

    @r = split(/\t/,$line2);
    $reads_hash{"$r[2]\t$r[3]"} = $line2;
  }

  #把seg,符合12 3和 1 23的情况找出来。
  #@sort_keys = sort { $a <=> $b } (keys %reads_hash);
  @sort_keys = sort (keys %reads_hash);
  for (my $i =0; $i < @sort_keys-2; $i++){
    ($chr1,$coor1)=split(/\t/,$sort_keys[$i]);
    ($chr2,$coor2)=split(/\t/,$sort_keys[$i+1]);
    ($chr3,$coor3)=split(/\t/,$sort_keys[$i+2]);
    if ($chr1 eq $chr2 and $chr2 eq $chr3){
      #12 3的情况
      if ($coor1+25 == $coor2 and abs($coor3-$coor2) < $max_gap and abs($coor3-$coor2) > $min_gap){
        #存储符合条件的三个reads
        $seg1 = $reads_hash{$sort_keys[$i]};
        $seg2 = $reads_hash{$sort_keys[$i+1]};
        $seg3 = $reads_hash{$sort_keys[$i+2]};
	#存储三个reads的tophat_no，供后续判断正负链，并找到缺失的那一条的信息
        ($read_name,$_,$tophat_no_1,$_) = split(/_/,$seg1);        
        ($_,$_,$tophat_no_2,$_) = split(/_/,$seg2);
        ($_,$_,$tophat_no_3,$_) = split(/_/,$seg3);
        $tophat_no_sum = $tophat_no_1+$tophat_no_2+$tophat_no_3;

        #1+26+_+76 正链
        if ($tophat_no_sum == 103){
#=cut
	#输出测试
        #print Out "read_name:$read_name\t$tophat_no_1:$sort_keys[$i]\t$tophat_no_2:$sort_keys[$i+1]\t$tophat_no_3:$sort_keys[$i+2]\n";

          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);

	  #寻找左边的ref_1
    	  $input_chr = $tophat_read_2_chr;
          $input_coo = $tophat_read_2_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg2的行号，取它的下面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -1;
            $start = 74;
          }else{
            $start += 24;
            $idx = $tdx{$input_chr} + $lines;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_1 = substr($twolines,$start,25);

	  #寻找右边的ref_2
    	  $input_chr = $tophat_read_3_chr;
          $input_coo = $tophat_read_3_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg3的行号，取它的前面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -2;
            $start = 74;
          }else{
            $start = 25+($start-1);
            $idx = $tdx{$input_chr} + $lines -1;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_2 = substr($twolines,$start,25);

	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"51"});chomp @junction_read;
          $junction_seq = $junction_read[9];

	  #把ref_1和ref_2拼接成25bp，判断splice form，然后再和junction_reads对比，小于等于两个mismatch就输出
	  $junction_point = -1;
	  for (my $i=0; $i<25; $i++){
		#考虑i>1的时候，考虑经典splice form
	    if ($i>1){
		#如果不符合要求，跳入下一个base
	      if ( not(substr($ref_1,$i,2) =~ /GT/i or substr($ref_1,$i,2) =~ /CT/i or substr($ref_1,$i,2) =~ /GC/i or substr($ref_1,$i,2) =~ /AT/i ) ){
	        next;
	      }
	      if ( not(substr($ref_2,$i-2,2) =~ /AG/i or substr($ref_2,$i-2,2) =~ /AC/i or substr($ref_2,$i-2,2) =~ /GC/i) or substr($ref_2,$i-2,2) =~ /AT/i){
	        next;
	      }
            }
	    $ref_seg_1 = substr($ref_1,0,$i);
	    $ref_seg_2 = substr($ref_2,$i,25-$i);
	    $ref_seg = $ref_seg_1.$ref_seg_2;

            $same_no = align_get_same_no($junction_seq,$ref_seg);
	    if ($same_no > 22){ 
	      $junction_point = $i;
	      #写输出，输出成一个read
	      $whole_seq = $seq_1.$seq_2.$junction_seq.$seq_3;
	      $whole_qua = $qua_1.$qua_2.$junction_read[10].$qua_3;
	      ($reads_real_name,$_) = split(/_/, $reads_tophat_name);
	      $cigar1 = 50+$junction_point;
	      $cigar2 = 50-$junction_point;
	      $N_no = $tophat_read_3_coo - $tophat_read_2_coo -50;
	      $out_read = "$reads_real_name\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	      print Out $out_read;
	      $splice_point_1 = $seg_coo_1 + $cigar1 - 1;#splice_L的最后一个碱基坐标
	      $splice_point_2 = $tophat_read_3_coo - 25 + $junction_point; #splice_R的第一个碱基坐标
	      print Out_Splice "${seg_chr_1}_$splice_point_1\t${seg_chr_1}_$splice_point_2\n";
	      #print Out "ref1:$ref_1\tref2:$ref_2\tref:$ref_seg\tjunc_read:$junction_seq\tpoint:$junction_point\n";
	    }
	  }
#=cut
        }
        #76+51+_+1 反链
        elsif ($tophat_no_sum == 128){
#=cut
	  	#输出测试
        #print Out "read_name:$read_name\t$tophat_no_1:$sort_keys[$i]\t$tophat_no_2:$sort_keys[$i+1]\t$tophat_no_3:$sort_keys[$i+2]\n";

          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);

	  #寻找左边的ref_1
    	  $input_chr = $tophat_read_2_chr;
          $input_coo = $tophat_read_2_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg2的行号，取它的下面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -1;
            $start = 74;
          }else{
            $start += 24;
            $idx = $tdx{$input_chr} + $lines;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_1 = substr($twolines,$start,25);

	  #寻找右边的ref_2
    	  $input_chr = $tophat_read_3_chr;
          $input_coo = $tophat_read_3_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg3的行号，取它的前面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -2;
            $start = 74;
          }else{
            $start = 25+($start-1);
            $idx = $tdx{$input_chr} + $lines -1;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_2 = substr($twolines,$start,25);

	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"26"});chomp @junction_read;
          $junction_seq = reverse_complement($junction_read[9]);

	  #把ref_1和ref_2拼接成25bp，判断splice form，然后再和junction_reads对比，小于等于两个mismatch就输出
	  $junction_point = -1;
	  for (my $i=0; $i<25; $i++){
		#考虑i>1的时候，考虑经典splice form
	    if ($i>1){
		#如果不符合要求，跳入下一个base
	      if ( not(substr($ref_1,$i,2) =~ /GT/i or substr($ref_1,$i,2) =~ /CT/i or substr($ref_1,$i,2) =~ /GC/i or substr($ref_1,$i,2) =~ /AT/i ) ){
	        next;
	      }
	      if ( not(substr($ref_2,$i-2,2) =~ /AG/i or substr($ref_2,$i-2,2) =~ /AC/i or substr($ref_2,$i-2,2) =~ /GC/i) or substr($ref_2,$i-2,2) =~ /AT/i){
	        next;
	      }
            }
	    $ref_seg_1 = substr($ref_1,0,$i);
	    $ref_seg_2 = substr($ref_2,$i,25-$i);
	    $ref_seg = $ref_seg_1.$ref_seg_2;

            $same_no = align_get_same_no($junction_seq,$ref_seg);
	    if ($same_no > 22){
	      $junction_point = $i;
	      #写输出，输出成一个read
	      $whole_seq = $seq_1.$seq_2.$junction_seq.$seq_3;
	      $whole_qua = $qua_1.$qua_2.reverse($junction_read[10]).$qua_3;
	      ($reads_real_name,$_) = split(/_/, $reads_tophat_name);
	      $cigar1 = 50+$junction_point;
	      $cigar2 = 50-$junction_point;
	      $N_no = $tophat_read_3_coo - $tophat_read_2_coo -50;
	      $out_read = "$reads_real_name\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	      print Out $out_read;
	      $splice_point_1 = $seg_coo_1 + $cigar1 - 1;#splice_L的最后一个碱基坐标
	      $splice_point_2 = $tophat_read_3_coo - 25 + $junction_point; #splice_R的第一个碱基坐标
	      print Out_Splice "${seg_chr_1}_$splice_point_1\t${seg_chr_1}_$splice_point_2\n";
	      #print Out "ref1:$ref_1\tref2:$ref_2\tref:$ref_seg\tjunc_read:$junction_seq\tpoint:$junction_point\n";
	    }
	  }#跑一边ref结尾
#=cut
	}#12 3反链结尾
      }#12 3结尾

      #1 23的情况
      elsif ($coor2+25 == $coor3 and abs($coor2-$coor1) < $max_gap and abs($coor2-$coor1) > $min_gap){
        #存储符合条件的三个reads
        $seg1 = $reads_hash{$sort_keys[$i]};
        $seg2 = $reads_hash{$sort_keys[$i+1]};
        $seg3 = $reads_hash{$sort_keys[$i+2]};
	#存储三个reads的tophat_no，供后续判断正负链，并找到缺失的那一条的信息
        ($read_name,$_,$tophat_no_1,$_) = split(/_/,$seg1);        
        ($_,$_,$tophat_no_2,$_) = split(/_/,$seg2);
        ($_,$_,$tophat_no_3,$_) = split(/_/,$seg3);
        $tophat_no_sum = $tophat_no_1+$tophat_no_2+$tophat_no_3;

        #1+_+51+76 正链
        if ($tophat_no_sum == 128){
#=cut
	#输出测试
        #print Out "read_name:$read_name\t$tophat_no_1:$sort_keys[$i]\t$tophat_no_2:$sort_keys[$i+1]\t$tophat_no_3:$sort_keys[$i+2]\n";

          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_1_chr,$tophat_read_1_coo,$_,$_,$_,$_,$_,$seq_1,$qua_1,$_) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);

	  #寻找左边的ref_1
    	  $input_chr = $tophat_read_1_chr;
          $input_coo = $tophat_read_1_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg2的行号，取它的下面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -1;
            $start = 74;
          }else{
            $start += 24;
            $idx = $tdx{$input_chr} + $lines;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_1 = substr($twolines,$start,25);

	  #寻找右边的ref_2
    	  $input_chr = $tophat_read_2_chr;
          $input_coo = $tophat_read_2_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg3的行号，取它的前面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -2;
            $start = 74;
          }else{
            $start = 25+($start-1);
            $idx = $tdx{$input_chr} + $lines -1;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_2 = substr($twolines,$start,25);

	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"26"});chomp @junction_read;
          $junction_seq = $junction_read[9];

	  #把ref_1和ref_2拼接成25bp，判断splice form，然后再和junction_reads对比，小于等于两个mismatch就输出
	  $junction_point = -1;
	  for (my $i=0; $i<25; $i++){
		#考虑i>1的时候，考虑经典splice form
	    if ($i>1){
		#如果不符合要求，跳入下一个base
	      if ( not(substr($ref_1,$i,2) =~ /GT/i or substr($ref_1,$i,2) =~ /CT/i or substr($ref_1,$i,2) =~ /GC/i or substr($ref_1,$i,2) =~ /AT/i ) ){
	        next;
	      }
	      if ( not(substr($ref_2,$i-2,2) =~ /AG/i or substr($ref_2,$i-2,2) =~ /AC/i or substr($ref_2,$i-2,2) =~ /GC/i) or substr($ref_2,$i-2,2) =~ /AT/i){
	        next;
	      }
            }
	    $ref_seg_1 = substr($ref_1,0,$i);
	    $ref_seg_2 = substr($ref_2,$i,25-$i);
	    $ref_seg = $ref_seg_1.$ref_seg_2;

            $same_no = align_get_same_no($junction_seq,$ref_seg);
	    if ($same_no > 22){ 
	      $junction_point = $i;
	      #写输出，输出成一个read
	      $whole_seq = $seq_1.$junction_seq.$seq_2.$seq_3;
	      $whole_qua = $qua_1.$junction_read[10].$qua_2.$qua_3;
	      ($reads_real_name,$_) = split(/_/, $reads_tophat_name);
	      $cigar1 = 25+$junction_point;
	      $cigar2 = 75-$junction_point;
	      $N_no = $tophat_read_2_coo - $tophat_read_1_coo -50;
	      $out_read = "$reads_real_name\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	      print Out $out_read;
	      $splice_point_1 = $seg_coo_1 + $cigar1 - 1;#splice_L的最后一个碱基坐标
	      $splice_point_2 = $tophat_read_2_coo - 25 + $junction_point; #splice_R的第一个碱基坐标
	      print Out_Splice "${seg_chr_1}_$splice_point_1\t${seg_chr_1}_$splice_point_2\n";
	      #print Out "ref1:$ref_1\tref2:$ref_2\tref:$ref_seg\tjunc_read:$junction_seq\tpoint:$junction_point\n";
	    }
	  }
#=cut
        }

        #76+_+26+1 反链
        elsif ($tophat_no_sum == 103){
	  	#输出测试
        #print Out "read_name:$read_name\t$tophat_no_1:$sort_keys[$i]\t$tophat_no_2:$sort_keys[$i+1]\t$tophat_no_3:$sort_keys[$i+2]\n";

          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
	  ($_,$_,$tophat_read_1_chr,$tophat_read_1_coo,$_,$_,$_,$_,$_,$seq_1,$qua_1,$_) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);

	  #寻找左边的ref_1
    	  $input_chr = $tophat_read_1_chr;
          $input_coo = $tophat_read_1_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg2的行号，取它的下面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -1;
            $start = 74;
          }else{
            $start += 24;
            $idx = $tdx{$input_chr} + $lines;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_1 = substr($twolines,$start,25);

	  #寻找右边的ref_2
    	  $input_chr = $tophat_read_2_chr;
          $input_coo = $tophat_read_2_coo;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg3的行号，取它的前面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -2;
            $start = 74;
          }else{
            $start = 25+($start-1);
            $idx = $tdx{$input_chr} + $lines -1;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg2后面的25个base
          $ref_2 = substr($twolines,$start,25);

	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"51"});chomp @junction_read;
          $junction_seq = reverse_complement($junction_read[9]);

	  #把ref_1和ref_2拼接成25bp，判断splice form，然后再和junction_reads对比，小于等于两个mismatch就输出
	  $junction_point = -1;
	  for (my $i=0; $i<25; $i++){
		#考虑i>1的时候，考虑经典splice form
	    if ($i>1){
		#如果不符合要求，跳入下一个base
	      if ( not(substr($ref_1,$i,2) =~ /GT/i or substr($ref_1,$i,2) =~ /CT/i or substr($ref_1,$i,2) =~ /GC/i or substr($ref_1,$i,2) =~ /AT/i ) ){
	        next;
	      }
	      if ( not(substr($ref_2,$i-2,2) =~ /AG/i or substr($ref_2,$i-2,2) =~ /AC/i or substr($ref_2,$i-2,2) =~ /GC/i) or substr($ref_2,$i-2,2) =~ /AT/i){
	        next;
	      }
            }
	    $ref_seg_1 = substr($ref_1,0,$i);
	    $ref_seg_2 = substr($ref_2,$i,25-$i);
	    $ref_seg = $ref_seg_1.$ref_seg_2;

            $same_no = align_get_same_no($junction_seq,$ref_seg);
	    if ($same_no > 22){
	      $junction_point = $i;
	      #写输出，输出成一个read
	      $whole_seq = $seq_1.$junction_seq.$seq_2.$seq_3;
	      $whole_qua = $qua_1.reverse($junction_read[10]).$qua_2.$qua_3;
	      ($reads_real_name,$_) = split(/_/, $reads_tophat_name);
	      $cigar1 = 25+$junction_point;
	      $cigar2 = 75-$junction_point;
	      $N_no = $tophat_read_2_coo - $tophat_read_1_coo -50;
	      $out_read = "$reads_real_name\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	      print Out $out_read;
	      $splice_point_1 = $seg_coo_1 + $cigar1 - 1;#splice_L的最后一个碱基坐标
	      $splice_point_2 = $tophat_read_2_coo - 25 + $junction_point; #splice_R的第一个碱基坐标
	      print Out_Splice "${seg_chr_1}_$splice_point_1\t${seg_chr_1}_$splice_point_2\n";
	      #print Out "ref1:$ref_1\tref2:$ref_2\tref:$ref_seg\tjunc_read:$junction_seq\tpoint:$junction_point\n";
	    }
	  }#跑一边ref结尾
	}#1_23反链结尾
      }
    }
  }
}

sub align_get_same_no
{
  @junction_seq_array = split(//,$_[0]);
  @reference_seq_array = split(//,$_[1]);
  $same_no = 0;
  for(my $i = 0; $i < 25; $i++){
    if($junction_seq_array[$i] =~ /$reference_seq_array[$i]/i){$same_no++;}
  }
  return $same_no;
}

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
