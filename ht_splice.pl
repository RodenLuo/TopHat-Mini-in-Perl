#!/usr/bin/env perl
use List::Util qw(max min sum);
die "usage: perl splice.pl <splitted_srt_name.sam> <output.sam> <junction.txt> <ref.tdx> <ref.fa> [max gap]\n" unless (@ARGV == 5 or @ARGV == 6);

$input = $ARGV[0];
$output = $ARGV[1];
$junction_file = $ARGV[2];
$ref_tdx = $ARGV[3];
$ref_fa = $ARGV[4];
# $ARGV[5] is used

open Reads,$input;
open Junction_file,$junction_file;
open Out,">$output";

#store all junction found by 12_3 and 1_23 before
%splice_L_hash;
%splice_R_hash;
while($line = <Junction_file>){
  chomp $line;
  ($L_point, $R_point) = split(/\t/, $line);
  $splice_L_hash{$L_point} = $R_point;
  $splice_R_hash{$R_point} = $L_point;
}

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

  #继续往下，只要是一个100bp read里面来的，就继续放到哈希里面去
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
    #取出i后四个segment的坐标信息
    ($chr1,$coor1)=split(/\t/,$sort_keys[$i]);
    ($chr2,$coor2)=split(/\t/,$sort_keys[$i+1]);
    ($chr3,$coor3)=split(/\t/,$sort_keys[$i+2]);
    ($chr4,$coor4)=split(/\t/,$sort_keys[$i+3]);
    if ($chr1 eq $chr2 and $chr2 eq $chr3){
      #12 3的情况
      if ($coor1+25 == $coor2 and abs($coor3-$coor2) < $max_gap and abs($coor3-$coor2) > $min_gap){
      }#12 3结尾

      #1 23的情况
      elsif ($coor2+25 == $coor3 and abs($coor2-$coor1) < $max_gap and abs($coor2-$coor1) > $min_gap){
      }

      #123的情况
      elsif ($coor1+25 == $coor2 and $coor2+25 == $coor3){
	#存储符合条件的三个reads
        $seg1 = $reads_hash{$sort_keys[$i]};
        $seg2 = $reads_hash{$sort_keys[$i+1]};
        $seg3 = $reads_hash{$sort_keys[$i+2]};
	#存储三个reads的tophat_no，供后续判断正负链，并找到缺失的那一条的信息
        ($read_name,$_,$tophat_no_1,$_) = split(/_/,$seg1);
        ($_,$_,$tophat_no_2,$_) = split(/_/,$seg2);
        ($_,$_,$tophat_no_3,$_) = split(/_/,$seg3);
        $tophat_no_sum = $tophat_no_1+$tophat_no_2+$tophat_no_3;
	
	#1+26+51 case
	if ($tophat_no_1 == 1 and $tophat_no_sum == 78){	 
          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);	 
	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"76"});chomp @junction_read;
          $junction_seq = $junction_read[9];

	  #寻找左边的ref_1
    	  $input_chr = $chr3;
          $input_coo = $coor3;
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
 
	  for (my $i = 0; $i <= 25; $i++){
	    $poteintial_junction = $coor3 + 24 +$i;
	    if( exists($splice_L_hash{"${chr3}_$poteintial_junction"}) ){
	      ($_, $ref_2_start) = split(/_/, $splice_L_hash{"${chr3}_$poteintial_junction"});
	      #print "coor3:$coor3\tref:$ref_2_start\n";
	      #寻找右边的ref_2
    	      $input_chr = $chr3;
              $input_coo = $ref_2_start;
              $lines = int($input_coo/50);
              $start = $input_coo % 50;
              #找到seg2的行号，取它的下面两行做ref
              if ($start == 0){
                $idx = $tdx{$input_chr} + $lines -1;
                $start = 49;
              }else{
                $start -= 1;
                $idx = $tdx{$input_chr} + $lines;
              }
              chomp $fa[$idx]; chomp $fa[$idx+1];
              $twolines = $fa[$idx].$fa[$idx+1];
	      #取seg2后面的25个base
              $ref_2 = substr($twolines,$start,25);
	      
	      #对比ref_seq和junction_seq,成功就输出
	      $ref_seg = substr($ref_1,0,$i).substr($ref_2,0,25-$i);
	      $same_no = align_get_same_no($junction_seq,$ref_seg);
	      if ($same_no > 22){
	        #写输出，输出成一个read
	        $whole_seq = $seq_1.$seq_2.$seq_3.$junction_seq;
	        $whole_qua = $qua_1.$qua_2.$qua_3.$junction_read[10];
	        $cigar1 = 75+$i;
	        $cigar2 = 25-$i;
	        $N_no = $ref_2_start - $coor3 - 25 -$i;
	        $out_read = "$read_no\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	        print Out $out_read;
	      } 
	    }
	  }
	}

	#26+51+76 case
	elsif ($tophat_no_1 == 26 and $tophat_no_sum == 153){
          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);	 
	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"1"});chomp @junction_read;
          $junction_seq = $junction_read[9];

	  #寻找26前面的ref_2，取它前面的25bp
    	  $input_chr = $chr1;
          $input_coo = $coor1;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg1的行号，取它的前面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -2;
            $start = 74;
          }else{
            $start = 25+($start-1);
            $idx = $tdx{$input_chr} + $lines -1;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg1前面的25个base
          $ref_2 = substr($twolines,$start,25);
 
	  for (my $i = 0; $i <= 25; $i++){
	    $poteintial_junction = $coor1 - 25 +$i;
	    if( exists($splice_R_hash{"${chr3}_$poteintial_junction"}) ){
	      ($_, $ref_1_start) = split(/_/, $splice_R_hash{"${chr1}_$poteintial_junction"});
	      #print "coor1:$coor1\tref:$ref_1_start\n";
	      #寻找左边的ref_1
    	      $input_chr = $chr1;
              $input_coo = $ref_1_start;
              $lines = int($input_coo/50);
              $start = $input_coo % 50;
	      #找到ref_1的行号，取它的前面两行做ref
              if ($start == 0){
                $idx = $tdx{$input_chr} + $lines -2;
                $start = 75;
              }else{
                $start = 26+($start-1);
                $idx = $tdx{$input_chr} + $lines -1;
              }
              chomp $fa[$idx]; chomp $fa[$idx+1];
              $twolines = $fa[$idx].$fa[$idx+1];
	      #取ref_1前面的25个base，含ref_1指向的那一个
              $ref_1 = substr($twolines,$start,25);
	      #print "ref_1:$ref_1\n";
	      #print "ref_2:$ref_2\n";
	      #print "i: $i\n";
	      #对比ref_seq和junction_seq,成功就输出
	      $ref_seg = substr($ref_1,25-$i,$i).substr($ref_2,$i,25-$i); #print "ref:$ref_seg\njuc:$junction_seq\n";
	      $same_no = align_get_same_no($junction_seq,$ref_seg);
	      if ($same_no > 22){
	        #写输出，输出成一个read
	        $whole_seq = $junction_seq.$seq_1.$seq_2.$seq_3;
	        $whole_qua = $junction_read[10].$qua_1.$qua_2.$qua_3;
	        $cigar1 = $i;
	        $cigar2 = 100-$i;
		$read_start = $ref_1_start - $i + 1;
	        $N_no = $coor1 - $read_start - 25;
	        $out_read = "$read_no\t$flag\t$seg_chr_1\t$read_start\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	        print Out $out_read;
	      } 
	    }
	  }
	}

	#76+51+26 case
	elsif ($tophat_no_1 == 76 and $tophat_no_sum == 153){
          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);	 
	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"1"});chomp @junction_read;
	  $junction_seq = reverse_complement($junction_read[9]);

	  #寻找左边的ref_1
    	  $input_chr = $chr3;
          $input_coo = $coor3;
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
 
	  for (my $i = 0; $i <= 25; $i++){
	    $poteintial_junction = $coor3 + 24 +$i;
	    if( exists($splice_L_hash{"${chr3}_$poteintial_junction"}) ){
	      ($_, $ref_2_start) = split(/_/, $splice_L_hash{"${chr3}_$poteintial_junction"});
	      #print "coor3:$coor3\tref:$ref_2_start\n";
	      #寻找右边的ref_2
    	      $input_chr = $chr3;
              $input_coo = $ref_2_start;
              $lines = int($input_coo/50);
              $start = $input_coo % 50;
              #找到seg2的行号，取它的下面两行做ref
              if ($start == 0){
                $idx = $tdx{$input_chr} + $lines -1;
                $start = 49;
              }else{
                $start -= 1;
                $idx = $tdx{$input_chr} + $lines;
              }
              chomp $fa[$idx]; chomp $fa[$idx+1];
              $twolines = $fa[$idx].$fa[$idx+1];
	      #取seg2后面的25个base
              $ref_2 = substr($twolines,$start,25);
	      
	      #对比ref_seq和junction_seq,成功就输出
	      $ref_seg = substr($ref_1,0,$i).substr($ref_2,0,25-$i);
	      $same_no = align_get_same_no($junction_seq,$ref_seg);
	      if ($same_no > 22){
	        #写输出，输出成一个read
	        $whole_seq = $seq_1.$seq_2.$seq_3.$junction_seq;
	        $whole_qua = $qua_1.$qua_2.$qua_3.reverse($junction_read[10]);
	        $cigar1 = 75+$i;
	        $cigar2 = 25-$i;
	        $N_no = $ref_2_start - $coor3 - 25 -$i;
	        $out_read = "$read_no\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	        print Out $out_read;
	      } 
	    }
	  }
	}

	#51+26+1 case
	elsif ($tophat_no_1 == 51 and $tophat_no_sum == 78){
          ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
          ($_,$_,$tophat_read_2_chr,$tophat_read_2_coo,$_,$_,$_,$_,$_,$seq_2,$qua_2,$_) = split(/\t/,$seg2);
          ($_,$_,$tophat_read_3_chr,$tophat_read_3_coo,$_,$_,$_,$_,$_,$seq_3,$qua_3,$_) = split(/\t/,$seg3);	 
	  #找到junction_read的信息
          @junction_read = ();
          @junction_read = split(/\t/,$reads_tophat{"76"});chomp @junction_read;
	  $junction_seq = reverse_complement($junction_read[9]);

	  #寻找26前面的ref_2，取它前面的25bp
    	  $input_chr = $chr1;
          $input_coo = $coor1;
          $lines = int($input_coo/50);
          $start = $input_coo % 50;
	  #找到seg1的行号，取它的前面两行做ref
          if ($start == 0){
            $idx = $tdx{$input_chr} + $lines -2;
            $start = 74;
          }else{
            $start = 25+($start-1);
            $idx = $tdx{$input_chr} + $lines -1;
          }
          chomp $fa[$idx]; chomp $fa[$idx+1];
          $twolines = $fa[$idx].$fa[$idx+1];
	  #取seg1前面的25个base
          $ref_2 = substr($twolines,$start,25);
 
	  for (my $i = 0; $i <= 25; $i++){
	    $poteintial_junction = $coor1 - 25 +$i;
	    if( exists($splice_R_hash{"${chr3}_$poteintial_junction"}) ){
	      ($_, $ref_1_start) = split(/_/, $splice_R_hash{"${chr1}_$poteintial_junction"});
	      #print "coor1:$coor1\tref:$ref_1_start\n";
	      #寻找左边的ref_1
    	      $input_chr = $chr1;
              $input_coo = $ref_1_start;
              $lines = int($input_coo/50);
              $start = $input_coo % 50;
	      #找到ref_1的行号，取它的前面两行做ref
              if ($start == 0){
                $idx = $tdx{$input_chr} + $lines -2;
                $start = 75;
              }else{
                $start = 26+($start-1);
                $idx = $tdx{$input_chr} + $lines -1;
              }
              chomp $fa[$idx]; chomp $fa[$idx+1];
              $twolines = $fa[$idx].$fa[$idx+1];
	      #取ref_1前面的25个base，含ref_1指向的那一个
              $ref_1 = substr($twolines,$start,25);
	      #print "ref_1:$ref_1\n";
	      #print "ref_2:$ref_2\n";
	      #print "i: $i\n";
	      #对比ref_seq和junction_seq,成功就输出
	      $ref_seg = substr($ref_1,25-$i,$i).substr($ref_2,$i,25-$i); #print "ref:$ref_seg\njuc:$junction_seq\n";
	      $same_no = align_get_same_no($junction_seq,$ref_seg);
	      if ($same_no > 22){
	        #写输出，输出成一个read
	        $whole_seq = $junction_seq.$seq_1.$seq_2.$seq_3;
	        $whole_qua = reverse($junction_read[10]).$qua_1.$qua_2.$qua_3;
	        $cigar1 = $i;
	        $cigar2 = 100-$i;
		$read_start = $ref_1_start - $i + 1;
	        $N_no = $coor1 - $read_start - 25;
	        $out_read = "$read_no\t$flag\t$seg_chr_1\t$read_start\t$quality\t${cigar1}M${N_no}N${cigar2}M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	        print Out $out_read;
	      } 
	    }
	  }
	}

      }

      #12_34的情况
      if ($chr4 eq $chr3){
	if ($coor1+25 == $coor2 and $coor3+25 == $coor4 and $min_gap < $coor3 - $coor2 and $coor3 - $coor2 < $max_gap){
	  #prepare for the output
	  $seg1 = $reads_hash{$sort_keys[$i]};
          $seg2 = $reads_hash{$sort_keys[$i+1]};
          $seg3 = $reads_hash{$sort_keys[$i+2]};
          $seg4 = $reads_hash{$sort_keys[$i+3]};
	  ($reads_tophat_name,$flag,$seg_chr_1,$seg_coo_1,$quality,$cigar,$_,$_,$_,$seq_1,$qua_1,$tail) = split(/\t/,$seg1);
	  ($reads_real_name,$_) = split(/_/, $reads_tophat_name);
	  @seg2 = split(/\t/, $seg2);
	  @seg3 = split(/\t/, $seg3);
	  @seg4 = split(/\t/, $seg4);
	  $N_no = $coor3 - $coor1 -50;
	  $whole_seq = $seq_1.$seg2[9].$seg3[9].$seg4[9];
	  $whole_qua = $qua_1.$seg2[10].$seg3[10].$seg4[10];
	  $out_read = "$reads_real_name\t$flag\t$seg_chr_1\t$seg_coo_1\t$quality\t50M${N_no}N50M\t*\t0\t0\t$whole_seq\t$whole_qua\t$tail\n";
	  print Out $out_read;
	}
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

