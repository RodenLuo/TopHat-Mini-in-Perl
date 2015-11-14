# TopHat-Mini-in-Perl
Implementation of the main function, splicing, of TopHat in perl

This is the course project for Bioinformatics(BI3204 2015.03-2015.07) at [SUSTC](http://www.sustc.edu.cn/). 


## Requirements
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), already tested on version 2.2.4

[samtools](http://www.htslib.org/), already tested on version 1.2

Human genome reference file, already tested on hg19.fa

## Introduction to TopHat pipeline
![Image of TopHat pipeline](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/tophat_pipeline.png)

>Reference link: [Trapnell C, Pachter L, Salzberg S L. TopHat: discovering splice junctions with RNA-Seq[J]. Bioinformatics, 2009, 25(9): 1105-1111.](http://bioinformatics.oxfordjournals.org/content/25/9/1105.full)

In this perl implementation, we split unaligned reads(100bp) into 4 segments(1,2,3,4), first we find junction position at segment2 and segment3, which is the case 1_34 and 12_4 can be aligned. And then use that junction postion to help find splicing point for segment1 and segment2. So this implementation does not need known junction site, which is usually recorded in a gtf file.

## Preprocessing from fastq
```bash
bowtie2 -x <ref.fa> -U clean.fastq -S aligned.sam --un unaligned.fastq -p 10 &>bowtie2_align_output.txt
# get the unaligned reads in fastq by adding --un option
# -p means the number of processors going to be used, change it if needed
```
```bash
perl split_reads.pl unaligned.fastq unaln_splitted.fastq
# split unaligned reads into 4 segments, and add identifiers following the reads' names
```

```bash
bowtie2 -x <ref.fa> -U unaln_splitted.fastq -S unaln_splitted.sam -a -p 10 &>bowtie2_unaln_splitted.txt
# -a output all alignments
# One can use -k mode: search for one or more alignments, report each
# Larger value or -a yields more alignments for the splitted reads, which will inturn lead to more splice possibilities but a lower speed for both bowtie2 and TopHat-Mini-in-Perl.
```

```bash
samtools sort -o unaln_splitted_srt.sam -O sam -T temp unaln_splitted.sam
samtools sort -n -o unaln_splitted_srt_name.sam -O sam -T temp unaln_splitted_srt.sam
# two steps sort, first sort by position, then by reads' names which includes the identifers
```

## Tophat-Mini-in-Perl

```bash
grep -n '>chr' <USCS ref.fa> > ref.tdx
# I used this method to build a special index, tdx.
# Sorry that when I rewrote it, I did not realize this step may cause a probelm for the reference which is not downloaded from UCSC.
```

```bash
perl splice.pl unaln_splitted_srt_name.sam 23.sam <ref.tdx> <ref.fa> [max gap]
# Find splice junction point at segments 2 and 3.
# This command will output 23.sam and 23.sam.splice_junction.txt. The second stores the position which is needed by the next step.
# Max gap is optional. [20000]

perl ht_splice.pl unaln_splitted_srt_name.sam 14.sam 23.sam.splice_junction.txt <ref.tdx> <ref.fa> [max gap]
# Find splice junction point at segments 1 and 4.
# This command will output 14.sam.
# Max gap is optional. [20000]
```

```bash
head -26 unaln_splitted_srt_name.sam >final.sam  
cat 14.sam >>final.sam
cat 23.sam >>final.sam
samtools sort -o final_srt.sam -O sam -T temp final.sam
## add header, combine 23.sam and 14.sam and output into final.sam
## get final sorted sam in fina_srt.sam
```

## See the results
