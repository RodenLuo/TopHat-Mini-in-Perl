# TopHat-Mini-in-Perl
Implementation of the main function, splicing, of TopHat in perl

This is the course project for Bioinformatics(BI3204 2015.03-2015.07) at [SUSTC](http://www.sustc.edu.cn/). 


## Requirements
Bowtie2, already tested on version 2.2.4

samtools, already tested on version 1.2

Human genome reference file, already tested on hg19.fa

## Introduction to TopHat pipeline
![Image of TopHat pipeline](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/tophat_pipeline.png)

>Reference link: [Trapnell C, Pachter L, Salzberg S L. TopHat: discovering splice junctions with RNA-Seq[J]. Bioinformatics, 2009, 25(9): 1105-1111.](http://bioinformatics.oxfordjournals.org/content/25/9/1105.full)

In this perl implementation, we split unaligned reads(100bp) into 4 segments(1,2,3,4), first we find junction position at segment2 and segment3, which is the case 1_34 and 12_4 can be aligned. And then use that junction postion to help find splicing point for segment1 and segment2. So this implementation does not need known junction site, which is usually recorded in a gtf file.

## Preprocessing from fastq
```bash
bowtie2 -x <ref> -U clean.fastq -S aligned.sam --un unaligned.fastq -p 10 &>bowtie2_align_output.txt

# get the unaligned reads in fastq by adding --un option
# -p means the number of processors going to be used, change it if needed
```
```bash
perl split_reads.pl unaligned.fastq unaln_splitted.fastq

# split unaligned reads into 4 segments, and add identifiers following the reads' name
```

```bash
bowtie2 -x <ref> -U unaln_splitted.fastq -S unaln_splitted.sam -a -p 10 &>bowtie2_unaln_splitted.txt

# -a output all alignments
# One can use -k mode: search for one or more alignments, report each
# Larger value or -a yields more alignments for the splitted reads, which will inturn lead to more splice possibilities but a lower speed for both bowtie2 and TopHat-Mini-in-Perl.
```

```bash
samtools sort -o unaln_splitted_srt.sam -O sam -T temp unaln_splitted.sam

samtools sort -n -o unaln_splitted_srt_name.sam -O sam -T temp unaln_splitted_srt.sam
```
