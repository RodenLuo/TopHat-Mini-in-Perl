<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [TopHat-Mini-in-Perl](#tophat-mini-in-perl)
  - [Requirements](#requirements)
  - [Introduction to TopHat pipeline](#introduction-to-tophat-pipeline)
  - [Preprocessing from fastq](#preprocessing-from-fastq)
  - [Tophat-Mini-in-Perl](#tophat-mini-in-perl)
  - [View the results](#view-the-results)
  - [Putative novel loci for WASH7P](#putative-novel-loci-for-wash7p)
      - [Other sites with anotation:](#other-sites-with-anotation)
      - [Other sites without anotations:](#other-sites-without-anotations)
      - [Evidence from BLAST:](#evidence-from-blast)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

# TopHat-Mini-in-Perl
Implementation of the main function, splicing, of TopHat in perl

This is the course project for Bioinformatics(BI3204 2015.03-2015.07) at [SUSTC](http://www.sustc.edu.cn/). 

# Table of Contents
  * [Requirements](#requirements)
  * [Introduction to TopHat pipeline](##introduction-to-tophat-pipeline)
  * [Preprocessing from fastq](#Preprocessing-from-fastq)
  * [Tophat-Mini-in-Perl](#Tophat-Mini-in-Perl)
  * [View the results](#View-the-results)
  * [Putative novel loci for WASH7P](#putative-novel-loci-for-WASH7P)
        * [Other sites with anotations:](#other-sites-with-anotations)
        * [Other sites without anotation:](#Other-sites-without-anotation)
        * [Evidence from BLAST:](#Evidence-from-BLAST)

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

##Output spliced reads only:
head -26 unaln_splitted_srt_name.sam >final.sam  
cat 14.sam >>final.sam
cat 23.sam >>final.sam
samtools sort -o final_srt.sam -O sam -T temp final.sam
## add header, combine 23.sam and 14.sam and output into final.sam
## get final sorted sam in fina_srt.sam

## OR output spliced reads and the aligned intact reads
cat aligned.sam >final.full.sam
cat 14.sam >>final.full.sam
cat 23.sam >>final.full.sam
samtools sort -o final_srt.sam -O sam -T temp final.sam
## combine 23.sam and 14.sam as well as aligned.sam, output into final.sam
## get sorted sam in fina_srt.sam
```

## View the results
```bash
samtools view -bS final_srt.sam > final_srt.bam  # Convert it to bam
samtools index final_srt.bam # Index it by samtools index
# View in IGV or other alignments viewer.
```
overall_figure

![overall_figure](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/overall_figure.png)

zoom_in

![zoom_in](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/zoom_in.png)

## Putative novel loci for WASH7P
View the same bam at different position let me find the putative novel gene loci

####Other sites with anotation:
chr2: WASH2P, DDX11L2

![chr2_WASH2P_DDX11L2](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chr2_WASH2P_DDX11L2.png)

chr9: WASH1, DDLX11L5

![chr9_WASH1_DDLX11L5](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chr9_WASH1_DDLX11L5.png)

chr12: NR_028269

![chr12_NR_028269.png](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chr12_NR_028269.png)

chr15: WASH3P, DDX11L9

![chr15_WASH3P_DDX11L9](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chr15_WASH3P_DDX11L9.png)

####Other sites without anotations:
chr16: ?, DDX11L10

![chr16_?_DDX11L10](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chr16_%3F_DDX11L10.png)

chrX:

![chrX](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chrX.png)

chrY:

![chrY](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/chrY.png)

####Evidence from BLAST:
Use WASH7P gene sequence([NR_024540.1.fa](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/NR_024540.1.fa)) to do BLAST, results:

![BLAST_result](https://github.com/RodenLuo/TopHat-Mini-in-Perl/blob/master/images/putative_novel_loci_for_WASH7P/BLAST_result.png)

>Reference link: http://blast.ncbi.nlm.nih.gov/Blast.cgi

So, the next step, we may need to explore in detail to see whether these three positions (on: chr16, chrX, chrY) is novel gene loci or not.
