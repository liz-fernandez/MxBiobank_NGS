---
layout: page
title: NGS – from fastq to variant annotation
subtitle: Read alignment
minutes: 5
---
> ## Learning Objectives {.objectives}
>
> * Align sequencing data to a reference genome
> * Understand how to interpret high-throughput sequencing alignment data
> * First approach to the SAM and BAM coordinate formats

We will use the fastq files that we used in the previous practice, as well as the reference genome of our organism, *Saccharomyces pombe*:

~~~ {.bash}
$ cd .. 
$ pwd 
~~~

~~~ {.output}
/Data/VariantCalling
~~~

Characteristics of the experiment:

Yeast genome: 12.5 Mbp; 16 chromosomes
Whole genome sequencing
Paired-end reads, 108bp, one library, 2 lanes

## Mapping the reads to the genome

Once we verify that the reads are in the correct format, we will align the reads and transcripts to the genome using BWA.

You can find the manual in the following[link](http://bio-bwa.sourceforge.net/bwa.shtml).


~~~ {.bash}
$ bwa
~~~

~~~ {.output}
Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
         pemerge       merge overlapping paired ends (EXPERIMENTAL)
         aln           gapped/ungapped alignment
         samse         generate alignment (single ended)
         sampe         generate alignment (paired ended)
         bwasw         BWA-SW for long queries

         shm           manage indices in shared memory
         fa2pac        convert FASTA to PAC format
         pac2bwt       generate BWT from PAC
         pac2bwtgen    alternative algorithm for generating BWT
         bwtupdate     update .bwt to the new format
         bwt2sa        generate SA from BWT and Occ

Note: To use BWA, you need to first index the genome with `bwa index'.
      There are three alignment algorithms in BWA: `mem', `bwasw', and
      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
      first. Please `man ./bwa.1' for the manual.
~~~

First we will generate a bwa index for the genome:

~~~ {.bash}
$ bwa index -a is Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa
~~~

~~~ {.output}
[bwa_index] Pack FASTA... 0.12 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 7.56 seconds elapse.
[bwa_index] Update BWT... 0.07 sec
[bwa_index] Pack forward-only FASTA... 0.07 sec
[bwa_index] Construct SA from BWT and Occ... 3.67 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index -a is Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa
[main] Real time: 12.089 sec; CPU: 11.490 sec
~~~

~~~ {.bash}
$ ls Saccharomyces*
~~~

~~~ {.output}
Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa  
Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa.amb
Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa.ann
Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa.bwt
Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa.pac
Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa.sa
~~~

All new files that do not end with .fa form the bwa index. 

For aligning the reads we will be using bwa mem, we check its manual:

~~~ {.bash}
$ bwa mem
~~~

Some commonly used options include:

-t	number of threads/processors to use – see PBS script at end of workbook
-p 	Assume the first input query file is interleaved paired-end FASTA/Q. See the command description for details. 
-a	Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
 
We will use default options:

~~~ {.bash}
$ bwa mem Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa lane1/s-7-1.fastq lane1/s-7-2.fastq > lane1.sam
~~~ 

We omit the output.

We explore the result, which is a file in SAM format.

~~~ {.bash}
$ head lane1.sam
~~~ 

~~~ {.output}
@SQ	SN:I	LN:230218
@SQ	SN:II	LN:813184
@SQ	SN:III	LN:316620
@SQ	SN:IV	LN:1531933
@SQ	SN:IX	LN:439888
@SQ	SN:Mito	LN:85779
@SQ	SN:V	LN:576874
@SQ	SN:VI	LN:270161
@SQ	SN:VII	LN:1090940
@SQ	SN:VIII	LN:562643
~~~ 

~~~ {.bash}
$ tail lane1.sam
~~~ 

~~~ {.output}
IL29_4505:7:120:19585:18832#2	77	*	0	0	*	*	0	0	GATCGGAAGAGCGGTTCAGCAGGAAATGCCGAGACCGATCTCCGATGTTTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAACAAACCACATAACTACATCTCCACG	CCBCBDCCCBBBBBBCBCBBBDCBCBAACBCBBBBAAAAAB=BB<?BABAAA.ABBAA?89BBCCBA;AA-B0B@AAAAB(--(%$&$-3$%$'*6<).635%'+*/*	AS:i:0	XS:i:0
IL29_4505:7:120:19585:18832#2	141	*	0	0	*	*	0	0	GATCGGAAGAGCGGCGGGGAGGGGAAGAGGGGAGATCTCGGGGGGCGCCGGATCATTAAAAAAAAAAAAAAAAAAAGGACGTAAACCGTGATATCTCCGTCCGCGTGA	@64<;?27134)2)0.&1'$8=9.14-'2&1(/*/,,2424*49'&1%0,&)/'-*0(43>>>>>>>>33&.*&)%$'&$&$(%$$$$*$%)%&)$%($),%$+$&$&	AS:i:0	XS:i:0
IL29_4505:7:120:19710:18114#2	77	*	0	0	*	*	0	0	GATCGGAAGAGCGGTTCAGCAGGGAATGCCGAGACCGATCTCCGATGTTTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAACACACAAACACATGACTGTGACGATG	AABBBBBBBABBBBABABB:BBBAABBBBBBBBAB?AABABBABBABAABAA5?ABA?BB@<B@@<9>>;)=:?4AA=9B,3(0$(-'(%-)0&%*.,'*$$,&&-13	AS:i:0	XS:i:0
IL29_4505:7:120:19710:18114#2	141	*	0	0	*	*	0	0	GATCGGAAGAGCGTCGGGTAGGGAAAGAGGGAAGATCTCGGGGGCCGCCGTATCATTAAAAAAAAAAAAAACACAACAATCAAACTCTACCACACCTGGCACACGCGA	?4+=@?3333,&7%-:,/0((;<-6/).5*2&16)(02*77&52&&1'*+*-//0333-5=71144%)+$$'%$%&&&(%%%-%*%-$$*&&*$$$$*$-$,%'%%'%	AS:i:0	XS:i:0
IL29_4505:7:120:19764:17570#2	77	*	0	0	*	*	0	0	GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCCGAAGTTTATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAACTAAGACACGATACCATCTCACCGCT	BBCB@?AAA?BABAA@C:CAAAB<AC?@7C<B?8BACABBB?BB<C?@?B?*9?==B?@??=ABAAABA?BAAAA0;97;(.%+<1*&*&(%&&$(',&(%+.,+5=%	AS:i:0	XS:i:0
IL29_4505:7:120:19764:17570#2	141	*	0	0	*	*	0	0	GATCGGAAGCGCGGCGTGTGGGGCCAGAGTGTAGAGCTGGGTGGTCGCCGCATACTTTAAAAAAGATAATGAGATCATGCGAAGTGTATCTGCAAGTGGATCTAATAC	<@+>=9(93);72'=9,,)+222&'(5*,.7'10&'3,'6*'1<%&9)*4%/&(%0.%*)5<>4%/$%%%)$('(+,%.*-&)/(%'(&%(*5&)'$&()31*,).-%	AS:i:0	XS:i:0
IL29_4505:7:120:19773:11846#2	77	*	0	0	*	*	0	0	CCGATGTTTATCTCGTATGCCGTCTTCTGCTTGAAAACACAAAAGAAACATACTCTATCAAGATACATGCTAGTAACGAAAACACACACATGACTCCACACTGTGATG	A@@*@@6@@@@>@@<:@;@@:A;A?@/@4@@6@@@:A+;4>&9,$$,'%%'1-(%6%/:,(*/%&&%'(7,55(04*''84%8,-.%+1/(&)&&)%&%&%%&*$%1.	AS:i:0	XS:i:0
IL29_4505:7:120:19773:11846#2	141	*	0	0	*	*	0	0	GATTTCATAAACACAGGCCTGAATGTAGACAGGTCACTAGAGATACCGTCACAAACAAAACATCAATACTAACAACAACTACTACAATCACCACACACATCACCACAC	*.(=4'0*'1'/;6/-(**(20((46'('2**(5-*,5(**'8764.*.((</'./'&3(.6,*(;-'*-*(8,+''*/(<0*(*'.'((-7'*3/'/)'(7)'5/'0	AS:i:0	XS:i:0
IL29_4505:7:120:19817:15837#2	77	*	0	0	*	*	0	0	ATTTTTAGTATACGTGAATAATGTATTTTTTAAATTGGTTATTGCGTAGGAACACAATATTTTAAACCTTATCGATTAATTTTACTTGACTGAGTACGAGCGCCAGCT	ABA8A@BB8BBBBBA?BBBB;ABAB@>9?B/;@8B@>A??B6:BA@AAB?BB/=A:?@;AB@BB8ABB@BA@B>/ABB>BB-?B/A3.<=.=>3<B/<B82'077.*9	AS:i:0	XS:i:0
IL29_4505:7:120:19817:15837#2	141	*	0	0	*	*	0	0	AAAAAAAAAAAAAACGGTGGACTAAAACACCCGGTCGTTGAAACAAAATTTGATAAATGTCTAAAGGCATTTACATGTCCTACATCCAATCATATGACTTTTGTTTCC	=2=,=,(-8179+4('(-))1&)-'(,/2**''&)*83(&:21))7)'6++%1)('3%'))%%'2%&(:2''1&*2&''&'8+040&3*2,<'2-%3&&'--'(*'%-	AS:i:0	XS:i:0
~~~

### The SAM format

Let's see a smaller example of this format.
Suppose we have the following alignment:

~~~ {.output}
Coor	12345678901234 5678901234567890123456789012345
ref	AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
+r001/1	      TTAGATAAAGGATA*CTG
+r002	     aaaAGATAA*GGATA
+r003	   gcctaAGCTAA
+r004	                 ATAGCT..............TCAGC
-r003	                        ttagctTAGGC
-r001/2	                                      CAGCGGCAT
~~~

The corresponding SAM format will be the following:

~~~ {.output}
@HD	VN:1.5	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r003	0	ref	9	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r004	0	ref	16	30	6M14N5M	*	0	0	ATAGCTTCAGC	*
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
~~~

The SAM format is a plain text format that allows us to save sequencing data
in ASCII format delimited by tabs. 

It is made up of two core sections:

*  The headers
*  The alignment

The section of the **header** starts with the character `@` followed by one of the codes
of two letters that denote the characteristics of the alignments in this file.
Each line is delimited by tabs and, in addition to the lines that begin with
`@CO`, each data field has the format` TAG:VALUE`, where `TAG` is a string
Two characters that define the format and content of `VALUE`.

The header is not indispensable but contains information about the version of the
file as well as if it is ordered or not. Therefore it is advisable to include it.

The **alignment section** contains the following information:

1. **QNAME** Name of the reference, QNAME (SAM) / Name of the read (BAM).
It is used to group alignments that are together, as in the case of alignments
of paired end reads or a read that appears in multiple alignments.
2. **FLAG** Information set describing the alignment. Provides the following information:
	* Are there multiple fragments?
	* Are all the fragments well aligned?
	* Is this fragment aligned?
	* Has the following fragment not been aligned?
	* Is the reference the reverse string?
	* Is the following fragment the reverse chain?
	* Is this the first fragment?
	* Is this the last fragment?
	* Is this a secondary alignment?
	* Did this reading fail quality filters?
	* Is this reading a PCR or optical duplicate?
3. **RNAME** Name of the reference sequence.
3. **POS** Left alignment position (base 1).
3. **MAP** Alignment quality.
3. **CIGAR** CIGAR chain.
3. **RNEXT** Name of the paired end (mate) or the next read.
3. **PNEXT** Position of the pair (mate) or the next read.
3. ** TLEN** Length of the alignment.
3. **SEQ** The test sequence of this alignment (in this case the sequence of the read).
3. **QUAL** The read quality.
3. **TAGs** Additional information.

> ## Cadenas CIGAR {.callout}
> The sequence aligned to the reference may have additional bases that are not in
> the reference or the read may be missing bases that are in the reference.
> The CIGAR string is a string that encodes each base and the characteristics of each in
> alignment.
>
> For example, the string CIGAR:
>
> ~~~ {.output}
> CIGAR: 3M1I3M1D5M
> ~~~
>
> indicates that the first 3 bases of the reading aligns with the reference (3M), the next base
> does not exist in the reference (1I), the following 3 bases align with the reference (3M), the
> next base does not exist in the reading (1D), and 5 more bases align with the reference (5M).

As you can see these files contain a lot of information that can be analyzed
using scripts that show alignment statistics. Programs like
[Picard](http://broadinstitute.github.io/picard/) perform
this type of analysis.

The compressed version of the SAM files is known as BAM (binary sam).
Let's convert the SAM file to BAM using samtools:

~~~ {.bash}
$ samtools view -b lane1.sam -o lane1.bam
~~~

We do not open the BAM files because, since they are in binary format, they are illegible.
However, we have to perform two last steps to visualize
the results:

~~~ {.bash}
$ samtools sort lane1.bam -o lane1_sorted.bam
$ samtools index lane1_sorted.bam
~~~

The first step orders the results by their coordinates and the second one creates indexes
to speed up the display using a browser.




















