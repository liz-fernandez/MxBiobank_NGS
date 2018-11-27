---
layout: page
title: NGS – from fastq to variant annotation
subtitle: Quality control of NGS data
minutes: 5
---

> ## Learning objectives {.objectives}
>
> *   Understand the FastQ format.
> *   Understand the concept of sequence quality. 
> *   Learn how to use software to measure and improve the quality of massive sequencing data.

## Downloading the data

The first thing we do when receiving our data is to download them. They are usually in a compressed format. We can decompress the data using the `unzip` command:

~~~ {.bash}
$ unzip 'Alignment&VariantCalling.zip'
~~~

This command decompresses a directory called `FastQC_Short`.
We enter that directory:

~~~ {.bash}
$ cd VariantCalling
~~~
 
Revise its contents:

~~~ {.bash}
$ ls
~~~

~~~ {.output}
lane1  lane2  mito.intervals  Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa
~~~

We will now look at the contents of the two "lane" directories:

~~~ {.bash}
$ ls lane*
~~~

~~~ {.output}
lane1:
s-7-1.fastq  s-7-2.fastq

lane2:
s-7-1.fastq  s-7-2.fastq
~~~

> ## The name of fastq files {.callout}
>
> In this case, the name of the files has been assigned in a
> arbitrary but often contain important information.
> For example, the Illumina files generally have the following format:
>
> ~~~ {.output}
> <nombre de la muestra>_<secuencia identificadora (barcode)>_L<línea (3 dígitos)>_R<número de lectura (read)>_<número del set (3 dígitos)>.fastq.gz
> ~~~
> Example: NA10831_ATCACG_L002_R1_001.fastq.gz
>
> It is ideal to respect these names to avoid losing information.
 
Let's review one of the files using the head command. It usually shows us the first
10 lines of a file but with the `-n` flag will show us the first 12:

~~~ {.bash}
$ head -n 12 ./lane1/s-7-1.fastq
~~~

~~~ {.output}
@IL29_4505:7:24:8932:6562#2/1
TAACGGTGGGTGAGTGGTAGTAAGTAGAGGGATGGATGGTGGTTCGGAGTGGTATGGTTGAATGGGACAGGGTAACGAGTGGAGAGTAGGGTAATGGAGGGTAAGTTC
+
CDDCDDABBBABABABB@BCACBDABCBBAB@BBCABBBABB?CBCCABABBABBBBABA?ACBAAAAA?BB;BCAABA7AA?B?A??AAA>?A:AA?AA?%?AA@=9
@IL29_4505:7:15:7929:11873#2/1
TTAACGTTTCAATATGGTAGGTAGAACAACAGTACAGTGAGTGGGACATGGTGGATGGTAAAAGAATGGTAGGGTAAGTGGCAGTGGGGTTGGATATGGGTAATTGGA
+
DDDDDDCDDCDDDDCCDBDDDACCCCCCCDCCBCDCCBCBCACCCBCBCCCACCACCCAAABABAABBC?ABBA=BBA?BAAAB?BABA7AAB?AAAAA@:A?AA?B-
@IL29_4505:7:109:18599:8300#2/1
CTTACCCTCCATTACCCTACCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAACCACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCT
+
D?BDCABBBABBDACCBBDCDCBBBBABBBCDABBBABBBBBCBA@ABCBCBB@BB0=BB4><B:BABAABBABA@AB=BB?ACBABA?AB;BBBABBBB'0=CB9<9
~~~

This file is in `fastq` format. This is the format that most modern sequencing platforms generate. This format is a modification of the sequence format
[Fasta](https://en.wikipedia.org/wiki/FASTA_format). You can recognize them because they end
with the extension `.fq` or `.fastq`.

FastQ files are plain ASCII text files that contain both the sequences
detected by the sequencer as well as information about the quality of each one
of the sequenced nucleotides.

> ## Do not change fastq files! {.callout}
>
> The fastq files that you receive after a sequencing experiment will be
> requested when you want to publish an article. Ideally, make a
> couple of copies on different devices as soon as you receive them and *do not modify them manually*.
> 
> It is also recommended that, if you have the information at hand, you create a small
> README file (plain text) that is in the same directory as your fastq files
> and describe how the experiment was performed, how many replicas there are of each condition, which
> indexes and controls were used, etc. This will help you greatly when trying to
> interpret your results.

We will take a single sequence as an example:

~~~ {.output}
@IL29_4505:7:24:8932:6562#2/1
TAACGGTGGGTGAGTGGTAGTAAGTAGAGGGATGGATGGTGGTTCGGAGTGGTATGGTTGAATGGGACAGGGTAACGAGTGGAGAGTAGGGTAATGGAGGGTAAGTTC
+
CDDCDDABBBABABABB@BCACBDABCBBAB@BBCABBBABB?CBCCABABBABBBBABA?ACBAAAAA?BB;BCAABA7AA?B?A??AAA>?A:AA?AA?%?AA@=9
~~~

We see that each sequence is represented by 4 lines.

1. **Sequence name** - This line begins with the symbol `@`, followed by the read name. 
2. The nucleotide **sequence**. 
3. **Segundo título** - Esta línea comienza con el símbolo `+`. Generally the information is the same as in the first line but it can also be blank as long as it begins with the symbol `+`.
4. **Quality information** - It contains an ASCII character string. Every
one of these characters corresponds to a nucleotide of the sequence and represents the quality of the same. The score of each nucleotide indicates the
level of confidence that you have in that base. High levels indicate that the base has been reported correctly while low levels suggest
that there is uncertainty about the actual sequence in this position.

The name of the sequence also contains important information. In
In this case, the Illumina data provide the following information:

~~~ {.output}
@<instrument>:<run>:<flow cell identifier>:<lane>:<square>:<x-coord>:<y-coord> <read>:<if the ead was filtered>:<control number>:<sequence index> length=<sequence size>
~~~

The level of quality represented in the fourth line is called Phred score. In its
original format a Phred score is simply the logarithm of the probabilities of
error:

~~~ {.output}
Phred score = - 10 * log10(error probability)
~~~

| Quality score | Error probability |
|:----------------------|:----------------------|
| Q40 0.0001 | (1 en 10,000) |
| Q30 0.001 | (1 en 1,000) |
| Q20 0.01 | (1 en 100) |
| Q10 0.1 | (1 en 10) |

Wait, if the Phred scores are numbers, why do the quality scores shown
in line four are a mixture of alphanumeric symbols?

~~~ {.output}
#1:=BDDDHDHB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
~~~

This is because these scores are encoded in ASCII (American Standard Code for Informational Interchange).
In short, ASCII is a coding system that equates a number to a character
alphanumeric. For example, the character 'A' is represented by the number 65 in the table
ASCII code, while '%' is represented by the number 37. This system allows us to
represent a total of 256 different characters.

Why use ASCII? Given the huge amount of data produced during massive parallel sequencing, it's all about reducing the data to the maximum. The ASCII system allows us to represent
two-digit numbers in a single bit (8 bits), which reduces the space they occupy
these files. If you are interested in the subject, you can read more about binary coding.

Finally, since this quality coding system was invented, which is already
used with Sanger-type sequencing, different companies have "slipped" the scales
of ASCII conversion because what is important to verify with the company or laboratory
in which version of this scale have been coded to use the correct conversion. Some
of the quality control tools that we will use infer the type of coding
used directly from the data provided.

## Verifying the quality of the sequences

We will verify the quality of the sequences using the program [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
This is one of the most used programs for this type of analysis,
since it provides a global vision of how the data is structured, as well as
being pretty fast.

FastQC reads a series of sequencing files and produces a report of
quality control for each one, which consists of a number of different modules,
each of which helps us identify different problems in our data.

Let's analyze our first sequencing file,
first we created a directory to store our results:

~~~ {.bash}
$ mkdir QUAL
~~~

And we perform our quality analysis using FastQC:

~~~ {.bash}
$ fastqc -O ./QUAL/ lane1/s-7-1.fastq
~~~

~~~ {.output}
Started analysis of s-7-1.fastq
Approx 5% complete for s-7-1.fastq
Approx 10% complete for s-7-1.fastq
Approx 15% complete for s-7-1.fastq
Approx 20% complete for s-7-1.fastq
Approx 25% complete for s-7-1.fastq
Approx 30% complete for s-7-1.fastq
Approx 35% complete for s-7-1.fastq
Approx 40% complete for s-7-1.fastq
Approx 45% complete for s-7-1.fastq
Approx 50% complete for s-7-1.fastq
Approx 55% complete for s-7-1.fastq
Approx 60% complete for s-7-1.fastq
Approx 65% complete for s-7-1.fastq
Approx 70% complete for s-7-1.fastq
Approx 75% complete for s-7-1.fastq
Approx 80% complete for s-7-1.fastq
Approx 85% complete for s-7-1.fastq
Approx 90% complete for s-7-1.fastq
Approx 95% complete for s-7-1.fastq
Analysis complete for s-7-1.fastq
~~~ 

Once finished, let's review the result:

~~~ {.bash}
$ cd QUAL
$ ls
~~~

~~~ {.output}
s-7-1_fastqc.html  s-7-1_fastqc.zip
~~~

The easiest way to explore these results is by opening the html file in
your browser. You can do it by double clicking on the file.

The results should be similar to [this](http://liz-fernandez.github.io/DATA/s-7-1_fastqc.html).

This page contains a lot of information broken down into the following sections:

1. **Basic Statistics** -  The basic statistics of the experiment.
2. **Per base sequence quality** - Box diagrams showing the quality of each base.
3. **Per tile sequence quality** - It contains the box (tile) from which each sequence comes. It only appears if the analyzes are made in an Illumina library that retains its original identifiers.
4. **Per sequence quality scores** - It allows to see if a subgroup of sequences has bad quality.
5. **Per base sequence content** - Shows the proportion of each base in each position. 
6. **Per sequence GC content** - Compares the observed GC distribution with a modeled GC distribution.
7. **Per base N content** - Indicates how many nucleotides could not be interpreted and are indicated with an N.
8. **Sequence Length Distribution** - Shows the sequence length distribution.
9. **Sequence Duplication Levels** - Indicate how many times each sequence is repeated. It is done only in the first 100,000 sequences.
10. **Overrepresented sequences** - Shows over represented sequences.
11. **Adapter Content** - Shows the content of adapters in the library.

Let's review the results once more by clicking on the [link](http://liz-fernandez.github.io/DATA/s-7-1_fastqc.html). 

In general we see that most of the criteria have a green tick of quality control pass, with the exception of 
per base sequence content, GC content per base and overrepresented sequences.

Although we did not provide this information, the program has inferred the coding
ASCII (Sanger / Illumina 1.9) as well as the number of sequences, their size and
its GC content.

When we observe the **quality per base**, we see that most of the bases are in the area
green or have good quality. Each base is represented by a box diagram that
provides a good idea about the distribution of quality scores in all
the sequences. It is evident that the quality decreases more markedly
at the end of the sequence. This is known as small errors accumulate while
the longer the sequence, this is one of the reasons why sequencers based
in nucleotide incorporation they have a limit on the maximum length of their readings.

The graph showing the **quality per tile** is almost completely blue showing that the average distribution of errors in the tile is very similar.

The graph of **quality score per sequence** shows that most of the sequences have high scores.

However, the **distribution of nucleotides per base** shows that, at the beginning of the sequence, the distribution is very disparate. 

We also observe that:

* The **GC distribution per sequence** is somewhat different to the hypothetical distribution.
* There are few **ambiguous bases (Ns)**.
* All sequences have the same **size**.
* Levels of **duplication of sequences** are very low.
* There are a few **over represented sequences**. Mostly Illumina PCR primers.
* There are no **adapters**.

When running FastQC, a compressed file called
`s-7-1_fastqc.zip`. 

We unpack this file using the following command:

~~~ {.bash}
$ unzip s-7-1_fastqc.zip
~~~

~~~ {.output}
Archive:  s-7-1_fastqc.zip
   creating: s-7-1_fastqc/
   creating: s-7-1_fastqc/Icons/
   creating: s-7-1_fastqc/Images/
  inflating: s-7-1_fastqc/Icons/fastqc_icon.png
  inflating: s-7-1_fastqc/Icons/warning.png
  inflating: s-7-1_fastqc/Icons/error.png
  inflating: s-7-1_fastqc/Icons/tick.png
  inflating: s-7-1_fastqc/summary.txt
  inflating: s-7-1_fastqc/Images/per_base_quality.png
  inflating: s-7-1_fastqc/Images/per_tile_quality.png
  inflating: s-7-1_fastqc/Images/per_sequence_quality.png
  inflating: s-7-1_fastqc/Images/per_base_sequence_content.png
  inflating: s-7-1_fastqc/Images/per_sequence_gc_content.png
  inflating: s-7-1_fastqc/Images/per_base_n_content.png
  inflating: s-7-1_fastqc/Images/sequence_length_distribution.png
  inflating: s-7-1_fastqc/Images/duplication_levels.png
  inflating: s-7-1_fastqc/Images/adapter_content.png
  inflating: s-7-1_fastqc/fastqc_report.html
  inflating: s-7-1_fastqc/fastqc_data.txt
  inflating: s-7-1_fastqc/fastqc.fo
~~~

If we enter the newly created directory, we can see the report in plain text format,
both the summary:

~~~ {.bash}
$ cd s-7-1_fastqc
$ more summary.txt
~~~

~~~ {.output}
PASS	Basic Statistics	s-7-1.fastq
PASS	Per base sequence quality	s-7-1.fastq
PASS	Per tile sequence quality	s-7-1.fastq
PASS	Per sequence quality scores	s-7-1.fastq
WARN	Per base sequence content	s-7-1.fastq
WARN	Per sequence GC content	s-7-1.fastq
PASS	Per base N content	s-7-1.fastq
PASS	Sequence Length Distribution	s-7-1.fastq
PASS	Sequence Duplication Levels	s-7-1.fastq
WARN	Overrepresented sequences	s-7-1.fastq
PASS	Adapter Content	s-7-1.fastq
~~~

and the full report file:

~~~ {.bash}
$ more fastqc_data.txt
~~~

We have omitted the content since it is very long.
Being able to review the results in flat files is extremely useful when
results are in a server that does not have a graphical interface.

FastQC programmers have compiled examples of quality sequencing data
diverse Let's analyze them:

* [Good data - Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
* [Bad data - Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)
* [Dimer contamination](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html)
* [Small RNAs with read through the adaptor ](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/small_rna_fastqc.html)

> ## How valid is our analysis? {.challenge}
>
> Since the data we used only less than 200,000 sequences, can we
> trust that the complete data will have similar quality? Why or why not?

> ## How do the rest of the files look? {.challenge}
>
> Obtain the fastQC report html files for all other three fastq files.
>

## Cleaning the sequences

There are a number of tools to clean sequences and different researchers
have different preferences. How with other sequencing problems there is a
active development of programs for this purpose. When the data is of decent quality there is no need to carry out this step for variant annotation.

Here are some examples of programs you can use in case your data is not great. 

- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [cutadapt](https://cutadapt.readthedocs.io/)





