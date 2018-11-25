---
layout: page
title: NGS – from fastq to variant annotation
subtitle: Quality control of NGS data
minutes: 5
---

> ## Learning objectives {.objectives}
>
> *   Entender el formato FastQ.
> *   Entender el concepto de calidad de una secuencia. 
> *   Aprender a utilizar software para medir y mejorar la calidad de datos de secuenciación masiva.

## Downloading the data

The first thing we do when receiving our data is to download them. They are usually in a compressed format called tar. We can decompress the data using that same command:

~~~ {.bash}
$ tar -xvf FastQC_Short.tar.gz
~~~

This command decompresses a directory called `FastQC_Short`.
We enter that directory:

~~~ {.bash}
$ cd FastQC_Short
~~~
 
Revise its content:

~~~ {.bash}
$ ls
~~~
~~~ {.output}
Partial_SRR2467141.fastq 
Partial_SRR2467142.fastq 
Partial_SRR2467143.fastq 
Partial_SRR2467144.fastq 
Partial_SRR2467145.fastq 
Partial_SRR2467146.fastq
Partial_SRR2467147.fastq
Partial_SRR2467148.fastq
Partial_SRR2467149.fastq
Partial_SRR2467150.fastq
Partial_SRR2467151.fastq
~~~

> ## El nombre de los archivos fastq {.callout}
>
> En este caso el nombre de los archivos ha sido asignado de manera
> arbitraria pero muchas veces contienen información importante. 
> Por ejemplo, los archivos de Illumina generalmente tienen el siguiente formato:
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
$ head -n 12 Partial_SRR2467141.fastq 
~~~

~~~ {.output}
@SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
NTTATTTTTTTCGTTCTTCTTGAAGAAGACGTTACCTACGGCGTATTTGCCCATCTCAGGTATGTCTAGATCAAGATCTAACTTGAATTCTCTTTTCATAA
+SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
#1:=BDDDHDHB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
@SRR2467141.2 SALLY:472:C6NCDACXX:8:1101:1326:1950 length=101
CACCCATTGACTGGCCAAATGCCCCATTATTTTGAGTATTGTTATTTCCAAATAAACTGTTACTATTACTGCCAGCGGCAGAAGTGAATCCACAGATCGGA
+SRR2467141.2 SALLY:472:C6NCDACXX:8:1101:1326:1950 length=101
??@DFFFFHHFH?GEHEHIDEHCCFEHIE@GIGCHG<DGGIHIGIGGIGHIIFIHFFGDHGGIIHIGIIIIGGEH@EADB>=AA>3@>CCCCCCBC@C###
@SRR2467141.3 SALLY:472:C6NCDACXX:8:1101:1477:1959 length=101
CATCTTTTCTTTAGGCACATCATCCTGATAAGTGTACTTACCAGGATATATACCATCGGTATTGATGTTATCGGCATCACATAAAACTAATTCACCAGAAA
+SRR2467141.3 SALLY:472:C6NCDACXX:8:1101:1477:1959 length=101
@CCFFFFFHHHHHJJJIIJJIJJJJJJEIJJJIIJJJIJIJJIJJIIIJJJIIIHIIIJDHIJJIGJJIJJJJJIHHHFFFFFEEEEEEDDEDDDDDDDDD
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
@SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
NTTATTTTTTTCGTTCTTCTTGAAGAAGACGTTACCTACGGCGTATTTGCCCATCTCAGGTATGTCTAGATCAAGATCTAACTTGAATTCTCTTTTCATAA
+SRR2467141.1 SALLY:472:C6NCDACXX:8:1101:1392:1873 length=101
#1:=BDDDHDHB?<?CFGGGC9@FF@GGGG>EEEDGDHGFHGE;AEFH>AC@D;@B>C>CCC@C>>DDCC3:>AA5>CC>>CCD>@CCDCDCCCCC@C@>C
~~~

We see that each sequence is represented by 4 lines.

1. **Nombre de la secuencia** - Esta línea comienza con el símbolo `@`, seguido por el nombre de la lectura. 
2. La **secuencia** nucleotídica. 
3. **Segundo título** - Esta línea comienza con el símbolo `+`. Generalmente la información es la misma que en la primera línea pero también puede estar en blanco siempre y cuando contenga el símbolo de `+`.
4. **Información de calidad** - Contiene una cadena de caracteres ASCII. Cada
uno de estos caracteres corresponde a un nucleótido de la secuencia y representa la calidad del mismo. El puntaje de cada nucleótido indica el 
nivel de confianza que se tiene en esa base. Niveles altos indican que la base se ha reportado correctamente mientras que niveles bajos sugieren 
que hay incertidumbre acerca de la secuencia real en esta posición. 

The name of the sequence also contains important information. In
In this case, the Illumina data provide the following information:

~~~ {.output}
@<instrumento>:<corrida>:<identificador de celda de flujo>:<línea>:<cuadro>:<x-coord>:<y-coord> <lectura>:<si fue filtrada>:<número de control>:<índice de la secuencia> length=<tamaño de la secuencia>
~~~

The level of quality represented in the fourth line is called Phred score. In its
original format a Phred score is simply the logarithm of the probabilities of
error:

~~~ {.output}
Phred score = - 10 * log10(error probability)
~~~

| Puntaje de calidad (quality score) | Probabilidad de error |
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
$ fastqc -O ./QUAL/ Partial_SRR2467141.fastq 
~~~

~~~ {.output}
Started analysis of Partial_SRR2467141.fastq
Approx 20% complete for Partial_SRR2467141.fastq
Approx 40% complete for Partial_SRR2467141.fastq
Approx 60% complete for Partial_SRR2467141.fastq
Approx 80% complete for Partial_SRR2467141.fastq
Analysis complete for Partial_SRR2467141.fastq
~~~ 

Once finished, let's review the result:

~~~ {.bash}
$ cd QUAL
$ ls
~~~
~~~ {.output}
Partial_SRR2467141_fastqc.html 
Partial_SRR2467141_fastqc.zip
~~~

The easiest way to explore these results is by opening the html file in
your browser. You can do it by double clicking on the file.

El resultado obtenido deberá ser similar a [este](http://liz-fernandez.github.io/transcriptome_analysis/Partial_SRR2467141_fastqc.html).

This page contains a lot of information broken down into the following sections:

1. **Basic Statistics** - Las estadística básicas del experimento.
2. **Per base sequence quality** - Diagramas de caja mostrando la calidad de cada base.
3. **Per tile sequence quality** - Contiene el cuadro (tile) del que proviene cada secuencia. Solo aparece si los analísis se realizan en una librería de Illumina que retiene sus identificadores originales.
4. **Per sequence quality scores** - Permite ver si un subgrupo de secuencias tiene mala calidad.
5. **Per base sequence content** - Muestra la proporción de cada base en cada posición. 
6. **Per sequence GC content** - Compara la distribución de GC observada con un distribución modelada de GC.
7. **Per base N content** - Indica cuantos nucléotidos no pudieron ser interpretados y se indican con una N.
8. **Sequence Length Distribution** - Muestra la distribución de la longitud de secuencias.
9. **Sequence Duplication Levels** - Indica cuantas veces se repite cada secuencia. Se realiza solo en las primeras 100,000 secuencias.
10. **Overrepresented sequences** - Muestra secuencias sobre representadas.
11. **Adapter Content** - Muestra el contenido de adaptadores en la librería.
12. **Kmer Content** - Muestra la distribución de kmers (subcadenas de tamaño específico) en cada posición de las lecturas. Solo se realiza en el 2% de la librería.

Revisemos los resultados una vez más dando click en el [link](http://liz-fernandez.github.io/transcriptome_analysis/Partial_SRR2467141_fastqc.html). 

En general vemos que la mayoría de los criterios tienen una palomita verde o certificado 
de pase de control de calidad, con la excepción de el contenido de GC por base.

A pesar de que no proporcionamos esta información, el programa ha inferido la codificación 
ASCII (Sanger / Illumina 1.9) así como contado el número de secuencias, su tamaño y 
su contenido de GC. 

Cuando observamos la **calidad por base**, vemos que la mayoría de las bases están en la zona 
verde o tienen buena calidad. Cada base esta representada por un diagrama de caja que 
proporciona una buena idea acerca de la distribución de los puntajes de calidad en todas
las secuencias. Es evidente que la calidad decrece levemente al inicio y más marcadamente
al final de la secuencia. Esto es conocido ya que pequeños errores se acumulan mientras 
más larga sea la secuencia, esta es una de las razones por la que los secuenciadores basados
en incorporación de nucleótidos tienen un límite en la longitud máxima de sus lecturas. 

La gráfica mostrando la **calidad por cuadro** esta completamente azul mostrando que la distribución promedio de errores en el cuadro (tile) es muy similar.

La gráfica de **puntaje de calidad por secuencia** muestra que la mayoría de las secuencias tienen puntajes altos.

Si embargo, la **distribución de nucleótidos por base** muestra que, al principio de la secuencia, la distribución es muy disparatada. Esto es uno de los errores comunes en RNA-Seq hecho por Illumina. Generalmente se soluciona cortando algunos de los nucleótidos en el 5'.

También observamos que: 

* La **distribución de GC por secuencia** es muy similar a la distribución hipotética. 
* Hay pocas **bases ambiguas (Ns)** .
* Todas las secuencias tienen el mismo **tamaño**.  
* Los niveles de **duplicación de secuencias** son bajos.
* No hay **secuencias sobre representadas**.
* Hay pocos **adaptadores**.
* No hay **kmers** sobre representados.

Al ejecutar FastQC también se generó un archivo comprimido llamado 
`Partial_SRR2467141_fastqc.zip`. 
Desempacamos este archivo usando el comando:

~~~ {.bash}
$ unzip Partial_SRR2467141_fastqc.zip
~~~
~~~ {.output}
Archive:  Partial_SRR2467141_fastqc.zip
   creating: Partial_SRR2467141_fastqc/
   creating: Partial_SRR2467141_fastqc/Icons/
   creating: Partial_SRR2467141_fastqc/Images/
  inflating: Partial_SRR2467141_fastqc/Icons/fastqc_icon.png
  inflating: Partial_SRR2467141_fastqc/Icons/warning.png
  inflating: Partial_SRR2467141_fastqc/Icons/error.png
  inflating: Partial_SRR2467141_fastqc/Icons/tick.png
  inflating: Partial_SRR2467141_fastqc/summary.txt
  inflating: Partial_SRR2467141_fastqc/Images/per_base_quality.png
  inflating: Partial_SRR2467141_fastqc/Images/per_tile_quality.png
  inflating: Partial_SRR2467141_fastqc/Images/per_sequence_quality.png
  inflating: Partial_SRR2467141_fastqc/Images/per_base_sequence_content.png
  inflating: Partial_SRR2467141_fastqc/Images/per_sequence_gc_content.png
  inflating: Partial_SRR2467141_fastqc/Images/per_base_n_content.png
  inflating: Partial_SRR2467141_fastqc/Images/sequence_length_distribution.png
  inflating: Partial_SRR2467141_fastqc/Images/duplication_levels.png
  inflating: Partial_SRR2467141_fastqc/Images/adapter_content.png
  inflating: Partial_SRR2467141_fastqc/fastqc_report.html
  inflating: Partial_SRR2467141_fastqc/fastqc_data.txt
  inflating: Partial_SRR2467141_fastqc/fastqc.fo
~~~

Si entramos a el directorio creado, podemos ver el reporte en formato de texto plano, 
tanto el resumen:

~~~ {.bash}
$ cd Partial_SRR2467141_fastqc
$ more summary.txt
~~~
~~~ {.output}
PASS    Basic Statistics        Partial_SRR2467141.fastq
PASS    Per base sequence quality       Partial_SRR2467141.fastq
PASS    Per tile sequence quality       Partial_SRR2467141.fastq
PASS    Per sequence quality scores     Partial_SRR2467141.fastq
FAIL    Per base sequence content       Partial_SRR2467141.fastq
PASS    Per sequence GC content Partial_SRR2467141.fastq
PASS    Per base N content      Partial_SRR2467141.fastq
PASS    Sequence Length Distribution    Partial_SRR2467141.fastq
PASS    Sequence Duplication Levels     Partial_SRR2467141.fastq
PASS    Overrepresented sequences       Partial_SRR2467141.fastq
PASS    Adapter Content Partial_SRR2467141.fastq
PASS    Kmer Content    Partial_SRR2467141.fastq
~~~

Como el archivo completo:

~~~ {.bash}
$ more fastqc_data.txt
~~~

Hemos omitido el contenido ya que es muy largo. 
El poder revisar los resultados en archivos planos es extremadamente útil cuando los 
resultados se encuentren en un servidor que no cuente tiene interface gráfica. 

Los programadores de FastQC han compilado ejemplos de datos de secuenciación de calidad
diversa. Analisemoslos:

* [Buenos datos - Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)
* [Malos datos - Illumina](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)
* [Contaminación por dímeros](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/RNA-Seq_fastqc.html)
* [RNAs pequeños con lecturas a través del adaptador](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/small_rna_fastqc.html)
* [Datos PacBio](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/pacbio_srr075104_fastqc.html)
* [Datos 454](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/454_SRR073599_fastqc.html)

> ## ¿Qué tan válido es nuestro análisis? {.challenge}
>
> Dado que los datos que descargamos solo tienen 5000 secuencias, ¿podemos 
> confiar en que los datos completos tendrán calidad similar? ¿por qué si o por qué no?

## Limpiando las secuencias

Existen un número de herramientas para limpiar secuencias y distintos investigadores
tienen distintas preferencias. Cómo con otros problemas de secuenciación existe un 
desarrollo activo de programas con este fin.

Examples of programs you can use in case your data is not great. 

- Trimmomatic
- 

> ## ¿Qué diferencias hay entre las secuencias crudas y limpias? {.challenge}
>
> Compara los archivos html entre las secuencias crudas y limpias. Resume:
> ¿Qué diferencias existen?
> ¿Cuál crees que sea la causa de esas diferencias?

Existen otras herramientas como `fastx-toolkit`,`scythe` y `sickle` que realizan procesos similares. 
Estás herramientas están en su versión de Biolinux si las quieren comparar.

> ## Tarea - Análisis de calidad {.challenge}
>
> Realiza un análisis de calidad y corte por Trimmomatic de cada uno de los archivos 
> fastq en el directorio que bajaron. 





