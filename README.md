# DTP_Transcriptomics_Oxford
DTP DTP course Bioinformatics - Transcriptomics


# Transcriptome Practical v.2025.10

This practical provides a complete, hands-on workflow for RNA-seq
transcriptome analysis, taking you from raw sequencing data all the way
to biological interpretation.

The dataset used in this exercise originates from a real-world research
project on the legume Chamaecrista mimosoides. For the purposes of this
practical, we focus on a single scaffold, which represents roughly 10%
of the full genome — sufficient to explore key transcriptomic principles
without the computational overhead of whole-genome analysis.

In the underlying experiment, Chamaecrista roots were grown under two
conditions:

-   Mock-treated roots (control)

-   Bacteria-treated roots, where symbiotic bacteria were applied to
    induce early infection responses

This setup mimics the initial stages of a symbiotic interaction between
legumes and nitrogen-fixing bacteria. Such interactions can lead to
nodulation — a remarkable biological process in which the plant forms
specialized root structures (nodules) that host nitrogen-fixing
bacteria. Within these nodules, the microbes convert atmospheric
nitrogen into a form the plant can use, in exchange for carbon
compounds.

Because nodulation involves complex signaling, infection, and organ
development pathways, it triggers widespread transcriptional changes —
making this system an excellent case study for differential gene
expression analysis.

Note: since this dataset represents only a subset of the full genome, we
naturally expect to detect a smaller number of genes and perhaps some
minor artifacts introduced during filtering or subsampling.
Nevertheless, it is based on genuine experimental data and provides a
realistic yet lightweight framework for learning transcriptome analysis
methods.

To complement the short-read Illumina RNA-seq data used here, the
original research project also employed long-read Nanopore sequencing.
By combining these two technologies, researchers can capture both
high-accuracy quantification (from Illumina) and full-length transcript
structures (from Nanopore). This hybrid approach enables more accurate
gene model curation and a deeper understanding of transcript diversity.

Follow the steps carefully, execute the commands, and reflect on the
interpretation questions along the way.

**Before You Begin**

This exercise assumes you have:

• A working WSL (Windows Subsystem for Linux) setup or similair

• Conda installed and a basic familiarity with using it

We’ll start by creating a conda environment and installing the required
tools (e.g. hisat2, minimap2, samtools, kallisto, etc.). You can name
the environment whatever you like — for example, rna\_practical.

## 1. Data and Directory Setup

For this practical, a dummy dataset has been prepared for you.

After downloading it, extract the ZIP file — this will set up all the
data needed for the exercises in this session.

Inside, you’ll find a subset of the Chamaecrista mimosoides genome
focused on a single scaffold:

• Genome sequence: Chammim\_Cha06\_scaffold\_6.fasta

• Annotation file: Chammim\_Cha06\_scaffold\_6.gtf

• (You may also see a file named Chammim\_Cha06\_scaffold\_6.fasta.fai —
this is just an index file used by mapping tools. You don’t need to load
it in IGV.)

Let’s start by getting familiar with the data.

Open IGV (Integrated Genomics Viewer) and load the genome and annotation
files:

1\. Load the genome:

• Go to Genomes → Load Genome from File…

Select Chammim\_Cha06\_scaffold\_6.fasta

2\. Load the annotation (and later the BAM or other tracks):

• Go to File → Load from File…

Select Chammim\_Cha06\_scaffold\_6.gtf

Once loaded, you should see the annotated genes displayed along the
scaffold.

If you zoom in far enough, you’ll even be able to view the individual
base pairs of the sequence.

#### Experimental Design

The toy dataset represents a small-scale version of a real RNA-seq
experiment investigating root nodulation in Chamaecrista.

It includes two treatments at two time points:

• Mock: untreated control roots

• Bacteria: roots inoculated with nitrogen-fixing bacteria

Each folder in the dataset corresponds to one biological replicate and
treatment condition.

Samples are labeled using this convention:

💡 Tip: make sure that you are in the working directory where your data
is stored you can go there by typing cd /path/to/where/you/want/tobe

Your directory structure for the short read RNA-seq data should rougly
follow this pattern:

Each folder in the dataset corresponds to one biological replicate and
treatment condition.

Samples are labeled using this convention:

-   Cxx = sample number

-   2dpi = days post infection

-   Mock = untreated roots or Bacteria = inoculated roots

💡 Tip: If any step fails or not complete on time or you encounter issues
running the commands, you can download the Results.zip file, which
contains all the output files generated throughout the practical. This
allows you to continue with the following steps without having to re-run
the failed analyses.

#### Quality Control of Raw Reads

Before we start mapping, let’s first assess the quality of the
sequencing reads.

Quality control helps identify problems such as low-quality bases,
adapter contamination, or uneven GC content, which can affect downstream
analyses.

We’ll use:

• FastQC — to assess the quality of raw FASTQ files

• MultiQC — to summarize all FastQC reports into a single overview

#### Run FastQC

To save time, we’ll only run FastQC on two representative sampels (R1
and/or R2) — for example, one from the Mock condition and one from the
Bacteria condition.

Navigate to the folder containing the FASTQ files and run:

This command will generate .html and .zip report files for each input
FASTQ.

You can open the .html files in your browser to inspect the results.

#### Summarize Reports with MultiQC

Even though we only ran FastQC on two. (or more) files, we can still use
MultiQC to compile the results into one easy-to-read summary:

**Questions:**

-   What is the average per-base sequence quality across samples?  
    Are there any adapter content issues?

-   If this would be your dataset what would you have trimmed or
    removed?

## 3. Mapping Reads with HISAT2 & Minimap2

In this part, we’ll align the RNA-seq reads from **Mock** and
**Bacteria-treated (Bac)** samples to the reference contig:
Chammim\_Cha06\_scaffold\_6.fasta.

Because we have several biological replicates for each treatment and
timepoint, we’ll first **combine reads** for each group. This reduces
the number of alignment files while still letting us compare overall
expression patterns between treatments. Here is an example of how to
combine the reads for **Mock 2 dpi** (day post infection). After doing
this successfully, please continue and generate the merged files for Bac
2 dpi.

Once complete, you should have **four merged FASTQ files**, two per
condition (R1/R2).

#### Build the HISAT2 index and preform the mapping

HISAT2 requires an index of the reference genome before mapping. This
index allows fast, memory-efficient alignment.

💡 Tip: The prefix (Chammim\_Cha06\_scaffold\_6\_index) is how HISAT2
will later find these files.

You’ll use the same prefix when mapping.

Now we’ll map each treatment group to the indexed reference. HISAT2 will
generate a **SAM** file (Sequence Alignment Map). To make the SAM files
smaller and easier to view in **IGV**, convert them to **sorted BAM**
files and create index files. You can use samtools flagstat to quickly
check mapping statistics:

**Questions**:

1.  What is the overall alignment rate?

2.  Compare across samples. Are the rates roughly similar, or do some
    align much better?

3.  (Remember: we are mapping only to one scaffold, not the entire
    genome.)

4.  Are there noticeable differences between Mock and Bac samples?

5.  Does one treatment show more or fewer mapped reads? What might that
    suggest about gene activity?

#### Long-read Nanopore reads mapping

Long Oxford Nanopore (ONT) cDNA reads are much longer than short
Illumina reads and often span entire transcripts. Because they cover
full-length isoforms, they are very useful to determine accurate splice
junctions and to identify alternative splice variants. For cDNA
(reverse-transcribed RNA) libraries, we use a spliced aligner tuned for
noisy long reads.

Map Nanopore cDNA reads with minimap2 (recommended preset for ONT cDNA):

Why these options?

-   -ax splice — spliced alignment mode (finds introns / exon
    junctions).

-   -uf -k14 — recommended for noisy long reads (u = long-reads, f = use
    forward-strand heuristics for cDNA; k14 lowers minimizer k to
    improve sensitivity).

-   -t 4 — use 4 threads (increase for faster runs if you have more
    CPUs).

Adjust -k/-uf if your reads are higher quality (e.g., use default -x
splice or -ax splice:hq for very high-quality (Q&gt;Q) reads).

1.  How many Nanopore reads mapped to scaffold 6? What is the overall
    mapping rate? (use flagstat)

2.  Pick one gene with multiple splice isoforms in IGV: how many
    distinct full-length isoforms are supported by single Nanopore
    reads?

3.  Are there splice junctions detected by minimap2 in the long reads
    that are not supported by short-read junctions? (possible novel
    isoforms)

4.  If you see many inconsistent or low-quality mappings, which minimap2
    options might you change to improve mapping sensitivity or
    specificity? (hint: presets like -ax splice:hq or changing -k).

## 4. Transcript Quantification (Kallisto)

Now that we’ve explored read mapping and visualization, let’s move on to
quantifying transcript abundance. We’ll use Kallisto, a fast and
accurate tool that performs pseudoalignment — instead of aligning each
read base-by-base, it finds which transcripts the reads could have come
from, and estimates transcript counts and expression levels.

Kallisto needs an index of the transcriptome FASTA file, which contains
the nucleotide sequences of all transcripts. This only needs to be built
once.

#### Run quantification for each sample

Each quantification run will take the paired-end FASTQ files for a
sample and output estimated transcript counts and TPM values
(Transcripts Per Million).

Here’s an example for C1\_2dpi\_Bac.

After you’ve successfully run this example, repeat for the other
samples.

## 5. Differential Expression Analysis (DESeq2) and visualization in R

This exercise is a bit more challenging, and can be overwhelming with
all the code involved. But let see how far we get. We will use R and
DESeq2 to test which genes change expression between the Bacteria (Bac)
treatment and Mock control. Kallisto gives us transcript-level estimates
(abundance.h5). We import those with tximport to obtain estimated counts
suitable for DESeq2.

From here we continue in Rstudio, open it and load the
‘analysis\_deseq2\_kallisto.R’ script.

NOTE: In this dataset C5 and C6 were harvested at a different time than
C1–C4; by giving C5 and C6 the same pair ID we allow DESeq2 to
explicitly account for that pairing and remove some unwanted
variability.

There are two analysis models that are helpful here:

-   Simple model: \~ condition — tests Bac vs Mock without accounting
    for other factors.

-   Paired model: \~ pair + condition — includes a pair (blocking)
    factor so DESeq2 can control for paired measurements or
    harvest/batch differences.

1.  Run the provided R script (or RStudio notebook) that imports
    Kallisto outputs and runs the two DESeq2 models.

2.  What is pseudo-mapping and how does it differ from normal mapping of
    reads?

3.  Inspect the PCA plot:

    -   Do samples separate by condition (Bac vs Mock) or by harvest
        (early vs late)?

    -   What does the observed clustering tell you about potential batch
        or harvest effects?

4.  Compare the number of significantly up-regulated and down-regulated
    genes between the simple and paired analyses (use thresholds: padj
    &lt; 0.05 and \|log2FC\| &gt; 1):

    -   Does pairing increase or decrease the number of significant
        genes?

    -   Why might that happen (hint: think about variance explained by
        pair)?

    -   Decide which model you trust more to call genes affected by
        treatment, given that C5 and C6 were harvested at a different
        time. Explain briefly in one sentence.
