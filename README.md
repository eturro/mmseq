# MMSEQ: Transcript and gene level expression analysis using multi-mapping RNA-seq reads

![mmseq collage](https://raw.github.com/eturro/mmseq/master/doc/mmseq-collage.png)

## What is MMSEQ?
The MMSEQ package contains a collection of statistical tools for analysing RNA-seq expression data. Expression levels are inferred for each transcript using the `mmseq` program by modelling mappings of reads or read pairs (fragments) to sets of transcripts. These transcripts can be based on reference, custom or haplotype-specific sequences. The latter allows haplotype-specific analysis, which is useful in studies of allelic imbalance. The posterior distributions of the expression parameters for groups of transcripts belonging to the same gene are aggregated to provide gene-level expression estimates. Other aggregations (e.g. of transcripts sharing the same UTRs) are also possible. Isoform usage (i.e., the proportion of a gene's expression due to each isoform) is also estimated. Uncertainty in expression levels is summarised as the standard deviation of the posterior mean of each expression parameter. When the uncertainty is large in all samples, a collapsing algorithm can be used for grouping transcripts into inferential units with reduced levels of uncertainty.

The package also includes a model-selection algorithm for differential analysis (implemented in `mmdiff`) that takes into account the posterior uncertainty in the expression parameters and can be used to select amongst an arbitrary number of models. The algorithm is regression based and thus it can accomodate complex experimental designs. The model selection algorithm can be applied at the level of transcripts or transcript aggregates such as genes and it can also be applied to detect differential isoform usage by modelling summaries of the posterior distributions of isoform usage proportions as the outcomes of the linear regression models.

## Citing MMSEQ
If you use the MMSEQ package please cite:

- Haplotype and isoform specific expression estimation using multi-mapping RNA-seq reads. Turro E, Su S-Y, Goncalves A, Coin L, Richardson S and Lewin A. Genome Biology, 2011 Feb; 12:R13. doi: [10.1186/gb-2011-12-2-r13](http://dx.doi.org/10.1186/gb-2011-12-2-r13).

If you use the `mmdiff` or `mmcollapse` programs, please also cite:

- Flexible analysis of RNA-seq data using mixed effects models. Turro E, Astle WJ and Tavar&eacute; S. *Submitted*.

## Key features
- Isoform-level expression analysis (works out-of-the-box with Ensembl cDNA and ncRNA files)
- Gene-level expression analysis, which is robust to changes in isoform usage
- Haplotype-specific analysis
- Multi-mapping of reads, including mapping to transcripts from different genes, is properly taken into account
- The insert size distribution is taken into account
- Sequence-specific biases can be taken into account
- Flexible differential analysis based on linear mixed models
- Polytomous model selection (i.e. selecting amongst numerous competing models) is straightforward
- Uncertainty in expression parameters is taken into account
- Collapsing of transcripts with high levels of uncertainty into inferential units which can be estimated with reduced uncertainty
- Multi-threaded C++ implementations

## Estimating expression levels

### Input files:
- FASTQ files containing reads from the experiment (e.g. `asample_1.fq` and `asample_2.fq` paired-end filefs)
- A FASTA file containing transcript sequences to align to ([find ready-made files](#reference-files)) (e.g. `Homo_sapiens.GRCh37.70.ref_transcripts.fa`)

### Step 1: Index the reference transcript sequences

    bowtie-build --offrate 3 Homo_sapiens.GRCh37.70.ref_transcripts.fa Homo_sapiens.GRCh37.70.ref_transcripts 

(It is advisable to use a lower-than-default value for --offrate (such as 2 or 3) as long as the resulting index fits in memory.)

### Step 2a: Trim out adapter sequences if necessary
If the insert size distribution overlaps the read length, trim back the reads to exclude adapter sequences. [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) works well:

<pre><code>trim_galore --phred64 -q 15 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC<span style="color:blue">TAACAAG</span>ATCTCGTATGCCGTCTTCTGCTTG \
  -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -s 3 -e 0.05 \
  --length 36 --trim1 --gzip --paired asample_1.fq.gz asample_2.fq.gz</code></pre>

(TruSeq index highlighted in blue.)

### Step 2b: Align reads with Bowtie 1 (not Bowtie 2)

    bowtie -a --best --strata -S -m 100 -X 500 --chunkmbs 256 -p 8 Homo_sapiens.GRCh37.70.ref_transcripts \
      -1 asample_1.fq -2 asample_2.fq | samtools view -F 0xC -bS - | samtools sort -n - asample.namesorted

- Always specify -a to ensure you get multi-mapping alignments
- Suppress alignments for reads that map to a huge number of transcripts with the -m option (e.g. -m 100)
- Adjust -I and -X according to the minimum and maximum insert size
- If the reference FASTA file doesn't use the vertebrate Ensembl naming convention, then also specify `--fullref`
- The read names must end with /1 or /2, not /3 or /4 (this can be corrected with `awk 'FNR % 4==1 { sub(/\/[34]$/, "/2") } { print }' secondreads.fq > secondreads-new.fq`).
- If there are multiple FASTQ files from the same library, feed them all together to Bowtie in one go (separate the FASTQ file names with commas)
- If you are getting many "Exhausted best-first chunk memory" warnings, try increasing `--chunkmbs` to 128 or 256.
- If the read names contain spaces, make sure the substring up to the first space in each read is unique, as Bowtie strips anything after a space in a read name
- The output BAM file must be sorted by read name. With paired-end data, only pairs where both reads have been aligned are used, so you might as well use the samtools `0xC` filtering flag above to reduce the size of the BAM file

### Step 3: Map reads to transcript sets

    bam2hits Homo_sapiens.GRCh37.70.ref_transcripts.fa asample.namesorted.bam > asample.hits

### Step 4: Obtain expression estimates

    mmseq asample.hits asample

Description of output:

- `asample.mmseq` contains a table with columns:
    1.  feature_id: name of transcript
    2.  log\_mu: **posterior mean of the log\_e expression parameter** (use this as your log expression measure)
    3.  sd: posterior standard deviation of the log\_e expression parameter
    4.  mcse: Monte Carlo standard error
    5.  iact: integrated autocorrelation time
    6.  effective\_length: effective transcript length
    7.  true\_length: length of transcript sequence
    8.  unique\_hits: number of reads uniquely mapping to the transcript
    9.  mean\_proportion: posterior mean isoform/gene proportion
    10.  mean\_probit\_proportion: posterior mean of the probit-transformed isoform/gene proportion
    11.  sd\_probit\_proportion: posterior standard deviation of the probit-transformed isoform/gene proportion
    12.  log\_mu\_em: log-scale transcript-level EM estimate
    13.  observed: whether or not a feature has hits
    14.  ntranscripts: number of isoforms for the gene of that transcript

- asample.identical.mmseq: as above but aggregated over transcripts sharing the same sequence (these estimates are usually far more precise than the corresponding individual estimates in the transcript-level table); log\_mu\_em and proportion summaries are not available for aggregates
- asample.gene.mmseq: as above but aggregated over genes and the effective\_length is an average of isoform effective lengths weighted by their expression
- Various other files (`output_base.*.trace_gibbs.gz`, `output_base.M` and `output_base.k`) containing more detailed output

MMSEQ operates on a sample-by-sample basis, and thus the expression estimates are roughly proportional to the RNA concentrations in each sample. Some scaling of the estimates may be required to make them comparable across biological replicates and conditions. The posterior standard deviations capture the uncertainty due to Poisson counting noise and due to the ambiguity in the mappings between reads and transcripts. The biological variance between samples can only be discerned with the use of biological replicates (see section on differential expression below).

## Differential expression analysis

## Collapsing transcripts

## Reference files
- ***Homo sapiens***: download transcriptome FASTA files containing cDNA and ncRNA transcript sequences (but excluding alternative haplotype/supercontig entries) for the following versions of Ensembl: [64](http://haemgen.haem.cam.ac.uk/eturro/hs_transcripts/Homo_sapiens.GRCh37.64.ref_transcripts.fa.gz), [68](http://haemgen.haem.cam.ac.uk/eturro/hs_transcripts/Homo_sapiens.GRCh37.68.ref_transcripts.fa.gz), [70](http://haemgen.haem.cam.ac.uk/eturro/hs_transcripts/Homo_sapiens.GRCh37.70.ref_transcripts.fa.gz).
- ***Mus musculus***: genome and transcriptome FASTA files based on the GRCm38 build and the [March 2013 SNP and indel calls](ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/) from the [Wellcome Trust Mouse Genomes Project](http://www.sanger.ac.uk/resources/mouse/genomes/) are [available here](http://haemgen.haem.cam.ac.uk/eturro/REL-1303-SNPs_Indels-GRCm38/) for the following strains: C57BL6, 129P2, 129S1, 129S5, AJ, AKRJ, BALBcJ, C3HHeJ, C57BL6NJ, CASTEiJ, CBAJ, DBA2J, FVBNJ, LPJ, NODShiLtJ, NZOHlLtJ, PWKPhJ, SPRETEiJ and WSBEiJ. For F1 data, you should append "\_STRAIN" to each transcript and each gene ID in the transcriptome FASTA headers (where "STRAIN" is the name of the strain) and concatenate the two relevant files into one hybrid FASTA. Then align your F1 reads to the hybrid reference as per the documentation above. For analysis with `mmdiff`, the `*.mmseq` files should be split into two, one file for each strain (use `head` to extract the headers and `grep STRAIN` to extract the rows for a particular strain).

## Installation and dependencies

## Software usage


