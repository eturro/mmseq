## Authors
Ernest Turro wrote all the code in this package except

  - hitsio library: J&uuml;rgen J&auml;nes
  - `readmmseq.R`: Angela Goncalves and Ernest Turro
  - `ensembl_gtf_to_gff.pl`: user 'genec' on seqanswers.com
  - `sokal.cc`/`hh`: Graeme Ambler
  - `fasta.c`/`h`: Bio5495/BME537 Computational Molecular Biology, Michigan State University

## Changelog
  - 2013-06-11: MMSEQ 1.0.4 released, fixed bugs in mmdiff apparent when M.n\_cols > 1
  - 2013-06-09: MMSEQ 1.0.3 released, fixed beta update (covariates common to all models can now be taken into account); added -de convenience option for mmdiff; also fixed harmless bug that could strip first character of sample names in mmdiff output header
  - 2013-04-26: MMSEQ 1.0.2 released, better prior for eta in `mmdiff` (Inverse-Gamma)
  - 2013-04-08: MMSEQ 1.0.1 released, including new binaries `mmcollapse` and `mmdiff` (publication forthcoming)
  - 2012-12-06: MMSEQ 1.0.0 beta released (more `mmseq` output: estimates for features even if they have no reads, output posterior standard deviations, whether or not features have been observed, etc. Improved normalisation procedure in `readmmseq.R` and `mmdiff` (using only features with unique hits in a minimum proportion of samples). Beta versions of `mmcollapse` and `mmdiff` programs. Inclusion of `t2g_hits` program.)
  - 2012-05-28: MMSEQ 0.11.2 released (output mean posterior isoform/gene expression proportions)
  - 2012-04-19: MMSEQ 0.11.1 released (output true sequence lengths in addition to effective transcript lengths)
  - 2012-04-16: MMSEQ 0.11.0 released (`bam2hits` updates: print nice ASCII insert size histogram to check it's roughly normal; `-i` option now uses `expected_isize` and `sd_of_isizes` (instead of `deviation_thres`); unless adjusted lengths are provided or calculated with `mseq` (options `-c` and `-u`), transcript lengths are adjusted using the insert size distribution - the prior fragment distribution given the transcript is assumed to be uniform up to the transcript length within the main support of the insert size distribution; proper repetitive transcript filter applied before simplified isize deviation filter; deal gracefully with truncated hits files)
  - 2012-03-26: MMSEQ 0.10.0 released (new hitsio library by J&uuml;rgen J&auml;nes results in ~7x reduction in hits file sizes - inspection of compressed binary hits files and conversion to plain text can be done with the `hitstools` utility; fixed default regex to work with new (&gt;v64) Ensembl fasta entry naming convention; more information in `*.mmseq` tables; output gzipped MCMC traces; improved `readmmseq.R` script)
  - 2011-12-16: MMSEQ 0.9.18 released (the `*.mmseq` files now contain all the transcripts or genes, regardless of whether they have counts; improved R script to read in the output (`readmmseq.R`); fixed potential integer overflow if `gibbs_iter` set too high; default `gibbs_iter` doubled to 2^14, though there can be a noticeable reduction in autocorrelation up to 2^17 iterations)
  - 2011-11-07: MMSEQ 0.9.17 released (misc. improvements and small bugfix in `bam2hits` code for identifying read mates)
  - 2011-10-21: MMSEQ 0.9.16 released (output mean(log) rather than log(mean) of expression parameter; output expression-weighted sum of transcript lengths in gene files; misc other changes)
  - 2011-10-17: MMSEQ 0.9.15 released (automatically set insert size deviation filter parameters; remove read length - 1 from effective transcript lengths; added `extract_transcripts` program to produce sequences from genome FASTA and Ensembl-like GTF files)
  - 2011-09-22: MMSEQ 0.9.14 released (workaround for bug in Mac `gcc` that causes `omp_get_max_threads()` to incorrectly return 1 if called from a parallel region - estimates obtained on Mac with previous versions of `mmseq` are probably invalid)
  - 2011-09-20: MMSEQ 0.9.13 released (changed default prior on mu to Gamma(0.1,0.1), which appears to reduce upwards bias for lowly expressed transcripts; all mmseq parameters can now be controlled from the command line; include `mmseq2counts.R` for DE analysis)
  - 2011-08-15: MMSEQ 0.9.12 released (bug fix for potential seg fault)
  - 2011-06-12: MMSEQ now available on Fedora 14/15 and EPEL 5/6 (`yum install mmseq`). Many thanks to Adam Huffman.
  - 2011-06-07: new version of polyHap2 released with support for random assignment of alleles to haplotypes
  - 2011-05-27: MMSEQ 0.9.11 released (fixed bug in insert size deviation filter in `bam2hits`)
  - 2011-05-04: MMSEQ 0.9.10c released (ensured `haploref.rb` works properly with custom regular expressions - same behaviour as `testregexp.rb` and `bam2hits`; fixed bug in identical transcript/gene amalgamation code. Prior to 0.9.10, incorrect values outputted in `*.amalgamated.mmseq` and `*.gene.mmseq` files if one of the transcripts had zero hits and no values outputted in `*.gene.mmseq` for genes with log Gibbs mean &lt; 0)
  - 2011-03-18: MMSEQ 0.9.9 released (better integration with 'mseq'; deal gracefully with Ns in counts file; harmless buffer overflow bugfix in bam2hits)
  - 2011-02-10: <a href="http://dx.doi.org/10.1186/gb-2011-12-2-r13">MMSEQ paper</a> published in Genome Biology
  - 2010-10-19: MMSEQ 0.9.8 released
