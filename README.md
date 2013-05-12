# MMSEQ: Transcript and gene level expression analysis using multi-mapping RNA-seq reads

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
