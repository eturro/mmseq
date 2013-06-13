#### Feeding MMSEQ estimated counts to DESeq or edgeR

The MMSEQ expression estimates are roughly in FPKM units (fragments per kilobase of transcript per million mapped reads or read pairs), which makes different samples broadly comparable. You may read in the output from multiple samples using the `mmseq.R` R script included in the `misc` directory (requires R &ge; 2.9). The `readmmseq` function defined in `mmseq.R` takes five optional arguments:

-  `mmseq_files`: vector of `.mmseq` files to read in
-  `sample_names`: vector of strings containing the sample names
-  `normalize`: logical scalar specifying whether to normalise the log expression values and log count equivalents for each sample by their median deviation from the mean (like DESeq)
-  `partition`: a string containing two suffixes corresponding to two different classes of haplotypes separated by a "|" character (e.g. "\_BL6|\_CAST"); each haplotype class is then treated as a separate sample
-  `common_set`: vector specifying the subset of features to use

The function returns a list in which columns have been collated from the MMSEQ files. In addition, the `counts` slot contains estimated counts for each feature.

To test for differential expression with [edgeR](http://dx.doi.org/10.1186/gb-2010-11-3-r25) or [DESeq](http://dx.doi.org/10.1186/gb-2010-11-10-r106) instead of using [mmdiff](https://github.com/eturro/mmseq#differential-expression-analysis), the estimated counts need to be used. E.g., to test for DE between two groups of two samples, run the following code in R from the directory containing the mmseq output files:

    source("/path/to/mmseq.R")
    library(edgeR)
    library(DESeq) # version >=1.5

    ms = readmmseq()
    tab = ms$counts[apply(ms$counts,1,function(x) any(round(x)>0)),] # only keep features with counts

    # edgeR
    d = DGEList(counts = tab, group = c("group1", "group1", "group2","group2"))
    d = estimateCommonDisp(d)
    de.edgeR = exactTest(d)

    # DESeq
    cds = newCountDataSet(round(tab), c("group1", "group1", "group2", "group2") )
    sizeFactors(cds) = rep(1, dim(cds)[2])
    cds = estimateDispersions(cds)
    de.DESeq = nbinomTest(cds,"group1", "group2")
