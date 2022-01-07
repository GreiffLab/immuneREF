.. _immuneREF_additional_analysis:

############################
Additional Analysis Tools
############################

.. toctree::
   :maxdepth: 2


Gene expression exploratory analysis
=====================================

immuneREF assumes that the input RNA-seq gene expression matrix is already (1) preprocessed with the appropriate methods for correcting possible technical noise and (2) normalized to make samples comparable. Briefly, RNA-seq technology biases include transcript length, GC content, PCR artifacts, uneven transcript read coverage, off-target transcript contamination and differences in transcript distribution. These potential biases can be corrected sequentially, or there are methods that perform in combination. For instance, the R Bioconductor package ``NOISeq`` includes exploratory plots to detect these errors and methods to reduce noise at each step. 


.. _qc_preproc:

1 - Count data quality control and pre-processing
===================================================

Once the gene expression matrix in read count format has been obtained from the raw fasta/fastq files, immuneREF assumes that the user has checked/corrected mainly (a) low count genes, (b) sequencing bias (gene length and GC content) and (c) different count distribution per sample. 

Genes with low counts are generally less reliable and increase data noise, making it harder to extract relevant information. The R package NOISeq incorporates procedures to filter out low counts genes in a given dataset with the function ``filtered.data``. We recommend CPM (method 1) when sample size in each condition is small, previously choosing a CPM threshold from the Sensitivity plot included NOISeq (``explo.plot`` function). If the number of replicates per condition is at least five, the Wilcoxon test (method 2) is more appropriate since it does not need to set any threshold.
NOISeq uses readData function to create a NOISeq object necessary for further functions. The dat function performs different calculations according to the type argument (see `NOISeq manual <https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf>`_ for adequate input format). The plots can be generated for each experimental condition using the argument factor. We provide a code example for low count detection and filtering:

a - Low count filter detection and correction
--------------------------------------------------

Genes with low counts are generally less reliable and increase data noise, making it harder to extract relevant information. The R package NOISeq incorporates procedures to filter out low counts genes in a given dataset with the function ``filtered.data``. We recommend CPM (method 1) when sample size in each condition is small, previously choosing a CPM threshold from the "Sensitivity plot" included NOISeq (``explo.plot`` function). If the number of replicates per condition is at least five, the Wilcoxon test (method 2) is more appropriate since it does not need to set any threshold.
NOISeq uses readData function to create a NOISeq object necessary for further functions. The dat function performs different calculations according to the type argument (see `NOISeq manual <https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf>`_ for adequate input format). The plots can be generated for each experimental condition using the argument factor. We provide a code example for low count detection and filtering:

.. code-block:: r

    library(NOISeq)
    mydata = readData(data = RNAseq, factors = myfactors) 
    mycounts = dat(mydata, factor = "Group", type = "countsbio")
    explo.plot(mycounts, toplot = 1, samples = NULL, plottype = "barplot")
    myfilt = filtered.data(RNAseq, factor = "Group", norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")


b -  Sequencing bias detection and correction
--------------------------------------------------

The length bias and the GC bias plots describe the relationship between the expression values and the feature length and the GC content, respectively. Two vectors,matrices containing transcript length and GC content for each gen have to be provided (see NOISeq manual for adequate format). As before, the plots can be generated for each experimental condition using the argument factor. Code example:

.. code-block:: r

    library(NOISeq)
    mydata = readData(data = myfilt, factors = myfactors, length = mylength, gc = mygc)
    mylengthbias = dat(mydata, factor = "Group", type = "lengthbias")
    explo.plot(mylenthbias, samples = NULL, toplot = "global")
    myGCbias = dat(mydata, factor = factors$group, type = "GCbias")
    explo.plot(myGCbias, samples = NULL, toplot = "global")


If length and GC bias are detected, there exist different normalization methods to correct them. Transcript length bias is corrected by dividing the expression matrix by a numeric vector containing the length of each feature. NOISeq allows internal length correction in the normalization step (point 2). There exist normalization methods external to NOISeq that include length gene and GC content correction, such as `CQN <https://bioconductor.org/packages/release/bioc/vignettes/cqn/inst/doc/cqn.pdf>`_ (conditional quantile normalization) or the R/Bioconductor package `EDASeq <https://www.bioconductor.org/packages/devel/bioc/vignettes/EDASeq/inst/doc/EDASeq.html>`_.



c - Different count distribution per sample
--------------------------------------------------

The count distribution for all samples can be visualized using a boxplot:
``explo.plot(mydata, samples = NULL, plottype = "boxplot")``
The normalization methods homogenize,  to a greater or lesser extent, inter-sample distributions to make samples comparable for subsequent analysis. The following paragraph defines some of the available normalization methods for RNA-seq data.



2 - Normalization
======================

NOISeq performs three types of normalizations (RPKM [1]_, UQUA [2]_ and TMM [3]_) that can be applied as follows:

.. code-block:: r
    myRPKM = rpkm(assayData(mydata)$exprs, long = mylength, k = 0, lc = 1) 
    myUQUA = uqua(assayData(mydata)$exprs, long = mylength, lc = 0.5, k = 0) 
    myTMM = tmm(assayData(mydata)$exprs, long = 1000, lc = 0)

The parameters long, k and lc are used for optional length bias correction (see `NOISeq manual <https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf>`_ for more details).


One of the lightest normalization methods is RPKM (Reads Per Kilobase Million) [1]_, that count the total reads in a sample and divide that number by 106 (“per million” scaling factor), then divide the read counts by the “per million” scale factor, which normalizes for sequencing depth (RPM), and finally divide the RPM values by the length of the gene in kilobases (RPKM).

Two commonly applied normalization methods are UQUA (Upper Quantile) [2]_ and TMM (Trimmed Mean of M values) [3]_, both scaling the samples to have the same upper quartile or median distribution, respectively.

Finally, two normalization methods have been mentioned in the previous section in case sequencing bias were detected and they allow correcting gene length, GC content or both. `CQN <https://bioconductor.org/packages/release/bioc/vignettes/cqn/inst/doc/cqn.pdf>`_ (conditional quantile normalization) applies full-quantile normalization across samples, being one of the more exaggerated normalization methods. The `EDASeq <https://www.bioconductor.org/packages/devel/bioc/vignettes/EDASeq/inst/doc/EDASeq.html>`_ R Bioconductor package includes loess normalization that transforms the data by regressing the counts on a gene feature (GC content, gene length) and subtracting the loess fit from the counts to remove the dependence.

To support quality control of the provided data, the immuneREF function ``print_omics_data()`` generates plots for the exploratory analysis. This allows an initial screening of sample quality to decide whether it is necessary to discard any of the samples (outliers) as well as previous check for differences in gene expression between experimental conditions. The implemented plots were PCA, boxplot and heatmap. All of them are useful for outlier detection (extreme points in the PCA score plot, boxes of the boxplot or branches in the dendrogram will point to them). Boxplots show inter-sample distributions, helping to check if the samples were comparable or the data need to be previously normalized. PCA and heatmap help to check whether gene expression differences exists in the dataset coloring by experimental condition.


.. figure:: /images/analysis_genexp.png


Finally, immuneREF also generates a histogram for the standard deviation (SD) of the gene expression in the provided dataset. This plot helps to select the most appropriate SD threshold for immuneREF. The immuneREF function ``print_sd_histogram()`` ask for the expression matrix as input and generates the following plot:


.. figure:: /images/histogram_sd.png



.. _reference_genexp_layer:

References 
==========

.. [1] Mapping and quantifying mammalian transcriptomes by RNA-Seq, Mortazavi et al., Nature Methods, 2008, https://www.nature.com/articles/nmeth.1226
.. [2] Evaluation of statistical methods for normalization and differential expression in mRNA-Seq experiments, Bullard et al., BMC Bioinformatics, 2010, https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94
.. [3] A scaling normalization method for differential expression analysis of RNA-seq data, Robinson et al., Genome Biology, 2010, https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25#:~:text=TMM%20normalization%20is%20a%20simple,statistical%20methods%20for%20DE%20analysis.

