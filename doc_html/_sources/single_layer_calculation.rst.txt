.. _single_layer_calculation:

###################################
Determining repertoire similarities 
###################################

.. toctree::
   :maxdepth: 2


.. _similarity_features:

Similarity per feature
======================

Having calculated the features of each repertoire (See: :ref:`single_feature_analysis`), the similarity value between each repertoire pair may be calculated. This is done on a per feature basis as laid out in the table below. These measures are the default settings of immuneREF. However, depending on application realm and repertoire characteristics, other methods might be preferred (JS-divergence, Morisita-Horn overlap, etc.) The user is free to deviate from these default settings.


.. _similarity_overview:

Table similarity calculation per feature
========================================

.. csv-table:: Similarity Calculation
   :header: "Layer", "Similarity calculation"
   :widths: 20, 40

   "Diversity", "Pearson correlation, [1-3]." 
   "VDJ Usage", "Pearson correlation coefficient was determined for each frequency vectors (Default V-, J-gene usage vectors) [3,4]."
   "Positional AA Frequency", "Pearson correlation of positional frequency distributions between repertoires."
   "Kmer occurrence", "Pearson correlation of gapped k-mer occurrence (k = 3, m ≤ 3) (Weber et al., 2019)."
   "Network Architecture", "Mean of four measures, (i) Pearson correlated cumulative degree distribution, (ii) absolute mean hub score difference (complement), (iii) complement of absolute difference in fraction of unconnected clusters and nodes and (iv) complement of absolute difference of percent of sequences in the largest connected component."
   "Repertoire overlap", "The clonal sequence overlap measure represents the similarity score between repertoires with respect to clonal convergence [3]."



.. _similarity_calculation:

Similarity calculation for all features
========================================

Having calculated the features :ref:`single_feature_analysis` of  each repertoire, the similarity value between each repertoire pair may be calculated. The easiest way is to use either ``calculate_similarities()`` or ``calculate_similarities_parallel()``. Alternatively, the user may choose to calculate each layer separately using the ``make_cormat()`` or ``make_cormat_parallel()`` functions:each repertoire the similarity score between each repertoire pair can be calculated. 


.. code-block:: r
  		
	# Calculate similarities across layers
	list_single_layers<-list()
		
	# Calculate diversity similarity 
	list_single_layers[["Diversity"]] <- make_cormat(
		repertoires_analyzed= repertoires_analyzed,	# list of analyzed repertoires containing results for all features
		weights_overall = c(1,0,0,0,0,0),			# weighting for each feature. since we calculate per feature set weight to 1 for diversity 
		correlation_method='pearson')				# correlation method is set.
	
	# Calculate AA frequency similarity 
	list_single_layers[["AAfreq"]] <- make_cormat(repertoires_analyzed, weights_overall = c(0,1,0,0,0,0))

	# Calculate Architecture similarity (uses pearson correlation)
	list_single_layers[["Architecture"]] <- make_cormat(repertoires_analyzed, weights_overall = c(0,0,1,0,0,0), correlation_method='pearson')

	# Add repertoire overlap as additional similarity layer
	list_single_layers[["repertoire_overlap"]] <- overlap_layer

  	# Calculate VDJ usage similarity (includes weighting of germlines weights_vdj )
	list_single_layers[["VDJ_usage"]] <- make_cormat(repertoires_analyzed, weights_VDJ=c("V"=1,"D"=0,"J"=1,"VJ"=0), weights_overall = c(0,0,0,0,1,0))
	  
	# Calculate k-mer occurrence similarity
	list_single_layers[["k-mers"]] <- make_cormat(repertoires_analyzed, weights_overall = c(0,0,0,0,0,1))
	  
	## Alternatively to overlap the user may choose to calculate Immunosignatures similarity   
	#list_single_layers[["Immunosignatures"]] <- make_cormat(repertoires_analyzed, weights_overall = c(0,0,0,1,0,0))
  		


.. _reference_chp_single_layer:

References 
==========

[1]  The TCR Repertoire Reconstitution in Multiple Sclerosis: Comparing One-Shot and Continuous Immunosuppressive Therapies, Amoriello et al., Frontiers in immunology, 2020, https://www.frontiersin.org/articles/10.3389/fimmu.2020.00559/full 


[2]  A bioinformatic framework for immune repertoire diversity profiling enables detection of immunological status, Greiff et al., Genome Medicine, 2015, https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-015-0169-8


[3]  Learning the High-Dimensional Immunogenomic Features That Predict Public and Private Antibody Repertoires, Greiff et al., Journal of Immunology, 99(8), 2017, http://www.jimmunol.org/content/199/8/2985


[4]  immuneSIM: tunable multi-feature simulation of B- and T-cell receptor repertoires for immunoinformatics benchmarking, Weber et al., Bioinformatics, 2020, https://academic.oup.com/bioinformatics/article/36/11/3594/5802461


