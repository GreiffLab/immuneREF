# immuneREF

Overview
========

immuneREF is an R package that enables the analysis of repertoire similarity on a one-to-one, one-to-many and many-to-many scale across repertoire features ranging from fully sequence- to fully frequency-dependent features. This results in a thorough characterization of repertoire datasets for applications ranging from quality control to the search for disease-associated repertoire characteristics.

Documentation: https://immuneREF.readthedocs.io

Publication: tba

![alt text](docs/source/images/immuneREF_Figures-01.jpg?raw=true)

Prerequisites
-------------

To be able to install immuneREF, the following prerequisites need to be fulfilled:

1.  R >= 3.4.0.
2.  Imports: kebabs, igraph, Biostrings, stringdist, vegan, doMC, foreach, dplyr, ggplot2, ggiraphExtra, grid, ComplexHeatmap

```r 
    # Check R version
        version[['version.string']]

    # Install required R packages hosted on CRAN
        install.packages(c("ggplot2","igraph","stringdist","vegan","doParallel","foreach","dplyr","grid"))

    # Install required R packages hosted on Bioconductor 
    
        #If R version ≥ "4.0""
        if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
        BiocManager::install(c("Biostrings","kebabs","ComplexHeatmap"))
    
        #If R version < "4.0"
        source("https://bioconductor.org/biocLite.R")
        biocLite(c("Biostrings","kebabs","ComplexHeatmap"))
```

Installing immuneREF
--------------------

The package can be installed via GitHub.

Installation via GitHub:
1.  Check if all the prerequisites are fulfilled/installed.
2.  Execute the following lines in R:

```r

    # Install the devtools package
    install.packages("devtools")
    
    # Load devtools and install immuneREF from github 
    library(devtools)
    install_github("GreiffLab/immuneREF")
    
    # Test if installation was successful
    library(immuneREF)
```    


Quickstart worfklow
===================

The quickstart workflow in immuneREF_quickstart.R shows the simplest application of 'immuneREF'. In it we run an analysis on a tutorial repertoire dataset consisting of four simulated immune repertoires (included in the package). At the end of the quickstart script, a heatmap visualizing the similarity landscape of the tutorial repertoires is generated. For a more detailed, step-by-step analysis we additionally provide a tutorial R script (immuneREF_tutorial.R)


Performing the analysis
-----------------------

immueREF_quickstart.R, provides a simple example of an immuneREF analysis that includes 3 steps: (i) Extraction of repertoire features, (ii) Calculation of feature-specific similarity between repertoire pairs and (iii) visualization of the results. 

```r
    library(immuneREF)
    
    # Calculate simularity networks for single features
    similarity_networks <- immuneREF_quickstart(repertoire_list = tutorial_repertoires)
    
    # Calculate network features and plot heatmap of repertoire similarities
    network_features <- analyze_similarity_network(similarity_networks[["Condensed"]])
    
    pheatmap::pheatmap(similarity_networks[["Condensed"]],scale='row')

```
