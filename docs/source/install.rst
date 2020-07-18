.. _install:

####################
Installing immuneREF
####################

Prerequisites
=============

To be able to run the code, the following prerequisites should be fulfilled:

1. R >= 3.4.0.
2. The following packages need to be installed:

    * ggplot2 
    * igraph
    * Biostrings
    * stringdist 
    * vegan
    * doMC
    * foreach
    * dplyr
    * grid
    * kebabs
    * ComplexHeatmap


.. code-block:: RST
 
    # Check R version
        version[['version.string']]

    # Install required R packages hosted on CRAN
        install.packages(c("ggplot2","igraph","stringdist","vegan","doMC","foreach","dplyr","grid"))

    # Install required R packages hosted on Bioconductor 
    
        #If R version ≥ "4.0""
        if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
        BiocManager::install(c("Biostrings","kebabs","ComplexHeatmap"))
    
        #If R version < "4.0"
        source("https://bioconductor.org/biocLite.R")
        biocLite(c("Biostrings","kebabs","ComplexHeatmap"))



Install immuneREF
=================

The package can be installed in R (via GitHub):

1.  Check if all the prerequisites are fulfilled/installed.
2.  Execute the following lines in R:

.. code-block:: RST

    # Install the devtools package
    install.packages("devtools")
    
    # Load devtools and install immuneREF from GitHub 
    library(devtools)
    install_github("GreiffLab/immuneREF")

    # Test if installation was successful by loading immuneREF
    library(immuneREF)

