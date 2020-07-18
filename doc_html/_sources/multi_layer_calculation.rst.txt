.. _multi_layer_calculation:

#####################
Multi-layer analysis 
#####################

.. toctree::
   :maxdepth: 2


.. _condense_layers:

Condensing multiple layers
==========================

The immuneREF package provides the ``condense_layer`` function to combine single layers into a multi-layer network. It allows the user to set weights for the different layers using the ``weight`` option, thereby enabling the analysis of repertoire datasets based on problem-specific feature combinations. 
Additionally, immuneREF has a method option, by which method the multi-layer network is calculated. The default method is 'standard', which takes a weighted mean across layers. Other methods, such as majority voting (that is conversion into a boolean network where edges are drawn if they are above a certain weihgt threshold in a majority of layers), will be added in the future. 

.. code-block:: r

	cormat <- condense_layers(list_single_layers = list_single_layers,
    			weights = c(1,1,1,1,1,1),
    			method = "standard")

