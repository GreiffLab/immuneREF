library(immuneREF)

# Calculate simularity networks for single features (Step 1 and 2)
similarity_networks <- immuneREF_quickstart(repertoire_list = tutorial_repertoires)

# Calculate network features and plot heatmap of repertoire similarities(Step 4)
network_features <- analyze_similarity_network(multilayer_network)

pheatmap::pheatmap(multilayer_network,scale='row')