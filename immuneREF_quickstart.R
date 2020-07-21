library(immuneREF)

# Calculate simularity networks for single features (Step 1 and 2)
similarity_networks <- immuneREF_quickstart(repertoire_list = tutorial_repertoires)

# Calculate network features and plot heatmap of repertoire similarities (condensed network) (Step 4)
network_features <- analyze_similarity_network(similarity_networks[["Condensed"]])

pheatmap::pheatmap(similarity_networks[["Condensed"]],scale='row')