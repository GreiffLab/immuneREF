#
# 01 Main script to preprocess repertoires, calculate characteristics and calculate correlation matrix
#
library(immuneREF)
library(foreach)
library(doMC)
registerDoMC(1)

# Set working directory
setwd("/Users/cedricweber/Desktop/tutorial_immuneREF")

# Load datasets (check for compatibility, 
# Compatible repertoires have: 
#       1) AIRR standard column naming. 
#       2) V(D)J naming compatible with included list_germline_genes dataset.
#       3) Clonal frequencies that sum up to 1.
list_simulated_repertoires <- tutorial_repertoires


# Set number of repertoires to be analyzed and names for reference dataframe
repertoire_names <- names(list_simulated_repertoires)
repertoire_lengths <- sapply(list_simulated_repertoires,nrow)
repertoire_species <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[1]])
repertoire_receptor <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[2]])
repertoire_chain <- sapply(strsplit(repertoire_names, "\\_"),function(x) x[[3]])

input_data_ref<-data.frame(sample_id = repertoire_names,
                           nb_sequences = repertoire_lengths,
                           species = repertoire_species,
                           receptor = repertoire_receptor,
                           chain = repertoire_chain,row.names=c(1:length(repertoire_names)))



# Calculate overlap layer on full datasets for Convergence layer
overlap_layer<-repertoire_overlap(list_simulated_repertoires,basis="CDR3_aa")



# Subsample for the ones that are not 10000
subsample_size<-10000

list_simulated_repertoires<-subset_input_repertoires(list_repertoires=list_simulated_repertoires,
                               subset_size=subsample_size,
                               random=FALSE)


# Calculate all features for each repertoire
repertoires_analyzed<-list()
for(i in 1:length(list_simulated_repertoires)){
  repertoires_analyzed[[repertoire_names[i]]]<-calc_characteristics(
                                repertoire_df=list_simulated_repertoires[[i]],
                                species=repertoire_species[i],
                                receptor=repertoire_receptor[i],
                                chain=repertoire_chain[i],
                                identifier_rep=repertoire_names[i])
}


# Calculate similarities between repertoire for each layer
list_single_layers<-calculate_similarities(
  repertoires_analyzed=repertoires_analyzed,
  overlap_layer=overlap_layer)




# Calculate condensed network (here equal weights for each layer)
cormat <- condense_layers(list_single_layers,
    weights = c(1,1,1,1,1,1),
    method = "standard")


###
# Draw heatmap of immuneREF layers
###

# Make list of all layers you want to plot heatmaps for
list_all_layers <- list_single_layers
list_all_layers[["Condensed"]] <- cormat


# Prepare list with heatmap annotations containing categories and colors
annotation_list<-list()
annotation_list[["categories"]]<-data.frame(Species=input_data_ref$species,
                                            Receptor = input_data_ref$receptor)

annotation_list[["colors"]]<-list(Species=c(mm='#ffffbf',hs='#fc8d59'),
                                  Receptor=c(ig='#91bfdb'))


# For each entry (immuneREF layer) plot a heatmap
# Note: Folder for which path is set in path_figure needs to exist already
print_heatmap_sims(list_similarity_matrices=list_all_layers,
                  annotation_list=annotation_list,
                  path_figure="figures")




###
# Draw Global Similarity Plots of immuneREF layers
###

# Define relevant subsets for splits
categories_list<-list()
categories_list[["categories"]]<-input_data_ref
categories_list[["color"]]<-c("white",'#91bfdb','#ffffbf')
categories_list[["subset"]]<-"species"

#Plot global similarity 
print_global_similarity(list_similarity_matrices=list_all_layers,
                      categories_list = categories_list,
                      path_figure="figures")




###
# Plot local similarity per category and identify max and min locally similar repertoires
##
max_min_reps<-print_local_similarity(list_similarity_matrices=list_all_layers,
                      categories_list = categories_list,
                      path_figure="figures")


###
# Radar plot to visualize similarity across all 6 layers
##

radar_list<-list()
radar_list[["category"]]<-c("Murine A","Murine B","Human A","Human B")
radar_list[["roi"]]<-c("mm_ig_h_2_0__0_0_0_A","mm_ig_h_4_0__0_0_0_A","hs_ig_h_2_0__0_0_0_A","hs_ig_h_4_0__0_0_0_A")
radar_list[["label"]]<-c("Murine A","Murine B","Human A","Human B")
radar_list[["colors"]]<-c("grey","blue",'red',"green")

comparison_list<-list(roi=radar_list[["roi"]],
  roi_names=radar_list[["category"]],
  ref="mm_ig_h_2_0__0_0_0_A",
  plot_names=radar_list[["label"]],
  colors=radar_list[["colors"]])

print_repertoire_radar(list_similarity_matrices=list_single_layers,
  to_compare=comparison_list,
  path_figure="figures",
  name_plot="tutorial")



####
# Classical repertoire analysis of maximally and minimally similar repertoires per category
####

# Print classic repertoires comparing max and min locally similar plots for:
# Simulated murine igh repertoires
mm_igh<-list()
mm_igh[["max"]]<-repertoires_analyzed[["mm_ig_h_2_0__0_0_0_A"]]
mm_igh[["min"]]<-repertoires_analyzed[["mm_ig_h_4_0__0_0_0_A"]]

print_repertoire_comparison(list_repertoires=mm_igh,name_plots="mm_igh",aa_freq_length=14,path_figure="figures")


# Simulated human igh repertoires
hs_igh<-list()
hs_igh[["max"]]<-repertoires_analyzed[["hs_ig_h_2_0__0_0_0_A"]]
hs_igh[["min"]]<-repertoires_analyzed[["hs_ig_h_4_0__0_0_0_A"]]

print_repertoire_comparison(list_repertoires=hs_igh,
  name_plots="hs_igh",
  aa_freq_length=17,
  path_figure="figures")


## Bonus: calculate network features of condensed immuneREF layer
network_features <- analyze_similarity_network(cormat)



