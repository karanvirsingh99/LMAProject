 ##################################################
  ## Project: IROC WHOLE SLIDE LMAs
  ## Script purpose: This script determines the optimal number of LMAs detected in the IROC whole slides
  ## Date: June 7th 2023
  ## Author: Karanvir Singh
  ##################################################
 
 # Function that takes in the output of getNeighborhoods (list of cell level data with neighborhood information)
 # and generates a matrix that is ready for clustering
 
 # @param neighborhoods_l is a list 
 
 getNeighborhoodSummary <- function(neighborhoods_l) {
   # Get summary statistics per neighborhood
   
   neighborhoods_l <- bind_rows(neighborhoods_l, .id="Image")
   summary_no_location <- neighborhoods_l %>% filter(!is.na(Neighborhood) & Phenotype!="Negative") %>% group_by(Image, Neighborhood, Phenotype) %>% summarise(sum = n()) %>%
     pivot_wider(names_from = "Phenotype", values_from = "sum")
   
   summary_no_location[is.na(summary_no_location)] <- 0
   
   summary_no_location <- summary_no_location %>% filter(Neighborhood != "Free_cell")
   
   # Convert to matrix
   communities_matrix <- summary_no_location %>% ungroup() %>% select(-Neighborhood, -Image) %>%
     as.matrix()
   
   # Convert to proportions
   communities_matrix_prop <- t(apply(communities_matrix,1, function(x) x/sum(x)))
 }
 
 # Import .RData object containing neighborhood detection results
 
load("/Users/karanvir/Documents/LMAProject/IROC_WHOLE_SLIDES/results/spiat/rdata/IROC_WS_15r_15n_SPIAT.RData")

 # Get matrix for clustering (proportion of cell types in each neighborhood)
 neighborhoods_m <- getNeighborhoodSummary(neighborhoods)
 
 # Run ConsensusClusterPlus to determine optimal number of clusters.
 # Using hierarchical clustering as the method.
 
 distances_bray <- vegan::vegdist(neighborhoods_m, type="bray")
 ConsensusClusterPlus::ConsensusClusterPlus(distances_bray, maxK = 8, clusterAlg="hc")
 
 # 
