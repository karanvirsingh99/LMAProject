# Check that the cell level data dataframe has these columns: Class, Centroid_x, Centroid_y, Image
# @param df a dataframe containing cell level data as exported from qupath

checkCellLevelData <- function(df){
  sum(c("Class", "Centroid_x", "Centroid_y", "Image") %in% colnames(df)) == 4
}


# This function converts a df with cell level data into a list where each entry is a spatial object
# @param cell_level_data a data.frame with columns Image, Class, Parent, Centroid_x, Centroid_y

images_as_list <- function(cell_level_data){
  imgDataList <- list() #empty list
  image_names <- unique(cell_level_data$Image) #get all image names
  
  # Trying to use data.table as it's faster..
  cell_level_dt <- data.table(cell_level_data)
  
  for (img in image_names) { 
    # Get all cells corresponding to image i and add them to imgDataList
    imgData <- cell_level_dt[Image == img, .(Image, Class, Parent, Centroid_x, Centroid_y)]
    # Add cell ID column
    imgData <- imgData[, ID := paste0("Cell_", .I)]
    
    # Extract spatial info
    coord_x = imgData$Centroid_x
    coord_y = imgData$Centroid_y
    phenotypes = imgData$Class
    dummy_intensity = rep(0, nrow(imgData))
    intensity_matrix = matrix(dummy_intensity, nrow=1, ncol=nrow(imgData))
    colnames(intensity_matrix) = imgData$ID
    
    # Use SPIAT function to create a spatial object
    imgDataList[[img]] = SPIAT::format_image_to_spe(intensity_matrix=intensity_matrix,
                                                    phenotypes = phenotypes, coord_x = coord_x,coord_y = coord_y)
  }
  return(imgDataList)
}

# Function which detects neighborhoods in a single image
# @param spe_object a spatial object, formatted as per SPIAT instructions
# @param radius the pairwise distance between cells under which they are considered to be interacting
# @neighborhood_size the minimum number of cells in a neighborhood to be detected
# @cell_types the types of cells to consider for neighborhood detection

detectNeighborhoods <- function(spe_object, radius, neighborhood_size, cell_types){
  
  error_catch <- tryCatch({
    identify_neighborhoods(spe_object = spe_object,
                           cell_types_of_interest = cell_types,
                           method="hierarchical",
                           radius=radius,
                           min_neighborhood_size = neighborhood_size,
                           feature_colname = "Phenotype")
    1
  }, error=function(e) 2)
  
  # error_catch == 1 when no error 
  if(error_catch == 1){
    
    communities <- identify_neighborhoods(spe_object = spe_object,
                                          cell_types_of_interest = cell_types,
                                          method="hierarchical",
                                          radius=radius,
                                          min_neighborhood_size = neighborhood_size,
                                          feature_colname = "Phenotype")
    
    communities_as_df <- colData(communities) %>% as.data.frame()
    communities_as_df <- cbind(communities_as_df, spatialCoords(communities))
    return(communities_as_df)
    
  } else(
    return(NULL) # return empty dataframe if no error
  )
} 


# Function that takes in cell level data and returns SPIAT output (cell level data with neighborhoods)
# @param cell_level_data a data.frame with cell level output from qupath

getNeighborhoods <- function(cell_level_data, radius, n_size, cell_types) {
  if(checkCellLevelData(cell_level_data)){
    cell_level_list <- images_as_list(cell_level_data)
    all_neighborhoods <- lapply(cell_level_list, detectNeighborhoods,
                                radius = radius,
                                neighborhood_size = n_size,
                                cell_types = cell_types)
  } else {
    print("Dataframe doesn't have Class, Image, Centroid_x, Centroid_y columns")
  }
}

# Function that takes in the output of getNeighborhoods (list of cell level data with neighborhood information)
# and clusters the data

# @param neighborhoods is a data.frame of cell level data

clusterNeighborhoods <- function(neighborhoods, k){
  
  # Get summary statistics per neighborhood
  
    summary_no_location <- neighborhoods %>% filter(!is.na(Neighborhood) & Phenotype!="Negative") %>% group_by(Image, Neighborhood, Phenotype) %>% summarise(sum = n()) %>%
      pivot_wider(names_from = "Phenotype", values_from = "sum")
    
    summary_no_location[is.na(summary_no_location)] <- 0
    
    summary_no_location <- summary_no_location %>% filter(Neighborhood != "Free_cell")
  
  # Convert to matrix
  communities_matrix <- summary_no_location %>% ungroup() %>% select(-Neighborhood, -Image) %>%
    as.matrix()
  
  # Convert to proportions
  communities_matrix_prop <- t(apply(communities_matrix,1, function(x) x/sum(x)))
  
  # Calculate dissimilarity and build dendrogram
  hier_bray <- vegan::vegdist(communities_matrix_prop, method = "bray")
  tree <- hclust(hier_bray, method="ward.D2")
  tree_dend <- as.dendrogram(tree)
  
  # Cut tree and cluster
  clust <- dendextend::cutree(tree_dend, k=k)
  neighborhood_df_with_LMA <- summary_no_location %>% ungroup() %>% mutate(LMA = clust)
}

# This function produces summary plots for clustered neighborhoods

getClusterPlots <- function(LMAtypes, cell_type_colors, breaks){
  # Generate plots ana summary statistics per LMA type
  
  # Plot Type 1: Barplots
  barPlotData <- LMAtypes %>% 
    mutate(COMMUNITY_ID = row_number()) %>%
    select(-c("Image", "Neighborhood")) %>% 
    pivot_longer(cols = -c("COMMUNITY_ID","LMA")) %>%
    group_by(COMMUNITY_ID) %>%
    mutate(prop = value/sum(value))
  
  p1 <- barPlotData %>%
    ggplot(aes(x=factor(COMMUNITY_ID), y=prop, fill=name))+
    geom_col()+
    facet_wrap(~LMA, scales="free")+
    scale_fill_manual(name="Cell Type",
                      breaks=breaks,
                      values = cell_type_colors)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    xlab(NULL)+
    ylab("# of cells")
  
  # Plot Type 2: Barplots of average composition
  averagePerLMAType <- barPlotData %>% group_by(LMA, name) %>%
    summarise(prop = mean(prop))
  
  p2 <- averagePerLMAType %>%
    ggplot(aes(x=factor(LMA), y=prop, fill=name))+
    geom_col()+
    scale_fill_manual(name="Cell Type",
                      breaks=breaks,
                      values = cell_type_colors)+
    xlab(NULL)+
    theme_cowplot()+
    ylab("% of cells")
  
  # Plot Type 3: Distribution of sizes for each LMA subtype
  sumCellsByLMA <- LMAtypes %>% mutate(totalNCells = rowSums(select(., starts_with("CD")))) %>%
    select(Image, LMA, totalNCells)
  
  p3 <- sumCellsByLMA %>% ggplot(aes(x=log10(totalNCells)))+
    geom_histogram(aes(y=..density..), colour="white")+
    geom_density(alpha=.2, fill="blue")+
    scale_x_continuous(breaks=c(1,2,3,4), limits = c(1,4))+
    facet_wrap(~LMA, ncol=1)+
    theme_cowplot()
  
  # Plot Type 4: How many images are contributing to each LMA subtype
  influentialImages <- LMAtypes  %>% group_by(Image, LMA) %>% summarise(nOfLMAs = n())
  
  p4 <- influentialImages %>% ggplot(aes(x=nOfLMAs))+
    geom_histogram()+
    xlab("Number of LMAs found in a single image")+
    ylab("Number of images")+
    theme_cowplot()+
    facet_wrap(~LMA, ncol=1)
  
  ggarrange(p1, p2, 
            common.legend=TRUE)
}