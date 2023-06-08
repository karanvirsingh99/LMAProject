## JUN28_2022 = Community detection, not differentiating between epithelium or stroma. Using SPIAT package

## Spatial analysis - B&T Panel - 20Apr22
library(data.table)
library(tidyverse)
library(tripack)
library(ggsci)
library(igraph)
library(scales)
library(ComplexHeatmap)
library(ggvoronoi)
library(cowplot)
library(SPIAT)

levels_B <- c("Negative_S", "Negative_E", "CD3pCD8n", "CD3pCD8p", "CD79ApCD20n", "CD20p", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p", "NA")

colors_B <- c("#999999","#984EA3",
              "#377EB8",
              "#E41A1C",
              "#4DAF4A",
              "#FF7F00",
              "#FFFF33",
              "#F781BF",
              "white")


# Load in cell level data. This df already has the corresponding pNums and updated image names
cell_level <- fread("E:/COEUR_mcIF/final_data/cell_level_with_area.csv")

# Annotation level data
annotation_level<- read.csv("E:/COEUR_mcIF/final_data/COEUR_B&T_final_counts_and_densities.csv")

#Cluster IDS
clusterIDs<- read.csv("E:/COEUR_mcIF/COEUR_ADAPT_COMBINED/New PDL1 classifier analysis (May 2022) final/09MAY22_ADAPT_BT_DENSITIES_WITH_CLUSTER_IDS_NOPDL1CK_NOPD1_5_CLUSTERS.csv")

#######

# Step 1: Import cell level data (imgDataList)
imgDataList <- readRDS("./final_data/cell_level_data_RDS_panel.rds")

#Step 2: Run community identification on each image
format_qupath_to_sce <- function(image){
  image_no_neg <- image %>% filter(!(Class %in% c("Negative_S", "Negative_E")))
  if(nrow(image_no_neg)!=0){
    # converting to sce
    sce <- SingleCellExperiment()
    dummy_intensity <- rep(0, nrow(image_no_neg))
    intensity_assay <- matrix(dummy_intensity, nrow=1, ncol=nrow(image_no_neg))
    colnames(intensity_assay) <- NULL
    rownames(intensity_assay) <- NULL
    sce <- SingleCellExperiment(assays = list(counts = intensity_assay))
    rownames(sce) <- "DAPI"
    colnames(sce) <- paste("Cell_", image_no_neg$ID, sep="")
    coldata_phenotype <- image_no_neg$Class
    coldata_combinedType <- image_no_neg$combinedType
    coldata_Xpos <- image_no_neg$Centroid_x
    coldata_Parent <- image_no_neg$Parent
    
    #Flipping y axis
    coldata_Ypos <- max(image_no_neg$Centroid_y) - image_no_neg$Centroid_y
    colData(sce)$Phenotype <- coldata_phenotype
    colData(sce)$Cell.X.Position <- coldata_Xpos
    colData(sce)$Cell.Y.Position <- coldata_Ypos
    colData(sce)$CombinedType <- coldata_combinedType
    colData(sce)$Location <- coldata_Parent
    return(sce)
  } else(return(NULL))
}

find_communities_list <- function(image, radius, neighborhood_size){
  
  image_as_sce <- format_qupath_to_sce(image)
  if(!is.null(image_as_sce)){
    error_catch <- tryCatch({
      identify_neighborhoods(image_as_sce, cell_types_of_interest = c("CD20p", "CD3pCD8p", "CD3pCD8n", "CD79ApCD20n", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p"),
                             method="hierarchical", radius=radius, feature_colname = "Phenotype",
                             min_neighborhood_size = neighborhood_size)
      1
    }, error=function(e) 2)
    
    if(error_catch == 1){
      
      communities <- identify_neighborhoods(image_as_sce, cell_types_of_interest = c("CD20p", "CD3pCD8p", "CD3pCD8n", "CD79ApCD20n", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p"),
                                            method="hierarchical", radius=radius, feature_colname = "Phenotype",
                                            min_neighborhood_size = neighborhood_size)
      
      
      
      communities_as_df <- colData(communities) %>% as.data.frame()
      
      return(communities_as_df)
      
    } else(
      return(NULL)
    )
  } else(return(NULL))
}

# Cell level data with community from SPIAT # THIS TAKES A WHILE!!!
cell_level_with_communities <- lapply(imgDataList, find_communities_list, radius=15, neighborhood_size = 15)

#Remove images without communities
cell_level_with_communities <- cell_level_with_communities %>% purrr::discard(is.null)

#Add image as a column in each list entry
for(i in 1:length(cell_level_with_communities)){
  cell_level_with_communities[[i]]$Image <- names(cell_level_with_communities)[i]
}

#Summarise each community in terms of numbers of each individual cell type
community_summary <- lapply(cell_level_with_communities, function(x) {
  summary_no_location <- x %>% group_by(Neighborhood, Phenotype) %>% summarise(sum = n()) %>%
    pivot_wider(names_from = "Phenotype", values_from = "sum")
  summary_no_location[is.na(summary_no_location)] <- 0
  
  x$Location <- factor(x$Location, levels=c("Stroma", "Tumor"))
  
  summary_location <- x %>% group_by(Neighborhood, Location, .drop=FALSE) %>% summarise(sum = n()) %>%
    pivot_wider(names_from = "Location", values_from="sum")
  summary_location[is.na(summary_location)] <- 0
  
  summary <- full_join(summary_no_location, summary_location)
  
  summary <- summary %>% dplyr::rename(Epithelium = "Tumor") %>%
    filter(Neighborhood != "Free_cell")
})

community_summary <- data.table::rbindlist(community_summary, fill=TRUE, idcol="Image")
community_summary[is.na(community_summary)] <- 0

#Add back pNum
community_summary$pNum <- cell_level$pNum[match(community_summary$Image, cell_level$Image)]

# Cluster ID - pNum matching
community_summary <- left_join(community_summary, clusterIDs %>% select(pNum, cluster_ID))

# remove images that did not make it to survival analysis / heatmap
community_summary <- community_summary %>% filter(!is.na(cluster_ID))
final_images <- community_summary$Image

#Filter cell level list
cell_level_with_communities <- cell_level_with_communities[names(cell_level_with_communities) %in% final_images]

# Make matrix for clustering
communities_matrix <- community_summary %>% select(-Neighborhood, -Image, -cluster_ID, -pNum, -Epithelium, -Stroma) %>%
  as.matrix()

# Histogram distribution of size
hist_data <- community_summary %>% mutate(sum = (CD3pCD8n + CD3pCD8p + CD3pCD8nFoxP3p + CD3pCD8pFoxP3p +
                                                   CD79ApCD20n + CD20p))

hist_data %>% ggplot(aes(x=sum))+
  geom_histogram()+
  theme_cowplot()+
  xlab("Total # of cells in community")+
  ylab("Number of communities")

hist_data %>% ggplot(aes(x=sum))+
  geom_histogram(binwidth = 5)+
  theme_cowplot()+
  xlim(c(0, 100))+
  xlab("Total # of cells in community")+
  ylab("Number of communities")

# Trying k=2 to k=10, checking how the communities are being clustered at each step
# in terms of overall cell frequencies and composition

# But first, I am going to compute the total number of cells and total number of positive cells in each image.
# So that I can see how much 'area' approximately each community is occupying

categories_of_interest  = c("CD3pCD8n", "CD3pCD8p", "CD79ApCD20n", "CD20p", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p")

total_cells_per_image <- cell_level %>% group_by(Image) %>% 
  summarise(total_n_cells = length(Class), total_positive = sum(Class %in% categories_of_interest))

community_summary <- left_join(community_summary, total_cells_per_image)

# number of total cells per community
community_summary <- community_summary %>% 
  mutate(sum = (CD3pCD8n + CD3pCD8p + CD3pCD8nFoxP3p + CD3pCD8pFoxP3p + CD79ApCD20n + CD20p))

# percent of total cells made up by each community, and percent of positive cells
community_summary <- community_summary %>% ungroup() %>% 
  mutate(perc_total = sum/total_n_cells,
         perc_positive = sum/total_positive)

#Hclust on raw counts
communities_matrix_prop <- t(apply(communities_matrix,1, function(x) x/sum(x)))

# matrix <- communities_matrix # raw counts
matrix <- communities_matrix_prop #proportions

# Calculate dissimilarity and build dendrogram
hier_bray <- vegan::vegdist(matrix, method = "bray")
tree <- hclust(hier_bray, method="ward.D2")
tree_dend <- as.dendrogram(tree)
plot(tree_dend, leaflab = "none")

number_of_clusters=7

png(paste(number_of_clusters,"_dendrogram",".png", sep = ""))
library(dendextend)
clust <- cutree(tree, k=number_of_clusters)
clust.cutree <- dendextend:::cutree(tree_dend, k=number_of_clusters, order_clusters_as_data = FALSE)
idx <- order(as.numeric(names(clust.cutree)))
clust.cutree <- clust.cutree[idx]
tbl <- table(clust, clust.cutree)
lbls <- apply(tbl,2,which.max)
dend1 <- color_branches(tree_dend, k = number_of_clusters, groupLabels = lbls)
dend1 %>% set("labels", NULL) %>% plot()
dev.off()

# Cut the tree at specified k
clusters<- cutree(tree, number_of_clusters)

# Heatmap to visualize composition of each cluster with viridis scale
col = circlize::colorRamp2(seq(0, 100, length=999),
                           viridis(999))

heatmap <- ComplexHeatmap::Heatmap(name="Number of cells",
                                   t(communities_matrix),
                                   cluster_columns = FALSE,
                                   cluster_rows = FALSE,
                                   column_split  = clust,
                                   col=col)

png(paste(number_of_clusters,"_heatmap", ".png", sep=""), res=300,
    width=7.5,
    height=2.5,
    units = "in")
draw(heatmap)
dev.off()

#Merge with original dataframe
community_summary$New_neighborhood <- paste("Cluster_", clust, sep="")
community_summary$New_neighborhood <- factor(community_summary$New_neighborhood)

levels(community_summary$New_neighborhood) <- c("T2", "T3", "P2", "P1", "P3", "T1", "B1")
community_summary$New_neighborhood <- ordered(community_summary$New_neighborhood, levels=c("T1", "T2", "T3", "P1", "P2", "P3", "B1"))

#### Let's look at each individual cluster ####

p1 <- community_summary %>% ggplot(aes(x=New_neighborhood, y=sum))+
  geom_boxplot()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cluster ID")+
  ylab("Number of cells")

## How much of the total cells do the cells occupy?

p2 <- community_summary %>% ggplot(aes(x=New_neighborhood, y=perc_total))+
  geom_boxplot()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cluster ID")+
  ylab("/tot. cells")

p3 <- community_summary %>% ggplot(aes(x=New_neighborhood, y=perc_positive))+
  geom_boxplot()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Cluster ID")+
  ylab("/pos. cells")

ggpubr::ggarrange(p1, p2, p3, nrow = 1)
ggsave(paste(number_of_clusters,"_boxplots",".png", sep=""),
       width=7,
       height=5,
       units="in")

#How many cells of each type in each cluster
bar_data <- community_summary %>%
  mutate(COMMUNITY_ID = row_number()) %>%
  select(-c("Image", "Neighborhood", "sum", "pNum", "cluster_ID", "Stroma", "Epithelium")) %>% pivot_longer(cols = 
                                                                                                              -c("COMMUNITY_ID", 
                                                                                                                 "New_neighborhood",
                                                                                                                 "total_n_cells",
                                                                                                                 "total_positive",
                                                                                                                 "perc_total",
                                                                                                                 "perc_positive"))

p4 <- bar_data %>% 
  ggplot(aes(x=factor(COMMUNITY_ID), y=value, fill=name))+
  geom_col()+
  facet_wrap(~New_neighborhood, scales="free")+
  scale_fill_manual(name = "Cell Type",
                    breaks = c("CD3pCD8n", "CD3pCD8p", "CD79ApCD20n", "CD20p", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p", "NA"),
                    values =  c("#377EB8",
                                "#E41A1C",
                                "#4DAF4A",
                                "#FF7F00",
                                "#FFFF33",
                                "#F781BF"))+
  theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab(NULL)+
  ylab("# of cells")

# Plot proportions instead

bar_data_prop <- bar_data %>% group_by(COMMUNITY_ID) %>% mutate(prop = value/sum(value))

p5 <- bar_data_prop %>% 
  ggplot(aes(x=factor(COMMUNITY_ID), y=prop, fill=name))+
  geom_col()+
  facet_wrap(~New_neighborhood, scales="free", ncol=1)+
  scale_fill_manual(name = "Cell Type",
                    breaks = c("CD3pCD8n", "CD3pCD8p", "CD79ApCD20n", "CD20p", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p", "NA"),
                    values =  c("#377EB8",
                                "#E41A1C",
                                "#4DAF4A",
                                "#FF7F00",
                                "#FFFF33",
                                "#F781BF"))+
  theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position= "bottom")+
  xlab(NULL)+
  ylab("% of total cells")

ggarrange(p4, p5, common.legend = TRUE, nrow=2)

ggsave(paste(number_of_clusters,"_bar_plots", ".png", sep=""),
       width=10,
       height=8,
       units="in")


# Plot proportion of stroma and epithelium in bar graph
how_much_epi <- community_summary %>%
  mutate(COMMUNITY_ID = row_number()) %>%
  select(COMMUNITY_ID, Image, Stroma, Epithelium, New_neighborhood)

how_much_epi$COMMUNITY_ID <- paste("C", how_much_epi$COMMUNITY_ID, sep="")

how_much_epi <- how_much_epi %>% pivot_longer(cols = c("Stroma", "Epithelium"), names_to = "Location") %>%
  group_by(COMMUNITY_ID) %>%
  mutate(value = value/sum(value))

communities_desc_epi <- how_much_epi %>% filter(Location=="Epithelium") %>% arrange(desc(value)) %>% .$COMMUNITY_ID

how_much_epi %>% ggplot(aes(x=ordered(COMMUNITY_ID, levels=communities_desc_epi), y=value, fill=Location))+
  geom_col()+
  facet_wrap(~New_neighborhood, scales = "free", ncol=1)+
  theme_cowplot()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position= "bottom")+
  xlab(NULL)+
  ylab("Proportion of cells in each community")+
  scale_fill_npg()

#Examples of communities with >0.75 epithelial cells from each cluster
images_with_epithelial_communities <- how_much_epi %>% filter(Location=="Epithelium") %>% filter(value > 0.75) %>%
  group_by(New_neighborhood) %>%
  sample_n(size=5, replace = TRUE)

#Add new cluster ID to cell level list
for(i in 1:length(cell_level_with_communities)){
  only_single_image <- community_summary[community_summary$Image %in% cell_level_with_communities[[i]]$Image,]
  
  cell_level_with_communities[[i]]$New_Neighborhood <- only_single_image$New_neighborhood[match(cell_level_with_communities[[i]]$Neighborhood,
                                                                                                only_single_image$Neighborhood)]
  
  cell_level_with_communities[[i]]$New_Neighborhood <- as.character(cell_level_with_communities[[i]]$New_Neighborhood)
  
  cell_level_with_communities[[i]]$New_Neighborhood[is.na(cell_level_with_communities[[i]]$New_Neighborhood)] <- "Free_cell"
}


# A function that produces three plots: the scatterplot with just the figure and markers,
# A plot with the neighborhoods (cluster output)
# And a heatmap showing the composition of each cluster

makePrettyPlots <- function(cell_level_data){
  
  new_cluster_labels <- unique(cell_level_data %>% filter(New_Neighborhood!= "Free_cell") %>% .$New_Neighborhood)
  og_cluster_labels <- unique(cell_level_data %>% filter(New_Neighborhood!= "Free_cell") %>% .$Neighborhood) #minus free cells
  label_location <- vector()
  for (Cluster in og_cluster_labels) {
    cells_in_Cluster <- cell_level_data[cell_level_data$Neighborhood == 
                                          Cluster, ]
    minX <- min(cells_in_Cluster$Cell.X.Position)
    maxX <- max(cells_in_Cluster$Cell.X.Position)
    minY <- min(cells_in_Cluster$Cell.Y.Position)
    maxY <- max(cells_in_Cluster$Cell.Y.Position)
    averageX <- (minX + maxX)/2
    averageY <- (minY + maxY)/2
    label_location <- rbind(label_location, c(Cluster, 
                                              averageX, averageY))
  }
  label_location <- as.data.frame(label_location)
  colnames(label_location) <- c("Neighborhood", "Xpos", "Ypos")
  label_location$Xpos <- as.numeric(label_location$Xpos)
  label_location$Ypos <- as.numeric(label_location$Ypos)
  
  label_location$New_Neighborhood <- cell_level_data$New_Neighborhood[match(label_location$Neighborhood,
                                                                            cell_level_data$Neighborhood)]
  
  cluster_colours <- dittoSeq::dittoColors()[1:length(og_cluster_labels)]
  
  p1 <- ggplot(cell_level_data, aes(x = Cell.X.Position, 
                                    y = Cell.Y.Position))
  p1 <- p1 + geom_point(aes(color = Neighborhood))
  p1 <- p1 + geom_text(data = label_location, aes(x = Xpos, 
                                                  y = Ypos, label = gsub("Cluster_", "",
                                                                         Neighborhood)))
  p1 <- p1 + scale_color_manual(breaks = c(as.character(unique(og_cluster_labels)), "Free_cell"), 
                                values = c(cluster_colours, "#808080"))
  p1 <- p1 + xlab("Cell.X.Position") + ylab("Cell.Y.Position") + 
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(), legend.position = "none")
  
  
  counts_per_group <- cell_level_data 
  counts_per_group$Neighborhood <- factor(counts_per_group$Neighborhood)
  counts_per_group$Phenotype <- factor(counts_per_group$Phenotype, levels=c("CD3pCD8n", "CD3pCD8p", "CD79ApCD20n", "CD20p", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p"))
  counts_per_group <- counts_per_group %>% group_by(Phenotype, Neighborhood, .drop = FALSE) %>%
    summarise(sum = n())
  
  counts_per_grp_m <- counts_per_group %>%
    pivot_wider(names_from = Neighborhood, values_from = sum) %>%
    select(-Free_cell) %>%
    column_to_rownames("Phenotype") %>%
    as.matrix()
  
  new_clusters <- cell_level_data$New_Neighborhood[match(colnames(counts_per_grp_m), cell_level_data$Neighborhood)]
  new_clusters <- gsub("Cluster_", "", new_clusters)
  
  colnames(counts_per_grp_m) <- gsub("Cluster_", "Community_", colnames(counts_per_grp_m))
  
  clust_col <- pal_npg()(7)
  names(clust_col) <- c(paste("Cluster_", seq(1,7), sep=""))
  
  column_annotation <- HeatmapAnnotation(Cells_Sum = anno_barplot(colSums(counts_per_grp_m)))
  
  proportions <- apply(counts_per_grp_m,2, function(x) x/sum(x))
  
  col_fun = circlize::colorRamp2(seq(0, 0.4, length=10), viridis(10))
  
  h1 <- ComplexHeatmap::Heatmap(name="Proportion of cell",
                                proportions,
                                top_annotation = column_annotation,
                                column_split = new_clusters,
                                col = col_fun,
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                heatmap_legend_param = list(
                                  legend_direction = "horizontal"))
  
  gridExtra::grid.arrange(grobs = list(p1, grid::grid.grabExpr(draw(h1, heatmap_legend_side = "top"))), nrow=1, ncol=2, vp=vp)
  
}


# Images in total = 454
# Patients in total = 346

# Distribution of patients with communities across 5 clusters

# read cluster data
logged_densities <- read.csv("E:/COEUR_mcIF/COEUR_ADAPT_COMBINED/New PDL1 classifier analysis (May 2022) final/09MAY22_ADAPT_BT_DENSITIES_WITH_CLUSTER_IDS_NOPDL1CK_NOPD1_5_CLUSTERS.csv")

logged_densities <- logged_densities %>% mutate(has_community = case_when(pNum %in% community_summary$pNum ~ "YES",
                                                                          TRUE ~ "NO"))

# How many patients have a detected community in each cluster?
has_community <- logged_densities %>% group_by(cluster_ID, has_community) %>% summarise(n = n())

has_community %>% ggplot(aes(cluster_ID, n, fill=has_community))+
  geom_col()+
  scale_fill_npg()+
  theme_cowplot()+
  xlab("Density cluster ID")+
  ylab("Number of cases")


# How does the distribution of the 7 community types differ across the 5 clusters?
community_summary %>% group_by(cluster_ID, New_neighborhood) %>% summarise(n = n()) %>%
  ggplot(aes(cluster_ID, n, fill=New_neighborhood))+
  geom_col()+
  scale_fill_npg(name="Type of community")+
  xlab("Density cluster ID (from heatmap)")+
  ylab("Number of communities")+
  theme_cowplot()


communities_tally <- community_summary %>%
  group_by(pNum, New_neighborhood) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = "New_neighborhood", values_from = "n") %>%
  select(pNum, T1, T2, T3, P1, P2, P3, B1)

communities_tally[is.na(communities_tally)] <- 0

cor(communities_tally[,-1]) %>%
  corrplot::corrplot(method = "number")

how_many_communities <- rowSums(communities_tally[,-1])

how_many_communities <- communities_tally %>% ungroup %>%  mutate(sum = rowSums(.[,-1]))

how_many_communities %>% ggplot(aes(x=sum))+
  geom_histogram(binwidth = 1)+
  theme_cowplot()+
  ylab("Number of patients")+
  xlab("Number of communities per patients")

how_many_communities %>% filter(sum==1)

how_many_communities$clusterID <- community_summary$cluster_ID[match(how_many_communities$pNum, community_summary$pNum)]

ggplot(how_many_communities, aes(x="", y=sum))+
  geom_boxplot()+
  facet_wrap(~clusterID, nrow=1)+
  theme_cowplot()+
  xlab("Density cluster")+
  ylab("Number of detected communities per patient")

community_summary %>% group_by(pNum, cluster_ID) %>% mutate(total_com_per_patient = n()) %>%
  filter(total_com_per_patient == 1) %>%
  group_by(pNum, cluster_ID, New_neighborhood) %>% summarise(n = n()) %>%
  ggplot(aes(x=New_neighborhood, n, fill=factor(cluster_ID)))+
  geom_col()+
  scale_fill_npg(name="Density Cluster")+
  xlab("Neighborhood ID")+
  ylab("Number of patients")+
  theme_cowplot()

community_summary %>% group_by(pNum, cluster_ID) %>% mutate(total_com_per_patient = n()) %>%
  filter(total_com_per_patient == 1) %>%
  group_by(pNum, cluster_ID, New_neighborhood) %>% summarise(n = n()) %>%
  ggplot(aes(cluster_ID, n, fill=New_neighborhood))+
  geom_col()+
  scale_fill_npg(name="Type of community")+
  xlab("Density cluster ID (from heatmap)")+
  ylab("Proportion of communities")+
  theme_cowplot()

more_than_one <- community_summary %>% group_by(pNum, cluster_ID) %>% mutate(total_com_per_patient = n()) %>%
  filter(total_com_per_patient > 1)

more_than_one <- more_than_one %>% group_by(pNum, New_neighborhood) %>% summarise(total = n()) %>%
  pivot_wider(names_from="New_neighborhood", values_from = total, values_fill = 0) %>%
  select(T1, T2, T3, P1, P2, P3, B1) %>%
  column_to_rownames("pNum")

more_than_one_binary <- decostand(more_than_one, method="pa")

more_than_one_binary[more_than_one_binary==1] <- "y"
more_than_one_binary[more_than_one_binary==0] <- "n"

more_than_one_binary <- more_than_one_binary %>% mutate(combination = paste(T1, T2, T3, P1, P2 , P3, B1, sep=""))

more_than_one_binary %>% group_by(combination) %>% summarise(total = n()) %>% arrange(desc(total))


# Combined with intraepithelial neighborhoods

all_neighborhoods <- full_join(communities_tally, ie_tally %>% ungroup() %>% select(-cluster_ID), by="pNum")

all_neighborhoods[is.na(all_neighborhoods)] <- 0

cor(all_neighborhoods[,-1]) %>%
  corrplot::corrplot(method = "number", p.mat=p_Values)

p_Values <- psych::corr.test(all_neighborhoods[,-1])$p




# Minimum distance
final_images <- cell_level$Image[cell_level$pNum %in% clusterIDs$pNum] %>% unique()

avg_min_distance_epi <- function(image){
  sce_object <- format_qupath_to_sce(imgDataList[[image]])
  if (is.null(sce_object)){
    return(NULL)
  }
  formatted_data <- data.frame(SummarizedExperiment::colData(sce_object))
  formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
  only_epithelium <- formatted_data %>% filter(Location=="Tumor")
  all_cell_cords <- only_epithelium[, c("Cell.X.Position", 
                                     "Cell.Y.Position")]
if (nrow(all_cell_cords) < 2) {
 return(NULL)
}
all_closest <- RANN::nn2(data = all_cell_cords, k = 2)
all_closest_dist <- all_closest$nn.dists[, 2]
average_min_distance <- mean(all_closest_dist)
return(average_min_distance)
}

distances_epi <- list()
for(image in final_images){
  distances_epi[[image]] <- avg_min_distance_epi(image)
}

distances_epi <- unlist(distances_epi)
distances_epi <- data.frame(distances_epi)

avg_min_distance_str <- function(image){
  sce_object <- format_qupath_to_sce(imgDataList[[image]])
  if (is.null(sce_object)){
    return(NULL)
  }
  formatted_data <- data.frame(SummarizedExperiment::colData(sce_object))
  formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
  only_epithelium <- formatted_data %>% filter(Location=="Stroma")
  all_cell_cords <- only_epithelium[, c("Cell.X.Position", 
                                        "Cell.Y.Position")]
  if (nrow(all_cell_cords) < 2) {
    return(NULL)
  }
  all_closest <- RANN::nn2(data = all_cell_cords, k = 2)
  all_closest_dist <- all_closest$nn.dists[, 2]
  average_min_distance <- mean(all_closest_dist)
  return(average_min_distance)
}

distances_str <- list()
for(image in final_images){
  distances_str[[image]] <- avg_min_distance_str(image)
}

distances_str <- unlist(distances_str)
distances_str <- as.data.frame(distances_str)

p <- ggplot(distances_epi, aes(x="", y=distances_epi))+
  geom_boxplot()+
  ylab("Minimum distance, microns")+
  xlab("Intraepithelial cells")+
  theme_cowplot()

pp <- ggplot(distances_str, aes(x="", y=distances_str))+
  geom_boxplot()+
  ylab(NULL)+
  xlab("Stromal cells")+
  scale_y_continuous(limits = c(0, 600))+
  theme_cowplot()

ggpubr::ggarrange(p, pp, align = "v")


per_patient_summary %>% ggplot(aes(x=New_neighborhood, y=perc_total))+ geom_boxplot() + 
  # facet_wrap(~New_neighborhood, nrow=1, scales='free')+
  theme_cowplot()+
  xlab("Community type")+
  ylab("Percent of total cells in the community")

per_patient_summary %>% ggplot(aes(x=factor(cluster_ID), y=perc_total))+ geom_boxplot() + 
  facet_wrap(~New_neighborhood, nrow=1, scales='free_x')+
  theme_cowplot()+
  xlab("Community type")+
  ylab("Percent of total cells in the community")

per_patient_summary <- community_summary %>% group_by(cluster_ID, pNum, New_neighborhood) %>%
  summarise(perc_positive = sum(perc_positive),
            perc_total = sum(perc_total))
