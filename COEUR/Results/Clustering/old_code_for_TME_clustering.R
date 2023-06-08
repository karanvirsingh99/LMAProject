# 26APR22 - with new adapt panel densities
# Take combined densities between adaptive and BT panel, and cluster on logged densities.

library(tidyverse)
library(RSKC)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(corrplot)
library(igraph)
library(ggsci)
library(scales)

# Next - I should make violin plots for the densities in each tissue and a quick RKSC clustering with ALL of the densities
# I need to figure out which phenotypes overlap..

COMBINED_PANEL_DENSITIES <- read.csv("E:/COEUR_mcIF/final_data/ADAPT_BT_COMBINED_INNER_JOIN_04MAY22.csv")

# This is the list of patients that are (1) pre-chemo (2) have outcomes data
final_pNum <- read.csv("./final_data/final_pNum.csv") %>% .[,1]


##### Cleaning up the data####
DENSITIES_BY_TISSUE <- COMBINED_PANEL_DENSITIES %>% select(-contains("area")) %>%
  select(-PDL1pCKp_S, -BT_str_prop, -BT_epi_prop, -ADAPT_epi_prop, -ADAPT_str_prop, -avg_epi_prop) %>%
  select(-CD3pCD8p_S,-CD3pCD8p_E)

DENSITIES_LOGGED <- DENSITIES_BY_TISSUE %>%
  mutate(across(-contains(c("prop", "area", "pNum")), ~log10(.+0.1))) #log all densities

# Which CD3 columns do I keep?

cormatrix <- cor(DENSITIES_LOGGED[,-1])
corrplot(cormatrix, method = 'number')

#### Making heatmap with only final patients (n=983) ####

# Making a logged densities matrix

densities_final_m <- DENSITIES_BY_TISSUE %>% 
  filter(pNum %in% final_pNum) %>%
  select(-contains(c("prop", "area")), -PDL1pCKp_E, -CD68nPDL1nPD1p_S, -CD68nPDL1nPD1p_E) %>% 
  column_to_rownames("pNum") %>% 
  mutate(across(everything(), ~log10(.+0.1))) %>% 
  as.matrix()

#optimal_n_final <- Clest(densities_final_m, 8, alpha=0.05, beta=0.05) #Trimming 10% of cases, significance level between clusters at least 0.05

#Optimal cluster N is 4 without or with PDL1+CK+
set.seed(1000)
#clusters <- RSKC(densities_final_m, 5, alpha=0.05)$labels

# Freeze output
clusters <- c(2, 3, 3, 2, 1, 3, 3, 3, 4, 2, 5, 3, 3, 1, 4, 2, 5, 4, 4, 1, 
              2, 4, 3, 1, 3, 4, 2, 5, 4, 1, 3, 4, 2, 1, 4, 4, 5, 4, 2, 1, 4, 
              5, 3, 2, 4, 2, 5, 3, 2, 4, 1, 4, 5, 2, 4, 1, 4, 5, 4, 1, 3, 2, 
              4, 3, 5, 1, 4, 3, 2, 3, 1, 5, 5, 5, 4, 5, 3, 5, 4, 4, 1, 1, 3, 
              1, 5, 1, 1, 3, 2, 3, 1, 5, 5, 3, 5, 2, 2, 4, 4, 1, 4, 1, 4, 4, 
              4, 5, 3, 4, 4, 3, 2, 1, 5, 2, 4, 2, 5, 3, 2, 2, 1, 2, 2, 4, 4, 
              2, 4, 1, 5, 5, 3, 5, 4, 2, 5, 1, 3, 5, 2, 2, 3, 5, 1, 4, 1, 4, 
              4, 4, 5, 5, 5, 5, 2, 4, 2, 4, 1, 5, 4, 3, 4, 3, 1, 4, 2, 5, 1, 
              3, 1, 5, 2, 3, 3, 3, 1, 3, 5, 5, 4, 4, 1, 4, 3, 2, 1, 5, 5, 4, 
              3, 3, 2, 4, 4, 4, 2, 5, 4, 4, 3, 3, 3, 4, 1, 2, 4, 2, 1, 2, 4, 
              3, 2, 2, 1, 5, 2, 4, 5, 5, 4, 1, 3, 3, 4, 5, 3, 4, 1, 5, 1, 2, 
              3, 2, 5, 4, 4, 5, 2, 3, 5, 3, 2, 3, 2, 1, 1, 5, 2, 5, 4, 2, 4, 
              3, 1, 5, 1, 4, 4, 4, 4, 3, 4, 4, 4, 2, 3, 2, 5, 3, 4, 4, 5, 2, 
              1, 3, 4, 5, 4, 4, 5, 4, 5, 4, 1, 5, 1, 2, 4, 3, 3, 1, 4, 5, 4, 
              2, 3, 2, 3, 2, 2, 5, 4, 1, 3, 3, 4, 1, 2, 3, 4, 3, 2, 4, 4, 5, 
              1, 1, 3, 1, 4, 2, 2, 3, 5, 2, 3, 4, 5, 2, 4, 3, 4, 4, 4, 5, 4, 
              1, 4, 2, 5, 3, 3, 4, 3, 3, 5, 4, 4, 2, 5, 3, 5, 1, 3, 4, 4, 1, 
              5, 3, 1, 1, 2, 2, 3, 2, 5, 4, 1, 4, 2, 1, 4, 2, 2, 4, 1, 4, 3, 
              1, 3, 1, 4, 4, 5, 4, 3, 1, 1, 3, 1, 2, 3, 4, 1, 5, 1, 2, 5, 4, 
              1, 4, 2, 4, 2, 2, 1, 4, 3, 1, 1, 3, 4, 2, 1, 2, 5, 3, 1, 1, 4, 
              1, 1, 1, 4, 3, 3, 2, 1, 2, 1, 4, 3, 4, 5, 1, 4, 5, 4, 1, 3, 5, 
              4, 5, 2, 4, 1, 5, 4, 1, 1, 5, 2, 3, 5, 1, 5, 4, 5, 4, 2, 3, 3, 
              1, 3, 3, 1, 2, 3, 5, 5, 3, 4, 2, 4, 3, 5, 4, 2, 5, 4, 5, 5, 1, 
              3, 1, 3, 3, 3, 2, 4, 5, 4, 3, 4, 2, 5, 4, 5, 4, 5, 5, 1, 1, 3, 
              1, 2, 5, 3, 4, 4, 4, 1, 4, 3, 3, 4, 2, 2, 3, 4, 2, 3, 5, 5, 5, 
              1, 2, 4, 5, 4, 5, 3, 1, 2, 2, 1, 2, 5, 2, 4, 4, 1, 3, 2, 1, 4, 
              5, 4, 4, 5, 1, 5, 4, 3, 3, 3, 5, 4, 4, 5, 2, 5, 4, 3, 3, 3, 2, 
              5, 5, 4, 3, 5, 1, 3, 4, 4, 1, 5, 2, 1, 2, 5, 1, 1, 4, 5, 5, 4, 
              2, 3, 1, 3, 1, 4, 1, 2, 4, 5, 1, 4, 1, 1, 1, 5, 2, 4, 5, 5, 4, 
              5, 5, 3, 5, 5, 5, 1, 3, 2, 5, 1, 1, 3, 5, 2, 4, 2, 3, 4, 2, 5, 
              5, 4, 5, 3, 1, 4, 5, 5, 3, 5, 5, 1, 1, 5, 1, 4, 3, 2, 3, 5, 1, 
              1, 2, 5, 5, 3, 5, 5, 3, 3, 3, 4, 3, 2, 3, 3, 5, 4, 2, 2, 1, 5, 
              3, 5, 2, 3, 3, 1, 5, 4, 2, 2, 2, 5, 3, 4, 3, 2, 2, 1, 2, 1, 5, 
              5, 1, 2, 2, 3, 2, 5, 4, 4, 1, 4, 2, 2, 2, 1, 4, 4, 5, 2, 5, 4, 
              3, 4, 5, 5, 1, 1, 4, 2, 5, 2, 5, 4, 4, 5, 4, 4, 1, 4, 2, 4, 1, 
              4, 4, 2, 3, 2, 4, 4, 5, 2, 2, 3, 2, 3, 2, 4, 4, 3, 4, 2, 5, 4, 
              3, 5, 5, 4, 1, 4, 5, 2, 1, 4, 3, 4, 1, 3, 1, 3, 5, 1, 3, 3, 1, 
              1, 5, 4, 5, 3, 5, 1, 5, 2, 4, 3, 4, 2, 2, 1, 3, 5, 4, 4, 4, 4, 
              2, 3, 4, 4, 5, 4, 5, 3, 1, 3, 2, 5, 3, 1, 4, 2, 2, 3, 4, 4, 5, 
              5, 5, 4, 5, 2, 1, 2, 4, 2, 4, 5, 5, 2, 5, 2, 2, 5, 2, 4, 4, 5, 
              5, 4, 3, 5, 3, 2, 5, 5, 1, 3, 5, 4, 5, 4, 4, 1, 3, 5, 4, 3, 4, 
              3, 5, 3, 2, 5, 2, 2, 3, 3, 4, 5, 2, 4, 4, 3, 4, 5, 4, 4, 5, 1, 
              4, 4, 1, 1, 3, 4, 3, 3, 4, 3, 2, 4, 3, 5, 4, 2, 4, 5, 3, 5, 5, 
              5, 1, 2, 2, 5, 1, 4, 3, 5, 2, 1, 4, 5, 2, 4, 5, 2, 5, 1, 3, 3, 
              1, 5, 5, 2, 5, 2, 3, 3, 4, 1, 3, 3, 3, 5, 4, 3, 5, 1, 1, 4, 1, 
              5, 1, 5, 1, 5, 4, 4, 5, 5, 1, 5, 2, 4, 1, 3, 5, 1, 5, 3, 3, 5, 
              2, 5, 4, 5, 4, 4, 1, 3, 3, 5, 1, 1, 1, 2, 1, 3)

# Change order of labels. 2 -> 1, 5 -> 2, 3 ->3, 4 -> 4, 1 -> 5
new_clusters <- recode(clusters, '3' = 1L, '1' = 2L, '2' = 3L, '4' = 4L, '5' = 5L)

#new densities DF with clusters
DENSITIES_LOGGED_FINAL <- DENSITIES_LOGGED %>% filter(pNum %in% final_pNum) %>% 
  mutate(cluster_ID = new_clusters)

# subset data to only include densities
H1 <- DENSITIES_LOGGED_FINAL %>%
  mutate(clusterID = new_clusters) %>%
  arrange(clusterID) %>%
  select(-contains(c("prop", "category", "cluster", "ratio")), -CD68nPDL1nPD1p_S, -CD68nPDL1nPD1p_E, -PDL1pCKp_E) %>% column_to_rownames("pNum")

#renaming the columns so that they're prettier
colnames(H1) <- gsub("p", "+", colnames(H1))
colnames(H1) <- gsub("n", "-", colnames(H1))
H1 <- rename(H1, c("CD79a+CD20+_E" = "CD20+_E",
                       "CD79a+CD20+_S" = "CD20+_S",
                   "CD79a+CD20-_E" = "CD79A+CD20-_E",
                        "CD79a+CD20-_S" = "CD79A+CD20-_S",
                   "CD3+CD8-FoxP3-_E" = "CD3+CD8-_E",
                   "CD3+CD8-FoxP3-_S" = "CD3+CD8-_S"))

#the order of the rows (grouped by cell type)
H1_row_order <- c(
  "CD3+CD8-FoxP3-_S", 
  "CD3+CD8-FoxP3-_E",
  "CD3+CD8-PD1-_S",
  "CD3+CD8-PD1-_E",
  "CD3+CD8-PD1+_S",
  "CD3+CD8-PD1+_E",
  "CD3+CD8+PD1-_S",
  "CD3+CD8+PD1-_E",
  "CD3+CD8+PD1+_S",
  "CD3+CD8+PD1+_E",
  "CD3+CD8+FoxP3+_S",
  "CD3+CD8+FoxP3+_E",
  "CD3+CD8-FoxP3+_S",
  "CD3+CD8-FoxP3+_E",
  "CD68+PDL1-_S",
  "CD68+PDL1-_E",
  "CD68+PDL1+_S",
  "CD68+PDL1+_E",
  "CD79a+CD20+_S",
  "CD79a+CD20+_E",
  "CD79a+CD20-_S",
  "CD79a+CD20-_E")

#split the rows by cell type
row_split <- c(rep("CD8-", 6),
               rep("CD8+", 6),
               rep("T reg", 2),
               rep("Macrophage", 4), 
               rep("B cell", 4))[match(colnames(H1), H1_row_order)]

#legend color breaks
col_fun <- colorRamp2(c(0,2,4), c("#2c7bb6", "#ffffbf", "#d7191c"))

#Cluster labels
clust_labels <- paste("C",1:5, " (", table(new_clusters), ")", sep="")

#Row names
row_names <- gsub("+", "<sup>+</sup>", colnames(H1), fixed = TRUE) %>%
  gsub("-", "<sup>-</sup>", ., fixed = TRUE) %>%
  gsub("_E", "", .) %>%
  gsub("_S", "", .)

lgd = Legend(labels = c("Stromal", "Intraepithelial"),
             title = "Cell Location", legend_gp = gpar(fill = c("#3C5488FF","#DC0000FF")),
             nrow = 1)

heatmap1 <- Heatmap(
  t(H1),
  # name = "hey",
  row_split = factor(
    row_split,
    levels = c("CD8+", "CD8-", "T reg", "Macrophage", "B cell")
  ),
  row_labels = gt_render(row_names),
  row_names_gp = gpar(col = c(rep(c("#3C5488FF","#DC0000FF"), times = 12), "#DC0000FF")[match(colnames(H1), H1_row_order)]),
  column_split = sort(new_clusters),
  row_order = H1_row_order,
  row_gap = unit(2.5, "mm"),
  cluster_row_slices = FALSE,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  top_annotation = HeatmapAnnotation(cluster = anno_block(
    gp = gpar(fill = pal_npg()(5)),
    labels = clust_labels,
    labels_gp = gpar(col ="white", fontsize = 10)
  )),
  heatmap_legend_param = list(
    title = gt_render("log<sub>10</sub>(cells/mm<sup>2</sup> + 0.1)"),
    legend_width = unit(6, "cm"),
    direction = "horizontal",
    at = c(-1, 0, 1, 2, 3, 4)
  ),
  col = col_fun,
  column_title = NULL,
  border = TRUE
)

draw(heatmap1, heatmap_legend_side = "bottom")
draw(lgd, x=unit(0.8, "npc"), y=unit(0.05, "npc"))

# Let's make some boxplots to compare the densities of each cell type in each cluster
DENSITIES_LONGER <- DENSITIES_LOGGED_FINAL %>% select(-avg_str_prop) %>%
  mutate(cluster_ID = factor(cluster_ID)) %>%
  pivot_longer(cols=-c("cluster_ID", "pNum"), names_sep = "_", names_to=c("density", "location"))
DENSITIES_LONGER$density <- gsub("p", "+", DENSITIES_LONGER$density, fixed=TRUE) %>%
  gsub("n", "-", ., fixed = TRUE) %>%
  gsub("CD20+", "CD79a+CD20+", ., fixed = TRUE) %>%
  replace(., .=="CD3+CD8-", "CD3+CD8-FoxP3-") %>%
  replace(., .=="CD79A+CD20-", "CD79a+CD20-") %>%
  replace(., .=="CD68-PDL1-PD1+", "CD3-CD8-PD1+") %>%
  replace(., .=="PDL1pCKp", "PDL1+CK+")

# add jitter..
DENSITIES_LONGER <- DENSITIES_LONGER %>% group_by(density, location, cluster_ID) %>%
  mutate(outlier = value > median(value) + 
           IQR(value)*1.5 | value < median(value) -
           IQR(value)*1.5) 

T_cells = c("CD3+CD8-FoxP3-", "CD3+CD8-PD1-",  "CD3+CD8-PD1+", "CD3+CD8+PD1-", "CD3+CD8+PD1+", "CD3+CD8-FoxP3+", "CD3+CD8+FoxP3+")
T_cells_subset = c("CD3+CD8-PD1-",  "CD3+CD8-PD1+", "CD3+CD8+PD1-", "CD3+CD8+PD1+", "CD3+CD8+FoxP3+")
B_cells = c("CD79a+CD20+", "CD79a+CD20-")
Mac = c("CD68+PDL1+", "CD68+PDL1-")
other = "PDL1+CK+"


stats <- DENSITIES_LONGER %>% group_by(density, location) %>%
  t_test(value ~ cluster_ID) %>%
  add_xy_position()

median_densities <- DENSITIES_LONGER %>% group_by(density, location, cluster_ID) %>%
  summarise(median_density = median(value)) %>%
  pivot_wider(names_from = cluster_ID, values_from = median_density, names_prefix="cluster_")

make_label <- function(value) {
  x <- as.character(value)
  bquote(x)
}

plot_labeller <- function(variable, value) {
  do.call(expression, lapply(levels(value), make_label))
}

T_cell_boxplots <- DENSITIES_LONGER %>%
  filter(density %in% T_cells) %>%
  mutate(density = ordered(density, levels=T_cells)) %>%
  ggplot(aes(x=cluster_ID, y=value))+
  geom_boxplot(outlier.shape = NA, position=position_dodge(2))+
  facet_grid(location~density, switch = "y", labeller = label_parsed)+
  geom_jitter(data = function(x) dplyr::filter(x, outlier), height=0, width = 0.1,
              alpha=0.5)+
  labs(y=NULL, x=NULL)+
  theme_cowplot(12)

T_cell_subset_boxplots <- DENSITIES_LONGER %>%
  filter(density %in% T_cells_subset) %>%
  mutate(density = ordered(density, levels=T_cells)) %>%
  ggplot(aes(x=cluster_ID, y=value))+
  geom_boxplot(outlier.shape = NA, position=position_dodge(2))+
  facet_grid(location~density, switch = "y", labeller = labeller(density = function(x) gsub("CD3+","",x, fixed = TRUE)))+
  geom_jitter(data = function(x) dplyr::filter(x, outlier), height=0, width = 0.1,
              alpha=0.5)+
  labs(y=NULL, x=NULL)+
  theme_cowplot(12)

B_cell_boxplots <- DENSITIES_LONGER %>%
  filter(density %in% B_cells) %>%
  ggplot(aes(x=cluster_ID, y=value))+
  geom_boxplot(outlier.shape = NA, position=position_dodge(2))+
  facet_grid(location~density, switch = "y")+
  geom_jitter(data = function(x) dplyr::filter(x, outlier), height=0, width = 0.1,
              alpha=0.5)+
  labs(y=NULL, x=NULL)+
  theme_cowplot(12)

Mac_boxplots <- DENSITIES_LONGER %>%
  filter(density %in% Mac) %>%
  ggplot(aes(x=cluster_ID, y=value))+
  geom_boxplot(outlier.shape = NA, position=position_dodge(2))+
  facet_grid(location~density, switch = "y")+
  geom_jitter(data = function(x) dplyr::filter(x, outlier), height=0, width = 0.1,
              alpha=0.5)+
  labs(y=NULL, x=NULL)+
  theme_cowplot(12)

other_boxplot <- DENSITIES_LONGER %>%
  filter(density %in% other) %>%
  ggplot(aes(x=cluster_ID, y=value))+
  geom_boxplot(outlier.shape = NA, position=position_dodge(2))+
  facet_grid(location~density, switch = "y")+
  geom_jitter(data = function(x) dplyr::filter(x, outlier), height=0, width = 0.1,
              alpha=0.5)+
  labs(y=NULL, x=NULL)+
  theme_cowplot(12)

yleft = richtext_grob("log<sub>*10*</sub>(cells/mm<sup>*2*</sup>+0.1)", rot=90)
xbottom = richtext_grob("Cluster")

gridExtra::grid.arrange(grobs = list(T_cell_boxplots, B_cell_boxplots, Mac_boxplots, other_boxplot),
                        layout_matrix = rbind(c(1,1,1,1,1),
                                              c(1,1,1,1,1),
                                              c(1,1,1,1,1),
                                              c(2,2,3,3,4),
                                              c(2,2,3,3,4)),
                        left=yleft,
                        bottom=xbottom)

gridExtra::grid.arrange(grobs = list(T_cell_subset_boxplots, B_cell_boxplots, Mac_boxplots),
                        layout_matrix = rbind(c(1,1,1,1),
                                              c(1,1,1,1),
                                              c(2,2,3,3),
                                              c(2,2,3,3)),
                        left=yleft,
                        bottom=xbottom)


# Which images make up each cluster?
images_per_cluster <- per_image_data %>% select(BT_pNum, BT_Image, ADAPT_Image) %>%
  right_join(DENSITIES_LOGGED_FINAL %>% select(pNum, cluster_ID), by=c("BT_pNum" = "pNum"))

images_per_cluster$BT_Image <- gsub(" - resolution #1", "", images_per_cluster$BT_Image)
images_per_cluster$BT_Image <- paste("C:\\\\Users\\\\Karanvir\\\\Documents\\\\Deeley\\\\COEUR_B&T_Images\\\\BR_COEUR_B&TMOTIF_5Aug21\\\\", images_per_cluster$BT_Image, sep="")

images_per_cluster$ADAPT_Image <- gsub(" - resolution #1", "", images_per_cluster$ADAPT_Image)
images_per_cluster$ADAPT_Image <- paste("E:\\\\COEUR_ADAPT_IMAGES\\\\BR_COEUR_AdaptiveResistance_3Aug21\\\\", images_per_cluster$ADAPT_Image, sep="")


# Pearson Correlation between different cell types per cluster (recerating Boudrette's paper analysis)

# First, reorder the df so that I have all stromal densities followed by all epithelial densities

cor_df <- DENSITIES_LOGGED_FINAL %>% 
  select(
    c(
      "CD3pCD8pPD1n_E",
      "CD3pCD8pPD1p_E",
      "CD3pCD8pFoxP3p_E",
      "CD3pCD8n_E",
      "CD3pCD8nPD1n_E",
      "CD3pCD8nPD1p_E",
      "CD3pCD8nFoxP3p_E",
      "CD68pPDL1n_E",
      "CD68pPDL1p_E",
      "CD20p_E",
      "CD79ApCD20n_E", "CD79ApCD20n_S", "CD20p_S", "CD68pPDL1p_S", 
      "CD68pPDL1n_S", "CD3pCD8nFoxP3p_S", "CD3pCD8nPD1p_S", "CD3pCD8nPD1n_S", 
      "CD3pCD8n_S", "CD3pCD8pFoxP3p_S", "CD3pCD8pPD1p_S", "CD3pCD8pPD1n_S",
      "cluster_ID",
    )
  )

colnames(cor_df) <- gsub("p", "+", colnames(cor_df)) %>%
  gsub("n", "-", .)

cor_df <- cor_df %>% rename("CD79a+CD20+_E" = "CD20+_E",
                            "CD79a+CD20+_S" = "CD20+_S",
                            "CD3+CD8-FoxP3-_E" = "CD3+CD8-_E",
                            "CD3+CD8-FoxP3-_S" = "CD3+CD8-_S")

correlation_cluster1 <- cor(cor_df %>% filter(cluster_ID == 1) %>% select(-cluster_ID))
correlation_cluster2 <- cor(cor_df %>% filter(cluster_ID == 2) %>% select(-cluster_ID))
correlation_cluster3 <- cor(cor_df %>% filter(cluster_ID == 3) %>% select(-cluster_ID))
correlation_cluster4 <- cor(cor_df %>% filter(cluster_ID == 4) %>% select(-cluster_ID))
correlation_cluster5 <- cor(cor_df %>% filter(cluster_ID == 5) %>% select(-cluster_ID))

# Setting a treshold of 0.3 for plotting purposes
correlation_cluster1[abs(correlation_cluster1) < 0.3] <- 0
correlation_cluster2[abs(correlation_cluster2) < 0.3] <- 0
correlation_cluster3[abs(correlation_cluster3) < 0.3] <- 0
correlation_cluster4[abs(correlation_cluster4) < 0.3] <- 0
correlation_cluster5[abs(correlation_cluster5) < 0.3] <- 0

# Average densities in each cluster, so that I can plot the size of the nodes accordingly
average_densities_cluster1 <- apply(cor_df %>% filter(cluster_ID == 1) %>% select(-cluster_ID), 2, mean)
average_densities_cluster2 <- apply(cor_df %>% filter(cluster_ID == 2) %>% select(-cluster_ID), 2, mean)
average_densities_cluster3 <- apply(cor_df %>% filter(cluster_ID == 3) %>% select(-cluster_ID), 2, mean)
average_densities_cluster4 <- apply(cor_df %>% filter(cluster_ID == 4) %>% select(-cluster_ID), 2, mean)
average_densities_cluster5 <- apply(cor_df %>% filter(cluster_ID == 5) %>% select(-cluster_ID), 2, mean)


library(igraph)

#Create networks
net1 <- graph_from_adjacency_matrix(correlation_cluster1,  weighted=T, mode="undirected", diag=F)
net2 <- graph_from_adjacency_matrix(correlation_cluster2,  weighted=T, mode="undirected", diag=F)
net3 <- graph_from_adjacency_matrix(correlation_cluster3,  weighted=T, mode="undirected", diag=F)
net4 <- graph_from_adjacency_matrix(correlation_cluster4,  weighted=T, mode="undirected", diag=F)
net5 <- graph_from_adjacency_matrix(correlation_cluster5,  weighted=T, mode="undirected", diag=F)

#Change colors so positive is red negative is blue
E(net1)$color[E(net1)$weight < 0] <- "blue"
E(net1)$color[E(net1)$weight > 0] <- "red"

E(net2)$color[E(net2)$weight < 0] <- "blue"
E(net2)$color[E(net2)$weight > 0] <- "red"

E(net3)$color[E(net3)$weight < 0] <- "blue"
E(net3)$color[E(net3)$weight > 0] <- "red"

E(net4)$color[E(net4)$weight < 0] <- "blue"
E(net4)$color[E(net4)$weight > 0] <- "red"

E(net5)$color[E(net5)$weight < 0] <- "blue"
E(net5)$color[E(net5)$weight > 0] <- "red"

# This pushes the labels outside the nodes
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}


lab.locs <- radian.rescale(x=1:24, direction=-1, start=0)

# I want to plot the edge width based on the correlation, but I want the widths to be comparable across the 5 plots.
linMap <- function(x, from, to){
  (x -(-0.985116))/(2.668254-(-0.985116)) * (to - from) + from}

#Plots!
#1
set.seed(1234)
p1 <- plot(net1,
           vertex.label.cex=0.8,
           layout = layout_in_circle,
           vertex.label.dist=1,
           vertex.label.degree=lab.locs,
           vertex.label.color="black",
           vertex.label.font=2,
           vertex.color = pal_npg()(5)[1],
           vertex.size=linMap(average_densities_cluster1, 1, 15),
           edge.width = rescale(E(net1)$weight, to=c(0, 5)),
           margin = -0.1,
           edge.color = E(net1)$color,
           main="Cluster 1")

set.seed(1234)
p2 <- plot(net2,
           layout = layout_in_circle,
           vertex.label.cex=0.8,
           vertex.label.dist=1,
           vertex.label.degree=lab.locs,
           vertex.label.color="black",
           vertex.label.font=2,
           vertex.color = pal_npg()(5)[2],
           vertex.size=linMap(average_densities_cluster2, 1, 15),
           edge.width = rescale(E(net2)$weight, to=c(0, 5)),
           margin = -0.1,
           edge.color = E(net2)$color,
           main="Cluster 2")

set.seed(1234)
p3 <- plot(net3,
           vertex.label.cex=0.8,
           layout = layout_in_circle,
           vertex.label.dist=1,
           vertex.label.degree=lab.locs,
           vertex.label.color="black",
           vertex.label.font=2,
           vertex.color = pal_npg()(5)[3],
           vertex.size=linMap(average_densities_cluster3, 1, 15),
           edge.width = rescale(E(net3)$weight, to=c(0, 5)),
           margin = -0.1,
           edge.color = E(net3)$color,
           main="Cluster 3")

set.seed(1234)
p4 <- plot(net4,
           layout = layout_in_circle,
           vertex.label.cex=0.8,
           vertex.label.dist=1,
           vertex.label.degree=lab.locs,
           vertex.label.color="black",
           vertex.label.font=2,
           vertex.color = pal_npg()(5)[4],
           vertex.size=linMap(average_densities_cluster4, 1, 15),
           edge.width = rescale(E(net4)$weight, to=c(0, 5)),
           margin = -0.1,
           edge.color = E(net4)$color,
           main="Cluster 4")

set.seed(1)
p5 <- plot(net5,
           layout = layout.fruchterman.reingold,
           vertex.label.cex=0.8,
           vertex.label.dist=1,
           vertex.label.degree=lab.locs,
           vertex.label.color="black",
           vertex.label.font=2,
           vertex.color = pal_npg()(5)[5],
           vertex.size=linMap(average_densities_cluster5, 1, 15),
           edge.width = rescale(E(net5)$weight, to=c(0, 5)),
           margin = -0.1,
           edge.color = E(net5)$color,
           main="Cluster 5")

library(gridExtra)

# Stepwise differences in correlation matrices between clusters

# Step 1 - How does correlation change between each cluster and the one before it?
delta_2_1 <- cor(cor_df %>% filter(cluster_ID == 2) %>% select(-cluster_ID)) - 
  cor(cor_df %>% filter(cluster_ID == 1) %>% select(-cluster_ID))

delta_3_2 <- cor(cor_df %>% filter(cluster_ID == 3) %>% select(-cluster_ID)) - 
  cor(cor_df %>% filter(cluster_ID == 2) %>% select(-cluster_ID))

delta_4_3 <- cor(cor_df %>% filter(cluster_ID == 4) %>% select(-cluster_ID)) - 
  cor(cor_df %>% filter(cluster_ID == 3) %>% select(-cluster_ID))

delta_5_4 <- cor(cor_df %>% filter(cluster_ID == 5) %>% select(-cluster_ID)) - 
  cor(cor_df %>% filter(cluster_ID == 4) %>% select(-cluster_ID))

delta_5_1 <- cor(cor_df %>% filter(cluster_ID == 5) %>% select(-cluster_ID)) - 
  cor(cor_df %>% filter(cluster_ID == 1) %>% select(-cluster_ID))

delta_2_1[abs(delta_2_1) < 0.2] <- 0 #threshold of 0.2
delta_3_2[abs(delta_3_2) < 0.2] <- 0 #threshold of 0.2
delta_4_3[abs(delta_4_3) < 0.2] <- 0 #threshold of 0.2
delta_5_4[abs(delta_5_4) < 0.2] <- 0 #threshold of 0.2
delta_5_1[abs(delta_5_1) < 0.2] <- 0 #threshold of 0.2

net_2_1 <- graph_from_adjacency_matrix(delta_2_1,  weighted=T, mode="undirected", diag=F)
net_3_2 <- graph_from_adjacency_matrix(delta_3_2,  weighted=T, mode="undirected", diag=F)
net_4_3 <- graph_from_adjacency_matrix(delta_4_3,  weighted=T, mode="undirected", diag=F)
net_5_4 <- graph_from_adjacency_matrix(delta_5_4,  weighted=T, mode="undirected", diag=F)
net_5_1 <- graph_from_adjacency_matrix(delta_5_1,  weighted=T, mode="undirected", diag=F)


E(net1)$color[E(net1)$weight < 0] <- "blue"
E(net1)$color[E(net1)$weight > 0] <- "red"

E(net_2_1)$color[E(net_2_1)$weight < 0] <- "blue"
E(net_2_1)$color[E(net_2_1)$weight > 0] <- "red"

E(net_3_2)$color[E(net_3_2)$weight < 0] <- "blue"
E(net_3_2)$color[E(net_3_2)$weight > 0] <- "red"

E(net_4_3)$color[E(net_4_3)$weight < 0] <- "blue"
E(net_4_3)$color[E(net_4_3)$weight > 0] <- "red"

E(net_5_4)$color[E(net_5_4)$weight < 0] <- "blue"
E(net_5_4)$color[E(net_5_4)$weight > 0] <- "red"

E(net_5_1)$color[E(net_5_1)$weight < 0] <- "blue"
E(net_5_1)$color[E(net_5_1)$weight > 0] <- "red"




p1 <- plot(net1,
           layout = layout_in_circle,
           vertex.label.cex=0.8,
           vertex.label.dist=1,
           vertex.label.degree=lab.locs,
           vertex.label.color="black",
           vertex.label.font=2,
           vertex.color = pal_npg()(5)[1],
           vertex.size=linMap(average_densities_cluster1, 1, 15),
           edge.width = rescale(E(net1)$weight, to=c(0, 10)),
           edge.color = E(net1)$color,
           margin = -0.1,
           main="Cluster 1")


p_2_1 <- plot(net_2_1,
              layout = layout.circle,
              vertex.label.cex=0.8,
              vertex.label.dist=1,
              vertex.label.degree=lab.locs,
              vertex.label.color="black",
              vertex.color = pal_npg()(5)[2],
              vertex.label.font=2,
              vertex.size=linMap(average_densities_cluster2, 1, 15),
              edge.width = rescale(abs(E(net_2_1)$weight), to=c(0, 10)),
              edge.color = E(net_2_1)$color,
              margin = -0.1,
              main="Cluster 2")

p_3_2 <- plot(net_3_2,
              layout = layout.circle,
              vertex.label.cex=0.8,
              vertex.label.dist=1,
              vertex.label.degree=lab.locs,
              vertex.label.color="black",
              vertex.label.font=2,
              vertex.color = pal_npg()(5)[3],
              vertex.size=linMap(average_densities_cluster3, 1, 15),
              edge.width = rescale(abs(E(net_3_2)$weight), to=c(0, 10)),
              edge.color = E(net_3_2)$color,
              margin = -0.1,
              main="Cluster 3")

p_4_3 <- plot(net_4_3,
              layout = layout.circle,
              vertex.label.color="black",
              vertex.label.cex=0.8,
              vertex.label.dist=1,
              vertex.label.degree=lab.locs,
              vertex.label.font=2,
              vertex.color = pal_npg()(5)[4],
              vertex.size=linMap(average_densities_cluster4, 1, 15),
              edge.width = rescale(abs(E(net_4_3)$weight), to=c(0, 10)),
              edge.color = E(net_4_3)$color,
              margin = -0.1,
              main="Cluster 4")

p_5_4 <- plot(net_5_4,
              layout = layout.circle,
              vertex.label.cex=0.8,
              vertex.label.font=2,
              vertex.label.dist=1,
              vertex.label.degree=lab.locs,
              vertex.label.color="black",
              vertex.color = pal_npg()(5)[5],
              vertex.size=linMap(average_densities_cluster5, 1, 15),
              edge.width = rescale(abs(E(net_5_4)$weight), to=c(0, 10)),
              edge.color = E(net_5_4)$color,
              margin = -0.1,
              main="Cluster 5")

p_5_1 <- plot(net_5_1,
              layout = layout.circle,
              vertex.label.cex=0.8,
              vertex.label.font=2,
              vertex.label.dist=1,
              vertex.label.degree=lab.locs,
              vertex.label.color="black",
              vertex.color = pal_npg()(5)[5],
              vertex.size=linMap(average_densities_cluster5, 1, 15),
              edge.width = rescale(abs(E(net_5_1)$weight), to=c(0, 10)),
              edge.color = E(net_5_1)$color,
              margin = -0.1,
              main="Cluster 5")

# Compare network in cluster 5 to the network in the whole dataset
all_cases_cor <- cor(cor_df %>% select(-cluster_ID))
average_densities <- apply(cor_df %>% select(-cluster_ID), 2, mean) #avg. density for size of nodes

all_cases_cor[abs(all_cases_cor) < 0.5] <- 0 #threshold of 0.2

net_all <- graph_from_adjacency_matrix(all_cases_cor,  weighted=T, mode="undirected", diag=F)

E(net_all)$color[E(net_all)$weight < 0] <- "blue"
E(net_all)$color[E(net_all)$weight > 0] <- "red"

p_all <- plot(net_all,
              layout = layout.circle,
              vertex.label.cex=0.8,
              vertex.label.font=2,
              vertex.label.dist=1,
              vertex.label.degree=lab.locs,
              vertex.label.color="black",
              vertex.color = pal_npg()(6)[6],
              vertex.size=linMap(average_densities, 1, 15),
              edge.width = rescale(abs(E(net_all)$weight), to=c(0, 10)),
              edge.color = E(net_all)$color,
              margin = -0.1,
              main="Correlation in entire dataset")

# Difference between network 5 and entire dataset
correlation_cluster5 <- cor(cor_df %>% filter(cluster_ID == 5) %>% select(-cluster_ID))
delta_5_all <- correlation_cluster5 - all_cases_cor

corrplot(delta_5_all, type = "full")

delta_5_all[abs(delta_5_all) < 0.2] <- 0 #threshold of 0.2

net_5_all <- graph_from_adjacency_matrix(delta_5_all,  weighted=T, mode="undirected", diag=F)

E(net_5_all)$color[E(net_5_all)$weight < 0] <- "blue"
E(net_5_all)$color[E(net_5_all)$weight > 0] <- "red"

p_net_5_all <- plot(net_5_all,
                    layout = layout.circle,
                    vertex.label.cex=0.8,
                    vertex.label.font=2,
                    vertex.label.dist=1,
                    vertex.label.degree=lab.locs,
                    vertex.label.color="black",
                    vertex.color = pal_npg()(5)[5],
                    vertex.size=linMap(average_densities_cluster5, 1, 15),
                    edge.width = rescale(abs(E(net_5_all)$weight), to=c(0, 10)),
                    edge.color = E(net_5_all)$color,
                    margin = -0.1,
                    main="Comparing cluster 5 to entire dataset")

# Last but not least, since the main survival differences are between clusters 1/2/3/4 and5
# I want to compare cluster 5 to all the other cases

correlation_cluster5 <- cor(cor_df %>% filter(cluster_ID == 5) %>% select(-cluster_ID))
delta_5_onehot <- correlation_cluster5 - cor(cor_df %>% filter(cluster_ID != 5) %>% select(-cluster_ID))

corrplot(delta_5_onehot)

delta_5_onehot[abs(delta_5_onehot) < 0.2] <- 0 #threshold of 0.2

net_5_onehot <- graph_from_adjacency_matrix(delta_5_onehot,  weighted=T, mode="undirected", diag=F)

E(net_5_onehot)$color[E(net_5_onehot)$weight < 0] <- "blue"
E(net_5_onehot)$color[E(net_5_onehot)$weight > 0] <- "red"

p_net_5_onehot <- plot(net_5_onehot,
                       layout = layout.circle,
                       vertex.label.cex=0.8,
                       vertex.label.font=2,
                       vertex.label.dist=1,
                       vertex.label.degree=lab.locs,
                       vertex.label.color="black",
                       vertex.color = pal_npg()(5)[5],
                       vertex.size=linMap(average_densities_cluster5, 1, 15),
                       edge.width = rescale(abs(E(net_5_onehot)$weight), to=c(0, 10)),
                       edge.color = E(net_5_onehot)$color,
                       margin = -0.1,
                       main="Comparing cluster 5 to clusters 1-2-3-4")


# B cell versus plasma cell densities in cluster 5
# Merge with B&T data
BT <- read.csv("E:/COEUR_mcIF/final_data/COEUR_B&T_final_counts_and_densities.csv") %>% select(1:22) %>% 
  select(-contains(c("Sector", "Row","Column","TMA", "Detections")))
colnames(BT) <- gsub("_t", "_E", colnames(BT)) %>%
  gsub("_s", "_S", .)

BT <- BT %>% group_by(pNum) %>% summarise(across(where(is.numeric), ~sum(.x, na.rm=TRUE)))

BT$total_B <- BT$CD20p_E + BT$CD20p_S
BT$total_plasma <- BT$CD79ApCD20n_E + BT$CD79ApCD20n_S

DENSITIES_LOGGED_FINAL <- left_join(DENSITIES_LOGGED_FINAL, BT %>% select(pNum, total_B, total_plasma))
DENSITIES_LOGGED_FINAL %>% filter(cluster_ID==5) %>% ggplot(aes(x=log10(total_B+0.1), y=log10(total_plasma+0.1)))+
  geom_point()+
  theme_cowplot()

DENSITIES_LOGGED_FINAL<- DENSITIES_LOGGED_FINAL %>% mutate(B_PLASMA_LOG_DIFF = log2(total_B+1) - log2(total_plasma+1))

DENSITIES_LOGGED_FINAL %>% filter(cluster_ID %in% c(3,4,5)) %>% mutate(pNum = factor(pNum), cluster_ID = factor(cluster_ID)) %>% ggplot(aes(x=reorder(pNum, B_PLASMA_LOG_DIFF), y=B_PLASMA_LOG_DIFF, fill=cluster_ID))+
  geom_col()+
  theme_cowplot()+
  theme(axis.ticks.y = element_blank())+
  scale_fill_manual(values=pal_npg()(5)[3:5])+
  scale_x_discrete(labels=NULL)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = -1)+
  labs(x=NULL, y="log2(B cell counts+1) - log2(plasma cell counts+ 1)")+
  coord_flip()


DENSITIES_LOGGED_FINAL %>% filter(cluster_ID %in% c(5)) %>% mutate(pNum = factor(pNum), cluster_ID = factor(cluster_ID)) %>% ggplot(aes(x=reorder(pNum, B_PLASMA_LOG_DIFF), y=B_PLASMA_LOG_DIFF, fill=cluster_ID))+
  geom_col()+
  theme_cowplot()+
  theme(axis.ticks.y = element_blank())+
  scale_fill_manual(values=pal_npg()(5)[5])+
  scale_x_discrete(labels=NULL)+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = -1)+
  labs(x=NULL, y="log2(B cell counts+1) - log2(plasma cell counts+ 1)")+
  coord_flip()

###########################################
############ NEGATIVE CELLS ###############
###########################################





