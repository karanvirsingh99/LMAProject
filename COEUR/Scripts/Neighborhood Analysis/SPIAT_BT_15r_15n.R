# Run SPIAT neighborhood detection on BT COEUR images

library(data.table)
library(tidyverse)
library(SPIAT)

# Import functions for neighborhood detection
source("/Users/karanvir/Documents/LMAProject/COEUR/Scripts/Neighborhood Analysis/SPIAT_functions.R")

# Import cell level data (BT panel)
cell_level_data <- vroom::vroom("/Users/karanvir/Documents/LMAProject/COEUR/Final_Data/BT_cell_data_fixed_images_with_pNum.csv")

# Phenotypes of interest
phenotypes <- c("CD3pCD8p", "CD3pCD8n", "CD3pCD8nFoxP3p", "CD3pCD8pFoxP3p", 
  "CD20p", "CD79ApCD20n")

# Run neighborhood detection with 15 cells in neighborhoods, max.distance is 15um
neighborhoods <- getNeighborhoods(cell_level_data = cell_level_data, radius=15, n_size=15, cell_types=phenotypes)

# Save as an .RData file for future analysis
save(neighborhoods, file="BT_15r_15n_SPIAT.RData")
