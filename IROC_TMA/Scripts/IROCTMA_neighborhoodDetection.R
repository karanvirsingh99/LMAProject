##################################################
## Project: LMA Project
## Script purpose: Run neighborhood detection on IROC TMA images, radius 15 and min_cell = 15
## Author: Karanvir Singh


library(data.table)
library(tidyverse)
library(SPIAT)

# Import functions for neighborhood detection
source("/Users/karanvir/Documents/LMAProject/COEUR/Scripts/Neighborhood Analysis/SPIAT_functions.R")

# Import cell level data (BT panel)
cell_level_data <- vroom::vroom("/Users/karanvir/Documents/LMAProject/IROC_TMA/Final_Data/23Mar23_IROC_TMA_CellLevelData.csv.gz")

# Phenotypes of interest
phenotypes <- c("CD3pCD8pPD1p", "CD68pPDL1n", "CD3pCD8pPD1n", 
                "CD68pPDL1p",  "PD1p", "CD20p", "CD79ApCD20n", 
                "CD3pCD8nPD1n", "CD3pCD8nPD1p")

# Run neighborhood detection with 15 cells in neighborhoods, max.distance is 15um
neighborhoods <- getNeighborhoods(cell_level_data = cell_level_data, radius=15, n_size=15, cell_types=phenotypes)

# Save as an .RData file for future analysis
save(neighborhoods, file="/Users/karanvir/Documents/LMAProject/IROC_TMA/Results/Neighborhood-analysis/RData/IROCTMA_15r_15n_SPIAT.RData")
