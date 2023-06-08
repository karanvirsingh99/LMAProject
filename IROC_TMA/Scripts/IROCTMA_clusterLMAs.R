##################################################
## Project: LMA Project
## Script purpose: Cluster IROC TMA LMAs into distinct types
## Author: Karanvir Singh

# Load spiat functions
source("/Users/karanvir/Documents/LMAProject/COEUR/Scripts/Neighborhood Analysis/SPIAT_functions.R")

# Load per-patient-densities and per-image counts
per_patient_densities <- vroom::vroom("/Users/karanvir/Documents/LMAProject/IROC_TMA/Results/IROC-tma-per-patient-densities.csv")
per_image_counts <- vroom::vroom("/Users/Karanvir/Documents/LMAProject/IROC_TMA/Results/IROC-tma-per-image-counts.csv")

# Load neighborhoods (SPIAT output) and bind rows
load("/Users/karanvir/Documents/LMAProject/IROC_TMA/Results/Neighborhood-analysis/RData/IROCTMA_15r_15n_SPIAT.RData")
neighborhoods <- bind_rows(neighborhoods, .id="Image")

# Add patient identifiers to neighborhoods
neighborhoods <- left_join(neighborhoods, per_image_counts %>% select(Image, ID))

#Cluster LMAs
IROC_LMAs <- clusterNeighborhoods(neighborhoods, 7)

# Save neighborhood summary spreadsheet
IROC_LMAs %>% write_tsv("/Users/karanvir/Documents/LMAProject/IROC_TMA/Results/Neighborhood-analysis/7LMAs_IROC_TMA.csv")
