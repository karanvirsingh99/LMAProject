##################################################
## Project: LMA Project
## Script purpose: To add patient identifiers to annotation level data to IROC TMA QuPath data and calculate densities
##                  as well as stroma percentages
## Author: Karanvir Singh


library(tidyverse)

# Import annotation level data
IROC_DENSITIES <- read.csv("/Users/karanvir/Documents/LMAProject/IROC_TMA/Final_Data/annotation_level_data_28Jul22.csv")

#Delete annotations that are not tissue annotations (remnants of me trying to take screenshots of the cores..)
IROC_DENSITIES <- IROC_DENSITIES %>% filter(Name != "PathAnnotationObject")

#Rename columns for clarity
colnames(IROC_DENSITIES) <- gsub('Num.', '', colnames(IROC_DENSITIES)) %>%
  gsub('\\..', '', .)

#Delete columns that are not of interest
IROC_DENSITIES <- IROC_DENSITIES %>% select(-Centroidm, -Centroidm, -PD1nbase., -Perimeterm, -PD1pbase., -PD1nPositive,
                                            -PD1pPositive, -PD1nD1pPositive)

#Change to long format
IROC_DENSITIES <- IROC_DENSITIES %>% pivot_longer(cols=c(7:13, 15:23), names_repair ="unique")

#Merge PD1p and PD1n phenotypes when they are not of interest/technical artifacts from spillover
IROC_DENSITIES <- IROC_DENSITIES %>% mutate(name = case_when(name == "PD1nNegative" ~ "Negative",
                                                             name %in% c("PD1nCD68pPDL1p", "PD1pCD68pPDL1p") ~ "CD68pPDL1p",
                                                             name %in% c("PD1nCD68pPDL1n", "PD1pCD68pPDL1n") ~ "CD68pPDL1n",
                                                             name %in% c("PD1nCKpPDL1p", "PD1pCKpPDL1p") ~ "CKpPDL1p",
                                                             name %in% c("PD1nCD79ApCD20n", "PD1pCD79ApCD20n") ~ "CD79ApCD20n",
                                                             name %in% c("PD1nCD79ApCD20p", "PD1pCD79ApCD20p") ~ "CD20p",
                                                             name == "PD1nCD3pCD8n" ~ "CD3pCD8nPD1n",
                                                             name == "PD1nCD3pCD8p" ~ "CD3pCD8pPD1n",
                                                             name == "PD1pCD3pCD8n" ~ "CD3pCD8nPD1p",
                                                             name == "PD1pCD3pCD8p" ~ "CD3pCD8pPD1p",
                                                             name == "PD1pNegative" ~ "PD1p"))

IROC_DENSITIES$value[is.na(IROC_DENSITIES$value)] <- 0

# Add the same phenotype counts together
IROC_DENSITIES <- IROC_DENSITIES %>% group_by(Image, Class, name, Aream) %>% summarise(value = sum(value))

# Change Tumor and Stroma to _E and _S
IROC_DENSITIES <- IROC_DENSITIES %>% mutate(name = case_when(Class == "Tumor" ~ paste(name, '_E', sep=""),
                                                             Class == "Stroma" ~ paste(name, '_S', sep="")))

# Change to wide format
IROC_DENSITIES_WIDE <- IROC_DENSITIES %>% pivot_wider(id_cols = "Image", names_from = name, values_from = value)

# Calculate areas and merge back to wide df
area_str <- IROC_DENSITIES %>% ungroup() %>% filter(Class=="Stroma") %>%  group_by(Image) %>% summarise(Area = unique(Aream))
area_epi <- IROC_DENSITIES %>% ungroup() %>% filter(Class=="Tumor") %>% group_by(Image) %>% summarise(Area = unique(Aream))

IROC_DENSITIES_WIDE$area_stroma <- (area_str$Area[match(IROC_DENSITIES_WIDE$Image, area_str$Image)])/1e6
IROC_DENSITIES_WIDE$area_epi <- (area_epi$Area[match(IROC_DENSITIES_WIDE$Image, area_epi$Image)])/1e6

# Add patient identifier

mid <- function(text, start, n) {
  substr(text, start, start + n - 1)
}

IROC_DENSITIES_WIDE$SEC_ROW_COL <- mid(IROC_DENSITIES_WIDE$Image, 23, 6)
IROC_DENSITIES_WIDE$SEC_ROW_COL <- gsub("]", "" ,IROC_DENSITIES_WIDE$SEC_ROW_COL)

# Key of core and ID pairs
key <- read.csv("/Users/karanvir/Documents/LMAProject/IROC_TMA/Final_Data/IROCPROS6_PNUM.csv")
key$SEC_ROW_COL <- paste(key$Sector, key$Row, key$Column, sep=",")

#combine
IROC_DENSITIES_WIDE$ID <- key$ID[match(IROC_DENSITIES_WIDE$SEC_ROW_COL, key$SEC_ROW_COL)]

# Filter out control cores and one endometrial case
IROC_DENSITIES_WIDE <- IROC_DENSITIES_WIDE %>% filter(!(ID %in% c("Liver", "Placenta", "FT", "198(en)")))

# Save per-image spreadsheet
IROC_DENSITIES_WIDE  %>% write_csv("IROC_TMA/Results/IROC-tma-per-image-counts.csv")
  
# Summarize on per-patient basis
per_patient_densities <- IROC_DENSITIES_WIDE %>% group_by(ID) %>% summarise(across(where(is.numeric), sum))

#Calculate densities
per_patient_densities <- per_patient_densities %>% mutate(across(ends_with("_S"), ~.x/area_stroma)) %>%
  mutate(across(ends_with("_E"), ~.x/area_epi)) 

# Calculate stroma proportion
per_patient_densities <- per_patient_densities %>% mutate(stroma_proportion = area_stroma/(area_stroma+area_epi))

# Save as new .csv
per_patient_densities %>% write_csv("IROC_TMA/Results/IROC-tma-per-patient-densities.csv")

# #log densities
# log_iroc <- per_patient_densities %>% select(-Negative_S, -Negative_E, -area_stroma, -area_epi) %>% 
#   mutate(across(where(is.numeric), ~log10(.+0.1)))
# 
# #move ID to rownames
# log_iroc <- log_iroc %>% select(-PD1p_S, -PD1p_E, -CKpPDL1p_S, -CKpPDL1p_E) %>%
#   column_to_rownames("ID")
