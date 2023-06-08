 ##################################################
  ## Project: IROC WHOLE SLIDE STAINING
  ## Script purpose: Import cell level and annotation level data for IROC whole slide staining, fix cell type names and 
  ##                 identify outliers. 
  ## Date: 07June2023
  ## Author: Karanvir Singh
  ##################################################
 

library(tidyverse)
 
# Import cell level data
 
cellLevel <- vroom::vroom("IROC_WHOLE_SLIDES/raw_data/cell_level_export_8jun23.tsv.gz")

# Split classification into two components

cellLevel$Classifier1 <- stringr::str_split_i(cellLevel$Class, pattern=": ", i=1)
cellLevel$Classifier2 <- stringr::str_split_i(cellLevel$Class, pattern=": ", i=2)

# I accidentally coded the PD1 as "+/-" instead of "p/n". Fixing that here
cellLevel$Classifier2 <- gsub("\\+", "p", cellLevel$Classifier2)
cellLevel$Classifier2 <- gsub("\\-", "n", cellLevel$Classifier2)

# Creating a composite phenotype column by adding PD1p/PD1n to T cell phenotypes,
# and leaving other cell types as is

Tcells <- c("CD3pCD8p", "CD3pCD8n")

cellLevel <- cellLevel %>% mutate(Class = case_when(Classifier1 %in% Tcells ~ paste0(Classifier1, Classifier2),
                                                    Classifier1=="Negative" & Classifier2=="PD1p" ~ "CD3nPD1p",
                                                      TRUE ~ Classifier1))

cellLevel %>% select(Image, Class, Parent, Centroid_x, Centroid_y, area) %>%
  write_tsv(gzfile("~/Documents/LMAProject/IROC_WHOLE_SLIDES/raw_data/08Jun23_clean_cell_level.tsv.gz"))
