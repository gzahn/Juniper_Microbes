#####################################################
#          Explore and visualize metadata           #
# This includes chemical and entomological data,    #
# as well as location and age                       #
# Author: Geoff Zahn                                #
#####################################################


# Packages
library(tidyverse)
library(readxl)
library(skimr)

# Load metadata
meta <- read_xlsx("./metadata/metadata.xlsx")

# First look
glimpse(meta)
skim(meta)

# Extra tidying and prep ####

# Fix column classes
meta$Living_Larvae <- as.numeric(meta$Living_Larvae)
meta$Raw_Exit_Holes <- as.numeric(meta$Raw_Exit_Holes)

# Find analytical chemistry columns
chem_cols <- grep(pattern = "CHEM_",names(meta))

# Add total chemical (relative area %)
meta$ChemTotal <- rowSums(meta[,chem_cols])

# Find entomology columns
names(meta)
ento_cols <- c("Bolt_Surface_Area_cm2","Raw_Exit_Holes_per_cm2","Raw_Exit_Holes","Living_Larvae")

# Export cleaned and prepped metadata ####
write_csv(meta,"./metadata/cleaned_metadata.csv")

