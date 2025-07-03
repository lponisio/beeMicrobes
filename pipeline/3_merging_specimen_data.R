######################
## Load in packages ##
######################
library(dplyr)
library(tidyr)
library(vctrs)

########################
## Load in dataframes ##
########################

## reset then change working directory to the 2025 Bee Microbial Meta Analysis folder (BMMA2025)
setwd("~/")
## run the lab_paths script to appropriately set path for your computer (must be set up for each lab member)
## once amended with your computer's information, lab_paths.R can be saved to your home directory
source("lab_paths.R")

## generate path from your folder to the BMMA folder
dir.meta.microbe <- file.path(local.path, "BMMA2025")
## change working directory
setwd(dir.meta.microbe)

## read in the 2022 Sky Island data, clean up columns and change some column classes
si.spec <- read.csv("data/SI2022_spec.csv", stringsAsFactors=FALSE) %>%
  select(-X) %>%
  mutate(WindEnd = as.numeric(WindEnd))
## read in the 2019 SunFlower data
## also rename the parasite columns to match those in the Sky Islands dataset
sun.spec <- read.csv("data/SF2019_spec.csv", stringsAsFactors=FALSE) %>%
  rename("ApicystisSpp" = "Apicystis",
         "AscosphaeraSpp" = "Ascosphaera",
         "AspergillusSpp" = "Aspergillus")

###############################################
## merge specimen data columns and pare down ##
###############################################

spec.combined <- vec_rbind(si.spec, sun.spec)
spec.full.trim <- spec.combined %>%
  filter(Order == "Hymenoptera",
         rowSums(is.na(select(., Apidae:NosemaBombi))) < 9)%>%
  select(where(~ any(!is.na(.))),
         -c(Method,) ) 
dim(spec.full.trim)
