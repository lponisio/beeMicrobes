## This script is designed to clean and prep bee trait data to then calculate
##  site-level bee trait diversity metrics.
## For the beeMicrobes project it has been modified from the Sky Islands script
##  to include data from the sky islands, sunflower and oregon survey projects.

# clear environment
rm(list=ls())

# load necessary libraries/packages
library(FD)
library(vegan)
library(tidyverse)
# set working directory to home
setwd("~/")
## run the lab_paths script to appropriately set path for your computer (must be set up for each lab member)
## once amended with your computer's information, lab_paths.R can be saved to your home directory
source("lab_paths.R")

## generate path from your folder to the BMMA folder
dir.beeMicrobes <- file.path(local.path, "beeMicrobes")
## change working directory
setwd(dir.beeMicrobes)

# Load in home-cooked functions
## calculate trait uniqueness and originality based on Coux et al. 2016 for each site-year-round
source("src/calcFuncUniqOrig.R")
source("src/misc.R")

# change working directory to the folder containing our data
setwd(file.path(local.path, "beeMicrobes_saved"))

# Start by loading in trait data from each project and collating it
# read in raw trait data
si.traits <- read.csv("./beeMicrobes_data/traits/SI_bee_traits.csv") %>% mutate(dataset = rep("SI", nrow(.)))
sf.traits <- read.csv("./beeMicrobes_data/traits/SF_bee_traits.csv") %>% mutate(dataset = rep("SF", nrow(.)))
or.traits <- read.csv("./beeMicrobes_data/traits/OR_bee_traits2023.csv") %>% mutate(dataset = rep("OR", nrow(.)))

# do an oversimplified bind rows to merge data
merged.traits.raw <- bind_rows(si.traits, sf.traits, or.traits)

# replace all blank cells with NA values for appropriate handling in R
merged.traits.raw <- merged.traits.raw %>%
  mutate(across(everything(), ~ {
    x <- as.character(.x)
    x[!grepl("[[:alnum:]]", x)] <- NA
    x
  }))

# Add column to flag duplicate species rows
merged.traits.raw$sp.dup <- merged.traits.raw$GenusSpecies %in% merged.traits.raw$GenusSpecies[duplicated(merged.traits.raw$GenusSpecies)]

# write trait data to csv for manual data cleaning in R (deduplicating, etc.)
# ONLY NEED TO RUN ONCE, SINCE FIX WAS MANUAL
#write.csv(merged.traits.raw, "../beeMicrobes_saved/beeMicrobes_data/traits/merged_bee_traits_WithDuplicates.csv", row.names = F)

### rectified duplicate species information manually in excel
# Load fixed trait information back into R
merge.traits.2 <- read.csv("./beeMicrobes_data/traits/merged_bee_traits_dupsFixed.csv") %>%
  mutate(GenusSpecies = trimws(GenusSpecies)
         )
merge.traits.2$sp.dup <- merge.traits.2$GenusSpecies %in% merge.traits.2$GenusSpecies[duplicated(merge.traits.2$GenusSpecies)]

write.csv(merge.traits.2, "~/Desktop/merged_bee_traitsFixed_final.csv", row.names = F)
### clean a few more elements of the data
# unify language in different columns (decide on language for one version of factor level
#   make all other versions the same: Rent & Renter  -> renter)
# deal with upper/lowercase issues (turn to lowercase)
traits.m <- read.csv("./beeMicrobes_data/traits/merged_bee_traitsFixed_final.edits.csv") %>%
  mutate(GenusSpecies = case_when(grepl("sp.", x = GenusSpecies) == T ~ paste0(GenusSpecies, "_", dataset),
                                  .default = GenusSpecies),
         across(-GenusSpecies, tolower),
         NestConstruction = case_when(grepl("excav", NestConstruction) == T ~ "excavator",
                                      grepl("rent", NestConstruction) == T ~ "renter",
                                      .default = NestConstruction),
         PrimaryNestMaterial = case_when(grepl("nomada", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("triepeolus", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("coelioxys", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("stelis", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("melecta", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("coelioxys", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("sphecodes", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         .default = PrimaryNestMaterial),
         NestPartitions = case_when(grepl("nomada", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("triepeolus", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("coelioxys", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("stelis", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("melecta", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("coelioxys", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         grepl("sphecodes", x = GenusSpecies, ignore.case = T) == T ~ "cleptoparasite",
                                         .default = NestPartitions),
         NestPartitions = case_when(NestPartitions == "wa" ~ "wax",
                                    NestPartitions == "m" ~ "mud",
                                    NestPartitions == "l" ~ "leaf",
                                    NestPartitions == "w" ~ "woody pith",
                                    NestPartitions == "p" ~ "cleptoparasite",
                                    NestPartitions == "c" ~ "cellophane",
                                    NestPartitions == "r" ~ "resin",
                                    NestPartitions == "l,r" ~ "cellophane",
                                    NestPartitions == "l,p" ~ "leaf, cleptoparasite",
                                    NestPartitions == "l,m" ~ "leaf, mud",
                                    NestPartitions == "r,g,m" ~ "resin, gum, mud",
                                    NestPartitions == "l,w,r" ~ "leaf, woody pith, resin",
                                    NestPartitions == "l,w,g" ~ "leaf, woody pith, gum",
                                    .default = NestPartitions),
         PrimaryNestMaterial = case_when(PrimaryNestMaterial == "#n/a" ~ NA,
                                         .default = PrimaryNestMaterial),
         NestLocation = case_when(is.na(NestLocation) == T ~ NestLoc,
                                  .default = NestLocation),
         NestLocation = case_when(NestLocation == "below ground" ~ "below-ground",
                                  NestLocation == "below" ~ "below-ground",
                                  NestLocation == "above ground" ~ "above-ground",
                                  NestLocation == "above/below ground" ~ "above/below-ground",
                                  NestLocation == "mixed" ~ "above/below-ground",
                                  NestLocation == "mix" ~ "above/below-ground",
                                  NestLocation == "above" ~ "above-ground",
                                  grepl("diadasia", x = GenusSpecies) == T ~ "below-ground",
                                  grepl("megachile prosopidis", x = GenusSpecies) == T ~ "above-ground",
                                  .default = NestLocation),
         NestConstruction = case_when(is.na(NestConstruction) == T ~ Excavate,
                                      .default = NestConstruction),
         NestConstruction = case_when(grepl("constr", NestConstruction, ignore.case = T) == T ~ "constructor",
                                      grepl("rent", NestConstruction, ignore.case = T) == T ~ "renter",
                                      grepl("ex", NestConstruction, ignore.case = T) == T ~ "excavator",
                                      .default = NestConstruction),
         Sociality = case_when(Sociality == "semisocial" ~ "social", # theres only one semisocial bee... lets just call it social
                               Sociality == "social parasite" ~ "brood parasite",
                               .default = Sociality),
         GenusSpecies = paste0(toupper(substr(GenusSpecies, 1, 1)),
                                       substr(GenusSpecies, 2, nchar(GenusSpecies))), # make first letter uppercase
         Genus = word(GenusSpecies, 1),
         GenusSpecies = str_remove(GenusSpecies, "_SI"),
         GenusSpecies = str_remove(GenusSpecies, "_SF"),
         MeanITD = as.numeric(MeanITD)
         ) %>%
  select(GenusSpecies, Genus, NestLocation, PrimaryNestMaterial,NestConstruction, NestPartitions,
         Sociality, MeanITD, Lecty, c(d:median.abund))

# generate a genus-mean of ITD, this can be used to fill in any missing values
genus.mean.itd <- traits.m %>%
  group_by(Genus) %>%
  summarize(GenusMeanITD = mean(MeanITD, na.rm=TRUE))

# fill missing values with the genus estimated mean
traits.m$MeanITD[is.na(traits.m$MeanITD)] <- 
  genus.mean.itd$GenusMeanITD[match(traits.m$Genus[is.na(traits.m$MeanITD)], genus.mean.itd$Genus)]

traits.filt <- traits.m %>%
  filter(GenusSpecies %in% spec.net$GenusSpecies) # filter to only the species for which we have microbe data

## To get the code running we'll assume missing lecty is polylectic for now
traits.filt$Lecty[is.na(traits.filt$Lecty) == T] <- "polylectic"

# assign binomials as row names
rownames(traits.filt) <- traits.filt$GenusSpecies

# Get rid of binomials (they are safely stored in the rownames)
traits.filt$GenusSpecies <- NULL

## We grouped the traits based on what they are related to. For
## example, we have 3 nesting-related traits, so weighted each one by
## 1/3. This will avoid overweighting any trait type.

bee.traits <- c("NestLocation", "PrimaryNestMaterial",
                "NestConstruction","NestPartitions","Sociality",
                "MeanITD")

bee.weights <- c(rep(1/4, 4), rep(1, 2))

## check on unique values... no NA values!
apply(traits.filt[, bee.traits], 2, unique)

traits.matrix <- traits.filt[, bee.traits]

## *****************************************************************
## species-level function diversity metrics
## *****************************************************************
## calculates data for two new columns we are adding to database,
## "originality" and "uniq"
bee.func <- calcFuncUniqOrig(traits.matrix,
                            traits.2.keep=bee.traits,
                            weights=bee.weights,
                            type="all spec")

# remake genus species column from rownames for functional traits data objects
bee.func$GenusSpecies <- rownames(bee.func)
traits.matrix$GenusSpecies <- rownames(traits.matrix)
# remove rownames
rownames(traits.matrix) <- NULL
rownames(bee.func) <- NULL

# merge functional uniqueness to traits matrix
traits.matrix <- merge(traits.matrix, bee.func, by="GenusSpecies")

# new list of bee traits
all.traits <- colnames(traits.matrix)[!colnames(traits.matrix) %in%
                                        c("GenusSpecies")]

## save prepped data
spec.net <- cbind(spec.net, traits.matrix[, all.traits][match(spec.net$GenusSpecies,
                                                       traits.matrix$GenusSpecies),])
spec.net$r.degree <- as.numeric(spec.net$r.degree)

## write.csv(traits, file="../data/traits.csv",
##           row.names=FALSE)

## save(spec.net, file="../data/spec_traits.Rdata")
## write.csv(spec.net, file="../data/spec_traits.csv", row.names=FALSE)

## *****************************************************************
## site-level function diversity metrics
## *****************************************************************
### NOTE: THE SITE-LEVEL FUNCTIONAL DIVERSITY R FUNCTION APPEARS TO GIVE THE SAME RESULTS AS THE SP-LEVEL FUNCTION.
###       FOR NOW, NOT USING THESE DATA.
# # make a singular aggregate site/date/sample-round grouping variable
# spec.net$SiteYearSr <- paste(spec.net$Site, spec.net$Date,
#                              sep=";")
# ## only working with bees
# bee.families <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae",
#                   "Megachilidae")
# unique(spec.net$SiteYearSr) %>% length()
# 
# # Generate community matrix for counts of each bee on each site/date/sampling round
# comms <- spec.net %>%
#   filter(Family %in% bee.families,# subset to bee families (no flies, etc.)
#          is.na(GenusSpecies) == F,
#          GenusSpecies != "") %>% # filter out missing genera (shouldn't be any, but just in case)
#   select(GenusSpecies, SiteYearSr) %>% # select only columns for binomial and site/year/sr aggregate
#   pivot_longer(cols = -SiteYearSr) %>% # ensure data is in longform before flipping to wide-form
#   group_by(SiteYearSr,value) %>% summarise(N=n()) %>% # generate counts of each species occurence in the dataset
#   pivot_wider(names_from = value, values_from=N) %>% # flip data to wide-form (community matrix format)
#   replace(is.na(.),0) %>% ungroup() %>% # ensure all NAs are counted as 0s
#   as.data.frame(.) %>% 
#   column_to_rownames(var = "SiteYearSr")
# 
# comms <- comms[, order(colnames(comms))]
# colSums(comms) %>% sort()
# 
# ## check that all the species in the community matrix are in traits
# colnames(comms)[!colnames(comms) %in% rownames(traits.matrix)]
# 
# # Move binomial into row names
# traits.matrix$GenusSpecies <- rownames(traits.matrix)
# traits.matrix$GenusSpecies <- NULL
# 
# # calculate bee functional uniqueness
# ##TODO: Need to troubleshoot "functionally singular" species issue
# ##  function not estimating FEVe, FDis, or FRic. Function only results in 2 columns:
# ##  "uniq" and "originality", which are necessary in then next step
# site.lvl.func <- calcFuncUniqOrigComm(traits = traits.matrix,
#                                       traits.2.keep=bee.traits,
#                                       weights=bee.weights,
#                                       type="all spec",
#                                       a = comms,
#                                       w.abun=TRUE
#                                       )
# 
# ## merge the site-year-round level functional diversity metrics
# site.fd <- data.frame(BeeFDis=bee.func$fd$FDis,
#                       BeeFEve=bee.func$fd$FEve,
#                       SiteYearSr=names(bee.func$fd$FDis))
# rownames(fd) <- NULL
# 
# spec.net <- merge(spec.net, fd, all.x=TRUE,
#                   by="SiteYearSr")
# dim(spec.net)
# 
# spec.net <- merge(spec.net, bee.func$by.comm.mets, all.x=TRUE)
# dim(spec.net)
# 
# 
# save(spec.net, file='../data/spec_net_fdiv.Rdata')
