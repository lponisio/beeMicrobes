## This script matches tree tips to bee samples in our dataset

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

dir.microbes <- file.path(local.path, "beeMicrobes")
setwd(dir.microbes)

# load packages
library(tidyverse)
library(gridExtra)
## trees
library(phytools)
library(ggtree)


source("../skyIslands/dataPrep/src/misc.R")
load('data/allNetSums.RData')


## Load tree from : Henriquez Piskulich, Patricia Andrea; Hugall,
## Andrew F.; Stuart-Fox, Devi (2023).  A supermatrix phylogeny of the
## world’s bees (Hymenoptera: Anthophila) [Dataset].
## Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

# Read in the downloaded phylogeny
big.bee.tree <- read.newick("./data/BEE_mat6b_fulltree.nwk")

# Clean tree tip labels for matching to our species list
clean_tips.df <- big.bee.tree$tip.label %>%
  as.data.frame() %>%
  rename("tip.labels" = ".") %>%
  mutate(tip.labels1 = str_split_fixed(big.bee.tree$tip.label, pattern = "_", n = 3)[,1],
         tip.labels2 = str_split_fixed(big.bee.tree$tip.label, pattern = "_", n = 3)[,2],
         GenSp_phylo = paste0(tip.labels1, " ", tip.labels2)
  )

clean_tips <- clean_tips.df$GenSp_phylo


# check which species don't have a representative tip in the phylogeny
no.tips <- unique(all.indiv.mets$GenusSpecies[!all.indiv.mets$GenusSpecies %in% clean_tips.df$GenSp_phylo])
no.tips %>% base::sort()

# Clean the species list and make a "map"
# GenusSpecies column is the original list, GS_clean includes changes to match the downloaded tree
## assigning species not in the big bee tree to their closest relatives
## This will generally be the most well-observed member of their subgenus on iNaturalist
##  if possible it will be a taxon within the western united states

species_df <- unique(trimws(all.indiv.mets$GenusSpecies)) %>%
  as.data.frame() %>%
  rename("GenusSpecies" = ".") %>%
  mutate(Genus = str_split_fixed(GenusSpecies, " ", n = 3)[,1],
         GenSp_phylo = str_remove(GenusSpecies, pattern = " cf."), # remove " cf."
         GenSp_phylo = case_when(GenusSpecies == "Megachile gemula fulvogemula" ~ "Megachile gemula",
                              GenusSpecies == "Megachile melanophaea rohweri" ~ "Megachile melanophaea",
                              GenusSpecies == "Lasioglossum pacatum/viridatulum" ~ "Lasioglossum pacatum",
                              GenusSpecies == "Xylocopa californica arizonensis" ~ "Xylocopa californica",
                              GenusSpecies == "Xylocopa tabaniformis orpifex" ~ "Xylocopa tabaniformis",
                              GenusSpecies == "Svastra obliqua expurgata" ~ "Svastra obliqua",
                              GenusSpecies == "Melissodes lupina" ~ "Melissodes lupinus",
                              GenusSpecies == "Melissodes communis alopex" ~ "Melissodes communis",
                              GenusSpecies == "Melissodes tepida timberlakei" ~ "Melissodes tepidus",
                              GenusSpecies == "Agapostemon obliquus" ~ "Agapostemon texanus",
                              GenusSpecies == "Andrena apacheorum" ~ "Andrena denticulata",
                              GenusSpecies == "Andrena birtwelli" ~ "Andrena fulva",
                              GenusSpecies == "Andrena pecosana" ~ "Andrena asteris",
                              GenusSpecies == "Andrena tegularis" ~ "Andrena asteris",
                              GenusSpecies == "Andrena segregans" ~ "Andrena astragali",
                              GenusSpecies == "Andrena hitei" ~ "Andrena milwaukeensis",
                              GenusSpecies == "Andrena livida" ~ "Andrena candida",
                              GenusSpecies == "Ceratina dallatorreana" ~ "Ceratina cyanea",
                              GenusSpecies == "Hylaeus conspicuus" ~ "Hylaeus communis",
                              GenusSpecies == "Protandrena porterae" ~ "Protandrena andrenoides",
                              GenusSpecies == "Syntrichalonia exquisita" ~ "Eucera actuosa",
                              GenusSpecies == "Lasioglossum diversopunctatum" ~ "Lasioglossum incompletum",
                              GenusSpecies == "Lasioglossum ruficorne" ~ "Lasioglossum foxii",
                              GenusSpecies == "Megachile comata" ~ "Megachile gemula",
                              GenusSpecies == "Megachile subexilis" ~ "Megachile campanulae",
                              GenusSpecies == "Panurginus cressoniellus" ~ "Panurginus crawfordi",
                              GenusSpecies == "Peponapis pruinosa" ~ "Eucera pruinosa",
                              GenusSpecies == "Pseudopanurgus albitarsis" ~ "Protandrena andrenoides",
                              GenusSpecies == "Pseudopanurgus atricornis" ~ "Protandrena andrenoides",
                              GenusSpecies == "Pseudopanurgus bakeri" ~ "Protandrena illinoiensis",
                              GenusSpecies == "Pseudopanurgus porterae" ~ "Protandrena andrenoides",
                              GenusSpecies == "Dufourea tinsleyi" ~ "Dufourea dilatipes",
                              GenusSpecies == "Svastra obliqua expurgata" ~ "Svastra obliqua",
                              Genus == "Triepeolus" ~ "Triepeolus utahensis",
                              Genus == "Sphecodes" ~ "Sphecodes ranunculi",
                              grepl("Hylaeus morpho", GenusSpecies) == T ~ "Hylaeus communis",
                              grepl("Osmia sp", GenusSpecies) == T ~ "Osmia pentstemonis",
                              grepl("Svastra sp", GenusSpecies) == T ~ "Svastra obliqua",
                              grepl("Ashmeadiella sp", GenusSpecies) == T ~ "Ashmeadiella opuntiae",
                              grepl("Lasioglossum sp", GenusSpecies) == T ~ "Lasioglossum sisymbrii",
                              .default = GenSp_phylo)
         ) %>%
  left_join(clean_tips.df %>% select(tip.labels, GenSp_phylo),
            by = "GenSp_phylo")

# check which species don't have a representative tip in the phylogeny
unique(species_df$GenSp_phylo[!species_df$GenSp_phylo %in% clean_tips.df$GenSp_phylo])
# none! We've addressed all of the missing taxa to the correct names or their nearest neighbors

species_vec_clean <- species_df$GenSp_phylo

# Match tree tips to our species list
keep <- clean_tips %in% species_vec_clean
sum(keep)

# Prune
pruned_tree <- keep.tip(big.bee.tree, big.bee.tree$tip.label[keep])

# visualize your tree make sure there are a reasonable number of tips
plot(pruned_tree)

# generate a phylogenetic relatedness matrix
co.var.mat <- ape::vcv.phylo(pruned_tree)

# phylogenetic distance matrix
bee.dist <- cophenetic.phylo(pruned_tree)

save(co.var.mat, pruned_tree, file="data/covarmatrix_community.Rdata")
save(bee.dist, pruned_tree, file = "data/bee_dist.RData")
write.csv(species_df, file = "data/phylo_names.csv", row.names = F) # the key to match tips and species

## ----------------------------------------------------------------------- ##
## Final cleaning of the bee traits dataframe before the phylogeny script  ##
## ----------------------------------------------------------------------- ##
# Read in a fresh version of the bee traits
# Some traits columns were removed upstream, presently it's just safer to re-merge the traits...
bee.traits <- c("NestLocation", "PrimaryNestMaterial", "NestConstruction",
                "NestPartitions", "Sociality", "MeanITD")
traits.clean <- read.csv("../beeMicrobes_saved/beeMicrobes_data/traits/merged_traits_clean.csv") %>%
  select(GenusSpecies, all_of(bee.traits))

# read in the naming "map"
# This links sample species ID's to their representative tip names on the phylogeny
phylo_names <- read.csv("./data/phylo_names.csv")

# List of bee families
bee.fams <- c("Apidae", "Andrenidae", "Halictidae", "Megachilidae", "Colletidae")

# left_join the original traits data to fill in NA values
all.indiv.mets <- all.indiv.mets %>%
  select(-any_of(bee.traits)) %>%
  left_join(traits.clean, by = "GenusSpecies") %>% # merge in the trait data
  left_join(phylo_names, by = "GenusSpecies") %>%
  mutate(Genus = word(GenusSpecies, 1),
         Apis = case_when(Genus == "Apis" ~ "Apis",
                          .default = "not Apis"), # generate the apis/not apis column
         PossibleParasite = 5,
         Family = case_when(Genus == "Hylaeus" ~ "Colletidae",
                            .default =  Family),
         Nest = paste(NestPartitions, NestLocation)
  ) %>%
  filter(Family %in% bee.fams, # filter any families not in the main 5 we collected
         is.na(Micro_originality) == F,
         tip.labels %in% rownames(co.var.mat) # filter to only bees on the phylogeny
  )

## add parasite richness, Crithidia spp. combined, NA if unscreened
para <- as.data.frame(lapply(all.indiv.mets[, c("AscosphaeraSpp", "ApicystisSpp", "NosemaCeranae",
                                                "NosemaBombi", "CrithidiaBombi", "CrithidiaExpoeki",
                                                "CrithidiaSpp")],
                             function(x) as.numeric(as.character(x))))
crith <- c("CrithidiaBombi", "CrithidiaExpoeki", "CrithidiaSpp")
para$Crithidia <- ifelse(rowSums(!is.na(para[, crith])) == 0, NA,
                         as.numeric(rowSums(para[, crith] > 0, na.rm = TRUE) > 0))
para <- para[, c("AscosphaeraSpp", "ApicystisSpp", "Crithidia", "NosemaCeranae", "NosemaBombi")]
all.indiv.mets$ParasiteRichness <- ifelse(rowSums(!is.na(para)) > 0,
                                          rowSums(para > 0, na.rm = TRUE), NA)
table(all.indiv.mets$ParasiteRichness, useNA = "ifany")

all.indiv.mets <- all.indiv.mets[!is.na(all.indiv.mets$ParasiteRichness),]

save(#all.sp.mets.binary,
  #all.sp.mets,
  all.indiv.mets.binary,
  all.indiv.mets,
  file = './data/allNetSums_final.RData')
