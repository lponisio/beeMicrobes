# The objective of this script is to collate existing specimen data together from across the 3 included studies
# The unaltered example script of how metrics were calculated and aggregated will remain here from the sky islands project.
# Toward the bottom of the script we will use existing dataframes that will be available from *DATA REPOSITORY X*, so
#   starting the script from the top is not necessary. All calculations and references to original raw data files are
#   available on github at the following urls
# Sky Islands: https://github.com/lponisio/skyIslands/blob/master/dataPrep/4_specimenPrep.R
# Sunflowers: https://github.com/lponisio/sunflower/blob/master/dataPrep/4_specimenPrep.R
# Oregon Survey: URL HERE (WILL BE ADDED WHEN READY)

library(igraph)
library(vegan)
library(fields)
library(fossil)
library(bipartite)
library(dplyr)
library(RSQLite)
library(tidyverse)
library(readr)
library(vctrs)
library(qiime2R)
library(phyloseq)
setwd("~/")
source("lab_paths.R")
local.path

## generate path from your folder to the beeMicrobes folder
beemicrobes <- file.path(local.path, "beeMicrobes_saved")
setwd(beemicrobes)

## read in the 2022 Sky Island data, clean up columns and change some column classes
si.spec <- read.csv("beeMicrobes_data/specimen_data/SI2022_spec.csv", stringsAsFactors=FALSE) %>%
  select(-X, -starts_with("Pan_")) %>%
  mutate(study_name = rep("SkyIslands", nrow(.)),
         WindEnd = as.numeric(WindEnd),
         across(where(is.character),
                ~ na_if(trimws(.), ""))) %>% # remove empty characters ("")
  rename("FloralAbundance" = "MeanFloralAbundance",
         "FloralRichness" = "MeanFloralRichness",
         "FloralDiversity" = "MeanFloralDiversity")

## read in the 2019 Sun Flower data
## also rename the parasite columns to match those in the Sky Islands dataset
sf.spec <- read.csv("beeMicrobes_data/specimen_data/SF2019_spec.csv", stringsAsFactors=FALSE) %>%
  mutate(study_name = rep("Sunflowers", nrow(.)),
         across(where(is.character),
                ~ na_if(trimws(.), ""))) %>% # remove empty characters ("")
  rename("ApicystisSpp" = "Apicystis",
         "AscosphaeraSpp" = "Ascosphaera",
         "AspergillusSpp" = "Aspergillus",
         "BeeDiversity" = "BeeDivShannon",
         "BeeRichness" = "Richness",
         "BeeAbundance" = "TotalAbundance",
         "FloralDiversity" = "FloralDiv"
         )

###############################################
## merge specimen data columns and pare down ##
###############################################

# combine the two datasets, don't worry how ugly they get when merged, we'll deal with that downstream
spec.combined <- vctrs::vec_rbind(si.spec, sf.spec) %>%
  mutate(UniqueID.fixed = case_when( Year == 2019 ~ str_remove(UniqueID, "^SF"),  # remove SF from 2019 sample names
                                     .default = UniqueID)
         ) %>%
  filter(is.na(GenusSpecies) == F)

## Read in QIIME2 Data and merge into phyloseq object
physeq.16s <- qza_to_phyloseq(
  features = "./beeMicrobes_pipeline_output/merged/16s/tablef4.qza",
  tree =     "./beeMicrobes_pipeline_output/merged/16s/rooted-tree16s.qza",
  taxonomy = "./beeMicrobes_pipeline_output/merged/16s/taxonomy16s.qza",
  #metadata = "./beeMicrobes_pipeline_output/merged/16s/HM2020_full_metadata_7-5-22.tsv"
  )

# generate a list of intersecting samples that appear in both specimen data and qiime2 samples
# For background, not all samples in the specimen data had gut microbes processed
# Additionally, some microbial samples were filtered out during bioinformatic processing
sample.metadata.intersection <- intersect(spec.combined$UniqueID.fixed, sample_names(physeq.16s))

# Use the intersecting list to prune the samples down to size so we can add the specimen data as phyloseq metadata
physeq.trim <- prune_samples(sample_names(physeq.16s) %in% sample.metadata.intersection, physeq.16s)

# use the list of sample names to filter the specimen data as the first step in cleaning it
spec.filtered <- spec.combined %>%
  filter(UniqueID.fixed %in% sample_names(physeq.trim)) %>% # filter the sample names down to those in the phyloseq object
  select(GenusSpecies, Site, Year, Date, SampleRound,Doy, Lat, Long, UniqueID.fixed, UniqueID, Sex, c(Apidae:NosemaBombi),
         Family, Genus, SubGenus, Species, PlantGenus, PlantSpecies, PlantGenusSpecies,
         Country, State, County, BeeAbundance, BeeRichness, BeeDiversity, FloralAbundance,
         FloralRichness, FloralDiversity)

# need to recalculate presence/absence columns for pathogens
# also bee species-level parasitism rates.
# maybe also calculate for higher taxa? (Genus, or even family)

## *****************************************************************
## Reorder and add spec data into phyloseq metadata
## *****************************************************************

# Set the rownames of spec to match the sample names/UniqueIDs
rownames(spec.filtered) <- spec.filtered$UniqueID.fixed

# reorder rows of spec data to match order of the phyloseq object
spec.net <- spec.filtered[match(sample_names(physeq.trim), rownames(spec.filtered)), ]

# add the spec data into the phyloseq object
sample_data(physeq.trim) <- sample_data(spec.reordered)

## *****************************************************************
## create relational database, add species IDs
## *****************************************************************

setwd(dir.bombus)
source('dataPrep/relational/1prep.R')

setwd(dir.bombus)
source('dataPrep/relational/2make.R')

setwd(dir.bombus)
source('dataPrep/relational/3join.R')

## *****************************************************************
## prep specimen data
## *****************************************************************

rm(list=ls())
setwd("~/")
source("lab_paths.R")
setwd(file.path(local.path, "skyIslands_saved"))

spec <-
    read.csv('data/relational/relational/traditional/specimens-complete.csv',
             stringsAsFactors=FALSE)

setwd("../skyIslands/dataPrep")
source("src/misc.R")
source("src/prepNets.R")
source("src/specialization.R")

spec$GenusSpecies <- fix.white.space(paste(spec$Genus,
                          spec$Species,
                          spec$SubSpecies))

spec$PlantGenusSpecies <-  fix.white.space(paste(spec$PlantGenus,
                                          spec$PlantSpecies,
                                          spec$PlantVar))

spec$Int <-  fix.white.space(paste(spec$GenusSpecies,
                                   spec$PlantGenusSpecies))

spec$Date <- as.Date(spec$Date, format='%Y-%m-%d')
spec$Doy <- as.numeric(strftime(spec$Date, format='%j'))
spec$Year <- as.numeric(format(spec$Date,'%Y'))

## all possible sample round, site, combos of net data, used for veg
## data subsequently
collections <- data.frame(unique(cbind(spec$Site, spec$Year,
                                       spec$SampleRound)))
colnames(collections) <- c("Site", "Year", "SampleRound")

## calculate orthoganol polynomials for doy
DoyPoly <- poly(spec$Doy, degree=2)
spec$DoyPoly1 <- DoyPoly[,'1']
spec$DoyPoly2 <- DoyPoly[,'2']

## also for Latitude
##spec$LatPoly <- poly(spec$Lat, degree=2)
##spec$LatPoly1 <- spec$LatPoly[,'1']
##spec$LatPoly2 <- spec$LatPoly[,'2']
##spec$LatPoly <- NULL

## check plant names
## library(Taxonstand)
## checked.plant.names <- TPL(id(spec$PlantGenusSpecies))

## write.csv(checked.plant.names,
## file="../../skyIslands_saved/data/checks/plant_names_check.csv")

## update plant names
checked.plant.names <-
    read.csv(file="../../skyIslands_saved/data/checks/normalized_plants.csv")

spec <- fixPlantNamesgBIF(spec, "PlantGenusSpecies",
                     checked.plant.names)

## DECISION POINT:  variabile identification fo Erigeron between years,
##  combine to Erigerson spp.

spec$PlantGenusSpecies[spec$PlantGenus == "Erigeron"] <-
    "Erigeron spp."

## *****************************************************************
## specimen-level parasite calculations
## *****************************************************************
dir.create("../data/", showWarnings = FALSE)
dir.create("../data/networks", showWarnings = FALSE)
dir.create("../data/splevel_network_metrics", showWarnings = FALSE)

crithidias <- c("CrithidiaExpoeki",
                "CrithidiaBombi", "CrithidiaSpp",
                "CrithidiaMellificae")
nosemas <- c("NosemaBombi", "NosemaCeranae")
parasites <- c( "AscosphaeraSpp",
               "ApicystisSpp", crithidias, nosemas)

spec[, parasites][is.na(spec[, parasites])] <- 0
spec[, parasites][spec[, parasites] == ""] <- 0
spec[, parasites] <- apply(spec[, parasites], 2,  as.numeric)

spec[spec$Apidae != 1 | is.na(spec$Apidae), parasites] <- NA

spec$ParasiteRichness <- rowSums(spec[, parasites],
                                 na.rm=TRUE)
spec$CrithidiaRichness <- rowSums(spec[, crithidias],
                                  na.rm=TRUE)
spec$NosemaRichness <- rowSums(spec[, nosemas],
                                 na.rm=TRUE)
spec$PossibleParasite <- apply(spec[, parasites], 1,
                               function(x) sum(!is.na(x)))
spec$ParasitePresence <- (spec$ParasiteRichness >= 1)*1
spec$CrithidiaPresence <- (spec$CrithidiaRichness >= 1)*1
spec$NosemaPresence <- (spec$NosemaRichness >= 1)*1

spec[spec$Apidae != 1  | is.na(spec$Apidae), "ParasiteRichness"] <- NA
spec[spec$Apidae != 1  | is.na(spec$Apidae), "ParasitePresence"] <- NA
spec[spec$Apidae != 1  | is.na(spec$Apidae), "CrithidiaRichness"] <-
    NA
spec[spec$Apidae != 1  | is.na(spec$Apidae), "NosemaRichness"] <-
  NA
spec[spec$Apidae != 1  | is.na(spec$Apidae), "CrithidiaPresence"] <- NA
spec[spec$Apidae != 1  | is.na(spec$Apidae), "NosemaPresence"] <- NA

check.spec <- spec[!is.na(spec$Apidae),]
check.spec <- check.spec[check.spec$GenusSpecies == "",]

write.csv(check.spec,
 file="../../skyIslands_saved/data/checks/screened_no_ID.csv",
 row.names=FALSE)


write.csv(spec, file="../data/spec_all_methods.csv", row.names=FALSE)

spec.net <- spec[spec$Method == "Net",]
spec.pan <- spec[spec$Method == "Pan",]
spec <- spec[spec$Method != "Vane",]

write.csv(spec, file="../data/spec_net_pan.csv", row.names=FALSE)

## net.only.columns <- c("PlantGenus", "PlantGenusSpecies",
##                       "PlantSpecies", "PlantSubSpecies", "PlantVar",
##                       "NetNumber", parasites, "Apidae",
##                       "AspergillusSpp", "PlantFamily",
##                       "ParasiteRichness", "PossibleParasite",
##                       "ParasitePresence")

## pan.only.columns <- c("PanColor", "PanLocation")

## spec.pan <- spec.pan[, !colnames(spec.pan) %in% net.only.columns]
## spec.net <- spec.net[, !colnames(spec.net) %in% pan.only.columns]

## write.csv(spec.net, file="../data/spec_net.csv", row.names=FALSE)
## write.csv(spec.pan, file="../data/spec_pan.csv", row.names=FALSE)

## ***********************************************************************
## site/species level insect data
## ***********************************************************************


bee.families <- c("Andrenidae", "Apidae", "Colletidae", "Halictidae",
                  "Megachilidae")

## DECISION POINT: drop non bees but keep the syrphids
# spec <- spec[spec$Family %in% c(bee.families, "Syrphidae"),]


calcSummaryStats <- function(spec.method, method){
    site.sp <- spec.method %>%
        group_by(Site, Year, SampleRound, GenusSpecies) %>%
        summarise(Abundance = length(GenusSpecies),
                  SpParasitismRate = mean(ParasitePresence, na.rm=TRUE),
                  SpCrithidiaParasitismRate = mean(CrithidiaPresence, na.rm=TRUE),
                  SpApicystisParasitismRate = mean(ApicystisSpp, na.rm=TRUE),
                  SpCrithidiaBombusParasitismRate = mean(CrithidiaPresence[Genus == "Bombus"],
                                                      na.rm=TRUE),
                  SpCrithidiaHBParasitismRate = mean (CrithidiaPresence
                                                   [GenusSpecies == "Apis mellifera"],
                                                   na.rm=TRUE),
                  SpParasitism = sum(ParasitePresence, na.rm=TRUE),
                  SpCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
                  SpApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
                  SpNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
                  SpNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
                  SpAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
                  SpCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
                  SpCrithidiaMellificae = sum(CrithidiaMellificae, na.rm=TRUE),
                  SpCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
                  SpCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
                  SpScreened = sum(!is.na(Apidae))
                  )

    site.sum <- spec.method %>%
        group_by(Site, Year, SampleRound) %>%
        summarise(PollAbundance = length(GenusSpecies),
                  BeeAbundance = length(GenusSpecies[Family %in%
                                                     bee.families]),
                  SyrphidAbundance =
                      length(GenusSpecies[Family == "Syrphidae"]),
                  HBAbundance = sum(GenusSpecies == "Apis mellifera"),
                  BombusAbundance = sum(Genus == "Bombus"),
                  NonBombusHBAbundance =
                      sum(Genus != "Bombus" &
                          Genus != "Apis" &
                          Family %in% bee.families),

                  PollRichness= length(unique(GenusSpecies)),
                  BeeRichness= length(unique(GenusSpecies[Family %in%
                                                          bee.families])),
                  SyrphidRichness= length(unique(
                    GenusSpecies[Family  == "Syrphidae"])),
                  BombusRichness= length(unique(
                    GenusSpecies[Genus == "Bombus"])),

                  PollDiversity=vegan:::diversity(table(GenusSpecies),
                                                  index="shannon"),
                  BeeDiversity=vegan:::diversity(table(
                    GenusSpecies[Family %in% bee.families]),
                                                 index="shannon"),
                  SyrphidDiversity=vegan:::diversity(table(
                    GenusSpecies[Family  == "Syrphidae"]),
                                                     index="shannon"),
                  BombusDiversity=vegan:::diversity(table(
                    GenusSpecies[Genus == "Bombus"]),
                                                    index="shannon"),

                  VisitedFloralRichness= length(unique(PlantGenusSpecies)),
                  VisitedFloralDiversity=vegan:::diversity(table(
                    PlantGenusSpecies),  index="shannon"),

                  SiteParasitismRate=mean(ParasitePresence, na.rm=TRUE),
                  MeanParasiteRichness=mean(ParasiteRichness, na.rm=TRUE),
                  HBSiteParasitismRate=mean(
                      ParasitePresence[GenusSpecies == "Apis mellifera"],
                      na.rm=TRUE),
                  BombusSiteParasitismRate=mean(
                      ParasitePresence[Genus == "Bombus"],
                    na.rm=TRUE),
                  SRDoyPoly1=mean(DoyPoly1),
                  SRDoyPoly2=mean(DoyPoly2),
                  SRDoy=mean(Doy),
                  CrithidiaParasitismRate=mean(CrithidiaPresence, na.rm=TRUE),
                  ApicystisParasitismRate=mean(ApicystisSpp, na.rm=TRUE),
                  CrithidiaBombusParasitismRate= mean(CrithidiaPresence[Genus == "Bombus"],
                                                      na.rm=TRUE),
                  CrithidiaHBParasitismRate= mean (CrithidiaPresence
                                                   [GenusSpecies == "Apis mellifera"],
                                                   na.rm=TRUE))

    site.sp.yr <- spec %>%
        group_by(Site, Year, GenusSpecies, Genus) %>%
        summarise(Abundance = length(GenusSpecies))

    site.sp.yr.round <- spec %>%
      group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%
      summarise(AbundanceSYR = mean(length(GenusSpecies)))

    bombus <- site.sp.yr[site.sp.yr$Genus == "Bombus",]
    bombus$Genus  <- NULL

    ## write species-level sumary data
    write.csv(bombus, file=sprintf('../data/bombus_year_site_%s.csv',
                                   method),
              row.names=FALSE)
    write.csv(site.sp.yr, file=sprintf('../data/sp_year_site_%s.csv',
                                       method),
              row.names=FALSE)
    write.csv(site.sp, file=sprintf('../data/spstats_%s.csv', method),
              row.names=FALSE)
    write.csv(site.sp.yr.round, file='../data/sp_year_site_round.csv',
              row.names=FALSE)
    return(site.sum)
}

net.site.sum <- calcSummaryStats(spec.net, "net")
pan.site.sum <- calcSummaryStats(spec.pan, "pan")
all.site.sum <- calcSummaryStats(spec, "all")

table(all.site.sum$Site)
table(net.site.sum$Site)
table(pan.site.sum$Site)

sum.cols <- c("PollAbundance", "BeeAbundance", "SyrphidAbundance",
              "HBAbundance", "BombusAbundance",
              "NonBombusHBAbundance", "PollRichness", "BeeRichness",
              "SyrphidRichness", "BombusRichness", "PollDiversity",
              "BeeDiversity", "SyrphidDiversity", "BombusDiversity",
              "VisitedFloralRichness", "VisitedFloralDiversity",
              "SiteParasitismRate", "MeanParasiteRichness",
              "HBSiteParasitismRate", "BombusSiteParasitismRate")

colnames(net.site.sum)[colnames(net.site.sum) %in% sum.cols] <-
  paste0("Net_", sum.cols)
colnames(pan.site.sum)[colnames(pan.site.sum) %in% sum.cols] <-
  paste0("Pan_", sum.cols)

## remove parasite columns from pan data because there will never be
## anything but NAs
pan.site.sum[, grepl("Parasit", colnames(pan.site.sum))] <- NULL

## drop day columns that might cause merge issues
dup.cols <- c("SRDoyPoly1", "SRDoyPoly2", "SRDoy")

pan.site.sum[, dup.cols] <- NULL
net.site.sum[, dup.cols] <- NULL

dim(all.site.sum)
site.sum <- merge(all.site.sum, net.site.sum, all.x=TRUE)
dim(site.sum)
site.sum <- merge(site.sum, pan.site.sum, all.x=TRUE)
dim(site.sum)


## add back site characteristics
sites <- unique(data.frame(Site=spec$Site,
                           Lat= spec$Lat,
                           Area=spec$Area,
                           Elev=spec$Elev,
                           Year=spec$Year))
site.sum <- merge(site.sum, sites)
site.sum$Year <- as.factor(site.sum$Year)

## write the site, year, sampling round summary data after merging
## with plant data


## ***********************************************************************
## For multi-study parasite summaries, a bit redundant but nice to
## have in a standardized format
## ***********************************************************************

par.site.sp <- spec.net %>%
    group_by(Site, Year, SampleRound, GenusSpecies, Genus) %>%
    summarise(SpAbundance = length(GenusSpecies),
              SpParasitism = sum(ParasitePresence, na.rm=TRUE),
              SpCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
              SpApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
              SpNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
              SpNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
              SpAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
              SpCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
              SpCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
              SpCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
              SpScreened = sum(!is.na(Apidae))
              )

par.site.gen <- spec.net %>%
    group_by(Site, Year, SampleRound, Genus) %>%
    summarise(GenusAbundance = length(GenusSpecies),
              GenusParasitism = sum(ParasitePresence, na.rm=TRUE),
              GenusCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
              GenusApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
              GenusNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
              GenusNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
              GenusAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
              GenusCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
              GenusCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
              GenusCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
              GenusScreened = sum(!is.na(Apidae))
              )

par.site.sum <- spec.net %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(SiteBeeAbundance = length(GenusSpecies[Family %in%
                                                     bee.families]),

              SiteBeeRichness= length(unique(GenusSpecies[Family %in%
                                                          bee.families])),

              SiteBeeDiversity=vegan:::diversity(table(
                                           GenusSpecies[Family %in% bee.families]),
                                           index="shannon"),

              SiteParasitism = sum(ParasitePresence, na.rm=TRUE),
              SiteCrithidiaPresence= sum(CrithidiaPresence, na.rm=TRUE),
              SiteApicystisSpp = sum(ApicystisSpp, na.rm=TRUE),
              SiteNosemaBombi= sum(NosemaBombi, na.rm=TRUE),
              SiteNosemaCeranae = sum(NosemaCeranae, na.rm=TRUE),
              SiteAscosphaeraSpp= sum(AscosphaeraSpp, na.rm=TRUE),
              SiteCrithidiaExpoeki = sum(CrithidiaExpoeki, na.rm=TRUE),
              SiteCrithidiaBombi = sum(CrithidiaBombi, na.rm=TRUE),
              SiteCrithidiaSpp = sum(CrithidiaSpp, na.rm=TRUE),
              SiteScreened = sum(!is.na(Apidae))
              )


all.sums <- left_join(par.site.sp, par.site.gen)

all.sums <- left_join(all.sums, par.site.sum)
all.sums$Project <- "SI"

save(all.sums, file="../data/SI_parasite_sums.Rdata")


## *******************************************************************
## create a giant plant-pollinator network to calculate specialization
## etc. across all SI
## ***************************************** **************************
spec.net.nets <- spec.net[!is.na(spec.net$PlantGenusSpecies),]

agg.spec <- aggregate(list(abund=spec.net.nets$GenusSpecies),
                      list(GenusSpecies=spec.net.nets$GenusSpecies,
                           PlantGenusSpecies=spec.net.nets$PlantGenusSpecies),
                      length)

nets.all <- samp2site.spp(agg.spec$PlantGenusSpecies,
                          agg.spec$GenusSpecies,
                          agg.spec$abund, FUN=sum)

all.traits <- specieslevel(nets.all)
## calculate rarified plant.pol degree
rare.plants.degree <- apply(nets.all, 1, chao1)
rare.pols.degree <- apply(nets.all, 2, chao1)

traits <- data.frame(GenusSpecies= unlist(sapply(all.traits,
                                                 rownames)),
                     do.call(rbind, all.traits))

traits$r.degree <-  rare.pols.degree[match(traits$GenusSpecies,
                                           names(rare.pols.degree))]
traits$r.degree[is.na(traits$r.degree)] <-
    rare.plants.degree[match(traits$GenusSpecies[is.na(traits$r.degree)],
                             names(rare.plants.degree))]

rownames(traits) <- NULL

write.csv(traits, file='../data/networks_traits.csv', row.names=FALSE)

## *******************************************************************
## create a giant pathogen-pollinator network to calculate
## specialization etc. across all SI
## *******************************************************************
agg.spec.sub <- spec[spec$Apidae == 1 & !is.na(spec$Apidae),]

agg.spec.para <- aggregate(agg.spec.sub[, parasites],
                      list(GenusSpecies=agg.spec.sub$GenusSpecies),
                                              sum, na.rm=TRUE)

para.gensp.counts <- table(agg.spec.sub$GenusSpecies)
para.gensp.counts

para.gen.counts <- table(agg.spec.sub$Genus)
para.gen.counts

## proportion of individuals screened
agg.spec.para[, parasites] <- agg.spec.para[, parasites]/
    para.gensp.counts[agg.spec.para$GenusSpecies]

nets.para <- agg.spec.para
rownames(nets.para) <- nets.para$GenusSpecies
nets.para$GenusSpecies <- NULL

all.traits.para <- specieslevel(nets.para)
## calculate rarified plant.pol degree

para.traits <- data.frame(GenusSpecies= unlist(sapply(all.traits.para,
                                                 rownames)),
                     do.call(rbind, all.traits.para))

rownames(para.traits) <- NULL

write.csv(traits, file='../data/parasitetraits.csv', row.names=FALSE)

## *******************************************************************
##  create site, SR, year level networks for plant-pollinators,
##  pollinator-parasites
## *******************************************************************

## plant-pollinator networks
bees.syr.yr.sr <- makeNets(spec.net.nets, net.type="YearSR",
                           poll.group="BeesSyrphids")
bees.syr.yr <- makeNets(spec.net.nets, net.type="YearSR", mean.by.year=TRUE,
         poll.group="BeesSyrphids")

bees.yr.sr <- makeNets(spec.net.nets[spec.net.nets$Family %in% bee.families,],
         net.type="YearSR",
         poll.group="Bees")

bees.yr <- makeNets(spec.net.nets[spec.net.nets$Family %in% bee.families,],
         net.type="YearSR", mean.by.year=TRUE, poll.group="Bees")

spec.sub <- agg.spec.sub %>%
  dplyr::select(UniqueID, GenusSpecies, Site, Year, SampleRound, AscosphaeraSpp,
         ApicystisSpp, CrithidiaExpoeki, CrithidiaBombi,
         CrithidiaMellificae,
         NosemaBombi, NosemaCeranae)

prep.para <- spec.sub %>%
  pivot_longer(cols=c("AscosphaeraSpp",
                      "ApicystisSpp", "CrithidiaExpoeki",
                      "CrithidiaBombi", "CrithidiaMellificae",
                      "NosemaBombi", "NosemaCeranae"),
               names_to = "Parasite", values_to = "count")
prep.para <- as.data.frame(prep.para)

par.bees.yr.sr <- makeNets(prep.para, net.type="YearSR",
                           species=c("Pollinator",
                                     "Parasite"),
                           lower.level="GenusSpecies",
                           higher.level="Parasite",
                           poll.group="Bees")


par.bees.yr <- makeNets(prep.para, net.type="YearSR",
                        species=c("Pollinator",
                                  "Parasite"),
                        lower.level="GenusSpecies",
                        higher.level="Parasite",
                        mean.by.year=TRUE,
                        poll.group="Bees")

## merge site summary metrics
print("specimen dim before network metric year-sr merge")
dim(spec.net)

bees.yr.sr <- bees.yr.sr[bees.yr.sr$speciesType == "higher.level",]
bees.yr.sr$SpSiteYear <- NULL

spec.net <- merge(spec.net, bees.yr.sr, all.x=TRUE)
dim(spec.net)
print("specimen dim after network metric year-sr merge")

spec.net$speciesType <- NULL

## *******************************************************************
##  Data checks
## *******************************************************************

print(paste("Pollinator species", length(unique(spec$GenusSpecies))))
print(paste("Bee species", length(unique(spec$GenusSpecies[
                                                spec$Family %in%
                                                bee.families]))))
print(paste("Syrphid species", length(unique(spec$GenusSpecies[
                                                spec$Family == "Syrphidae"]))))
print(paste("Plant species", length(unique(spec$PlantGenusSpecies))))
print(paste("Pollinator genera", length(unique(spec$Genus))))
print(paste("Interactions", length(unique(spec$Int))))
print(paste("Specimens", nrow(spec)))

## table(spec$GenusSpecies)
## tab <- table(spec$GenusSpecies, spec$Site)
## tab2 <- table(spec$PlantGenusSpecies, spec$Site)
## table(spec$PlantGenusSpecies)
## table(spec$PlantGenusSpecies, spec$Site)
## table(spec$PlantGenusSpecies, spec$Year)
## table(spec$PlantGenusSpecies, spec$Family)

## spec <- spec[spec$Sex == "f",]
## tab <- table(spec$GenusSpecies, spec$Site, spec$Year)

## all.occ <- apply(tab, 1, sum)
## sp.too.few <- names(all.occ)[all.occ < 5]

## sp.possible <- names(all.occ)[all.occ > 5]
## cleptos <-  c("Sphecodes sp. a", "Sphecodes sp. b", "Triepeolus sp. a")

## sp.possible <- sp.possible[!sp.possible %in% cleptos]
## sp.sub <- spec[spec$GenusSpecies %in% sp.possible,]

## tab.sub <- sp.sub  %>%
##     group_by(GenusSpecies, Site, Year, SampleRound) %>%
##     summarise(n = n())

##  tab.sub$n30p <- ceiling(tab.sub$n*0.3)
##  tab.sub$n20p <- ceiling(tab.sub$n*0.2)

## sum(tab.sub$n30p)
## sum(tab.sub$n20p)

## tab.sub.cutoff <- tab.sub

## tab.sub.cutoff$n[tab.sub.cutoff$n > 5] <- 5
## sum(tab.sub.cutoff$n)

## *******************************************************************
##  Veg and bloom cleaning
## *******************************************************************

bloom <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/bloom-complete.csv",
             stringsAsFactors=FALSE)
veg <-
    read.csv("../../skyIslands_saved/data/relational/relational/traditional/veg-complete.csv",
             stringsAsFactors=FALSE)

veg$PlantGenusSpecies <-  fix.white.space(paste(veg$PlantGenus,
                                          veg$PlantSpecies,
                                          veg$PlantVar))

## people in the field put NA as a placehold where there were no
## flowers to express data was collected, just nothing was
## there. Not needed for site level calculations.
veg <- veg[veg$PlantGenusSpecies != "",]

bloom$PlantGenusSpecies <-  fix.white.space(paste(bloom$PlantGenus,
                                          bloom$PlantSpecies,
                                          bloom$PlantVar))
bloom <- bloom[bloom$PlantGenusSpecies != "",]

## check plant names, run july 2024
scientificName <- unique(c(id(bloom$PlantGenusSpecies), id(veg$PlantGenusSpecies)))
write.csv(scientificName,
          file="../../skyIslands_saved/data/checks/all_plants.csv",
          row.names=FALSE)

## use gbif backbone  https://www.gbif.org/tools/species-lookup

## update plant names
checked.plant.names <-
    read.csv(file="../../skyIslands_saved/data/checks/normalized_plants.csv")

veg <- fixPlantNamesgBIF(veg, "PlantGenusSpecies",
                     checked.plant.names)

bloom <- fixPlantNamesgBIF(bloom, "PlantGenusSpecies",
                     checked.plant.names)

id(veg$PlantGenusSpecies)
id(bloom$PlantGenusSpecies)

## fix dates
veg$Date <- as.Date(veg$Date,  "%Y-%m-%d")
veg$Year <- format(veg$Date,  "%Y")

bloom$Date <- as.Date(bloom$Date,  "%Y-%m-%d")
bloom$Year <- format(bloom$Date,  "%Y")

## 2017-2022 quad options
quads.2017.2022 <- unique(veg$Quadrat)
quads.2012 <- unique(veg$Quadrat[veg$Year == 2012])
## quads not done in 2012
not.in.2012 <- quads.2017.2022[!quads.2017.2022 %in% quads.2012]

## prep a matrix of all the quads for site richness etc.
combos <- collections[rep(seq_len(nrow(collections)),
                          each = length(quads.2017.2022)), ]
combos$Quadrat <- quads.2017.2022
## drop the quads not done in 2012
combos <- combos[!(combos$Year == 2012 & combos$Quad %in%
                   not.in.2012),]

## blooming flowers
veg.blooming <- veg[veg$NumBlooms != "" &
                    veg$NumBlooms != "0",]

## didn't count the exact number of flowers in 2017-2018, use
## midpoints of bins
veg.blooming$NumBloomsNum <- veg.blooming$NumBlooms
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "<10"] <- 5
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "10-100"] <- 50
veg.blooming$NumBloomsNum[veg.blooming$NumBlooms == "100-1000"] <- 100

veg.blooming$NumBloomsNum <- as.numeric(veg.blooming$NumBloomsNum)

## *******************************************************************
##  veg site summaries as mean of quads
## *******************************************************************

## by quadrat for a site average
veg.bloom.quad.sp <- veg.blooming %>%
    group_by(Quadrat, Site, Year, PlantGenusSpecies, SampleRound) %>%
    summarise(FloralAbundance = sum(NumBloomsNum),
              FloweringPlantAbundance= sum(PlantCount, na.rm=TRUE))

veg.bloom.quad.sp$KeyQuad <- paste(veg.bloom.quad.sp$Site,
                              veg.bloom.quad.sp$Quadrat,
                              veg.bloom.quad.sp$Year,
                              veg.bloom.quad.sp$SampleRound)

veg.div.quad <- tapply(veg.bloom.quad.sp$FloralAbundance,
                  veg.bloom.quad.sp$KeyQuad,
                  vegan::diversity)

veg.plant.div.quad <- tapply(veg.bloom.quad.sp$FloweringPlantAbundance,
                  veg.bloom.quad.sp$KeyQuad,
                  vegan::diversity)


veg.bloom.sum.quad <- veg.bloom.quad.sp %>%
    group_by(Quadrat, Site, Year, SampleRound) %>%
    summarise(FloralAbundance = sum(FloralAbundance),
              FloralRichness= length(unique(PlantGenusSpecies)),
              FloweringPlantAbundance = sum(FloweringPlantAbundance)
              )

veg.name <-   names(veg.div.quad)
names(veg.div.quad) <- NULL

## add div data into site summary data
veg.bloom.sum.quad$FloralDiversity <-  veg.div.quad[match(
    paste(veg.bloom.sum.quad$Site,
          veg.bloom.sum.quad$Quadrat,
          veg.bloom.sum.quad$Year,
          veg.bloom.sum.quad$SampleRound),
    veg.name)]

veg.bloom.sum.quad$FloweringPlantDiversity <-
    veg.plant.div.quad[match(
        paste(veg.bloom.sum.quad$Site,
              veg.bloom.sum.quad$Quadrat,
              veg.bloom.sum.quad$Year,
              veg.bloom.sum.quad$SampleRound),
        veg.name)]

## merge together quad combos and quad level data
combos <- merge(combos, veg.bloom.sum.quad, all.x=TRUE)

cols.to.fill <-  c("FloralAbundance", "FloralRichness",
                   "FloralDiversity", "FloweringPlantAbundance",
                   "FloweringPlantDiversity")
## set NAs to zero
combos[, cols.to.fill][is.na(combos[, cols.to.fill])] <- 0

veg.bloom.mean <- combos %>%
    group_by(Site, Year, SampleRound) %>%
    summarise(MeanFloralAbundance = mean(FloralAbundance),
              MeanFloralRichness= mean(FloralRichness),
              MeanFloralDiversity =mean(FloralDiversity),
              MeanFloweringPlantDiversity =mean(FloweringPlantDiversity),
              MeanFloweringPlantAbundance =mean(FloweringPlantAbundance)
              )

veg.bloom.mean[order(veg.bloom.mean$MeanFloralRichness), ]

write.csv(veg.bloom.mean, file="../data/veg.csv", row.names=FALSE)

## merge with specimen summary data
site.sum <- merge(site.sum, veg.bloom.mean, all.x=TRUE)
write.csv(site.sum, file='../data/sitestats.csv', row.names=FALSE)

## floral richness across the entire meadow across sampling rounds
veg.year.sum <- veg.blooming %>%
    group_by(Site, Year) %>%
    summarise(TotalFloralRichness= length(unique(PlantGenusSpecies)))

write.csv(veg.year.sum, file="../data/veg_species_richness.csv",
          row.names=FALSE)

## *******************************************************************
## merging site and spec data
## *******************************************************************

## Merging specimen data with site characteristic data.
print("Before merge with site characteristics")
print(dim(spec.net))
spec.net <- merge(spec.net, site.sum, all.xy=TRUE)
print("After merge with site characteristics")
print(dim(spec.net))

## Load the csv with the parasitism rate by site/SR/year/spp. Create a
## key that matches the spec.net dataset and then merge the species
## level information to the full dataset.
spec.indiv <- read_csv("../data/spstats_net.csv")
spec.indiv$SiteSRYearSpp <- paste(spec.indiv$Site, spec.indiv$Year, spec.indiv$SampleRound,
                                spec.indiv$GenusSpecies, sep = "_")
spec.indiv <- subset(spec.indiv, select = -c(Site, Year, SampleRound, GenusSpecies))

spec.net$SiteSRYearSpp <- paste(spec.net$Site, spec.net$Year, spec.net$SampleRound,
                                spec.net$GenusSpecies, sep = "_")

## Merging species data with individual level data.
print("Before merge with ind level data")
print(dim(spec.net))
spec.net <- merge(spec.net, spec.indiv, all.x=TRUE)
spec.net$SiteSRYearSpp <- NULL
print("After merge with ind level data")
print(dim(spec.net))

## bee traits
bee.traits <-
    read.csv("../../skyIslands_saved/data/raw/bee_traits.csv")
bee.traits$GenusSpecies <- fix.white.space(bee.traits$GenusSpecies)
bee.traits <- bee.traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]

## network traits
net.traits <- read.csv("../data/networks_traits.csv")
net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]

## merge network traits to specimen data
print("Before merge with network traits to species bee.traits")
print(dim(bee.traits))
all.traits <- merge(bee.traits, net.traits,
                    by="GenusSpecies", all.x=TRUE)
print("After merge with network bee.traits to species bee.traits")
print(dim(all.traits))

print("Before merge with traits")
print(dim(spec.net))
spec.net <- merge(spec.net, all.traits, all.x=TRUE)
print("After merge with traits")
print(dim(spec.net))

spec.net$GenusSpecies[is.na(spec.net$GenusSpecies)] <- ""
save(spec.net, file="../data/spec_net.Rdata")


## not using site-level bloom data currently, has not been checked
## over as of Nov 2022

## *******************************************************************
## checking veg data between years
## *******************************************************************

## to.sample <- c("CH", "HM", "JC", "MM", "PL", "RP", "SC", "SM")

## veg.prep <- veg[veg$Site %in% to.sample,]
## veg.split <- split(veg.prep, veg.prep$Site)

## spec.veg <- spec[spec$Site %in% to.sample,]
## spec.split <- split(spec.veg, spec.veg$Site)

## checkveg <- function(x){
##     print(unique(x$Site))
##     years <- tapply(x$PlantGenusSpecies, x$Year, unique)
##     if(length(years) == 2){
##         badmatch12 <- years[[1]][!years[[1]] %in% years[[2]]]
##         badmatch21 <- years[[2]][!years[[2]] %in% years[[1]]]
##         return(list(badmatch12=badmatch12,
##                     badmatch21=badmatch21))
##     } else if(length(years) == 3){
##         badmatch12 <- years[[1]][!years[[1]] %in% years[[2]]]
##         badmatch21 <- years[[2]][!years[[2]] %in% years[[1]]]

##         badmatch13 <- years[[1]][!years[[1]] %in% years[[3]]]
##         badmatch31 <- years[[3]][!years[[3]] %in% years[[1]]]

##         badmatch23 <- years[[2]][!years[[2]] %in% years[[3]]]
##         badmatch32 <- years[[3]][!years[[3]] %in% years[[2]]]

##         return(list(badmatch12=badmatch12,
##                     badmatch21=badmatch21,
##                     badmatch13= badmatch13,
##                     badmatch31=badmatch31,
##                     badmatch23=badmatch23,
##                     badmatch32=badmatch32))
##     } else{
##         return(NULL)
##     }
## }

## badmatches.veg <- lapply(veg.split, checkveg)

## badmatches.spec <- lapply(spec.split, checkveg)


## all.plant.sp <- rbind(data.frame(Site=veg.prep$Site,
##                                  PlantGenusSpecies=
##                                      veg.prep$PlantGenusSpecies),
##                       data.frame(Site=spec.veg$Site,
##                                  PlantGenusSpecies=
##                                      spec.veg$PlantGenusSpecies))

## all.plant.sp <- unique(all.plant.sp)
## all.plant.sp <- sort(all.plant.sp)

## write.csv(all.plant.sp, row.names=FALSE,
##           file="../data/vegbySite.csv")

