#Network building for 16s and rbcl

# Load in software libraries
library(bipartite)
library(parallel)
library(dplyr)
library(readxl)


# Set the working directory

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

dir.microbes <- file.path(local.path, "beeMicrobes")
setwd(dir.microbes)

# load in-house functions for network generation
source('./src/nets_init.R')
source('./src/calcFuncUniqOrig.R')

ndim.func <- 10

## ********************************************************
## pollinator abundance by site to multiple interactions by
## ********************************************************
spec.net <- spec.net %>% select(-UniqueID) %>%
  rename("UniqueID" = "UniqueID.fixed") %>%
  mutate(SiteYear = paste0(Site, Year))
pol.abund <- spec.net %>%
  group_by(Site, SampleRound, GenusSpecies) %>%
  summarise(Abundance = n()) %>% ungroup

## prep for making a vector of bee abund at each site
pol.abund.site <- split(pol.abund, pol.abund$Site)


## RBCL temporarily commented out until the full RBCL data arrive
# ## ********************************************************
# ## RBCL
# ## ********************************************************

# ## list of rbcl data columns
# rbcl_names <- grep('RBCL',names(spec.net))

# # note check out RBCL:NA

# ## split out specimens that were sequenced
# seq.df <- spec.net[apply(spec.net[, rbcl_names], 1, function(x) !all(is.na(x))),]

# ## makes community matrix/network of plants by individuals
# indivNet_rbcl <- illumSplit(seq.df, "Site", rbcl_names)
# ## Format: each element of this list is a large matrix for that site


# pol.sp.rbcl <- lapply(indivNet_rbcl,
#                 function(x){
#                 spec.net$GenusSpecies[match(colnames(x), spec.net$UniqueID)]
# })


# pol.abund.site.rbcl <- pol.abund.site[names(indivNet_rbcl)]

# ## species level network
# spNet_rbcl <-  mapply(makeSpNet,
#     indiv.nets = indivNet_rbcl,
#     sp.names = pol.sp.rbcl,
#     sp.abund=pol.abund.site.rbcl,
#     SIMPLIFY = FALSE)


# ## save these matrices for use elsewhere
# save(indivNet_rbcl, spNet_rbcl, pol.sp.rbcl,
#      file=file.path(save.dir, 'rbclNets.RData'))

# ## ********************************************************

# ## calculate summary statistics for each "species" (individual) in the
# ## networks

# ## set m to the min number of traits. The Function's default is to use
# ## each trait as a dimensions, which get's stuck due to the high
# ## number of zeros, and "traits"
# ## see https://stat.ethz.ch/pipermail/r-sig-ecology/2016-January/005264.html

# nindiv.rbcl <- sapply(indivNet_rbcl, ncol) > 1

# nsp.rbcl <- sapply(spNet_rbcl, ncol) > 1

# indivNetSums_rbcl  <-  netSums(indivNet_rbcl[nindiv.rbcl],
#                                binary=FALSE,
#                                spec=spec.net, FUN=mclapply,
#                                m=ndim.func)

# indivNetSums_rbclBinary  <-  netSums(indivNet_rbcl[nindiv.rbcl],
#                                      binary=TRUE,
#                                      spec=spec.net,
#                                      FUN=mclapply,
#                                      m=ndim.func)

# spNetSums_rbcl  <-  netSums(spNet_rbcl[nsp.rbcl],
#                             binary=FALSE,
#                             spec=spec.net,
#                             type="sp")

# spNetSums_rbclBinary  <-  netSums(spNet_rbcl[nsp.rbcl],
#                                   binary=TRUE,
#                                   spec=spec.net,
#                                   type="sp")


# save(indivNetSums_rbcl, indivNetSums_rbclBinary,
#      spNetSums_rbcl,
#      spNetSums_rbclBinary,
#      file='./data/NetSums_rbcl.RData')

## ********************************************************
## 16s
## ********************************************************
#list of column names for 16s
names_16s  <-  grep('^16s:', names(spec.net), value=TRUE)

#split out specimens that were sequenced
seq.micro <- spec.net[!is.na(spec.net$Site) &
                        rowSums(spec.net[, names_16s], na.rm=TRUE) > 0, ]

## individual and species networks for a grouping column
buildMicroNets <- function(group){
  indivNet <- illumSplit(seq.micro, group, names_16s)
  pol.sp <- lapply(indivNet, function(x){
    spec.net$GenusSpecies[match(colnames(x), spec.net$UniqueID)]
  })
  pol.abund <- spec.net %>% count(.data[[group]], GenusSpecies, name = "Abundance")
  pol.abund <- split(pol.abund, pol.abund[[group]])[names(indivNet)]
  spNet <- mapply(makeSpNet,
                  indiv.nets = indivNet,
                  sp.names = pol.sp,
                  sp.abund = pol.abund,
                  SIMPLIFY = FALSE)
  list(indivNet = indivNet, spNet = spNet, pol.sp = pol.sp)
}

nets.site <- buildMicroNets("Site")
nets.siteyear <- buildMicroNets("SiteYear")

indivNet_micro <- nets.site$indivNet
spNet_micro <- nets.site$spNet
pol.sp.micro <- nets.site$pol.sp

indivNet_microSY <- nets.siteyear$indivNet
spNet_microSY <- nets.siteyear$spNet
pol.sp.microSY <- nets.siteyear$pol.sp

#save these matrices for use elsewhere
save(indivNet_micro, spNet_micro, pol.sp.micro,
     indivNet_microSY, spNet_microSY, pol.sp.microSY,
     file = './data/microNets.RData')

## calculate network summary statistics
## site level
indivNetSums_micro  <-  netSums(indivNet_micro, spec=spec.net,
                                m=ndim.func, scale.out=FALSE)

indivNetSums_microBinary  <-  netSums(indivNet_micro, binary=TRUE,
                                      spec=spec.net, m=ndim.func, scale.out=FALSE)

spNetSums_micro  <-  netSums(spNet_micro, spec=spec.net, type="sp",
                             m=ndim.func, scale.out=FALSE)

spNetSums_microBinary  <-  netSums(spNet_micro, binary=TRUE,
                                   spec=spec.net, type="sp", m=ndim.func, scale.out=FALSE)

## site-year level
indivNetSums_microSY  <-  netSums(indivNet_microSY, spec=spec.net,
                                  group.col="SiteYear", m=ndim.func, scale.out=FALSE)

indivNetSums_microSYBinary  <-  netSums(indivNet_microSY, binary=TRUE,
                                        spec=spec.net, group.col="SiteYear",
                                        m=ndim.func, scale.out=FALSE)

spNetSums_microSY  <-  netSums(spNet_microSY, spec=spec.net, type="sp",
                               group.col="SiteYear", m=ndim.func, scale.out=FALSE)

spNetSums_microSYBinary  <-  netSums(spNet_microSY, binary=TRUE,
                                     spec=spec.net, type="sp",
                                     group.col="SiteYear", m=ndim.func, scale.out=FALSE)

save(indivNetSums_micro, indivNetSums_microBinary,
     spNetSums_micro, spNetSums_microBinary,
     indivNetSums_microSY, indivNetSums_microSYBinary,
     spNetSums_microSY, spNetSums_microSYBinary,
     file = './data/NetSums_micro.RData')


## ***********************************************************************
## Parasite networks
## ***********************************************************************
### NOTE: All parasite metrics cause issues and do not have enough data to
###   generate useful or accurate results. Thus they are commented out and
###   not run as of the 3/16/2026 run of this script.
# 
# ## get columns of parasite presence to feed into illumSplit
# 
# names.Para <- c("AscosphaeraSpp", "ApicystisSpp",
#                 "CrithidiaBombi", "CrithidiaExpoeki",
#                 "NosemaCeranae", "NosemaBombi", "CrithidiaSpp")
# 
# spec.net[, names.Para] <- apply(spec.net[, names.Para], 2, as.numeric)
# seq.para <- spec.net[apply(spec.net[, names.Para], 1, function(x) !all(is.na(x))),]
# 
# indivNet_para  <-  illumSplit(seq.para,"Site",names.Para)
# 
# 
# pol.sp.para <- lapply(indivNet_para,
#                 function(x){
#                     spec.net$GenusSpecies[match(colnames(x),
#                                             spec.net$UniqueID)]
# })
# 
# ## species level
# pol.abund.site.para <- pol.abund.site[names(indivNet_para)]
# 
# ## species level network
# spNet_para <-  mapply(makeSpNet,
#     indiv.nets = indivNet_para,
#     sp.names = pol.sp.para,
#     sp.abund=pol.abund.site.para,
#     SIMPLIFY = FALSE)
# 
# 
# #save networks themselves
# save(indivNet_para, spNet_para, pol.abund.site.para,
#      file =file.path(save.dir, 'paraNets.RData'))
# 
# #calculate network summary statistics
# 
# indivNetSums_para  <-  netSums3(indivNet_para, spec=spec.net)
# 
# spNetSums_para  <-  netSums3(spNet_para, spec=spec.net,
#                             type="sp", FUN=lapply)
# 
# 
# save(indivNetSums_para, spNetSums_para,
#      file = file.path('./data/NetSums_para.RData'))


## ***********************************************************************
## combined dataframe of all the network metrics
## ***************************************************************
load(file='./data/NetSums_micro.RData')
## RBCL temporarily commented out until the full RBCL data arrive
# load(file='./data/NetSums_rbcl.RData')
#load(file='./data/NetSums_para.RData')

## prefix metric columns with the data type
makeNewColnames <- function(x, name){
    id.cols <- c("Site", "SiteYear", "GenusSpecies", "UniqueID")
    new <- !colnames(x) %in% id.cols
    colnames(x)[new] <- paste(name, colnames(x)[new], sep="_")
    return(colnames(x))
}

# commenting out all parasite-related code
#colnames(indivNetSums_para) <- makeNewColnames(indivNetSums_para, "Parasite")
## RBCL temporarily commented out until the full RBCL data arrive
# colnames(indivNetSums_rbcl) <- makeNewColnames(indivNetSums_rbcl, "RBCL")
# colnames(indivNetSums_rbclBinary) <- makeNewColnames(indivNetSums_rbclBinary, "RBCL")
# colnames(spNetSums_rbcl) <- makeNewColnames(spNetSums_rbcl, "RBCL")
colnames(indivNetSums_micro) <- makeNewColnames(indivNetSums_micro, "Micro")
colnames(indivNetSums_microBinary) <- makeNewColnames(indivNetSums_microBinary, "Micro")
colnames(spNetSums_micro) <- makeNewColnames(spNetSums_micro, "Micro")
colnames(spNetSums_microBinary) <- makeNewColnames(spNetSums_microBinary, "Micro")
colnames(indivNetSums_microSY) <- makeNewColnames(indivNetSums_microSY, "MicroSY")
colnames(indivNetSums_microSYBinary) <- makeNewColnames(indivNetSums_microSYBinary, "MicroSY")
colnames(spNetSums_microSY) <- makeNewColnames(spNetSums_microSY, "MicroSY")
colnames(spNetSums_microSYBinary) <- makeNewColnames(spNetSums_microSYBinary, "MicroSY")

## individual metrics at both levels, one row per bee with microbe data
combineIndivMets <- function(site.mets, siteyear.mets){
  site.mets <- site.mets[, !colnames(site.mets) %in% c("Site", "GenusSpecies")]
  siteyear.mets <- siteyear.mets[, !colnames(siteyear.mets) %in% c("SiteYear", "GenusSpecies")]
  mets <- merge(site.mets, siteyear.mets, by="UniqueID", all=TRUE)
  merge(spec.net, mets, by="UniqueID")
}

all.indiv.mets <- combineIndivMets(indivNetSums_micro, indivNetSums_microSY)
all.indiv.mets.binary <- combineIndivMets(indivNetSums_microBinary,
                                          indivNetSums_microSYBinary)

## RBCL temporarily commented out until the full RBCL data arrive
# all.indiv.mets <- merge(all.indiv.mets,
#                         indivNetSums_rbcl[, !colnames(indivNetSums_rbcl) %in% c("Site", "GenusSpecies")],
#                         by="UniqueID", all.x=TRUE)
# all.indiv.mets.binary <- merge(all.indiv.mets.binary,
#                                indivNetSums_rbclBinary[, !colnames(indivNetSums_rbclBinary) %in% c("Site", "GenusSpecies")],
#                                by="UniqueID", all.x=TRUE)

## combine species level metrics
# all.sp.mets <- merge(spNetSums_micro, spNetSums_rbcl, all=TRUE)
# all.sp.mets <- merge(all.sp.mets, spNetSums_para, all=TRUE)
# all.sp.mets.binary <- merge(spNetSums_microBinary,
#                             spNetSums_rbclBinary,
#                             all=TRUE)
# all.sp.mets.binary <- merge(all.sp.mets.binary, spNetSums_para, all=TRUE)


## species level merge
#sp.by.site <- read.csv("../skyIslands/data/spstats.csv")

#all.sp.mets <- merge(all.sp.mets, sp.by.site)
#all.sp.mets.binary <- merge(all.sp.mets.binary, sp.by.site)

save(#all.sp.mets.binary,
     #all.sp.mets,
     all.indiv.mets.binary,
     all.indiv.mets,
     file = './data/allNetSums.RData')
