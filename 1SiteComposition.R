setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

rm(list=ls())

#prepare matrices for bees and parasites (site by species)
source("src/commPrep.R")

#create an object called Parasite that includes the 7 we screened
#create spec, all specimenms that tested positive for Apidae gene control
#create spec.wild, everything but honeybees
#2018 data only for now
source("src/initialize.R")
source("src/runMRM.R")

## ***************************************************************
##  distance matrices for parasites, bees, plants, microbiome, RBCL by
##  site
## WILD BEES ONLY
##  ****************************************************************
## 16s
microbes <- colnames(spec.wild)[grepl("16s", colnames(spec.wild))]

## make community stucture site by species
microbes.comm <- makeStructFromComm(spec.wild, microbes)
dist.microbes <- as.matrix(vegan::vegdist(microbes.comm$comm,
                                   "gower"))
## 16s site dist matrix by phylogeny
microbes.phylo.dist <-  as.matrix(unifrac(microbes.comm$comm,
                                          tree.16s))

# RBCL
RBCL <- colnames(spec.wild)[grepl("RBCL", colnames(spec.wild))]

RBCL.comm <- makeStructFromComm(spec.wild, RBCL)
## rbcl phylo distance
RBCL.phylo.dist <-  as.matrix(unifrac(RBCL.comm$comm, tree.rbcl))

## bees
bee.comm <- makeCommStruct(spec.wild, "GenusSpecies")
dist.bee <- as.matrix(vegdist(bee.comm$comm,
                              "gower"))

## parasites
parasite.comm <- makeStructFromComm(spec.wild, parasites)
dist.parasite <- as.matrix(vegdist(parasite.comm$comm,
                                   "gower"))

## geographic distance
geo <- unique(spec.wild[,c("Site", "Lat", "Long")])
dist.geo <- rdist.earth(geo[,c("Long", "Lat")])
rownames(dist.geo) <- colnames(dist.geo)  <- geo$Site

## **********************************************************
##  MRM
## **********************************************************
## multiple regression on distance matrices

## make sure the sites line up
dist.bee <- dist.bee[rownames(dist.geo), rownames(dist.geo)]
dist.parasite <- dist.parasite[rownames(dist.geo), rownames(dist.geo)]
microbes.phylo.dist <- microbes.phylo.dist[rownames(dist.geo),
                                           rownames(dist.geo)]
RBCL.phylo.dist <- RBCL.phylo.dist[rownames(dist.geo),
                                   rownames(dist.geo)]

## 16s phylo, with RBCL phylo dist as x var
MRM(as.dist(microbes.phylo.dist) ~ as.dist(dist.bee) +
        as.dist(dist.parasite) +
        as.dist(RBCL.phylo.dist) + as.dist(dist.geo),  nperm=10^4)

## parasites with rbcl phylo
MRM(as.dist(dist.parasite) ~ as.dist(dist.bee) +
      as.dist(RBCL.phylo.dist) + as.dist(dist.geo),
    nperm=10^4)
