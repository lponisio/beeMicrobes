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
microbes <- colnames(spec)[grepl("16s", colnames(spec))]

## make community stucture site by species
microbes.comm <- makeStructFromComm(spec, microbes)
dist.microbes <- as.matrix(vegan::vegdist(microbes.comm$comm,
                                   "gower"))
## 16s site dist matrix by phylogeny
dist.phylo.microbes <-  as.matrix(unifrac(microbes.comm$comm,
                                          tree.16s))

# RBCL
RBCL <- colnames(spec)[grepl("RBCL", colnames(spec))]

RBCL.comm <- makeStructFromComm(spec, RBCL)
dist.RBCL <- as.matrix(vegdist(RBCL.comm$comm,
                               "gower"))
## rbcl phylo distance
dist.phylo.rbcl <-  as.matrix(unifrac(RBCL.comm$comm, tree.rbcl))

## wild bees, excluding honey bees
bee.comm <- makeCommStruct(spec.wild, "GenusSpecies")
dist.bee <- as.matrix(vegdist(bee.comm$comm,
                              "gower"))

## parasites
parasite.comm <- makeStructFromComm(spec, parasites)
dist.parasite <- as.matrix(vegdist(parasite.comm$comm,
                                   "gower"))

## geographic distance
geo <- unique(spec[,c("Site", "Lat", "Long")])
dist.geo <- rdist.earth(geo[,c("Long", "Lat")])
rownames(dist.geo) <- colnames(dist.geo)  <- geo$Site

## **********************************************************
##  MRMs ## multiple regression on distance matrices
## **********************************************************

## make sure the sites line up within the matrices
dist.bee <- dist.bee[rownames(dist.geo), rownames(dist.geo)]
dist.parasite <- dist.parasite[rownames(dist.geo), rownames(dist.geo)]
dist.phylo.microbes <- dist.phylo.microbes[rownames(dist.geo),
                                           rownames(dist.geo)]
dist.phylo.rbcl <- dist.phylo.rbcl[rownames(dist.geo),
                                   rownames(dist.geo)]

## 16s phylo, with RBCL phylo dist as x var
mrm.microbes <- MRM(as.dist(dist.phylo.microbes) ~ as.dist(dist.bee) +
        as.dist(dist.parasite) +
        as.dist(dist.phylo.rbcl) + as.dist(dist.geo),  nperm=10^5)
print(mrm.microbes)

## parasites with rbcl phylo
mrm.parasites <- MRM(as.dist(dist.parasite) ~ as.dist(dist.bee) +
      as.dist(dist.phylo.rbcl) + as.dist(dist.geo),
    nperm=10^5)

print(mrm.parasites)

## plotting
site.types <- spec$SiteType[match(rownames(dist.geo), spec$Site)]

combos <- outer(site.types, site.types, FUN=paste)

all.dists <- data.frame(parasites=dist.parasite[lower.tri(dist.parasite)],
                        bee=dist.bee[lower.tri(dist.bee)],
                        rbcl=dist.phylo.rbcl[lower.tri(dist.phylo.rbcl)],
                        microbes=dist.phylo.microbes[lower.tri(dist.phylo.microbes)],
                        geo=dist.geo[lower.tri(dist.geo)],
                        site.type=combos[lower.tri(combos)])

par <- ggplot(all.dists, aes(x=rbcl,y=parasites)) + geom_point() +
  xlab("Pollen phylogenetic distance") +
    ylab("Pathobiome dissimilarity")

microbe <- ggplot(all.dists, aes(x=rbcl,y=microbes)) + geom_point() +
      xlab("Pollen phylogenetic distance") +
    ylab("Microbiome  dissimilarity")

all <- grid.arrange(par, microbe)

ggsave("figures/MRMs.pdf", all, height=8, width=4)
