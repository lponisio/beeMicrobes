setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

rm(list=ls())
source("src/initialize.R")

## plot comm dist by group,
source("src/runMRM.R")
source("src/mod_phylo_funs.R")
source("src/makeIndivComm.R")

## ***************************************************************
## 16s
##  **************************************************************
microbes <- colnames(spec.wild)[grepl("16s",
                                      colnames(spec.wild))]


comm.microbes.indiv <- makeIndivComm(spec.wild, microbes)

## species
dist.microbes <- as.matrix(vegan::vegdist(comm.microbes.indiv*100,
                                "altGower"))
## phylo make distance matrix using merged tree, prune tree to match
## community dataset
prune.tree.16s <- prune.sample(comm.microbes.indiv, tree.16s)

## TAKES A LONGGGGG TIME
dist.phylo.microbes <- mod.unifrac(comm.microbes.indiv*100,
                                          tree.16s)
dist.phylo.microbes <- as.matrix(dist.phylo.microbes)

## drop the no apidae specimens
dist.phylo.microbes <- dist.phylo.microbes[
    !rownames(dist.phylo.microbes) %in% no.apidae,
    !colnames(dist.phylo.microbes) %in% no.apidae]

dist.microbes <- dist.microbes[
    !rownames(dist.microbes) %in% no.apidae,
    !colnames(dist.microbes) %in% no.apidae]

save(dist.phylo.microbes, dist.microbes,
     file="saved/distmats/indiv_16s.Rdata")

##  ****************************************************************
## RBCL
##  ****************************************************************
rbcl <- colnames(spec.wild)[grepl("RBCL", colnames(spec.wild))]
comm.rbcl.indiv <- makeIndivComm(spec.wild, rbcl)

## species
dist.rbcl <- as.matrix(vegan::vegdist(comm.rbcl.indiv*100,
                                "altGower"))
## phylo: prune tree to match community dataset
prune.tree.rbcl <- prune.sample(comm.rbcl.indiv, tree.rbcl)

dist.phylo.rbcl <- mod.unifrac(comm.rbcl.indiv*100,
                                          tree.rbcl)

dist.phylo.rbcl <- as.matrix(dist.phylo.rbcl)

save(dist.phylo.rbcl, dist.rbcl,
     file="saved/distmats/indiv_rbcl.Rdata")

##  ****************************************************************
## parasites
##  ****************************************************************
## include Apidae as a dummy species to avoid the issue of having
## individuals dropped if they did not have any parasites.

parasite.comm <- spec.wild[, c("UniqueID", "Apidae", parasites)]
parasite.comm <- parasite.comm[parasite.comm$Apidae == 1,]
parasite.comm[, c("Apidae", parasites)] <-
    apply(parasite.comm[, c("Apidae", parasites)], 2, as.numeric)
rownames(parasite.comm) <- parasite.comm$UniqueID
parasite.comm$UniqueID <- NULL
parasite.comm <- as.matrix(parasite.comm)

dist.parasite <- as.matrix(vegan::vegdist(parasite.comm,
                                          "jaccard"))

save(dist.parasite, file="saved/distmats/indiv_parasite.Rdata")

## **********************************************************
##  MRMs  multiple regression on distance matrices
## **********************************************************

## make sure the sites line up
in.all <- rownames(dist.phylo.microbes)[
    rownames(dist.phylo.microbes) %in% rownames(dist.parasite)]

in.all <- in.all[in.all %in% rownames(dist.phylo.rbcl)]

dist.parasite <- dist.parasite[in.all, in.all]
dist.rbcl <- dist.rbcl[in.all, in.all]
dist.phylo.rbcl <- dist.phylo.rbcl[in.all, in.all]
dist.phylo.microbes <- dist.phylo.microbes[in.all, in.all]
dist.microbes <- dist.microbes[in.all, in.all]

## **********************************************************
## 16s phylo vs. parasite, rbcl species dissimilarity,
MRM(as.dist(dist.phylo.microbes) ~  +
        as.dist(dist.parasite) +
        as.dist(dist.rbcl),  nperm=10^4)

## 16s phylo vs. parasite, rbcl phylo dissimilarity
MRM(as.dist(dist.phylo.microbes) ~  +
        as.dist(dist.parasite) +
        as.dist(dist.phylo.rbcl),  nperm=10^4)

##**********************************************************
## parasite vs. 16s phylo, rbcl phylogenetic dissimilariy
MRM(as.dist(dist.parasite) ~  +
        as.dist(dist.phylo.microbes) +
        as.dist(dist.phylo.rbcl),  nperm=10^4)

## parasite vs. 16s phylo, rbcl species dissimilariy
MRM(as.dist(dist.parasite) ~  +
        as.dist(dist.phylo.microbes) +
        as.dist(dist.rbcl),  nperm=10^4)

## **********************************************************
##  species specific MRMs
## **********************************************************

## run these models for each species individually, but only if >10
## individuals for each species

species <- split(spec, spec$GenusSpecies)
sp.ids <- lapply(species, function(x) x$UniqueID)
sp.ids <- sp.ids[sapply(sp.ids, length) >= 10]

## 16s
mrm.16sPhylo.by.sp <- lapply(sp.ids, runMantelBeeSpecies,
                    dist.microbe=dist.phylo.microbes,
                    dist.plant=dist.rbcl,
                    dist.parasite=dist.parasite)

mrm.16sPhylo.by.sp <- mrm.16sPhylo.by.sp[!sapply(mrm.16sPhylo.by.sp,
                                                 function(x)
                                                     is.na(x[1]))]

save(mrm.16sPhylo.by.sp, mrm.16s.by.sp,
     file="saved/distmats/mrm16sPhylo.Rdata")

