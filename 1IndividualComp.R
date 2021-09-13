## setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

## this code takes individual-level bee data to see if the composition of microbiome is related to composition of pathobiome
## uses MRM, an extension of mantel tests, which allows for multiple matrics. We included RBCL and geo to account for diet breadth
## for 16s, use phylogenetic distances to matrix distance matrix
## for RBCL, use species taxonomic distances to construct distance matrix

##setwd("~/Dropbox/beeMicrobes")
rm(list=ls())
source("src/initialize.R")

## plot comm dist by group,
source("src/runMRM.R")
source("src/mod_phylo_funs.R")
source("src/makeIndivComm.R")
source("src/commPrep.R")

## ***************************************************************
## 16s
##  **************************************************************
microbes <- colnames(spec)[grepl("16s",
                                      colnames(spec))]


comm.microbes.indiv <- makeIndivComm(spec, microbes)

## species
dist.microbes <- as.matrix(vegan::vegdist(comm.microbes.indiv,
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


group.info <- regmatches(prune.tree.16s$tip.label,regexpr("D\\_4\\_\\_[A-Za-z]+",
                                            prune.tree.16s$tip.label))
group.info <- gsub("D\\_4\\_\\_", "", group.info)

tip.lab <- regmatches(prune.tree.16s$tip.label,regexpr("D\\_5\\_\\_[A-Za-z]+",
                                            prune.tree.16s$tip.label))
tip.lab <- gsub("D\\_5\\_\\_", "", tip.lab)


prune.tree.16s <- groupOTU(prune.tree.16s, group.info)
ggtree(prune.tree.16s, aes(color=group), layout='circular') +
    geom_tiplab(size=1, aes(angle=angle))


prune.tree.16s$tip.label <- tip.lab

names(group.info) <- prune.tree.16s$tip.label

td <- data.frame(node = nodeid(prune.tree.16s, names(group.info)),
                 trait = group.info)

td$node <- as.numeric(td$node)
tree.d <- full_join(prune.tree.16s, td, by = 'node')


p1 <- ggtree(tree.d, aes(color=trait), layout = 'circular',
        ladderize = FALSE, size=2)  +
    geom_tiplab(hjust = -.1) +
    xlim(0, 1.2) +
    theme(legend.position = c(.05, .85))


prune.tree.16s <- groupClade(prune.tree.16s, groupInfo)
ggtree(prune.tree.16s, aes(color=groupInfo), layout='circular') + geom_tiplab(size=1, aes(angle=angle))



ggtree(prune.tree.16s, layout='circular') + geom_tiplab(size=1, aes(angle=angle))

##  ****************************************************************
## RBCL
##  ****************************************************************
rbcl <- colnames(spec)[grepl("RBCL", colnames(spec))]
comm.rbcl.indiv <- makeIndivComm(spec, rbcl)

## species
dist.rbcl <- as.matrix(vegan::vegdist(comm.rbcl.indiv,
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

parasite.comm <- spec[, c("UniqueID", "Apidae", parasites)]
parasite.comm <- parasite.comm[parasite.comm$Apidae == 1,]
parasite.comm[, c("Apidae", parasites)] <-
    apply(parasite.comm[, c("Apidae", parasites)], 2, as.numeric)
rownames(parasite.comm) <- parasite.comm$UniqueID
parasite.comm$UniqueID <- NULL
parasite.comm <- as.matrix(parasite.comm)

dist.parasite <- as.matrix(vegan::vegdist(parasite.comm,
                                          "jaccard"))

save(dist.parasite, file="saved/distmats/indiv_parasite.Rdata")

##  ****************************************************************
## geo and bee commmunities
##  ****************************************************************
indiv <- rownames(dist.parasite)
geo <- spec[,c("Site", "Lat", "Long", "UniqueID")][spec$UniqueID %in%
                                                   rownames(dist.parasite),]
dist.geo <- rdist.earth(geo[,c("Long", "Lat")])
rownames(dist.geo) <- colnames(dist.geo)  <- geo$UniqueID


## wild bees, excluding honey bees
bee.comm <- makeCommStruct(spec.wild, "GenusSpecies")
dist.bee <- as.matrix(vegdist(bee.comm$comm,
                              "gower"))

## **********************************************************
##  MRMs  multiple regression on distance matrices
## **********************************************************
load(file="saved/distmats/indiv_parasite.Rdata")
load(file="saved/distmats/indiv_16s.Rdata")
load(file="saved/distmats/indiv_rbcl.Rdata")

## make sure the sites line up
in.all <- rownames(dist.phylo.microbes)[
    rownames(dist.phylo.microbes) %in% rownames(dist.parasite)]

in.all <- in.all[in.all %in% rownames(dist.phylo.rbcl)]

dist.parasite <- dist.parasite[in.all, in.all]
dist.rbcl <- dist.rbcl[in.all, in.all]
dist.phylo.rbcl <- dist.phylo.rbcl[in.all, in.all]
dist.phylo.microbes <- dist.phylo.microbes[in.all, in.all]
dist.geo <- dist.geo[in.all, in.all]

sites.all <- spec$Site[match(in.all, spec$UniqueID)]

dist.bee <- dist.bee[sites.all, sites.all]

## **********************************************************
## 16s phylo vs. parasite, rbcl species dissimilarity,

MRM(as.dist(dist.phylo.microbes) ~  +
        as.dist(dist.parasite) +
        as.dist(dist.geo) +
        as.dist(dist.bee)  +
        as.dist(dist.rbcl),  nperm=10^4)

## 16s phylo vs. parasite, rbcl phylo dissimilarity
MRM(as.dist(dist.phylo.microbes) ~  +
        as.dist(dist.parasite) +
        as.dist(dist.geo) +
        as.dist(dist.bee) +
        as.dist(dist.phylo.rbcl),  nperm=10^4)

## same result with phylo and rbcl taxon dist matrices

##**********************************************************
## parasite vs. 16s phylo, rbcl species dissimilariy
MRM(as.dist(dist.parasite) ~  +
        as.dist(dist.phylo.microbes) +
        as.dist(dist.geo) +
        as.dist(dist.bee)  +
        as.dist(dist.rbcl),  nperm=10^4)

## parasite vs. 16s phylo, rbcl phylogenetic dissimilariy
MRM(as.dist(dist.parasite) ~  +
        as.dist(dist.phylo.microbes) +
        as.dist(dist.geo) +
        as.dist(dist.bee) +
        as.dist(dist.phylo.rbcl),  nperm=10^4)

## same result with phylo and rbcl taxon dist matrices

## plotting

## combine into one dataset
all.dists <- data.frame(parasites=dist.parasite[lower.tri(dist.parasite)],
                        rbcl.phylo=dist.phylo.rbcl[lower.tri(dist.phylo.rbcl)],
                        rbcl=dist.rbcl[lower.tri(dist.rbcl)],
                        microbes=dist.phylo.microbes[lower.tri(dist.phylo.microbes)],
                        geo=dist.geo[lower.tri(dist.geo)])

sites <- spec$Site[match(rownames(dist.geo), spec$UniqueID)]
combos <- outer(sites, sites, FUN=paste)
all.dists$Sites <-  combos[lower.tri(combos)]

species <- spec$GenusSpecies[match(rownames(dist.geo), spec$UniqueID)]
combos.sp <- outer(species, species, FUN=paste)
all.dists$Species <-  combos.sp[lower.tri(combos.sp)]

par <- ggplot(all.dists, aes(y=parasites,x=microbes)) + geom_point() +
  ylab("Pathobiome community dissimilarity") +
    xlab("Microbiome dissimilarity")

microbe <- ggplot(all.dists, aes(y=rbcl,x=microbes)) + geom_point() +
      ylab("Pollen community dissimilarity") +
    xlab("Microbiome  dissimilarity")

all <- grid.arrange(par, microbe)

ggsave("figures/MRMsIndiv.pdf", all, height=8, width=4)

## ## **********************************************************
## ##  species specific MRMs
## ## **********************************************************

## ## run these models for each species individually, but only if >10
## ## individuals for each species

## species <- split(spec, spec$GenusSpecies)
## sp.ids <- lapply(species, function(x) x$UniqueID)
## sp.ids <- sp.ids[sapply(sp.ids, length) >= 10]

## ## 16s
## mrm.16sPhylo.by.sp <- lapply(sp.ids, runMantelBeeSpecies,
##                     dist.microbe=dist.phylo.microbes,
##                     dist.plant=dist.rbcl,
##                     dist.parasite=dist.parasite)

## mrm.16sPhylo.by.sp <- mrm.16sPhylo.by.sp[!sapply(mrm.16sPhylo.by.sp,
##                                                  function(x)
##                                                      is.na(x[1]))]

## save(mrm.16sPhylo.by.sp, mrm.16s.by.sp,
##      file="saved/distmats/mrm16sPhylo.Rdata")

