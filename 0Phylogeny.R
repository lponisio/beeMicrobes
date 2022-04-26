## setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

rm(list=ls())
source("src/initialize.R")

## load the compiled tree
dated.tree <- read.tree("data/phylogenies/hedgerowsPhylo.new")

## cleaning so the tips match the species in the dataset

dated.tree$tip.label <- gsub("_", " ", dated.tree$tip.label,
                             fixed=TRUE)
dated.tree$tip.label <- gsub("-Dialictus-", " ", dated.tree$tip.label,
                             fixed=TRUE)
dated.tree$tip.label <- gsub("-Evylaeus-", " ", dated.tree$tip.label,
                             fixed=TRUE)
dated.tree$tip.label <- gsub("   ", " ", dated.tree$tip.label,
                             fixed=TRUE)

spec$GenusSpecies[spec$GenusSpecies == "Lasioglossum tegulariforme"] <- "Lasioglossum tegulare group"

spec$GenusSpecies[spec$GenusSpecies == "Lasioglossum sp. e"] <- "Lasioglossum granosum"

spec$GenusSpecies[spec$GenusSpecies == "Lasioglossum sp. d"] <-
    "Lasioglossum diversopunctatum"

spec$GenusSpecies[spec$GenusSpecies == "Ashmeadiella sp."] <- "Ashmeadiella bucconis denticulata"

spec$GenusSpecies[spec$GenusSpecies == "Hylaeus morphoA"] <-
    "Hylaeus bisinuatus"
spec$GenusSpecies[spec$GenusSpecies == "Hylaeus morphoB"] <-
    "Hylaeus calvus"
spec$GenusSpecies[spec$GenusSpecies == "Hylaeus morphoC"] <-
    "Hylaeus conspicuus"
spec$GenusSpecies[spec$GenusSpecies == "Hylaeus morphoD"] <-
    "Hylaeus episcopalis"
spec$GenusSpecies[spec$GenusSpecies == "Triepeolus sp. a"] <-
    "Triepeolus concavus"
spec$GenusSpecies[spec$GenusSpecies == "Agapostemon obliquus"] <-
    "Agapostemon texanus"


unique(spec$GenusSpecies[!spec$GenusSpecies %in% dated.tree$tip.label])

## pruning the tree so tips match
phylo <-
    drop.tip(dated.tree, dated.tree$tip.labe[!dated.tree$tip.label %in%
                                            spec$GenusSpecies])

co.var.mat  <- ape::vcv.phylo(phylo)

save(co.var.mat, phylo, file="data/covarmatrix_taxon.Rdata")


## **********************************************************************
## micheal's tree
## **********************************************************************

## load the compiled tree
dated.tree <- read.tree("/Volumes/bombus/Dropbox (University of Oregon)/sunflower_saved/data/tree/calibees-206t-95p-spruceup-iqtree-swsc-gtrg.tre")

## cleaning so the tips match the species in the dataset

dated.tree$tip.label <- gsub("_BLX\\d+", "", dated.tree$tip.label)

dated.tree$tip.label <- gsub("_", " ", dated.tree$tip.label)

## dated.tree <- drop.tip(dated.tree, dated.tree$tip.label[duplicated(dated.tree$tip.label)])

dated.tree$tip.label[dated.tree$tip.label == "Lasioglossum spE"] <-
    "Lasioglossum sp. e"

dated.tree$tip.label[dated.tree$tip.label == "Sphecodes sp"] <-
    "Sphecodes spp."

no.tips <- unique(spec$GenusSpecies[!spec$GenusSpecies %in%
                                    dated.tree$tip.label])

no.tips <- no.tips[no.tips != "Apis mellifera"]

genus.no.tips <- sapply(strsplit(no.tips, " "),
                       function (x) x[[1]])


## ## pruning the tree so tips match
## phylo <-
##     drop.tip(dated.tree, dated.tree$tip.labe[!dated.tree$tip.label %in%
##                                             spec$GenusSpecies])

co.var.mat  <- ape::vcv.phylo(dated.tree)


write.csv(co.var.mat, file=
                          "/Volumes/bombus/Dropbox (University of Oregon)/sunflower_saved/data/tree/hedgerow_covar_mat.csv")

## multi.tips <- colnames(co.var.mat)[duplicated(colnames(co.var.mat))]


## take the average for species with multiple specimens sequenced
tips <- unique(colnames(co.var.mat))

out.var <- vector(mode="list", length=length(tips))

for(i in 1:length(unique(tips))){
    print(tips[i])
    co.var.mat.sub <- co.var.mat[rownames(co.var.mat) == tips[i], ]
    print(length(co.var.mat.sub))
    if(length(co.var.mat.sub) != 206){
        out.var[[i]] <- apply(co.var.mat.sub, 2, mean)
    } else {
        out.var[[i]] <- co.var.mat.sub
    }
}

out.var2 <- do.call(rbind, out.var)

rownames(out.var2) <-
    colnames(out.var2)[!duplicated(colnames(out.var2))]

out.var3 <- out.var2[, !duplicated(colnames(out.var2))]

write.csv(out.var3, file=
                          "/Volumes/bombus/Dropbox (University of Oregon)/sunflower_saved/data/tree/hedgerow_covar_mat_no_dups.csv")

problem.genera <- c("Agapostemon", "Hylaeus", "Lasioglossum",
                    "Xylocopa",
                    "Triepeolus", "Ashmeadiella")

new.tips <- unique(colnames(out.var3))

genus.co.var <- sapply(strsplit(new.tips, " "),
                       function (x) x[[1]])

## new.tips.problem.gen <- new.tips[genus.co.var %in% problem.genera]

empty.problem.gen.cols <- matrix(rep(NA,
                                nrow(out.var3)*length(no.tips)),
                            ncol=length(no.tips))

rownames(empty.problem.gen.cols) <- rownames(out.var3)
colnames(empty.problem.gen.cols) <- no.tips

out.var4 <- cbind(out.var3, empty.problem.gen.cols)

empty.problem.gen.rows <- matrix(rep(NA,
                                ncol(out.var4)*length(no.tips)),
                            nrow=length(no.tips))

colnames(empty.problem.gen.rows) <- colnames(out.var4)
rownames(empty.problem.gen.rows) <- no.tips

out.var5 <- rbind(out.var4, empty.problem.gen.rows)


for(genus in problem.genera){
    print(genus)
    co.var.mat.sub <- out.var3[genus.co.var  ==
                               genus, ]
    if(is.matrix(co.var.mat.sub)){
        mean.mat <- apply(co.var.mat.sub, 2, mean, rm.na=TRUE)
    } else {
        mean.mat <- co.var.mat.sub
    }
    ## out.var4  <- cbind(out.var4, mean.mat)
    ## colnames(out.var4)[colnames(out.var3) == "mean.mat"] <- genus

    mean.genus <- rep(NA, length(problem.genera))
    names(mean.genus) <- problem.genera
    for(j in problem.genera){
        if(is.matrix(co.var.mat.sub)){
            co.var.mat.sub.sub <-
                co.var.mat.sub[, genus.co.var  ==  j]
        } else{
            co.var.mat.sub.sub <-
                co.var.mat.sub[genus.co.var  ==  j]
        }
        mean.genus[j] <- mean(co.var.mat.sub.sub)
    }
    reped.mean.genus <- mean.genus[genus.no.tips]
    names(reped.mean.genus) <- no.tips

    these.tips <- no.tips[genus.no.tips == genus]

    for(k in these.tips){
        out.var5[k, ] <- c(mean.mat, reped.mean.genus)
        out.var5[, k ] <- c(mean.mat,  reped.mean.genus)
    }


}

write.csv(out.var5, file=
                          "/Volumes/bombus/Dropbox (University of Oregon)/sunflower_saved/data/tree/hedgerow_covar_mat_problem_genera.csv")

co.var.mat <- out.var5

save(co.var.mat, dated.tree, file="data/covarmatrix.Rdata")

unique(spec$GenusSpecies[!spec$GenusSpecies %in%
                         rownames(out.var5)])
