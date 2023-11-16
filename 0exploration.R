


traits <- c( "NestLoc", "Excavate",
            "Sociality", "MeanITD", "PollenCarry", "r.degree", "NestPartitions")


gensp <- unique(all.indiv.mets[, c("GenusSpecies", traits)])

rownames(gensp) <- gensp$GenusSpecies
gensp$GenusSpecies <- NULL

gensp <- gensp[!is.na(gensp$r.degree),]

gensp$MeanITD <- scale(gensp$MeanITD)

gensp$r.degree <- scale(gensp$r.degree)

trait_famd <- FAMD(gensp,
                  graph=FALSE)

fviz_famd_ind(trait_famd,col.ind = "cos2",
             gradient.cols = c("blue", "orange", "red"),
             repel = TRUE)
