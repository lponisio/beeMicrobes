# This script prepares 16s and RBCL community diversity metrics and
#   adjusts phylogenetic trees and their tip labels
## NOTE: AS OF 10/22/26 RBCL DATA IS CAUSING ISSUES AND IS NOT RUN

# clear environment
rm(list=ls())

################################
# Load packages and functions
################################

if (!requireNamespace("qiime2R", quietly = TRUE)){
  if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
  devtools::install_github("jbisanz/qiime2R")
}
library(qiime2R)
library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(ape)
library(picante)

# set working directory to the beeMicrobes box folder
setwd("~/")
source("lab_paths.R")
local.path
dir.microbes <- file.path(local.path, "beeMicrobes"); setwd(dir.microbes)

## run misc functions script
source("./src/misc.R")

################################
# Read in 16s data
################################

# set working directory
setwd("../beeMicrobes_saved/")

## create filepaths to quickly read in qiime2 artifacts
qza.16s.path  <- "beeMicrobes_pipeline_output/merged/16s/final"

## Access formatted 16s phylogeny via phyloseq
physeq16s <- qza_to_phyloseq(
  features = file.path(qza.16s.path, "core_metrics16s/rarefied_table.qza"),
  tree =     "./beeMicrobes_pipeline_output/merged/16s/rooted-tree16s.qza",
  taxonomy = "./beeMicrobes_pipeline_output/merged/16s/taxonomy16s.qza",
  metadata = "./beeMicrobes_pipeline_output/maps/16s/beeMicrobes_combined_map_no-ctrls.txt"
  )

# check number of samples, etc to make sure everything is there.
physeq16s

# load in taxonomy with feature ID column (the matching key)
feature.2.tax.16s <- read.table("beeMicrobes_pipeline_output/merged/16s/taxonomy16s.tsv",
                                sep="\t", header=TRUE, quote="", comment.char="")

feature.2.tax.16s$Taxon <- paste("16s", feature.2.tax.16s$Taxon, sep=':')

## feature ID -> taxon label
lab <- setNames(feature.2.tax.16s$Taxon, feature.2.tax.16s$Feature.ID)
lab <- lab[taxa_names(physeq16s)]

## drop unassigned/bacteria-only features from table and tree
keep <- !is.na(lab) & lab != "16s:Unassigned" & lab != "16s:d__Bacteria"
ps16s <- prune_taxa(names(lab)[keep], physeq16s)

## convert to a phylo class which is more useful downstream
tree.16s <- phy_tree(ps16s, errorIfNULL=TRUE)
tree.16s$tip.label <- unname(lab[tree.16s$tip.label])

# peek at the tree to catch any glaring issues
plot(tree.16s, show.tip.label = FALSE)
save(physeq16s, ps16s, tree.16s, file= "../beeMicrobes/data/physeq16s.Rdata")

## ***********************************************************************
## 16s networks
## ***********************************************************************

## samples x taxa, ASVs with the same label summed
asv <- as(otu_table(ps16s), "matrix")
if (taxa_are_rows(ps16s)) asv <- t(asv)
indiv.comm.16s <- bipartite::empty(t(rowsum(t(asv), group = unname(lab[colnames(asv)]))))

## saving this out so we can make the tree visualizations
save(indiv.comm.16s, file="../beeMicrobes/data/indiv.comm16s.Rdata")

# Convert communities to relative abundance
finalASVtable <- indiv.comm.16s/rowSums(indiv.comm.16s)

## table and tree labels should match
stopifnot(all(colnames(finalASVtable) %in% tree.16s$tip.label))

save(finalASVtable, file= "../beeMicrobes/data/presAbsTable.Rdata")

merged.comm.16s <- as.data.frame(finalASVtable)

## ***********************************************************************
## 16s tree for PD
## ***********************************************************************

## one tip per taxon: keep the most abundant ASV
rep.ids <- names(sort(taxa_sums(ps16s), decreasing = TRUE))
rep.ids <- rep.ids[!duplicated(lab[rep.ids])]
tree.16s <- drop.tip(phy_tree(ps16s), setdiff(taxa_names(ps16s), rep.ids))
tree.16s$tip.label <- unname(lab[tree.16s$tip.label])

## RBCL commented out until the finalized dataset arrives
# ## ***********************************************************************
# ## RBCL - Read-in data
# ## ***********************************************************************

# # This section runs the same general processes as above, but combines data
# #   directly from the original sunflowers and sky islands projects. These data
# #   are split up across several sequencing runs, requiring somewhat repetitive
# #   processing to prepare the RBCL data to be combined with the bee specimen data.

# ### Sunflowers data

# # change working directory to the sunflowers saved folder
# sf.wd.path <- file.path(local.path, "sunflower_saved"); setwd(sf.wd.path)

# taxonomyRBCLR0 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza")
# taxonomyRBCLR0 <- taxonomyRBCLR0$data

# taxonomyRBCLR1 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
# taxonomyRBCLR1 <- taxonomyRBCLR1$data

# taxonomyRBCLR2 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza")
# taxonomyRBCLR2 <- taxonomyRBCLR2$data

# taxonomyRBCLR3 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza")
# taxonomyRBCLR3 <- taxonomyRBCLR3$data

# taxonomyRBCLR4 <- read_qza("ffar_pipeline/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza")
# taxonomyRBCLR4 <- taxonomyRBCLR4$data

# feature.2.tax.rbcl <- read.table("ffar_pipeline/merged/RBCL/taxonomyRBCL.txt",
#                                  sep="\t", header=TRUE)

# ## ## R knows to ignore the row in the txt with a # (see defaults read.table)

# ## add 16s/rbcl to make the columns easy to find when they are added
# ## to the specimen data

# feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", feature.2.tax.rbcl$Taxon)
# feature.2.tax.rbcl$Taxon <- paste("RBCL", feature.2.tax.rbcl$Taxon,
#                                   sep=':')

# ## ***********************************************************************
# ## rbcl networks
# ## ***********************************************************************

# indiv.comm.rbclR0 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR0,
#                                       feature.2.tax.rbcl)))
# indiv.comm.rbclR0 <- indiv.comm.rbclR0/rowSums(indiv.comm.rbclR0)

# indiv.comm.rbclR1 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR1,
#                                       feature.2.tax.rbcl)))
# indiv.comm.rbclR1 <- indiv.comm.rbclR1/rowSums(indiv.comm.rbclR1)

# indiv.comm.rbclR2 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR2,
#                                       feature.2.tax.rbcl)))
# indiv.comm.rbclR2 <- indiv.comm.rbclR2/rowSums(indiv.comm.rbclR2)

# indiv.comm.rbclR3 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR3,
#                                       feature.2.tax.rbcl)))
# indiv.comm.rbclR3 <- indiv.comm.rbclR3/rowSums(indiv.comm.rbclR3)

# indiv.comm.rbclR4 <-
#   bipartite::empty(catchDups(makeComm(taxonomyRBCLR4,
#                                       feature.2.tax.rbcl)))
# indiv.comm.rbclR4 <- indiv.comm.rbclR4/rowSums(indiv.comm.rbclR4)

# ### Sky Islands RBCL processing

# # change working directory to the sunflowers saved folder
# si.wd.path <- file.path(local.path, "skyislands_saved"); setwd(si.wd.path)

# # R0
# si.taxonomyRBCLR0 <-
#   read_qza("SI_pipeline/merged/RBCL/final/core_metricsRBCLR0/rarefied_table.qza")
# si.taxonomyRBCLR0 <- si.taxonomyRBCLR0$data

# # R1
# si.taxonomyRBCLR1 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR1/rarefied_table.qza")
# si.taxonomyRBCLR1 <- si.taxonomyRBCLR1$data

# # R2
# si.taxonomyRBCLR2 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR2/rarefied_table.qza")
# si.taxonomyRBCLR2 <- si.taxonomyRBCLR2$data

# # R3
# si.taxonomyRBCLR3 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR3/rarefied_table.qza")
# si.taxonomyRBCLR3 <- si.taxonomyRBCLR3$data

# # R4
# si.taxonomyRBCLR4 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR4/rarefied_table.qza")
# si.taxonomyRBCLR4 <- si.taxonomyRBCLR4$data

# # R5
# si.taxonomyRBCLR5 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR5/rarefied_table.qza")
# si.taxonomyRBCLR5 <- si.taxonomyRBCLR5$data

# # R6
# si.taxonomyRBCLR6 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR6/rarefied_table.qza")
# si.taxonomyRBCLR6 <- si.taxonomyRBCLR6$data

# # R7
# si.taxonomyRBCLR7 <-
#   read_qza("SI_pipeline/R2023/merged/RBCL/final/core_metricsRBCLR7/rarefied_table.qza")
# si.taxonomyRBCLR7 <- si.taxonomyRBCLR7$data

# si.feature.2.tax.rbcl <- read.table("SI_pipeline/R2023/merged/RBCL/taxonomyRBCL.txt",
#                                     sep="\t", header=TRUE)

# ## ## R knows to ignore the row in the txt with a # (see defaults read.table)

# ## add 16s/rbcl to make the columns easy to find when they are added
# ## to the specimen data

# si.feature.2.tax.rbcl$Taxon  <- gsub(" ", "_", si.feature.2.tax.rbcl$Taxon)
# si.feature.2.tax.rbcl$Taxon <- paste("RBCL", si.feature.2.tax.rbcl$Taxon, sep=':')

# ## ***********************************************************************
# ## rbcl networks
# ## ***********************************************************************

# #R0
# si.indiv.comm.rbclR0 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR0,
#                                       si.feature.2.tax.rbcl)))
# # changed to "Feature.ID" to match 'makeComm' function

# si.indiv.comm.rbclR0 <- si.indiv.comm.rbclR0/rowSums(si.indiv.comm.rbclR0)

# ## R1
# si.indiv.comm.rbclR1 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR1,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR1 <- si.indiv.comm.rbclR1/rowSums(si.indiv.comm.rbclR1)

# ## R2
# si.indiv.comm.rbclR2 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR2,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR2 <- si.indiv.comm.rbclR2/rowSums(si.indiv.comm.rbclR2)

# ## R3
# si.indiv.comm.rbclR3 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR3,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR3 <- si.indiv.comm.rbclR3/rowSums(si.indiv.comm.rbclR3)

# ## R4
# si.indiv.comm.rbclR4 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR4,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR4 <- si.indiv.comm.rbclR4/rowSums(si.indiv.comm.rbclR4)

# ## R5
# si.indiv.comm.rbclR5 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR5,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR5 <- si.indiv.comm.rbclR5/rowSums(si.indiv.comm.rbclR5)

# ## R6
# si.indiv.comm.rbclR6 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR6,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR6 <- si.indiv.comm.rbclR6/rowSums(si.indiv.comm.rbclR6)

# ## R7
# si.indiv.comm.rbclR7 <-
#   bipartite::empty(catchDups(makeComm(si.taxonomyRBCLR7,
#                                       si.feature.2.tax.rbcl)))

# si.indiv.comm.rbclR7 <- si.indiv.comm.rbclR7/rowSums(si.indiv.comm.rbclR7)


# ### Combine RBCL data from Sunflower and Sky Islands projects

# # label each rowname with an project-specific prefix
# bees.rbcl <- c(rownames(indiv.comm.rbclR0),
#                paste0("SF", rownames(indiv.comm.rbclR1)),
#                paste0("SF", rownames(indiv.comm.rbclR2)),
#                paste0("SF", rownames(indiv.comm.rbclR3)),
#                paste0("SF", rownames(indiv.comm.rbclR4)),
#                rownames(si.indiv.comm.rbclR0),
#                rownames(si.indiv.comm.rbclR1),
#                rownames(si.indiv.comm.rbclR2),
#                rownames(si.indiv.comm.rbclR3), 
#                rownames(si.indiv.comm.rbclR4),
#                rownames(si.indiv.comm.rbclR5),
#                rownames(si.indiv.comm.rbclR6),
#                rownames(si.indiv.comm.rbclR7)
#                )

# # generate a list of rbcl communities
# comms <- list(sfR0 = indiv.comm.rbclR0, 
#               sfR1 = indiv.comm.rbclR1,
#               sfR2 = indiv.comm.rbclR2, 
#               sfR3 = indiv.comm.rbclR3,
#               sfR4 = indiv.comm.rbclR4,
#               siR0 = si.indiv.comm.rbclR0, 
#               siR1 = si.indiv.comm.rbclR1,
#               siR2 = si.indiv.comm.rbclR2, 
#               siR3 = si.indiv.comm.rbclR3,
#               siR4 = si.indiv.comm.rbclR4, 
#               siR5 = si.indiv.comm.rbclR5,
#               siR6 = si.indiv.comm.rbclR6, 
#               siR7 = si.indiv.comm.rbclR7
#               )

# # bind rows of each community together
# # bind_rows generates appropriate NA values, these may need to be filled in w/zeros
# merged.comm.rbcl <- bind_rows(lapply(comms, as.data.frame))

# rownames(merged.comm.rbcl) <- bees.rbcl

## ***********************************************************************
## Make mega dataset
## ***********************************************************************
# Key data objects include:
# tree16s = the cleaned, merged ASV-level 16s phylogeny
# merged.comm.16s = the cleaned, merged microbial community matrix
# spec.net = specimen data, this gets redefined as the...
#   merged "spec.net" data from the sunflower and sky islands datasets...
#   Thus, this data will come from two csv files that are read-in below.

# Change working directory again
setwd(file.path(local.path,"beeMicrobes_saved"))

# Merge together specimen data from each project
spec.net <- read.csv("./beeMicrobes_data/specimen_data/spec.net_traits.csv")

# ## SF prefix check: rows matching UniqueID.fixed per run, as is vs. with "SF"
# sapply(comms, function(x)
#   c(n = nrow(x),
#     as.is = sum(rownames(x) %in% spec.net$UniqueID.fixed),
#     add.SF = sum(paste0("SF", rownames(x)) %in% spec.net$UniqueID.fixed)))
## same for 16s, as is vs. "SF" removed
c(n = nrow(merged.comm.16s),
  as.is = sum(rownames(merged.comm.16s) %in% spec.net$UniqueID.fixed),
  drop.SF = sum(gsub("^SF", "", rownames(merged.comm.16s)) %in% spec.net$UniqueID.fixed))
# ## rbcl rows that won't merge
# bees.rbcl[!bees.rbcl %in% spec.net$UniqueID.fixed]

#16s
indiv.comm.16s <- as.data.frame(merged.comm.16s) # convert community matrix to dataframe
bact <- colnames(indiv.comm.16s) # list of column names for matching
indiv.comm.16s$UniqueID <- rownames(indiv.comm.16s)

# #RBCL
# indiv.comm.rbcl <- as.data.frame(merged.comm.rbcl)  # convert community matrix to dataframe
# pollen <- colnames(indiv.comm.rbcl) # list of column names for matching
# indiv.comm.rbcl$UniqueID <- rownames(indiv.comm.rbcl)

# Merge
spec.net <-cbind(spec.net,
                 indiv.comm.16s[, bact][match(spec.net$UniqueID.fixed, indiv.comm.16s$UniqueID),])

# spec.net <-cbind(spec.net,
#                  indiv.comm.rbcl[, pollen][match(spec.net$UniqueID.fixed, indiv.comm.rbcl$UniqueID),])

## FULL MICROBE DATASET

## pull out 16s columns
microbes <- colnames(spec.net)[grepl("16s:", colnames(spec.net))] 
microbes <- microbes[microbes %in% tree.16s$tip.label]

## check which only have NAs in these columns (not screened) and drop them
screened.microbes <- apply(spec.net, 1, function(x) all(is.na(x[microbes])))
spec.microbes <- spec.net[!screened.microbes, ]

#3 rows have 0 for all microbes, need to drop
spec.microbes <- spec.microbes[rowSums(spec.microbes[,microbes])!=0,]

## QUESTION: should include root = TRUE? if false gives warning 3x
## warning: Rooted tree and include.root=TRUE argument required to calculate PD of single-species communities. Single species community assigned PD value of NA.

## phylogenetic distance function, modified from picante
PD <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- picante::prune.sample(t(this.bee), tree.16s)
  picante::pd(t(this.bee), this.tree, include.root = TRUE)
  #browser()
})

## phylogenetic species evenness for each bee
PSE <- apply(spec.microbes[,microbes], 1, function(x){
  this.bee <- x[x > 0]
  this.tree <- picante::prune.sample(t(this.bee), tree.16s)
  picante::pse(t(this.bee), this.tree)
  #browser()
})

## psd() phylogenetic species diversity metrics/phylogenetic species richness (PSR)
## ses.pd standardized effect size of phylo diversity in communities 

# Bind together PD list-class object to combine with spec.microbes dataframe
PD <- do.call(rbind, PD)
spec.microbes <- cbind(spec.microbes, PD)

# Bind together PSE list-class object to combine with spec.microbes dataframe
PSE <- do.call(rbind, PSE) %>%
  select(PSEs)
spec.microbes <- cbind(spec.microbes, PSE)

## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)

# add column for whether the sample has been screened for microbes or not
spec.net$ScreenedMicrobes <- ifelse(spec.net$Site %in% unique(spec.microbes$Site), 1, 0)

## change microbe NAs to 0
spec.net <- spec.net %>%
  mutate(PD = replace_na(PD, 0),
         SR = replace_na(SR, 0),
         PSEs = replace_na(PSEs, 0))

# Replace NA values with zeros in the microbe columns
spec.net[,microbes][is.na(spec.net[,microbes])] <- 0

## OBLIGATE AND TRANSIENT MICROBES DATASETS

## PD and SR for a subset of microbe columns, one row per bee
calcPD <- function(comm, taxa, tree){
  out <- apply(comm[, taxa, drop=FALSE], 1, function(x){
    if (!any(x > 0)) return(data.frame(PD = 0, SR = 0))  # no taxa in this subset
    tryCatch({
      this.bee <- x[x > 0]
      this.tree <- picante::prune.sample(t(this.bee), tree)
      pd_value <- picante::pd(t(this.bee), this.tree, include.root = TRUE)
      data.frame(PD = pd_value[[1]], SR = pd_value[[2]])
    }, error = function(e) {
      if (grepl("Tree has no branch lengths", e$message)) {
        return(data.frame(PD = 0, SR = 0))
      } else {
        stop(e)
      }
    })
  })
  do.call(rbind, out)
}

## obligate bee microbes: all (Zheng and Moran + Voulgari Kokota et al),
## social (Zheng and Moran), solitary (Voulgari Kokota et al)
obligate.sets <- list(
  all      = "Lactobacillaceae|Bifidobacteriaceae|Neisseriaceae|Orbaceae|Bartonella|Acetobacteraceae|Bacillaceae|Burkholderiaceae|Clostridiaceae|Comamonadaceae|Enterobacteriaceae|Lachnospiraceae|Methylobacteriaceae|Moraxellaceae|Sphingomonadaceae|Oxalobacteraceae|Erwiniaceae|Anaplasmataceae|Hafniaceae",
  social   = "Lactobacillaceae|Bifidobacteriaceae|Neisseriaceae|Orbaceae|Bartonella|Acetobacteraceae|Erwiniaceae|Anaplasmataceae|Hafniaceae",
  solitary = "Acetobacteraceae|Bacillaceae|Burkholderiaceae|Clostridiaceae|Comamonadaceae|Enterobacteriaceae|Lachnospiraceae|Lactobacillaceae|Methylobacteriaceae|Moraxellaceae|Sphingomonadaceae|Oxalobacteraceae|Erwiniaceae|Anaplasmataceae|Hafniaceae"
)

## obligate and transient PD/SR for each set
for (set in names(obligate.sets)){
  is.obligate <- grepl(obligate.sets[[set]], microbes)
  sfx <- if (set == "all") "" else paste0(".", set)   # PD.obligate, PD.obligate.social, ...
  for (type in c("obligate", "transient")){
    taxa <- if (type == "obligate") microbes[is.obligate] else microbes[!is.obligate]
    res <- calcPD(spec.microbes, taxa, tree.16s)
    colnames(res) <- paste0(c("PD.", "SR."), type, sfx)
    spec.microbes <- cbind(spec.microbes, res)
  }
}

## Merge back onto specimen data
spec.net <- merge(spec.net, spec.microbes, all.x=TRUE)

## change microbe NAs to 0
pd.cols <- grep("^(PD|SR)\\.(obligate|transient)", colnames(spec.net), value=TRUE)
spec.net <- spec.net %>%
  mutate(across(all_of(pd.cols), ~replace_na(.x, 0)))
spec.net[,microbes][is.na(spec.net[,microbes])] <- 0

# ## RBCL columns from the SI spec
# load("../beeMicrobes_saved/beeMicrobes_data/rds_files/SI_files/spec_RBCL_16s.Rdata")
# 
# ## drop bees and RBCL columns with no data
# pollen <- colnames(spec)[startsWith(colnames(spec), "RBCL:")]
# rbcl <- spec[, pollen, drop = FALSE]
# rbcl[is.na(rbcl)] <- 0
# keep.rows <- rowSums(rbcl) > 0
# keep.cols <- colSums(rbcl) > 0
# indiv.comm.rbcl <- rbcl[keep.rows, keep.cols, drop = FALSE]
# pollen <- colnames(indiv.comm.rbcl)
# indiv.comm.rbcl$UniqueID <- spec$UniqueID[keep.rows]
# any(duplicated(indiv.comm.rbcl$UniqueID))
# 
# ## rbcl bees matched in spec.net
# c(matched = sum(indiv.comm.rbcl$UniqueID %in% spec.net$UniqueID.fixed),
#   n       = nrow(indiv.comm.rbcl))
# 
# ## Merge
# spec.net <- cbind(spec.net,
#                   indiv.comm.rbcl[, pollen][match(spec.net$UniqueID.fixed, indiv.comm.rbcl$UniqueID), ])

## drop duplicate sampleID. LP: This should not happen?
## spec.net <- spec.net[!duplicated(spec.net$UniqueID),]
any(duplicated(spec.net$UniqueID.fixed))

# remove rows with no microbes
spec.net <- spec.net[rowSums(spec.net[, microbes], na.rm = TRUE) > 0, ]

save(spec.net, file= "../beeMicrobes/data/spec_RBCL_16s.Rdata")

write.csv(spec.net, file= "../beeMicrobes/data/spec_RBCL_16s.csv",
          row.names=FALSE)

save(tree.16s,
     file="../beeMicrobes/data/trees.Rdata")
