library(vegan)
library(bipartite)
library(tidyr)
library(dplyr)
library(ecodist)
library(fields)
library(picante)

load('data/spec_RBCL_16s.Rdata')

source("src/misc.R")

parasites <- c("Apicystis", "Ascosphaera", "CrithidiaSpp",
               "CrithidiaBombi", "CrithidiaExpoeki",
               "NosemaCeranae", "NosemaBombi" )

no.apidae <- spec$UniqueID[spec$Apidae != 1
                           & spec$Gut ==1]

print(paste("not dropping w/o apidae", no.apidae))

spec.wild <- spec[spec$Family != "Syrphidae",]
spec.wild <- spec[spec$GenusSpecies != "Apis mellifera",]

load('data/trees.Rdata')

dir.create("figures", showWarnings = FALSE)
dir.create("figures/diagnostics", showWarnings = FALSE)
dir.create("figures/mods", showWarnings = FALSE)
dir.create("saved/tables", showWarnings = FALSE)
dir.create("saved/distmats", showWarnings = FALSE)

