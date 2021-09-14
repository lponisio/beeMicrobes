## setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

rm(list=ls())
source("src/init_bayes.R")
load('data/allNetSums.RData')

all.indiv.mets$PossibleParasite <- 5

## run models to check VIF, other assessment
yvars <- c("Parasite_originality",
        "Micro_originality",
        "Micro_partner.diversity")

xvars <-      c("scale(RBCL_degree)",
               "scale(originality)", "scale(r.degree)",
                "(1|GenusSpecies)", "(1|Genus)", "(1|Site)")

x.form <- paste(paste(xvars,  collapse="+"))

func.formulas <-lapply(yvars, function(y) {
  as.formula(paste(y, "~", x.form))
})

names(func.formulas) <- yvars

func.mod  <-  lapply(func.formulas, lmer,
                     data = all.indiv.mets,
                     control = lmerControl(optimizer="bobyqa"))

## not quite the right model, p success changes
func.mod.bin <- glmer(cbind(Parasite_degree, PossibleParasite) ~
               scale(RBCL_degree) +
               scale(originality) + scale(r.degree) + (1|Genus) +
                (1|GenusSpecies) + (1|Site),
                      family=  "binomial",
               data = all.indiv.mets,
                control = glmerControl(optimizer="bobyqa")
               )

func.mod <- c(func.mod, func.mod.bin)
yvars <- c(yvars, "Parasite.degree")
names(func.mod) <- yvars

## model assessment

## all vif  < 2
lapply(func.mod, vif)
lapply(func.mod, r.squaredGLMM)
lapply(func.mod, check_heteroscedasticity)
