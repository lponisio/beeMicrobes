## Bee trait distinctness and diet breadth may drive microbiome and
## pathobiome dynamics We build formulas for the following response
## variables: microbiome distincess, microbiome diversity, parasite
## distincess, and parasite richness

## With the following explanatory variables: individual-level diet
## breadth (RBCl degree), species-level diet breadth (r.degree), and
## trait distinctness

## Note that what we call distinctiness in the manuscript is
## calculated as "originality" in Coux et al. 2016

rm(list=ls())
source("src/init_bayes.R")
load('data/allNetSums.RData')
load('data/covarmatrix_community.RData')

all.indiv.mets$Family[all.indiv.mets$Genus == "Hylaeus"] <- "Colletidae"
bee.fams <- c("Apidae", "Halictidae", "Megachilidae", "Colletidae")
all.indiv.mets <- all.indiv.mets[all.indiv.mets$Family %in% bee.fams,]

all.indiv.mets$Apis <- "not Apis"
all.indiv.mets$Apis[all.indiv.mets$Genus == "Apis"] <-
  "Apis"

all.indiv.mets$PossibleParasite <- 5
all.indiv.mets$GenusSpecies2 <- all.indiv.mets$GenusSpecies
all.indiv.mets$onSF <- all.indiv.mets$PlantGenusSpecies
all.indiv.mets$onSF[all.indiv.mets$onSF != "Helianthus annuus"] <-
  "not sf"
all.indiv.mets$Nest <- paste(all.indiv.mets$NestPartitions,
                             all.indiv.mets$NestLoc)
all.indiv.mets <- all.indiv.mets[all.indiv.mets$GenusSpecies %in%
                                 rownames(co.var.mat),]

all.indiv.mets$Parasite_originality <-
  all.indiv.mets$Parasite_originality^(1/3)
all.indiv.mets$Micro_originality <-
  all.indiv.mets$Micro_originality^(1/3)

## change to the number of cores you would like to run on
ncores <- 1

## *************************************************************
## set up response and explanatory variable sets
## *************************************************************
yvars <- c("Parasite_originality",
           "ParasiteRichness | trials(PossibleParasite)",
           "Micro_originality",
           "Micro_partner.diversity")

## for model with trait uniqueness (combination of all traits)
xvars.1 <-      c("scale(RBCL_degree)",
                "scale(originality)", "scale(r.degree)",
                "(1|gr(GenusSpecies, cov = co.var.mat))",
                "onSF",
                "(1|Site)")

## three sets of traits that can be put together in a model without
## being strongly colinear
xvars.traits2 <-      c("scale(RBCL_degree)",
                       "scale(r.degree)",
                       "Sociality",
                       "onSF",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)")

xvars.traits3 <-      c("scale(RBCL_degree)",
                       "scale(r.degree)",
                       "scale(MeanITD)",
                       "onSF",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)")

xvars.traits4 <-      c("scale(RBCL_degree)",
                       "scale(r.degree)",
                       "NestLoc",
                       "onSF",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)")

xvars.5 <-      c("Genus",
                   "onSF",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)")

xvars.6 <-      c("Apis",
                  "onSF",
                  "(1|gr(GenusSpecies, cov = co.var.mat))",
                  "(1|Site)")


makeforms <- function(xvars, yvars){
  x.form <- paste(paste(xvars,  collapse="+"))

  func.formulas <-lapply(yvars, function(y) {
    as.formula(paste(y, "~", x.form))
  })
  names(func.formulas) <- yvars
  return(func.formulas)
}

mod.set1 <- makeforms(xvars.1, yvars)
mod.set2 <- makeforms(xvars.traits2, yvars)
mod.set3 <- makeforms(xvars.traits3, yvars)
mod.set4 <- makeforms(xvars.traits4, yvars)
mod.set5 <- makeforms(xvars.5, yvars)
mod.set6 <- makeforms(xvars.6, yvars)

sub.indiv.mets <- all.indiv.mets[!is.na(all.indiv.mets$r.degree) &
                                 !is.na(all.indiv.mets$RBCL_degree),]

## *************************************************************
## set up response and explanatory variable sets
## *************************************************************

## function for running brms models with different xvar sets and y variables
runmodels <- function(mod.set, yvar,
                      indiv.mets,
                      phylo.mat,
                      runcores=ncores,
                      niter=10^4,
                      nchains=1,
                      name.file="set1",
                      mod.gaussian=TRUE,
                      family=gaussian(),
                      calcr2 = TRUE,
                      ...){

  fit  <- brm(mod.set[[yvar]],
              data = indiv.mets,
              data2 = list(co.var.mat = phylo.mat),
              cores=runcores,
              family = family,
              iter = niter,
              chains = nchains,
              init= 0,
              control = list(adapt_delta = 0.999,
                             max_treedepth = 11),
              ...)
  ## r2
  if(calcr2){
    br2 <- bayes_R2(fit)
    loor2 <- loo_R2(fit)
  } else{
    br2 <- NA
    loor2 <- NA
  }

  save(fit, br2, loor2, indiv.mets,
       file=sprintf("saved/%s_%s.Rdata", yvar, name.file))

  write.ms.table(fit, sprintf("%s_%s", yvar, name.file))

  ## phylogenic signal
  try(phyloHyp(fit, sprintf("%s_%s", yvar, name.file),
           mod.gaussian=mod.gaussian), silent=TRUE)

  return(list(fit=fit, br2=br2, loor2=loor2))
}

## ************************************************************************
## parasite originality
## ************************************************************************

all.model.sets <- list(mod.set1, mod.set2, mod.set3,
                       mod.set4, mod.set5, mod.set6)
names(all.model.sets) <- paste0("set", 1:6)

runAllModels <- function(yvar.index, family.fun = student(),
                         all.model.sets, save.path="saved/mod-comparisons"){  
  model.sets <- lapply(names(all.model.sets), function(this.set){
    runmodels(mod.set=all.model.sets[[this.set]],
              yvar=yvars[yvar.index],
              indiv.mets=sub.indiv.mets,
              phylo.mat=co.var.mat,
              name.file=this.set,
              family=family.fun)

  })
  save(model.sets, file=file.path(save.path, sprintf("%s.Rdata",
                                                     yvars[yvar.index])))
  return(model.sets)
}

plotppCheck <- function(models){
  for(i in 1:length(models)){
    plot(pp_check(models[[i]]$fit, ndraws=100))
  }
}

extractR2 <- function(x){  
    x$br2
}

## parasite originality
par.orig.models <- runAllModels(yvar.index=1, family.fun = gaussian(),
                                all.model.sets=all.model.sets)
plotppCheck(par.orig.models)

par.orig.r2s <- do.call(rbind, lapply(par.orig.models, extractR2))
par.orig.r2s

## microbe originality
micro.orig.models <- runAllModels(yvar.index=3, family.fun = gaussian(),
                                  all.model.sets=all.model.sets)
plotppCheck(micro.orig.models)
micro.orig.r2s <- do.call(rbind, lapply(micro.orig.models, extractR2))


## parasite diversity
par.div.models <- runAllModels(yvar.index=2, family.fun = beta_binomial(),
                               all.model.sets=all.model.sets)
plotppCheck(par.div.models)
par.div.r2s <- do.call(rbind, lapply(par.div.models, extractR2))

## microbe diversity
micro.div.models <- runAllModels(yvar.index=4, family.fun = gaussian(),
                                 all.model.sets=all.model.sets)
plotppCheck(micro.div.models)
micro.div.r2s <- do.call(rbind, lapply(micro.div.models, extractR2))

## ************************************************************************
## parasite originality model selection
## ************************************************************************

load(file="saved/mod-comparisons/Parasite_originality.Rdata")

par.orig.waic <- waic(par.orig.models[[1]]$fit,
     par.orig.models[[2]]$fit, par.orig.models[[3]]$fit,
     par.orig.models[[4]]$fit,
     par.orig.models[[5]]$fit, par.orig.models[[6]]$fit)
par.orig.waic

par.orig.loo <- loo(par.orig.models[[1]]$fit,
     par.orig.models[[2]]$fit, par.orig.models[[3]]$fit,
     par.orig.models[[4]]$fit,
     par.orig.models[[5]]$fit, par.orig.models[[6]]$fit)
par.orig.loo

## expected log predictive density (elpd_loo)
## the estimated effective number of parameters (p_loo)
## the Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO; looic).

## lowests waic and loo is mod set 6, the next lowest is 3 (diff 1.1,
## se 1.1) and 1 (diff 1.4, se 1.2)

## ************************************************************************
## microbe originality model selection
## ************************************************************************

load(file="saved/mod-comparisons/Micro_originality.Rdata")

micro.orig.waic <- waic(micro.orig.models[[1]]$fit,
     micro.orig.models[[2]]$fit, micro.orig.models[[3]]$fit,
     micro.orig.models[[4]]$fit,
     micro.orig.models[[5]]$fit, micro.orig.models[[6]]$fit)
micro.orig.waic

micro.orig.loo <- loo(micro.orig.models[[1]]$fit,
     micro.orig.models[[2]]$fit, micro.orig.models[[3]]$fit,
     micro.orig.models[[4]]$fit,
     micro.orig.models[[5]]$fit, micro.orig.models[[6]]$fit)
micro.orig.loo

## model set 3 is best fit, model 1 is next (diff 0.7, se 1.0). Then
## model set 2 (diff 0.9, 1.3)  Model set 1,2,3 are therefore the top models.

## ************************************************************************
## parasite richness model selection
## ************************************************************************

load(file="saved/mod-comparisons/ParasiteRichness | trials(PossibleParasite).Rdata")


par.div.waic <- waic(par.div.models[[1]]$fit,
     par.div.models[[2]]$fit, par.div.models[[3]]$fit,
     par.div.models[[4]]$fit,
     par.div.models[[5]]$fit, par.div.models[[6]]$fit)
par.div.waic

par.div.loo <- loo(par.div.models[[1]]$fit,
     par.div.models[[2]]$fit, par.div.models[[3]]$fit,
     par.div.models[[4]]$fit,
     par.div.models[[5]]$fit, par.div.models[[6]]$fit)
par.div.loo

## model 6 is the best fit, but 3,5,1,2 all have diffs less than the
## se of the diff.
## ************************************************************************
## microbe diversity model selection
## ************************************************************************

load(file="saved/mod-comparisons/Micro_partner.diversity.Rdata")

micro.div.waic <- waic(micro.div.models[[1]]$fit,
     micro.div.models[[2]]$fit, micro.div.models[[3]]$fit,
     micro.div.models[[4]]$fit,
     micro.div.models[[5]]$fit, micro.div.models[[6]]$fit)
micro.div.waic

micro.div.loo <- loo(micro.div.models[[1]]$fit,
     micro.div.models[[2]]$fit, micro.div.models[[3]]$fit,
     micro.div.models[[4]]$fit,
     micro.div.models[[5]]$fit, micro.div.models[[6]]$fit)
micro.div.loo


## model 6 is the best fit, but -1.0 different from 3 (se 1.9, so
## greater than the difference between the models). Model 2 & 1 are also
## onlt 1.2 different from model 6 (se 1.4). So models 6,1,2,3 are the top
## models. 
