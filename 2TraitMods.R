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

## *************************************************************
## change to the number of cores you would like to run on
ncores <- 1
## *************************************************************

all.indiv.mets$PossibleParasite <- 5

all.indiv.mets$GenusSpecies2 <- all.indiv.mets$GenusSpecies

all.indiv.mets$onSF <- all.indiv.mets$PlantGenusSpecies
all.indiv.mets$onSF[all.indiv.mets$onSF != "Helianthus annuus"] <-
  "not sf"


all.indiv.mets$Nest <- paste(all.indiv.mets$NestPartitions, all.indiv.mets$NestLoc)

yvars <- c("Parasite_originality",
           "ParasiteRichness | vint(PossibleParasite)",
           "Micro_originality",
           "Micro_partner.diversity")

## for model with trait uniqueness (combination of all traits)
xvars.1 <-      c("scale(RBCL_degree)",
                "scale(originality)", "scale(r.degree)",
                "(1|gr(GenusSpecies, cov = co.var.mat))",
                "onSF",
                "Genus",
                "(1|Site)",
                "(1|GenusSpecies2)")


## three sets of traits that can be put together in a model without
## being strongly colinear
xvars.traits2 <-      c("scale(RBCL_degree)",
                       "scale(r.degree)",
                       "Sociality",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)",
                       "(1|GenusSpecies2)")

xvars.traits3 <-      c("scale(RBCL_degree)",
                       "scale(r.degree)",
                       "scale(MeanITD)",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)",
                       "(1|GenusSpecies2)")

xvars.traits4 <-      c("scale(RBCL_degree)",
                       "scale(r.degree)",
                       "Nest",
                       "(1|gr(GenusSpecies, cov = co.var.mat))",
                       "(1|Site)",
                       "(1|GenusSpecies2)")

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

all.indiv.mets <- all.indiv.mets[all.indiv.mets$GenusSpecies %in%
                                 rownames(co.var.mat),]

## function for running brms models with different xvar sets and y variables
runmodels <- function(mod.set, yvar,
                      indiv.mets,
                      phylo.mat,
                      runcores=ncores,
                      niter=10^4,
                      nchains=1,
                      name.file="set1",
                      mod.gaussian=TRUE){
  
  fit  <- brm(mod.set[[yvar]],
              data = indiv.mets,
              data2 = list(co.var.mat = phylo.mat),
              cores=runcores,
              family = gaussian(),
              iter = niter,
              chains = nchains,
              init= 0,
              control = list(adapt_delta = 0.99))
  ## r2
  br2 <- bayes_R2(fit)
  loor2 <- loo_R2(fit)

  save(fit, br2, loor2, indiv.mets,
       file=sprintf("saved/%s_%s.Rdata", yvar, name.file))

  write.ms.table(fit, sprintf("%s_%s", yvar, name.file))

  ## phylogenic signal
  phyloHyp(fit, sprintf("%s_%s", yvar, name.file), mod.gaussian=mod.gaussian)

  return(list(fit=fit, br2=br2, loor2=loor2))
}


## ************************************************************************
## parasite originality
## ************************************************************************
par.orig.set1 <- runmodels(mod.set=mod.set1, yvar=yvars[1],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set1")
## posterior prdictive checks
plot(pp_check(par.orig.set1$fit, ndraws=100))


par.orig.set2 <- runmodels(mod.set=mod.set2, yvar=yvars[1],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set2")
## posterior prdictive checks
plot(pp_check(par.orig.set2$fit, ndraws=100))

par.orig.set3 <- runmodels(mod.set=mod.set3, yvar=yvars[1],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set3")
## posterior prdictive checks
plot(pp_check(par.orig.set3$fit, ndraws=100))

par.orig.set4 <- runmodels(mod.set=mod.set4, yvar=yvars[1],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set4")
## posterior prdictive checks
plot(pp_check(par.orig.set4$fit, ndraws=100))

## model selection
waic(par.orig.set1$fit,
     par.orig.set2$fit, par.orig.set3$fit, par.orig.set4$fit) 

loo(par.orig.set1$fit,
     par.orig.set2$fit, par.orig.set3$fit, par.orig.set4$fit) 

## asses r2s
par.orig.set1$br2
par.orig.set2$br2
par.orig.set3$br2
par.orig.set4$br2

## ************************************************************************
## microbe originality
## ************************************************************************
micro.orig.set1 <- runmodels(mod.set=mod.set1, yvar=yvars[3],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set1")

## posterior prdictive checks
plot(pp_check(micro.orig.set1$fit, ndraws=100))

micro.orig.set2 <- runmodels(mod.set=mod.set2, yvar=yvars[3],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set2")
## posterior prdictive checks
plot(pp_check(micro.orig.set2$fit, ndraws=100))

micro.orig.set3 <- runmodels(mod.set=mod.set3, yvar=yvars[3],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set3")
## posterior prdictive checks
plot(pp_check(micro.orig.set3$fit, ndraws=100))

micro.orig.set4 <- runmodels(mod.set=mod.set4, yvar=yvars[3],
                           indiv.mets=all.indiv.mets,
                           phylo.mat=co.var.mat,
                           name.file="set4")
## posterior prdictive checks
plot(pp_check(micro.orig.set4$fit, ndraws=100))

## model selection
waic(micro.orig.set1$fit,
     micro.orig.set2$fit, micro.orig.set3$fit, micro.orig.set4$fit) 

loo(micro.orig.set1$fit,
     micro.orig.set2$fit, micro.orig.set3$fit, micro.orig.set4$fit) 

## assess r2s
micro.orig.set1$br2
micro.orig.set2$br2
micro.orig.set3$br2
micro.orig.set4$br2

## ************************************************************************
## parasite richness
## ************************************************************************

beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"
stanvars <- stanvar(scode = stan_funs, block = "functions")

par.div.set1 <- runmodels(mod.set=mod.set1, yvar=yvars[2],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set1",
                          mod.gaussian=FALSE)

## posterior prdictive checks
plot(pp_check(par.div.set1$fit, ndraws=100))


par.div.set2 <- runmodels(mod.set=mod.set2, yvar=yvars[2],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set2",
                          mod.gaussian=FALSE)
## posterior prdictive checks
plot(pp_check(par.div.set2$fit, ndraws=100))

par.div.set3 <- runmodels(mod.set=mod.set3, yvar=yvars[2],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set3",
                          mod.gaussian=FALSE)
## posterior prdictive checks
plot(pp_check(par.div.set3$fit, ndraws=100))

par.div.set4 <- runmodels(mod.set=mod.set4, yvar=yvars[2],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set4",
                          mod.gaussian=FALSE)
## posterior prdictive checks
plot(pp_check(par.div.set4$fit, ndraws=100))

## model selection
waic(par.div.set1$fit,
     par.div.set2$fit, par.div.set3$fit, par.div.set4$fit) 

loo(par.div.set1$fit,
     par.div.set2$fit, par.div.set3$fit, par.div.set4$fit) 

## assess r2s
par.div.set1$br2
par.div.set2$br2
par.div.set3$br2
par.div.set4$br2

## ************************************************************************
## microbe diversity
## ************************************************************************

micro.div.set1 <- runmodels(mod.set=mod.set1, yvar=yvars[4],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set1",
                          mod.gaussian=FALSE)

## posterior prdictive checks
plot(pp_check(micro.div.set1$fit, ndraws=100))


micro.div.set2 <- runmodels(mod.set=mod.set2, yvar=yvars[4],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set2",
                          mod.gaussian=FALSE)
## posterior prdictive checks
plot(pp_check(micro.div.set2$fit, ndraws=100))

micro.div.set3 <- runmodels(mod.set=mod.set3, yvar=yvars[4],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set3",
                          mod.gaussian=FALSE)
## posterior prdictive checks
plot(pp_check(micro.div.set3$fit, ndraws=100))

micro.div.set4 <- runmodels(mod.set=mod.set4, yvar=yvars[4],
                          indiv.mets=all.indiv.mets,
                          phylo.mat=co.var.mat,
                          name.file="set4",
                          mod.gaussian=FALSE)
## posterior prdictive checks
plot(pp_check(micro.div.set4$fit, ndraws=100))

## model selection
waic(micro.div.set1$fit,
     micro.div.set2$fit, micro.div.set3$fit, micro.div.set4$fit) 

loo(micro.div.set1$fit,
     micro.div.set2$fit, micro.div.set3$fit, micro.div.set4$fit) 

## assess r2s
micro.div.set1$br2
micro.div.set2$br2
micro.div.set3$br2
micro.div.set4$br2
