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

bee.fams <- c("Apidae", "Halictidae", "Megachilidae")

all.indiv.mets <- all.indiv.mets[all.indiv.mets$Family %in% bee.fams,]

## *************************************************************
## change to the number of cores you would like to run on
ncores <- 10
## *************************************************************

all.indiv.mets$PossibleParasite <- 5

all.indiv.mets$GenusSpecies2 <- all.indiv.mets$GenusSpecies

yvars <- c("Parasite_originality",
           "ParasiteRichness | vint(PossibleParasite)",
           "Micro_originality",
           "Micro_partner.diversity")

xvars <-      c("scale(RBCL_degree)",
                "scale(originality)", "scale(r.degree)",
                "(1|gr(GenusSpecies, cov = co.var.mat))",
                "(1|Site)", "(1|GenusSpecies2)")

x.form <- paste(paste(xvars,  collapse="+"))

func.formulas <-lapply(yvars, function(y) {
    as.formula(paste(y, "~", x.form))
})

names(func.formulas) <- yvars

all.indiv.mets <- all.indiv.mets[all.indiv.mets$GenusSpecies %in%
                                 rownames(co.var.mat),]

## parasite originality
mod.par.orig  <- brm(func.formulas[[1]],
                     data = all.indiv.mets,
                     data2 = list(co.var.mat = co.var.mat),
                     cores=ncores,
                     family = gaussian(),
                     iter = 10^5,
                     chains = 3,
                     inits=0,
                     control = list(adapt_delta = 0.99))

save(mod.par.orig,
     file="saved/parasiteOrigFit.Rdata")
phyloHyp(mod.par.orig, "parOrig")

write.ms.table(mod.par.orig, "parOrig")


## parasite richness

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

mod.par.div  <- brm(func.formulas[[2]],
                    family = beta_binomial2, stanvars = stanvars,
                    data = all.indiv.mets,
                    data2 = list(co.var.mat = co.var.mat),
                    cores=1,
                    iter = 10^5,
                    chains = 3,
                    inits=0,
                    control = list(adapt_delta = 0.99))

save(mod.par.div,
     file="saved/parasiteDivFit.Rdata")
phyloHyp(mod.par.div, "parDiv",  mod.gaussian=FALSE)

write.ms.table(mod.par.div, "parDiv")

### micro originality
mod.micro.orig  <- brm(func.formulas[[3]],  data = all.indiv.mets,
                       cores=ncores,
                       data2 = list(co.var.mat = co.var.mat),
                       iter = 10^5,
                       chains = 3,
                       inits=0,
                       control = list(adapt_delta = 0.99))

save(mod.micro.orig,
     file="saved/microbeOrigFit.Rdata")
phyloHyp(mod.micro.orig, "microOrig")

write.ms.table(mod.micro.orig, "microOrig")

## micro diversity
mod.micro.div  <- brm(func.formulas[[3]],  data = all.indiv.mets,
                      cores=ncores,
                      data2 = list(co.var.mat = co.var.mat),
                      iter = 10^5,
                      chains = 3,
                      inits=0,
                      control = list(adapt_delta = 0.99))

save(mod.micro.div,
     file="saved/microbeDivFit.Rdata")
phyloHyp(mod.micro.div, "microDiv",  mod.gaussian=FALSE)

write.ms.table(mod.micro.div, "microDiv")

