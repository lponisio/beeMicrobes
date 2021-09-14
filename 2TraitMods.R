## setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

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

## change to the number of cores you would like to run on

ncores <- 2

all.indiv.mets$PossibleParasite <- 5

yvars <- c("Parasite_originality",
           "ParasiteRichness | vint(PossibleParasite)",
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

## parasite originality
mod.par.orig  <- brm(func.formulas[[1]],  data = all.indiv.mets,
                     cores=ncores,
                     iter = 10^5,
                     chains = 3,
                     inits=0,
                     control = list(adapt_delta = 0.99))


save(mod.par.orig, all.indiv.mets,
     file="saved/parasiteOrigFit.Rdata")


mcmc_trace(mod.par.orig)
ggsave("figures/diagnostics/parOrigBayesDiag.pdf",
       height=11, width=8.5)

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
                    cores=ncores,
                    iter = 10^5,
                    chains = 3,
                    inits=0,
                    control = list(adapt_delta = 0.99))

save(mod.par.div, all.indiv.mets,
     file="saved/parasiteDivFit.Rdata")


mcmc_trace(mod.par.div)
ggsave("figures/diagnostics/parDivBayesDiag.pdf",
       height=11, width=8.5)
write.ms.table(mod.par.div, "parDiv")

### micro originality
mod.micro.orig  <- brm(func.formulas[[3]],  data = all.indiv.mets,
                       cores=ncores,
                       iter = 10^3,
                       chains = 3,
                       control = list(adapt_delta = 0.99))


save(mod.micro.orig, all.indiv.mets,
     file="saved/microbeOrigFit.Rdata")


mcmc_trace(mod.micro.orig)
ggsave("figures/bayesMods/microOrigBayesDiag.pdf",
       height=11, width=8.5)

write.ms.table(mod.micro.orig, "microOrig")

## micro diversity
mod.micro.div  <- brm(func.formulas[[3]],  data = all.indiv.mets,
                      cores=ncores,
                      iter = 10^3,
                      chains = 3,
                      control = list(adapt_delta = 0.99))

mcmc_trace(mod.micro.div)
ggsave("figures/diagnostics/microDivBayesDiag.pdf",
       height=11, width=8.5)

write.ms.table(mod.micro.div, "micrDiv")


save(mod.micro.div, all.indiv.mets,
     file="saved/microbeDivFit.Rdata")

