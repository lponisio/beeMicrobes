## setwd('/Volumes/bombus/Dropbox (University of Oregon)/beeMicrobes')

## Bee trait distinctness and diet breadth may drive microbiome and pathobiome dynamics
## We build formulas for the following response variables: microbiome distincess, microbiome diversity, parasite distincess, and parasite richness
## With the following explanatory variables: individual-level diet breadth (RBCl degree), species-level diet breadth (r.degree), and trait distinctness
## Note that what we call distinctiness in the manuscript is calculated as "originality" in Coux et al. 2016

rm(list=ls())
source("src/init_bayes.R")
load('data/allNetSums.RData')

## change to the number of cores you would liek to run on

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
                    control = list(adapt_delta = 0.99))

mcmc_trace(mod.par.div)
ggsave("figures/diagnostics/parDivBayesDiag.pdf",
       height=11, width=8.5)
write.ms.table(mod.par.div, "parDiv")

### micro originality
mod.micro.orig  <- brm(func.formulas[[2]],  data = all.indiv.mets,
                       cores=ncores,
                       iter = 10^5,
                       chains = 3,
                       control = list(adapt_delta = 0.99))

mcmc_trace(mod.micro.orig)
ggsave("figures/bayesMods/microOrigBayesDiag.pdf",
       height=11, width=8.5)

write.ms.table(mod.micro.orig, "microOrig")

## micro diversity
mod.micro.div  <- brm(func.formulas[[3]],  data = all.indiv.mets,
                      cores=ncores,
                      iter = 10^5,
                      chains = 3,
                      control = list(adapt_delta = 0.99))

mcmc_trace(mod.micro.div)
ggsave("figures/diagnostics/microDivBayesDiag.pdf",
       height=11, width=8.5)

write.ms.table(mod.micro.div, "micrDiv")

###########################################################################
## plotting
###########################################################################
sp.cols <- viridis(length(unique(all.indiv.mets$GenusSpecies)))

sp.means <- all.indiv.mets  %>%
    group_by(GenusSpecies, originality, Site, r.degree) %>%
    summarize(Parasite_originality = median(Parasite_originality,
                                            na.rm=TRUE),
              Micro_originality =  median(Micro_originality,
                                          na.rm=TRUE),
              Micro_partner.diversity =  median(Micro_partner.diversity,
                                            na.rm=TRUE),
              )

## originality parasites
dd.orig <- expand.grid(originality=seq(
                           from= min(all.indiv.mets$originality,
                                     na.rm=TRUE),
                           to= max(all.indiv.mets$originality,
                                   na.rm=TRUE),
                           length=35),
                       RBCL_degree=mean(all.indiv.mets$RBCL_degree,
                                        na.rm=TRUE),
                       r.degree=mean(all.indiv.mets$r.degree,
                                     na.rm=TRUE))

fe_only <- dd.orig  %>%
    add_fitted_draws(mod.par.orig,
                     re_formula = NA,
                     scale = "response", n = 1e3)

fe_only_mean <- fe_only %>%
    group_by(originality) %>%
    summarize(.value = mean(.value))


parasite1 <- ggplot(fe_only,
       aes(x = originality, y = .value)) +
    stat_interval(alpha = 0.5) +
    geom_line(data = fe_only_mean,
              color = "red", lwd = 2) +
    labs(y= "Pathobiome distinctness",
         x = "Trait distinctness") +
    scale_color_viridis(discrete=TRUE) +
    geom_point(data = sp.means,
               mapping = aes(x = originality, y = Parasite_originality,
                             color =  GenusSpecies)) +
    guides(color = FALSE) +
    ylim(quantile(sp.means$Parasite_originality, c(0, .95),
       na.rm=TRUE)) +
    theme(axis.title.x = element_text(size=14, vjust=-2),
          axis.title.y = element_text( size=14, vjust=3),
          plot.margin = unit(c(1,0,1,1), "cm"),
          text = element_text(size=18)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

ggsave("figures/bayesMods/parOrigBayes.pdf", height=4, width=6)


###########################################################################
## originality microbiome

## rbcl
dd.rbcl <- expand.grid(RBCL_degree=seq(
                           from= min(all.indiv.mets$RBCL_degree,
                                     na.rm=TRUE),
                           to= max(all.indiv.mets$RBCL_degree,
                                   na.rm=TRUE),
                           length=35),
                       originality=mean(all.indiv.mets$originality,
                                        na.rm=TRUE),
                       r.degree=mean(all.indiv.mets$r.degree,
                                     na.rm=TRUE))


fe_only <- dd.rbcl  %>%
    add_fitted_draws(mod.micro.orig,
                     re_formula = NA,
                     scale = "response", n = 1e3)

fe_only_mean <- fe_only %>%
    group_by(RBCL_degree) %>%
    summarize(.value = mean(.value))


micro1 <- ggplot(fe_only,
       aes(x = RBCL_degree, y = .value)) +
    stat_interval(alpha = 0.5) +
    geom_line(data = fe_only_mean,
              color = "red", lwd = 2) +
    labs(y= "Microbiome distinctness",
         x = "Individual-level diet breadth") +
    scale_color_viridis(discrete=TRUE) +
    geom_point(data = all.indiv.mets,
               mapping = aes(x = RBCL_degree, y = Micro_originality,
                             color =  GenusSpecies)) +
    guides(color = FALSE) +
    ylim(quantile(sp.means$Micro_originality, c(0, .95),
       na.rm=TRUE)) +
    theme(axis.title.x = element_text(size=14, vjust=-2),
          axis.title.y = element_text( size=14, vjust=3),
          plot.margin = unit(c(1,0,1,1), "cm"),
          text = element_text(size=18)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))


ggsave("figures/bayesMods/MicroOrigRBCLBayes.pdf", height=4, width=6)


### degree
dd.degree <- expand.grid(r.degree=seq(
                             from= min(all.indiv.mets$r.degree,
                                       na.rm=TRUE),
                             to= max(all.indiv.mets$r.degree,
                                     na.rm=TRUE),
                             length=35),
                         originality=mean(all.indiv.mets$originality,
                                          na.rm=TRUE),
                         RBCL_degree=mean(all.indiv.mets$RBCL_degree,
                                          na.rm=TRUE))


fe_only <- dd.degree  %>%
    add_fitted_draws(mod.micro.orig,
                     re_formula = NA,
                     scale = "response", n = 1e3)

fe_only_mean <- fe_only %>%
    group_by(r.degree) %>%
    summarize(.value = mean(.value))


micro2 <- ggplot(fe_only,
       aes(x = r.degree, y = .value)) +
    stat_interval(alpha = 0.5) +
    geom_line(data = fe_only_mean,
              color = "red", lwd = 2) +
    labs(y= "Microbiome distinctness",
         x = "Species-level diet breadth") +
    scale_color_viridis(discrete=TRUE) +
    geom_point(data = all.indiv.mets,
               mapping = aes(x = r.degree, y = Micro_originality,
                             color =  GenusSpecies)) +
    guides(color = FALSE) +
    ylim(quantile(sp.means$Micro_originality, c(0, .95),
       na.rm=TRUE)) +
    theme(axis.title.x = element_text(size=14, vjust=-2),
          axis.title.y = element_text( size=14, vjust=3),
          plot.margin = unit(c(1,0,1,1), "cm"),
          text = element_text(size=18)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

ggsave("figures/bayesMods/MicroOrigDegreeBayes.pdf", height=4, width=6)

###########################################################################
## partner div microbiome

dd.div <- expand.grid(RBCL_degree=seq(
                          from= min(all.indiv.mets$RBCL_degree,
                                    na.rm=TRUE),
                          to= max(all.indiv.mets$RBCL_degree,
                                  na.rm=TRUE),
                          length=35),
                      originality=mean(all.indiv.mets$originality,
                                       na.rm=TRUE),
                      r.degree=mean(all.indiv.mets$r.degree,
                                    na.rm=TRUE))

fe_only <- dd.div  %>%
    add_fitted_draws(mod.micro.div,
                     re_formula = NA,
                     scale = "response", n = 1e3)

fe_only_mean <- fe_only %>%
    group_by(RBCL_degree) %>%
    summarize(.value = mean(.value))


micro3 <- ggplot(fe_only,
       aes(x = RBCL_degree, y = .value)) +
    stat_interval(alpha = 0.5) +
    geom_line(data = fe_only_mean,
              color = "red", lwd = 2) +
    labs(y= "Microbiome diversity",
         x = "Individual-level diet breadth") +
       scale_color_viridis(discrete=TRUE) +
    geom_point(data = all.indiv.mets,
               mapping = aes(x = RBCL_degree, y = Micro_partner.diversity,
                             color =  GenusSpecies)) +
    guides(color = FALSE) +
    ylim(quantile(sp.means$Micro_partner.diversity, c(0, .95),
       na.rm=TRUE)) +
    theme(axis.title.x = element_text(size=14, vjust=-2),
          axis.title.y = element_text( size=14, vjust=3),
          plot.margin = unit(c(1,0,1,1), "cm"),
          text = element_text(size=18)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))


ggsave("figures/bayesMods/MicroDivRBCLBayes.pdf", height=4, width=6)


all.micro <- grid.arrange(micro1, micro3, nrow = 2)

ggsave("figures/bayesMods/MicroBayes.pdf", all.micro, height=8, width=6)
