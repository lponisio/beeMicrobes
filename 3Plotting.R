rm(list=ls())
source("src/init_bayes.R")
load('data/allNetSums.RData')

library(gridExtra)
library(viridis)

load(file="saved/parasiteOrigFit.Rdata")
load(file="saved/parasiteDivFit.Rdata")
load(file="saved/microbeOrigFit.Rdata")
load(file="saved/MicrobeDivFit.Rdata")


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

