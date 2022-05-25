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

p1.parasite  <- mod.par.orig %>%
    spread_draws(b_Intercept,
                 b_scaleoriginality) %>%
    mutate(originality =
               list(seq(min(all.indiv.mets$originality, na.rm=TRUE),
                        max(all.indiv.mets$originality, na.rm=TRUE),
                        0.01))) %>%
    unnest(originality) %>%
    mutate(pred = (b_Intercept +
               b_scaleoriginality*originality))   %>%
               group_by(originality) %>%
               summarise(pred_m = mean(pred, na.rm = TRUE),
                         pred_low_95 = quantile(pred, prob = 0.025),
                         pred_high_95 = quantile(pred, prob = 0.975),
                         pred_low_90 = quantile(pred, prob = 0.05),
                         pred_high_90 = quantile(pred, prob = 0.95),
                         pred_low_85 = quantile(pred, prob = 0.075),
                         pred_high_85 = quantile(pred, prob = 0.925)) %>%
               ggplot(aes(x = originality, y = pred_m)) +
               geom_line() +
               geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                           fill="goldenrod") +
               geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                           fill="dodgerblue") +
               geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                           fill="darkolivegreen") +
               ylab("Pathobiome distinctness") +
               xlab("Host trait distinctness") +
               ylim(quantile(sp.means$Parasite_originality, c(0, .95),
                             na.rm=TRUE)) +
               xlim(range(all.indiv.mets$originality))  +
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     text = element_text(size=16),
                     legend.position = "none") +
               geom_point(data = sp.means,
                          mapping = aes(x = originality, y = Parasite_originality,
                                        color =  GenusSpecies)) +
    scale_color_viridis(discrete=TRUE)

ggsave("figures/parasite_originality.pdf",
       height=4, width=5)

p0.parasite  <- mod.par.orig %>%
    spread_draws(b_Intercept,
                 b_scaleRBCL_degree) %>%
    mutate(RBCL_degree =
               list(seq(min(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        max(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        0.01))) %>%
    unnest(RBCL_degree) %>%
    mutate(pred = (b_Intercept +
               b_scaleRBCL_degree*RBCL_degree))   %>%
               group_by(RBCL_degree) %>%
               summarise(pred_m = mean(pred, na.rm = TRUE),
                         pred_low_95 = quantile(pred, prob = 0.025),
                         pred_high_95 = quantile(pred, prob = 0.975),
                         pred_low_90 = quantile(pred, prob = 0.05),
                         pred_high_90 = quantile(pred, prob = 0.95),
                         pred_low_85 = quantile(pred, prob = 0.075),
                         pred_high_85 = quantile(pred, prob = 0.925)) %>%
               ggplot(aes(x = RBCL_degree, y = pred_m)) +
               geom_line() +
               geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                           fill="goldenrod") +
               geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                           fill="dodgerblue") +
               geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                           fill="darkolivegreen") +
               ylab("Pathobiome distinctness") +
               xlab("Individual diet breadth") +
               ## ylim(quantile(sp.means$Parasite_originality, c(0, .95),
               ##               na.rm=TRUE)) +
               xlim(range(all.indiv.mets$RBCL_degree))  +
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     text = element_text(size=16),
                     legend.position = "none") +
               geom_point(data = all.indiv.mets,
                          mapping = aes(x = RBCL_degree, y = Parasite_originality,
                                        color =  GenusSpecies)) +
    scale_color_viridis(discrete=TRUE)

ggsave("figures/parasite_RBCL_degree.pdf",
       height=4, width=5)



###########################################################################
## originality microbiome


p1.micro  <- mod.micro.orig %>%
    spread_draws(b_Intercept,
                 b_scaleRBCL_degree) %>%
    mutate(RBCL_degree =
               list(seq(min(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        max(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        0.01))) %>%
    unnest(RBCL_degree) %>%
    mutate(pred = (b_Intercept +
               b_scaleRBCL_degree*RBCL_degree))   %>%
               group_by(RBCL_degree) %>%
               summarise(pred_m = mean(pred),
                         pred_low_95 = quantile(pred, prob = 0.025),
                         pred_high_95 = quantile(pred, prob = 0.975),
                         pred_low_90 = quantile(pred, prob = 0.05),
                         pred_high_90 = quantile(pred, prob = 0.95),
                         pred_low_85 = quantile(pred, prob = 0.075),
                         pred_high_85 = quantile(pred, prob = 0.925)) %>%
               ggplot(aes(x = RBCL_degree, y = pred_m)) +
               geom_line() +
               geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                           fill="goldenrod") +
               geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                           fill="dodgerblue") +
               geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                           fill="darkolivegreen") +
               ylab("Microbiome distinctness") +
    xlab("Individual-level diet breadth") +
    ylim(0,6) +
    xlim(0,20) +
    ## ylim(quantile(sp.means$Micro_originality, c(0, .95),
    ##               na.rm=TRUE)) +
    ## ylim(quantile(all.indiv.mets$Micro_originality, c(0, .95),
    ##                          na.rm=TRUE)) +
    ##            xlim(range(all.indiv.mets$RBCL_degree, na.rm=TRUE))  +
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     text = element_text(size=16),
                     legend.position = "none") +
               geom_point(data = all.indiv.mets,
                          mapping = aes(x = RBCL_degree, y = Micro_originality,
                                        color =  GenusSpecies)) +
    scale_color_viridis(discrete=TRUE)

ggsave("figures/mirco_RBCL_degree.pdf",
       height=4, width=5)

##
p0.micro  <- mod.micro.orig %>%
    spread_draws(b_Intercept,
                 b_scaleoriginality) %>%
    mutate(originality =
               list(seq(min(all.indiv.mets$originality, na.rm=TRUE),
                        max(all.indiv.mets$originality, na.rm=TRUE),
                        0.01))) %>%
    unnest(originality) %>%
    mutate(pred = (b_Intercept +
               b_scaleoriginality*originality))   %>%
               group_by(originality) %>%
               summarise(pred_m = mean(pred),
                         pred_low_95 = quantile(pred, prob = 0.025),
                         pred_high_95 = quantile(pred, prob = 0.975),
                         pred_low_90 = quantile(pred, prob = 0.05),
                         pred_high_90 = quantile(pred, prob = 0.95),
                         pred_low_85 = quantile(pred, prob = 0.075),
                         pred_high_85 = quantile(pred, prob = 0.925)) %>%
               ggplot(aes(x = originality, y = pred_m)) +
               geom_line() +
               geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                           fill="goldenrod") +
               geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                           fill="dodgerblue") +
               geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                           fill="darkolivegreen") +
               ylab("Microbiome distinctness") +
    xlab("Host trait distinctness") +
    ## ylim(quantile(sp.means$Micro_originality, c(0, .95),
    ##               na.rm=TRUE)) +
    ## ylim(quantile(all.indiv.mets$Micro_originality, c(0, .95),
    ##                          na.rm=TRUE)) +
    ##            xlim(range(all.indiv.mets$originality, na.rm=TRUE))  +
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     text = element_text(size=16),
                     legend.position = "none") +
               geom_point(data = all.indiv.mets,
                          mapping = aes(x = originality, y = Micro_originality,
                                        color =  GenusSpecies)) +
    scale_color_viridis(discrete=TRUE)

ggsave("figures/mirco_originality.pdf",
       height=4, width=5)




## micro diversity

p2.micro  <- mod.micro.div %>%
    spread_draws(b_Intercept,
                 b_scaleRBCL_degree) %>%
    mutate(RBCL_degree =
               list(seq(min(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        max(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        0.01))) %>%
    unnest(RBCL_degree) %>%
    mutate(pred = (b_Intercept +
               b_scaleRBCL_degree*RBCL_degree))   %>%
               group_by(RBCL_degree) %>%
               summarise(pred_m = mean(pred),
                         pred_low_95 = quantile(pred, prob = 0.025),
                         pred_high_95 = quantile(pred, prob = 0.975),
                         pred_low_90 = quantile(pred, prob = 0.05),
                         pred_high_90 = quantile(pred, prob = 0.95),
                         pred_low_85 = quantile(pred, prob = 0.075),
                         pred_high_85 = quantile(pred, prob = 0.925)) %>%
               ggplot(aes(x = RBCL_degree, y = pred_m)) +
               geom_line() +
               geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                           fill="goldenrod") +
               geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                           fill="dodgerblue") +
               geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                           fill="darkolivegreen") +
               ylab("Microbiome diversity") +
    xlab("Individual-level diet breadth") +
    ylim(0,6) +
    xlim(0,20) +
    ## ylim(quantile(sp.means$Micro_originality, c(0, .95),
    ##               na.rm=TRUE)) +
    ## ylim(quantile(all.indiv.mets$Micro_originality, c(0, .95),
    ##                          na.rm=TRUE)) +
    ##            xlim(range(all.indiv.mets$RBCL_degree, na.rm=TRUE))  +
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.title.x = element_text(size=16),
                     axis.title.y = element_text(size=16),
                     text = element_text(size=16),
                     legend.position = "none") +
               geom_point(data = all.indiv.mets,
                          mapping = aes(x = RBCL_degree, y = Micro_partner.diversity,
                                        color =  GenusSpecies)) +
    scale_color_viridis(discrete=TRUE)

ggsave("figures/mirco_div_RBCL_degree.pdf",
       height=4, width=5)


all.micro <- grid.arrange(p1.micro, p2.micro, nrow = 2)

ggsave("figures/all_microBayes.pdf", all.micro, height=8, width=6)


