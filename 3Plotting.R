rm(list=ls())
source("src/init_bayes.R")
source("src/misc.R")
load('data/allNetSums.RData')

library(dplyr)
library(gridExtra)
library(viridis)
library(tidyr)
library(lemon)

load(file="saved/parasiteOrigFit.Rdata")
load(file="saved/parasiteDivFit.Rdata")
load(file="saved/microbeOrigFit.Rdata")
load(file="saved/MicrobeDivFit.Rdata")

bee.fams <- c("Apidae", "Halictidae", "Megachilidae")
all.indiv.mets <- all.indiv.mets[all.indiv.mets$Family %in% bee.fams,]

sp.cols <- viridis(length(unique(all.indiv.mets$GenusSpecies)))

sp.means <- all.indiv.mets  %>%
    group_by(GenusSpecies, originality, Site, r.degree) %>%
    summarize(Parasite_originality = median(Parasite_originality,
                                            na.rm=TRUE),
              Micro_originality =  median(Micro_originality,
                                          na.rm=TRUE),
              Micro_partner.diversity =  median(Micro_partner.diversity,
                                                na.rm=TRUE),
              ParasiteRichness =  median(ParasiteRichness,
                                         na.rm=TRUE),
              )


###########################################################################
## parasite distinctness
###########################################################################

## parasite distinctness, trait orig


p0.parasite  <- mod.par.orig %>%
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
                             color =  GenusSpecies),
               alpha=0.7) +
    scale_color_viridis(discrete=TRUE)


grid.p0.parasites <- shared_legend(p0.parasite, ncol=1, position='bottom')

ggsave(grid.p0.parasites, file= "figures/parasite_orig_originality.pdf",
       height=6, width=7)

print("here 0")

## parasite distinctness, indiv degree
p1.parasite  <- mod.par.orig %>%
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
                             color =  GenusSpecies),
               alpha=0.7) +
    scale_color_viridis(discrete=TRUE)


grid.p1.parasites <- shared_legend(p1.parasite, ncol=1, position='bottom')

ggsave(grid.p1.parasites, file= "figures/parasite_orig_RBCL_degree.pdf",
       height=6, width=7)

print("here 1")

###########################################################################
## parasite richness
###########################################################################

## parasite richness, indiv degree
p2.parasite  <- mod.par.div %>%
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
    ylab("Pathobiome richness") +
    xlab("Individual diet breadth") +
    xlim(range(all.indiv.mets$RBCL_degree, na.rm=TRUE))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none") +
    geom_point(data = all.indiv.mets[!is.na(all.indiv.mets$RBCL_degree),],
               mapping = aes(x = RBCL_degree, y = ParasiteRichness,
                             color =  GenusSpecies),
               alpha=0.7) +
    scale_color_viridis(discrete=TRUE)

grid.p2.parasites <- shared_legend(p2.parasite, ncol=1, position='bottom')

ggsave(grid.p2.parasites, file= "figures/parasite_rich_RBCL_degree.pdf",
       height=6, width=7)

print("here 2")

## parasite richness, species degree
p3.parasite  <- mod.par.div %>%
    spread_draws(b_Intercept,
                 b_scaler.degree) %>%
    mutate(r.degree =
               list(seq(min(all.indiv.mets$r.degree, na.rm=TRUE),
                        max(all.indiv.mets$r.degree, na.rm=TRUE),
                        0.1))) %>%
    unnest(r.degree) %>%
    mutate(pred = (b_Intercept +
                   b_scaler.degree*r.degree))   %>%
    group_by(r.degree) %>%
    summarise(pred_m = mean(pred, na.rm = TRUE),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = r.degree, y = pred_m)) +
    geom_line() +
    geom_ribbon(aes(ymin = pred_low_95, ymax = pred_high_95), alpha=0.2,
                fill="goldenrod") +
    geom_ribbon(aes(ymin = pred_low_90, ymax = pred_high_90), alpha=0.2,
                fill="dodgerblue") +
    geom_ribbon(aes(ymin = pred_low_85, ymax = pred_high_85), alpha=0.2,
                fill="darkolivegreen") +
    ylab("Pathobiome richness") +
    xlab("Species diet breadth") +
    xlim(range(all.indiv.mets$r.degree))  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none") +
    geom_point(data = all.indiv.mets,
               mapping = aes(x = r.degree, y = ParasiteRichness,
                             color =  GenusSpecies),
                             alpha=0.7) +
    scale_color_viridis(discrete=TRUE)


grid.p3.parasites <- shared_legend(p3.parasite, ncol=1, position='bottom')

ggsave(grid.p3.parasites, file= "figures/parasite_rich_sp_degree.pdf",
       height=6, width=7)

print("here 3")

###########################################################################
## microbiome distinctness
###########################################################################

## 16s distinctness, trait orig
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
                             color =  GenusSpecies),
                             alpha=0.7) +
    scale_color_viridis(discrete=TRUE)

grid.p0.micro <- shared_legend(p0.micro, ncol=1, position='bottom')

ggsave(grid.p0.micro, file= "figures/mirco_orig_originality.pdf",
       height=6, width=7)

print("here 0")

## 16s distinctness, indiv diet
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
                             color =  GenusSpecies),
                             alpha=0.7) +
    scale_color_viridis(discrete=TRUE)

grid.p1.micro <- shared_legend(p1.micro, ncol=1, position='bottom')

ggsave(grid.p1.micro, file= "figures/mirco_orig_RBCL_degree.pdf",
       height=6, width=7)

print("here 1")

###########################################################################
## microbiome richness
###########################################################################

## micro richness, indiv diet

p2.micro  <- mod.micro.div %>%
    spread_draws(b_Intercept,
                 b_scaleRBCL_degree) %>%
    mutate(RBCL_degree =
               list(seq(min(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        max(all.indiv.mets$RBCL_degree, na.rm=TRUE),
                        0.1))) %>%
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
                             color =  GenusSpecies),
               alpha=0.7) +
    scale_color_viridis(discrete=TRUE)

grid.p2.micro <- shared_legend(p2.micro, ncol=1, position='bottom')

ggsave(grid.p2.micro, file= "figures/mirco_div_RBCL_degree.pdf",
       height=6, width=7)

print("here 2")


all.micro <- shared_legend(p1.micro, p2.micro, ncol=1,
                           nrow = 2, position='bottom')

ggsave("figures/all_microBayes.pdf", all.micro, height=8, width=7)

print("here all")


## micro richness, species diet breadth
p3.micro  <- mod.micro.div %>%
    spread_draws(b_Intercept,
                 b_scaler.degree) %>%
    mutate(r.degree =
               list(seq(min(all.indiv.mets$r.degree, na.rm=TRUE),
                        max(all.indiv.mets$r.degree, na.rm=TRUE),
                        0.1))) %>%
    unnest(r.degree) %>%
    mutate(pred = (b_Intercept +
                   b_scaler.degree*r.degree))   %>%
    group_by(r.degree) %>%
    summarise(pred_m = mean(pred),
              pred_low_95 = quantile(pred, prob = 0.025),
              pred_high_95 = quantile(pred, prob = 0.975),
              pred_low_90 = quantile(pred, prob = 0.05),
              pred_high_90 = quantile(pred, prob = 0.95),
              pred_low_85 = quantile(pred, prob = 0.075),
              pred_high_85 = quantile(pred, prob = 0.925)) %>%
    ggplot(aes(x = r.degree, y = pred_m)) +
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
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none") +
    geom_point(data = all.indiv.mets,
               mapping = aes(x = r.degree, y = Micro_partner.diversity,
                             color =  GenusSpecies), alpha=0.7) +
    scale_color_viridis(discrete=TRUE)

grid.p3.micro <- shared_legend(p3.micro, ncol=1, position='bottom')

ggsave(grid.p3.micro, file= "figures/mirco_div_sp_degree.pdf",
       height=6, width=7)

print("here 3")

