# This script runs differential abundance analyses to detect microbes especially
#   differentially abundant in infected or uninfected bees

##-----------------------------------------##
## Load packages and set working directory ##
##-----------------------------------------##
rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

dir.microbes <- file.path(local.path, "beeMicrobes")
setwd(dir.microbes)

library(readxl)
library(tidyverse)
library(phyloseq) #bioconductor
library(qiime2R) #github
library(microbiome) #bioconductor
library(ANCOMBC) # bioconductor
#remotes::install_version("rbiom", version = "1.0.3") # need older version of rbiom
#BiocManager::install("treeclimbR")
library(treeclimbR)
#BiocManager::install("mia")
library(mia) # bioconductor

# Load in bee/microbe samples
load('data/allNetSums_final.RData')

##---------- ##
## Data prep ##
##---------- ##

# one object for all pathogen names
path_vars  <- c("ApicystisSpp", "AscosphaeraSpp", "CrithidiaBombi", "CrithidiaExpoeki",
                "CrithidiaSpp", "NosemaCeranae", "NosemaBombi")

#... and all trait names
trait_vars <- c("NestLocation", "PrimaryNestMaterial", "NestConstruction",
                "NestPartitions", "Sociality", "MeanITD", "Family", "Apis")

# select relevant columns from specimen data
meta_dat <- all.indiv.mets %>%
  select(UniqueID, Site, SiteYear, GenusSpecies, Genus,
         all_of(c(trait_vars, path_vars))
         ) %>%
  mutate(AnyPathogen = as.integer(if_any(all_of(path_vars), ~ .x == 1))
         ); rownames(meta_dat) <- meta_dat$UniqueID

# Load up qiime output into phyloseq
qiime_dir <- "../beeMicrobes_saved/beeMicrobes_pipeline_output"

# build phyloseq object
bee_ps <- qza_to_phyloseq(
  features = file.path(qiime_dir, "merged/16s/tablef4.qza"),
  tree     = file.path(qiime_dir, "merged/16s/rooted-tree16s.qza"),
  taxonomy = file.path(qiime_dir, "merged/16s/taxonomy16s.qza"),
  metadata = file.path(qiime_dir, "maps/16s/beeMicrobes_combined_map_no-ctrls.txt"))

# filter phyloseq
# Filter to only our final sample set
bee_ps <- prune_samples(sample_names(bee_ps) %in% meta_dat$UniqueID, bee_ps)
#replace phyloseq boilerplate metadata with the trimmed down metadata
sample_data(bee_ps) <- sample_data(meta_dat)
# prune out any empty samples
bee_ps <- prune_samples(sample_sums(bee_ps) > 0, bee_ps)

# redefine path_vars with AnyPathogen
path_vars  <- c("AnyPathogen", "ApicystisSpp", "AscosphaeraSpp", "CrithidiaBombi", "CrithidiaExpoeki",
                "CrithidiaSpp", "NosemaCeranae", "NosemaBombi")

for (v in path_vars) {
  sample_data(bee_ps)[[v]] <- factor(sample_data(bee_ps)[[v]],
                                     levels = c(0, 1),
                                     labels = c("Uninfected", "Infected"))
}

##------------------- ##
## ANCOM-BC2 Analysis ##
##------------------- ##
# just use site random effect, confounded with year in ffar

# For any pathogen
daa.anypath <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "AnyPathogen",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "AnyPathogen", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

# For any ApicystisSpp
daa.apsp <- ancombc2(data = bee_ps,
                   tax_level = "Family",
                   fix_formula = "ApicystisSpp",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "ApicystisSpp", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

# For AscosphaeraSpp
daa.assp <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "AscosphaeraSpp",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "AscosphaeraSpp", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

# For CrithidiaBombi
daa.crbo <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "CrithidiaBombi",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "CrithidiaBombi", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

# For CrithidiaExpoeki
daa.crex <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "CrithidiaExpoeki",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "CrithidiaExpoeki", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

# For CrithidiaSpp
daa.crsp <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "CrithidiaSpp",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "CrithidiaSpp", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())


# For NosemaCeranae
daa.noce <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "NosemaCeranae",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "NosemaCeranae", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

# For NosemaBombi
daa.nobo <- ancombc2(data = bee_ps,
                   tax_level = "Genus",
                   fix_formula = "NosemaBombi",
                   rand_formula = "(1|Site) + (1|GenusSpecies)",
                   p_adj_method = "fdr", pseudo_sens = TRUE,
                   prv_cut = 0.05, lib_cut = 1000, s0_perc = 0.05,
                   group = "NosemaBombi", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = FALSE, pairwise = FALSE,
                   dunnet = FALSE, trend = FALSE,
                   iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl())

##----------------------- ##
## Plot ANCOM-BC2 results ##
##----------------------- ##
# prepare the data
res_prim <- daa.anypath$res

term <- "AnyPathogenInfected"

df_path <- res_prim %>%
  dplyr::select(taxon,
                lfc    = paste0("lfc_",term),
                se     = paste0("se_",term),
                diff   = paste0("diff_",term),
                robust = paste0("diff_robust_", term))

df_fig_path <- df_path %>%
  filter(diff == 1) %>%
  arrange(desc(lfc)) %>%
  mutate(direction = ifelse(lfc > 0, "Higher in infected", "Lower in infected"),
         color  = ifelse(robust, "aquamarine3", "black"))

# make taxon into a factor in the order of positive to negative LFC
df_fig_path$taxon  <- factor(df_fig_path$taxon, levels = df_fig_path$taxon)
df_fig_path$direction <- factor(df_fig_path$direction,
                             levels = c("Higher in infected", "Lower in infected"))

fig_path <- df_fig_path %>%
  ggplot(aes(x = taxon, y = lfc, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7, color = "black",
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se),
                width = 0.2, position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log fold change",
       title = "Log fold changes: pathogen-infected vs. uninfected bees",
       caption = "Error bars: ± 1 SE. Taxon labels in aquamarine passed the pseudo-count sensitivity analysis.") +
  scale_fill_discrete(name = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_path$color))
fig_path



##-------------------- ##
## TreeclimbR Analysis ##
##-------------------- ##

path_vars  <- c("AnyPathogen", "ApicystisSpp", "AscosphaeraSpp", "CrithidiaBombi", "CrithidiaExpoeki",
                "CrithidiaSpp", "NosemaCeranae", "NosemaBombi")

# format data for treeclimbR (it uses tree summarized experiment objects)
tse <- convertFromPhyloseq(bee_ps)

# 1. Add aggregated counts for every internal node
all_nodes <- showNode(tree = rowTree(tse), only.leaf = FALSE) # put nodes into object
tse_agg <- aggTSE(x = tse, rowLevel = all_nodes, rowFun = sum) # sum reads at each node

# 2. Run DA on all leaves and nodes
# re-run from "design <-..." with each pathogen
md <- colData(tse_agg) # get metadata from tse object
design <- model.matrix(~ Genus + Site + ApicystisSpp, data = md) # define experiment, and which parasite
da_res <- runDA(tse_agg, assay = "counts", option = "glmQL",
                design = design, contrast = NULL,  # NULL tests the last coefficient
                filter_min_count = 10)
da_tbl <- nodeResult(da_res, n = Inf, type = "DA")

# 3. Propose candidate resolutions, 4. pick the best one
cand <- getCand(tree = rowTree(tse_agg), score_data = da_tbl,
                node_column = "node", p_column = "PValue",
                sign_column = "logFC", message = FALSE)
best <- evalCand(tree = rowTree(tse_agg), levels = cand$candidate_list,
                 score_data = da_tbl, node_column = "node",
                 p_column = "PValue", sign_column = "logFC")
topNodes(best)   # final significant nodes/leaves with direction

write_rds(da_res, file = "~/Desktop/beeMicrobees_treeclimbR_NosemaBombi_11-19-26.rds")
da_res.apsp <- readRDS("~/Desktop/beeMicrobees_treeclimbR_ApicystisSpp_11-18-26.rds")
