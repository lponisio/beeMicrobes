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
library(ALDEx3)
#library(microbiome) #bioconductor
#library(ANCOMBC)
#remotes::install_version("rbiom", version = "1.0.3") # need older version of rbiom
#BiocManager::install("treeclimbR")
#library(treeclimbR)
#BiocManager::install("mia")
#library(mia) # bioconductor

# Load in bee/microbe samples
load('data/allNetSums_final.RData')
#source("~/Desktop/node_taxonomy_labels.R")

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

# make each pathogen column into a 0/1 factor
for (v in path_vars) {
  sample_data(bee_ps)[[v]] <- factor(sample_data(bee_ps)[[v]],
                                     levels = c(0, 1),
                                     labels = c("Uninfected", "Infected"))
}

##-------------------- ##
## ALDEX3 Prep and Run ##
##-------------------- ##

###--- Function to lapply ALDEx3 Runs ---###

run_aldex3 <- function(var, Y, md, nsample = 250, gamma = 0.5) {
  message(format(Sys.time(), "%H:%M:%S"), "  Running ", var,
          " (", match(var, path_vars), " of ", length(path_vars), ")")
  keep <- !is.na(md[[var]])                        # drop bees missing this pathogen
  f <- as.formula(paste("~", var, "+ (1 | Genus) + (1 | Site)"))
  set.seed(42)
  aldex(Y[, keep], f,
        data    = md[keep, ],
        method  = "lme4",
        nsample = nsample,
        scale   = clr.sm,
        gamma   = gamma)
}

###--- All bee species ---###

# Keep only samples with infection status, species and site recorded
ald.md <- data.frame(sample_data(bee_ps))
ald.ps <- prune_samples(complete.cases(ald.md[, c("NosemaCeranae", "Genus", "Site")]),
                    bee_ps)
ald.md <- data.frame(sample_data(ald.ps))

# Raw counts, taxa as rows
Y <- as(otu_table(ald.ps), "matrix")
if (!taxa_are_rows(ald.ps)) Y <- t(Y)

# Pool rare taxa (present in < 5% of samples) into one "other" row,
# following ALDEx3 examples
rare <- rowMeans(Y > 0) < 0.05
Y <- rbind(Y[!rare, ], other = colSums(Y[rare, , drop = FALSE]))

all.bee.fits <- setNames(lapply(path_vars, run_aldex3, Y = Y, md = ald.md), path_vars)

# Run ALDEx3: infection effect with random intercepts for Genus (fewer factor levels) and site
set.seed(42)
fit <- aldex(Y,
             ~ NosemaCeranae + (1 | Genus) + (1 | Site),
             data    = ald.md,
             method  = "lme4",   # exact mixed model; "blmm" is a faster approximation
             nsample = 250,      # Monte Carlo draws; raise (e.g. 1000+) for final results
             scale   = clr.sm,   # CLR scale model (like ALDEx2)
             gamma   = 0.5)      # scale uncertainty; gamma = 0 means none

# Results
res <- summary(fit)
head(res)                                         # check the column names
res[res$parameter == "NosemaCeranaeInfected", ]    # infection effect per taxon

### Plot ALDEx results
par(mfrow=c(2,2))
aldex.plot(fit, contrast="NosemaCeranaeInfected", plot='vol', main="Volcano plot")
aldex.plot(fit, contrast="NosemaCeranaeInfected", plot='eff', main="Effect plot")
aldex.plot(fit, contrast="NosemaCeranaeInfected", plot='MA', main="MA plot")
aldex.plot(fit, contrast="NosemaCeranaeInfected", plot="water", 
           water.col=c("orange", "purple"), main="Waterfall plot")


###--- Bombus Only ---###
# Keep only samples with infection status, species and site recorded
ald.md <- data.frame(sample_data(bee_ps))
ald.ps <- prune_samples(complete.cases(ald.md[, c("NosemaCeranae", "Genus", "Site")]),
                        bee_ps)
bombus_ps <- subset_samples(bee_ps, Genus == "Bombus")
?prune_samples

ald.md <- data.frame(sample_data(ald.ps))


###--- Mellisodes Only ---###


###--- Apis Only ---###



# ##--------------- ##
# ## Run ANCOM-BC2  ##
# ##--------------- ##
# 
# # For any ApicystisSpp
# daa.apsp <- ancombc2(data = bee_ps,
#                      tax_level = "Genus",
#                      fix_formula = "ApicystisSpp",
#                      rand_formula = "(1|Site) + (1|Genus)",
#                      p_adj_method = "fdr", pseudo_sens = TRUE,
#                      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
#                      group = "ApicystisSpp", struc_zero = TRUE, neg_lb = TRUE,
#                      alpha = 0.05, n_cl = 2, verbose = TRUE,
#                      global = FALSE, pairwise = FALSE,
#                      dunnet = FALSE, trend = FALSE,
#                      iter_control = list(tol = 1e-2, max_iter = 50, verbose = TRUE),
#                      em_control = list(tol = 1e-5, max_iter = 100),
#                      lme_control = lme4::lmerControl())
# 
# ##----------------------- ##
# ## Plot ANCOM-BC2 results ##
# ##----------------------- ##
# # prepare the data
# res_prim <- daa.apsp$res
# 
# term <- "ApicystisSppInfected"
# 
# df_path <- res_prim %>%
#   dplyr::select(taxon,
#                 lfc    = paste0("lfc_",term),
#                 se     = paste0("se_",term),
#                 diff   = paste0("diff_",term),
#                 robust = paste0("diff_robust_", term))
# 
# df_fig_path <- df_path %>%
#   filter(diff == 1) %>%
#   arrange(desc(lfc)) %>%
#   mutate(direction = ifelse(lfc > 0, "Higher in infected", "Lower in infected"),
#          color  = ifelse(robust, "aquamarine3", "black"))
# 
# # make taxon into a factor in the order of positive to negative LFC
# df_fig_path$taxon  <- factor(df_fig_path$taxon, levels = df_fig_path$taxon)
# df_fig_path$direction <- factor(df_fig_path$direction,
#                              levels = c("Higher in infected", "Lower in infected"))
# 
# fig_path <- df_fig_path %>%
#   ggplot(aes(x = taxon, y = lfc, fill = direction)) +
#   geom_bar(stat = "identity", width = 0.7, color = "black",
#            position = position_dodge(width = 0.4)) +
#   geom_errorbar(aes(ymin = lfc - se, ymax = lfc + se),
#                 width = 0.2, position = position_dodge(0.05), color = "black") +
#   labs(x = NULL, y = "Log fold change",
#        title = "Log fold changes: pathogen-infected vs. uninfected bees",
#        caption = "Error bars: ± 1 SE. Taxon labels in aquamarine passed the pseudo-count sensitivity analysis.") +
#   scale_fill_discrete(name = NULL) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.minor.y = element_blank(),
#         axis.text.x = element_text(angle = 60, hjust = 1,
#                                    color = df_fig_path$color))
# fig_path

# ##-------------------- ##
# ## TreeclimbR Analysis ##
# ##-------------------- ##
# 
# path_vars  <- c("AnyPathogen", "ApicystisSpp", "AscosphaeraSpp", "CrithidiaBombi", "CrithidiaExpoeki",
#                 "CrithidiaSpp", "NosemaCeranae", "NosemaBombi")
# 
# # format data for treeclimbR (it uses TreeSummarizedExperiment objects -> tse)
# tse <- convertFromPhyloseq(bee_ps)
# 
# # 1. Add aggregated counts for every internal node
# all_nodes <- showNode(tree = rowTree(tse), only.leaf = FALSE) # put nodes into object
# tse_agg <- aggTSE(x = tse, rowLevel = all_nodes, rowFun = sum) # sum reads at each node
# 
# # 2. Run DA on all leaves and nodes
# # re-run from "design <-..." with each pathogen
# md <- colData(tse_agg) # get metadata from tse object
# design <- model.matrix(~ Genus + Site + ApicystisSpp, data = md) # define experiment, and which parasite
# da_res <- runDA(tse_agg, assay = "counts", option = "glmQL",
#                 design = design, contrast = NULL,  # NULL tests the last coefficient
#                 filter_min_count = 10)
# da_tbl <- nodeResult(da_res.apsp, n = Inf, type = "DA")
# 
# # 3. Propose candidate resolutions, 4. pick the best one
# cand <- getCand(tree = rowTree(tse_agg), score_data = da_tbl,
#                 node_column = "node", p_column = "PValue",
#                 sign_column = "logFC", message = FALSE)
# best <- evalCand(tree = rowTree(tse_agg), levels = cand$candidate_list,
#                  score_data = da_tbl, node_column = "node",
#                  p_column = "PValue", sign_column = "logFC")
# topNodes(best)   # final significant nodes/leaves with direction
# 
# #write_rds(da_res, file = "~/Desktop/beeMicrobees_treeclimbR_NosemaBombi_11-19-26.rds")
# 
# da_res.apsp <- readRDS("~/Desktop/beeMicrobees_treeclimbR_ApicystisSpp_11-18-26.rds")

##----------------- ##
## TreeclimbR Plots ##
##----------------- ##

### Prep data for plotting

## Label nodes with taxonomy

# extract the tree, taxonomy table (rows are asvs, columns are tax levels), and top flagged nodes
tree <- rowTree(tse) # same tree as rowTree(tse_agg)
tax  <- as.data.frame(rowData(tse)) # ASV-level taxonomy
sel_nodes <- topNodes(best, n = nrow(best$output), p_value = 0.05)$node
node_labels <- label_nodes_by_taxonomy(nodes = sel_nodes, tree = tree, tax = tax)

## prep the dataframe for plotting log fold change values
df_treeclimb_lfc <- da_tbl %>%
  dplyr::filter(node %in% sel_nodes) %>%
  dplyr::left_join(node_labels, by = "node") %>%
  # edgeR's QL F-test has 1 df here (infected vs not), so F = t^2 and SE = |logFC| / sqrt(F)
  dplyr::mutate(se = abs(logFC) / sqrt(F), # calculate approximate standard error
                direction = ifelse(logFC > 0, "Higher in infected",
                                "Lower in infected")) %>%
  dplyr::arrange(desc(logFC)) # arrange by lfc values, high to low

df_treeclimb_lfc$label  <- factor(df_treeclimb_lfc$label, levels = df_treeclimb_lfc$label)
df_treeclimb_lfc$direction <- factor(df_treeclimb_lfc$direction,
                                     levels = c("Higher in infected", "Lower in infected"))

# Results table for export
res_table_treeclimb <- df_treeclimb_lfc %>%
  dplyr::select(node, label, shared_rank, shared_name, n_asv,
                logFC, se, PValue, FDR)

### Log fold change plots ###

wf.plot.apsp <- ggplot(df_treeclimb_lfc,
                       aes(x = label, y = logFC, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = logFC - se, ymax = logFC + se),
                width = 0.2, color = "black") +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(x = NULL,
       y = "log2 fold change\n(infected vs. uninfected with Apicystis)") +
  scale_fill_discrete(name = NULL, drop = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
wf.plot.apsp

### Phylogeny plots ###

