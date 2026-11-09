# This script runs random forest analyses to determine importance of microbe-parasite associations
#   and compares importance of traits and taxonomy relative to microbes in predicting infection status.

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
library(ranger)
library(randomForestSRC)
library(pROC)
library(cluster)
library(phyloseq)
library(devtools)
library(qiime2R)
#install_github("hk1785/PIRF", force=T)
library(PIRF)

set.seed(1)

# Load in bee/microbe samples
load('data/allNetSums_final.RData')
# load in a few helper functions
source("./src/rf_helpers.R")
source("./src/pirf_plots.R")

## -------------------------------------- ##
## Prep data for Random Forests Analysis  ##
## -------------------------------------- ##

# Get vectors of trait names and pathogen/parasite names
trait_vars <- c("NestLocation", "PrimaryNestMaterial", "NestConstruction",
                "NestPartitions", "Sociality", "MeanITD", "Family", "Apis")
path_vars  <- c("ApicystisSpp", "CrithidiaBombi", "CrithidiaExpoeki",
                "CrithidiaSpp", "NosemaCeranae", "NosemaBombi")

# ensure unique IDs are character and nest location has no spaces in character values
all.indiv.mets <- all.indiv.mets %>%
  mutate(UniqueID = as.character(UniqueID),
         NestLocation = sub("-", " ", NestLocation))

# trim the dataset down to just the trait, pathogen and metadata columns
# lump(): merges rare trait categories into "other"
meta <- all.indiv.mets %>%
  select(UniqueID, Site, GenusSpecies, all_of(c(trait_vars, path_vars))) %>%
  mutate(AnyPathogen = as.integer(if_any(all_of(path_vars), ~ .x == 1)),
         across(all_of(setdiff(trait_vars, "MeanITD")), lump),
         fold = as.integer(factor(Site)) %% 5) # assign bees to one of 5 "folds" for cross validation

# Pivot microbes to longform and filter to microbes that occur in at least 1% of bees
micro_long <- all.indiv.mets %>%
  select(UniqueID, starts_with("16s:")) %>%
  pivot_longer(-UniqueID, names_to = "taxon", values_to = "rel") %>%
  group_by(taxon) %>%
  filter(mean(rel > 0) >= 0.01) %>% # prevalence filter - try with and without prevalence filter
  ungroup()

# Make a list of shorthand taxon names for interpretability
taxa <- micro_long %>%
  distinct(taxon) %>%
  mutate(name = make.unique(make.names(sub(".*[a-z]__", "", taxon))))

# tiny value (half the smallest relative abundance) is added to each value
# needed for calculating the centered log ratio of each microbe
pseudocount <- min(micro_long$rel[micro_long$rel > 0]) / 2

# add in short-hand taxon names and CLR transform data
micro_long <- micro_long %>%
  left_join(taxa, by = "taxon") %>%
  group_by(UniqueID) %>%
  mutate(clr = log(rel + pseudocount) - mean(log(rel + pseudocount))) %>% # the CLR transformation
  ungroup()

# pivot the CLR-transformed data wider with the short-hand taxon names
micro <- micro_long %>%
  select(UniqueID, name, clr) %>%
  pivot_wider(names_from = name, values_from = clr) %>%
  select(-UniqueID)

# make separate trait and pathogen/parasite data frames
traits    <- meta %>% select(all_of(trait_vars)) # trait values
pathogens <- meta %>% select(all_of(path_vars), AnyPathogen) # 0/1 bee x pathogen dataframe
bee.sp    <- meta %>% transmute(GenusSpecies = factor(GenusSpecies)) # trait values

# ==== 1. pathogens: microbes, traits, microbes + traits models ====
# make blocks of predictors for each of the 3 models
blocks <- list(micro = micro,
               #traits = traits,
               bee.sp = bee.sp,
               #micro_traits = bind_cols(micro, traits),
               micro_sp = bind_cols(micro, bee.sp)#,
               #sp.traits = bind_cols(traits, bee.sp)
               )

# cv_auc(): site-fold CV AUC with balanced class weights (uses cv_pred() + bal_wts())
#   ^ allows us to hold out the site of each bee being classified to avoid site-level correlations
#   ^ weights AUC by number of bees in the class (infected/not) to avoid high accuracy for only
#       uninfected bees or whatever the rare class is.
# run random forests for each of the three models using the ranger function called by cv_pred
auc_tab <- expand_grid(pathogen = names(pathogens), block = names(blocks)) %>%
  rowwise() %>%
  mutate(auc = cv_auc(pathogens[[pathogen]], blocks[[block]], meta$fold)) %>%
  ungroup()

auc_tab %>%
  pivot_wider(names_from = block, values_from = auc)

## ------------------------- ##
##  Null model for microbes  ##
## ------------------------- ##

null_auc_sp <- replicate(99, {
  shuffled <- ave(pathogens$AnyPathogen, meta$GenusSpecies,
                  FUN = \(grp) grp[sample.int(length(grp))]
                  ) # shuffle microbes within each species
  cv_auc(shuffled, blocks$micro_sp, meta$fold) # 
})

real_auc_sp <- cv_auc(pathogens$AnyPathogen, blocks$micro_sp, meta$fold)


# p_perm is proportion of permutations where the null model AUC was greater than the real AUC
tibble(real_auc      = real_auc_sp,
       null_auc_mean = mean(null_auc_sp),
       null_auc_max  = max(null_auc_sp),
       p_perm        = (sum(null_auc_sp >= real_auc_sp) + 1) / (length(null_auc_sp) + 1)
       )

tibble(null_auc = null_auc_sp) %>%
  ggplot(aes(null_auc)) +
  geom_histogram(bins = 20, fill = "grey70") +
  geom_vline(xintercept = real_auc_sp, colour = "red") +
  labs(x = "AUC with labels shuffled within species", y = "Count")

# ==== 2. rfsrc importance with CIs (microbes + traits) ====
# Fit random forest models, generate confidence intervals around importance of top predictive features

# path_imp(): rfsrc fit + subsampled importance CIs

# 2A. make dataframe from the anypathogen data and the microbes+traits block
importance_data <- as.data.frame(bind_cols(y = factor(pathogens$AnyPathogen),
                                           blocks$micro_traits))

# 2B. Fit random forest model
any.path.rf <- rfsrc(y ~ ., data = importance_data,
                     ntree = 500,
                     rfq = TRUE,
                     perf.type = "gmean",
                     importance = "permute")

# 2C. refit on 50 random subsets of bees (runs pretty slowly)
# uses the forest from 2B to test if variable importance is constant across bees
any_sub <- subsample(any.path.rf, B = 50, verbose = FALSE)

# 2D. get importance from subsamples
importance_ci <- extract.subsample(any_sub)$var.jk.sel.Z %>%
  as_tibble(rownames = "variable") %>% # convert rownames to a column for inspection
  arrange(desc(mean)) # arrange by mean

# 2E. look at the results
importance_ci %>% filter(signif) # variables whose CI excludes 0
plot.subsample(any_sub)

top_vars <- importance_ci %>% filter(signif) %>% slice_head(n = 4) %>% pull(variable)
plot.variable(any.path.rf, xvar.names = "Family",
              partial = TRUE,
              target = "1")

# 2F. Use path_imp (runs 2A-2D in a function) to run other pathogens with other variable sets

## Any Pathogen
# Microbes only
imp_any.mic <- path_imp(pathogens$AnyPathogen, blocks$micro)
imp_any.mic$imp_ci %>% filter(signif) # look at the top importance microbes
plot.subsample(imp_any.mic$sub_fit)

top_vars <- imp_any.mic$imp_ci %>% filter(signif) %>% slice_head(n = 4) %>% pull(variable)
plot.variable(imp_any.mic$fit, xvar.names = top_vars, partial = TRUE, target = "1")   # direction

# Traits only
imp_any.trt <- path_imp(pathogens$AnyPathogen, blocks$traits)
imp_any.trt$imp_ci %>% filter(signif) # look at the top importance microbes
plot.subsample(imp_any.trt$sub_fit)

top_vars <- imp_any.trt$imp_ci %>% filter(signif) %>% slice_head(n = 4) %>% pull(variable)
plot.variable(imp_any.trt$fit, xvar.names = top_vars, partial = TRUE, target = "1")   # direction

# ==== 3. PIRF (microbes only, ASV level) ====
qiime_dir <- "../beeMicrobes_saved/beeMicrobes_pipeline_output"

# build phyloseq object
bee_ps <- qza_to_phyloseq(
  features = file.path(qiime_dir, "merged/16s/final/core_metrics16s/rarefied_table.qza"),
  tree     = file.path(qiime_dir, "merged/16s/rooted-tree16s.qza"),
  taxonomy = file.path(qiime_dir, "merged/16s/taxonomy16s.qza"),
  metadata = file.path(qiime_dir, "maps/16s/beeMicrobes_combined_map_no-ctrls.txt"))

mean(meta$UniqueID %in% sample_names(bee_ps))                  # should be ~1
# if not, relabel from the map column that holds UniqueID:
# sample_names(bee_ps) <- as.character(sample_data(bee_ps)$UniqueID)

# Filter to only our final sample set
bee_ps <- prune_samples(sample_names(bee_ps) %in% meta$UniqueID, bee_ps)
ntaxa(bee_ps)
bee_ps <- tax_glom(bee_ps, taxrank = "Genus", NArm = FALSE)
ntaxa(bee_ps)         
bee_ps <- transform_sample_counts(bee_ps, \(abund) abund / sum(abund))

# Filter low-prevalence taxa
bee_ps <- filter_taxa(bee_ps, \(abund) mean(abund > 0) >= 0.02, prune = TRUE) 

meta_ps  <- meta %>% filter(UniqueID %in% sample_names(bee_ps)) # filter metadata by ps object samples
asv_tree <- phy_tree(bee_ps) # extract microbe phylogeny

uninformative <- "uncultured|unidentified|unknown|metagenome|Incertae_Sedis|^$"
asv_tax <- tax_table(bee_ps)[asv_tree$tip.label, ] %>%
  as.data.frame() %>%
  as_tibble(rownames = "asv") %>%
  mutate(across(-asv, \(lab) if_else(grepl(uninformative, lab, ignore.case = TRUE), NA_character_, lab)),
         id = paste0("ASV", row_number()),                    # safe names for ranger formula
         finest_taxon = coalesce(Genus, Family, Order, Class, Phylum, Kingdom),
         finest_rank  = case_when(!is.na(Genus)  ~ "genus",
                                  !is.na(Family) ~ "family",
                                  !is.na(Order)  ~ "order",
                                  !is.na(Class)  ~ "class",
                                  !is.na(Phylum) ~ "phylum",
                                  .default       = "kingdom")) %>%
  select(id, asv, finest_taxon, finest_rank, genus = Genus)

# PIRF needs tip order, so rearrange appropriately
asv_tab <- as.data.frame(t(otu_table(bee_ps)))[meta_ps$UniqueID, asv_tax$asv]
names(asv_tab) <- asv_tree$tip.label <- asv_tax$id

# run PIRF to associate microbes with having any gut parasite
pirf_any <- PIRF(asv_tab, meta_ps$AnyPathogen, asv_tree,
                 num.trees = 1000, num.threads = 4)

# PIRF returns 3 fits (mtry = 10% / sqrt / log2 of ASVs); keep the one with lowest OOB error
best_any  <- best_mtry(pirf_any)

## -- Check model performance --##
sapply(pirf_any, `[[`, "oob.err") # Brier score per mtry setting (lower = better)
best_any$mtry

roc_any <- roc(meta_ps$AnyPathogen, best_any$fit$predictions, quiet = TRUE)   # OOB predictions
auc(roc_any)
plot(roc_any)

## -- Rank ASV importance --##
pirf_imp <- tibble(id = names(best_any$sel.prob), sel_prob = best_any$sel.prob) %>%
  mutate(cluster = best_any$clust[id]) %>%
  left_join(asv_tax, by = "id") %>%
  arrange(desc(sel_prob))
pirf_imp %>% print(n = 20)

setdiff(best_any$rev.features, "") # ASVs PIRF dropped (zero selection probability)

## -- Importance of phylogenetic clusters --##
pirf_imp %>%
  slice_head(n = 20) %>%
  ggplot(aes(sel_prob, reorder(paste(id, finest_taxon), sel_prob), fill = factor(cluster))) +
  geom_col() +
  labs(x = "Selection probability", y = NULL, fill = "Phylo cluster")

## -- Infected VS Uninfected directionality --##
top_ids <- pirf_imp %>% slice_head(n = 10) %>% pull(id)

asv_tab %>%
  select(all_of(top_ids)) %>%
  mutate(AnyPathogen = meta_ps$AnyPathogen) %>%
  pivot_longer(-AnyPathogen, names_to = "id", values_to = "rel") %>%
  group_by(id, AnyPathogen) %>%
  summarise(prev = mean(rel > 0), mean_rel = mean(rel), .groups = "drop") %>%
  pivot_wider(names_from = AnyPathogen, values_from = c(prev, mean_rel)) %>%
  left_join(asv_tax, by = "id") %>%
  arrange(match(id, top_ids))

# this applies PIRF to each pathogen and the "any pathogen" scenario
pirf <- setNames(lapply(names(pathogens), \(pathogen)
                        PIRF(asv_tab, meta_ps[[pathogen]],
                             asv_tree,
                             num.trees = 1000,
                             num.threads = 4)),
                 names(pathogens)) # 0/1 y -> regression (Brier)


tibble(pathogen = names(pirf)) %>%
  rowwise() %>%
  mutate(oob_brier = best_mtry(pirf[[pathogen]])$oob.err,
         auc = as.numeric(auc(meta_ps[[pathogen]],
                              best_mtry(pirf[[pathogen]])$fit$predictions, quiet = TRUE))
         ) %>%
  ungroup() # compare to auc_tab micro

pirf_any <- best_mtry(pirf$AnyPathogen)
pirf_apic <- best_mtry(pirf$ApicystisSpp)
pirf_nosb <- best_mtry(pirf$NosemaBombi)
pirf_nosc <- best_mtry(pirf$NosemaCeranae)
pirf_critsp <- best_mtry(pirf$CrithidiaSpp)
pirf_critex <- best_mtry(pirf$CrithidiaExpoeki)
pirf_critbo <- best_mtry(pirf$CrithidiaBombi)

tibble(id = names(pirf_apic$sel.prob), sel_prob = pirf_apic$sel.prob) %>%
  mutate(cluster = pirf_any$clust[id]) %>%
  left_join(asv_tax, by = "id") %>%
  arrange(desc(sel_prob))

pd_any <- pirf_pd(pirf_any, asv_tab, asv_tax, n_top = 6)
plot_pirf_pd(pd_any, asv_tab, meta_ps$AnyPathogen,
             y_lab = "Predicted probability of any pathogen") + theme_bw()

pd_apic <- pirf_pd(pirf_apic, asv_tab, asv_tax, n_top = 6)
plot_pirf_pd(pd_apic, asv_tab, meta_ps$ApicystisSpp,
             y_lab = "Predicted probability of Apicystis") + theme_bw()

pd_nosb <- pirf_pd(pirf_nosb, asv_tab, asv_tax, n_top = 6)
plot_pirf_pd(pd_nosb, asv_tab, meta_ps$NosemaBombi,
             y_lab = "Predicted probability of Nosema Bombi") + theme_bw()

pd_nosc <- pirf_pd(pirf_nosc, asv_tab, asv_tax, n_top = 6)
plot_pirf_pd(pd_nosc, asv_tab, meta_ps$NosemaCeranae,
             y_lab = "Predicted probability of Nosema Ceranae") + theme_bw()

## ----------------- ##
##  Learning Curves  ##
## ----------------- ##

train_fracs <- seq(0.1, 0.9, by = 0.1)    # 10% to 90% training data
n_rep       <- 5                          # random splits per fraction
curve_list  <- list()

for (pathogen in names(pathogens)) {
  resp <- factor(pathogens[[pathogen]])
  
  for (frac in train_fracs) {
    for (rep_i in 1:n_rep) {
      
      # stratified split: same infection rate in training and test bees
      train_rows <- tibble(row = seq_along(resp), resp) |>
        group_by(resp) |>
        slice_sample(prop = frac) |>
        pull(row)
      test_rows <- setdiff(seq_along(resp), train_rows)
      
      # train on the training bees, with balanced class weights
      fit <- ranger(x = micro[train_rows, ], y = resp[train_rows],
                    case.weights = bal_wts(resp[train_rows]),
                    probability = TRUE, num.trees = 200)
      
      # predict the test bees and score
      prob_1 <- predict(fit, micro[test_rows, ])$predictions[, "1"]
      pred   <- ifelse(prob_1 > 0.5, "1", "0")
      obs    <- as.character(resp[test_rows])
      
      sens <- mean(pred[obs == "1"] == "1")    # infected bees correctly flagged
      spec <- mean(pred[obs == "0"] == "0")    # uninfected bees correctly cleared
      
      curve_list[[length(curve_list) + 1]] <-
        tibble(pathogen, frac, rep = rep_i, bal_acc = (sens + spec) / 2)
    }
  }
}

curves <- bind_rows(curve_list)

# mean and 95% CI across the random splits
curve_summary <- curves |>
  group_by(pathogen, frac) |>
  summarise(mean_acc = 100 * mean(bal_acc),
            se       = 100 * sd(bal_acc) / sqrt(n()),
            .groups  = "drop")

ggplot(curve_summary, aes(100 * frac, mean_acc, colour = pathogen, fill = pathogen)) +
  geom_ribbon(aes(ymin = mean_acc - 1.96 * se, ymax = mean_acc + 1.96 * se),
              alpha = 0.3, colour = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +
  labs(x = "Training data %", y = "Balanced accuracy %",
       colour = "Pathogen", fill = "Pathogen") +
  theme_bw()

# ==== 6. differential abundance check (traits as covariates) ====
library(Maaslin2)
rel_wide <- micro_long %>%
  select(UniqueID, name, rel) %>%
  pivot_wider(names_from = name, values_from = rel) %>%
  column_to_rownames("UniqueID")
meta_df <- meta %>% column_to_rownames("UniqueID")

for (pathogen in names(pathogens))
  Maaslin2(rel_wide, meta_df, output = paste0("maaslin_", pathogen),
           fixed_effects = c(pathogen, "Sex", "Sociality", "MeanITD"),
           random_effects = "Site", reference = "Sociality,solitary",
           normalization = "NONE", transform = "LOG", min_prevalence = 0,
           plot_heatmap = FALSE, plot_scatter = FALSE)
