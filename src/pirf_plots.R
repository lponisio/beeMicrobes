# Partial dependence plots for PIRF fits
# needs: dplyr, tidyr, ggplot2 (loaded in main script)

# partial dependence for the top ASVs of one PIRF fit (one element of best_mtry() output):
# for each ASV, set its abundance to each grid value for every bee, predict with the
# PIRF forest, and average -> predicted P(detection) across that ASV's abundance range.
# grid runs from 0 to the 95th percentile of nonzero values (avoids stretching to outliers).
# slow step: n_top x n_grid predictions on the full data
pirf_pd <- function(pirf_fit, asv_tab, asv_tax, n_top = 6, n_grid = 15) {
  top_ids <- head(names(sort(pirf_fit$sel.prob, decreasing = TRUE)), n_top)

  pd_list <- list()
  for (asv_id in top_ids) {
    asv_vals  <- asv_tab[[asv_id]]
    grid_vals <- seq(0, quantile(asv_vals[asv_vals > 0], 0.95), length.out = n_grid)

    for (val in grid_vals) {
      tab_mod <- asv_tab
      tab_mod[[asv_id]] <- val                                  # every bee gets this value
      pd_list[[length(pd_list) + 1]] <- tibble(
        id   = asv_id,
        rel  = val,
        pred = mean(predict(pirf_fit$fit, data = tab_mod)$predictions))
    }
  }

  bind_rows(pd_list) |>
    left_join(asv_tax, by = "id") |>
    mutate(label = paste0(finest_taxon, " (", id, ")"),
           label = factor(label, levels = unique(label)))      # panels ordered by selection prob.
}

# plot pirf_pd() output: one panel per ASV.
# dashed line = overall detection rate (resp = the 0/1 outcome the model was fit to);
# ticks along the bottom = observed abundances, showing where the data actually are
plot_pirf_pd <- function(pd_df, asv_tab, resp, y_lab = "Predicted probability of detection") {
  obs_df <- asv_tab |>
    select(all_of(unique(pd_df$id))) |>
    pivot_longer(everything(), names_to = "id", values_to = "rel") |>
    left_join(pd_df |> group_by(id, label) |> summarise(max_rel = max(rel), .groups = "drop"),
              by = "id") |>
    filter(rel <= max_rel)                                      # ticks within plotted range

  ggplot(pd_df, aes(rel, pred)) +
    geom_hline(aes(yintercept = mean(resp)),
               linetype = "dashed", colour = "grey60") +
    geom_line() +
    geom_rug(data = obs_df, aes(x = rel), inherit.aes = FALSE, sides = "b", alpha = 0.1) +
    facet_wrap(~ label, scales = "free_x") +
    labs(x = "Relative abundance", y = y_lab) +
    theme(legend.position = "bottom")
}
