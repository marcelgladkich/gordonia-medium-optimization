# plot_linearity_heatmap.R

# Part 4 of the analysis: two questions about the measurement system.
#   1.  Linearity: does the plate reader OD halve when the culture is diluted 1:2?
#       (needs a dilution in the layout, e.g. plate 2; skipped otherwise)
#   2.  Position effect: do wells with the same content give the same OD,
#       or does it depend on where the well sits on the plate? (heatmap)
# Input: long table from prepare_plate_long.R   Figures are saved as png.

library(tidyverse)

# ---- helper: relative cell density from the condition label -----------------
# "S01 unverd." -> 1, "S03 1:4" -> 1/4, anything else -> NA
dilution_from_label <- function(x) {
  case_when(
    str_detect(tolower(x), "unverd|undil") ~ 1,
    str_detect(x, "1:\\d+") ~ 1 / as.numeric(str_extract(str_extract(x, "1:\\d+"),
                                                         "\\d+$")),
    TRUE ~ NA_real_
  )
}

# ---- 1. linearity ---------------------------------------------------------------
# Start of the run (first cycle): the cell number is still the pipetted one.
# Expected line = straight line through the origin, fitted only to the
# dilution levels in 'anchor_levels' (low OD, where the reader is assumed linear).
linearity_data <- function(plate_long, anchor_levels) {
  start <- plate_long |>
    filter(type == "sample", cycle == min(cycle)) |>
    mutate(rel_density = dilution_from_label(condition)) |>
    filter(!is.na(rel_density))
  
  anchor <- start |> filter(rel_density %in% anchor_levels)
  if (nrow(anchor) == 0) stop("None of the anchor levels found in the layout")
  slope <- sum(anchor$rel_density * anchor$od_net) / sum(anchor$rel_density^2)
  
  start <- start |>
    mutate(expected = slope * rel_density,
           pct_of_expected = 100 * od_net / expected)
  list(start = start, slope = slope)
}

level_label <- function(rel) ifelse(rel == 1, "unverd.", paste0("1:", round(1 / rel)))

plot_linearity <- function(lin, anchor_levels) {
  start <- lin$start |> filter(od_net > 0)       # log axis needs values > 0
  means <- start |> group_by(rel_density) |>
    summarise(od_net = mean(od_net), .groups = "drop")
  xs <- 10^seq(log10(min(start$rel_density)), 0, length.out = 50)
  ref <- tibble(rel_density = xs, od_net = lin$slope * xs)
  
  ggplot(start, aes(rel_density, od_net)) +
    geom_line(data = ref, linetype = "dashed", linewidth = 0.8) +
    geom_point(alpha = 0.35, size = 1.8, colour = "steelblue") +
    geom_point(data = means, size = 3.2, colour = "firebrick") +
    scale_x_log10(breaks = sort(unique(start$rel_density)),
                  labels = level_label) +
    scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
    labs(x = "Verd\u00fcnnung der Kultur", y = "OD zu Beginn (netto)",
         caption = paste0("gestrichelt: lineare Erwartung, angepasst an die Stufen ",
                          paste(level_label(sort(anchor_levels, decreasing = TRUE)),
                                collapse = ", "),
                          "\nblau: einzelne Wells | rot: Mittelwert pro Stufe")) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_deviation <- function(lin, min_net_od = 0.01) {
  lev <- lin$start |>
    group_by(rel_density) |>
    summarise(mean_net = mean(od_net), pct = mean(pct_of_expected), .groups = "drop") |>
    filter(mean_net >= min_net_od)             # below this the blank error dominates
  ggplot(lev, aes(mean_net, pct)) +
    geom_hline(yintercept = 100, linetype = "dashed") +
    geom_line(colour = "grey50") +
    geom_point(size = 3, colour = "firebrick") +
    geom_text(aes(label = level_label(rel_density)), vjust = -1, size = 3.8) +
    scale_x_log10() +
    scale_y_continuous(limits = c(0, 120)) +
    labs(x = "gemessene OD zu Beginn (netto, log)",
         y = "gemessen in % der linearen Erwartung",
         caption = paste0("nur Stufen mit Netto-OD >= ", min_net_od)) +
    theme_bw(base_size = 14)
}

# ---- 2. position effect -----------------------------------------------------
# Deviation of each sample well from the mean of its own condition (in %).
# Wells of conditions with a very low OD are grey (relative error too large).
plot_position_heatmap <- function(plate_long, cycle_use = NULL, min_net_od = 0.02) {
  if (is.null(cycle_use)) cycle_use <- max(plate_long$cycle)
  
  d <- plate_long |> filter(cycle == cycle_use)
  cond <- d |> filter(type == "sample") |>
    group_by(condition) |>
    summarise(cond_mean = mean(od_net), .groups = "drop")
  d <- d |>
    left_join(cond, by = "condition") |>
    mutate(rel = if_else(type == "sample" & cond_mean >= min_net_od,
                         100 * (od_net / cond_mean - 1), NA_real_),
           row = factor(row, levels = rev(LETTERS[1:8])))
  
  ggplot(d, aes(factor(col), row, fill = rel)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = round(rel)), size = 3.2, na.rm = TRUE) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-30, 30), oob = scales::squish,
                         na.value = "grey85",
                         name = "Abweichung vom\nMittel der \nBedingung / %") +
    scale_x_discrete(position = "top") +
    coord_equal() +
    labs(x = NULL, y = NULL,
         title = paste0("Positionseffekt nach ", round(unique(d$time_min)), "min"),
         caption = paste0("grau: Blank, leet oder Netto-OD < ", min_net_od)) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank())
}

# ---- settings for this measurement (edit here) ------------------------------
long_csv <- "data/processed/14102026-Ecoli_Platte1-2_long.csv"
anchor_levels <- 1 / c(32, 64, 128, 256)    # dilutions assumed to be in linear range
make_heatmap <- TRUE                        # TRUE only if all sample wells should be identical

# ---- run --------------------------------------------------------------------
plate_long <- read_csv(long_csv, show_col_types = FALSE)
plate_id <- sub("_long\\.csv$", "", basename(long_csv))
fig_dir <- "results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# linearity (only if the layout contains a dilution series)
n_levels <- plate_long |>
  filter(type == "sample") |>
  distinct(condition) |>
  mutate(rel = dilution_from_label(condition)) |>
  filter(!is.na(rel)) |> nrow()

if (n_levels >= 4) {
  lin <- linearity_data(plate_long, anchor_levels)
  cat("Slope of the linear expectation (net OD for undiluted culture):",
      round(lin$slope, 3), "\n")
  print(lin$start |>
          group_by(Stufe = level_label(rel_density)) |>
          summarise(Start_OD_netto = round(mean(od_net), 4),
                    CV_prozent     = round(100 * sd(od_net) / mean(od_net), 1),
                    pct_der_Erwartung = round(mean(pct_of_expected), 0),
                    .groups = "drop") |>
          arrange(desc(Start_OD_netto)), n = 20)
  ggsave(file.path(fig_dir, paste0(plate_id, "_linearity_loglog.png")),
         plot_linearity(lin, anchor_levels), width = 8, height = 5.5, dpi = 200)
  ggsave(file.path(fig_dir, paste0(plate_id, "_linearity_deviation.png")),
         plot_deviation(lin), width = 8, height = 5.5, dpi = 200)
} else {
  message("Linearity part skipped: no dilution series found in the layout labels")
}

# position effect heatmap (last cycle)
if (make_heatmap) {
  n_cond <- plate_long |> filter(type == "sample") |> distinct(condition) |> nrow()
  if (n_cond > 1) {
    message("Note: ", n_cond, " different conditions on this plate - deviations are ",
            "relative to the mean of each condition, low-OD conditions are noisy")
  }
  ggsave(file.path(fig_dir, paste0(plate_id, "_heatmap_lastcycle.png")),
         plot_position_heatmap(plate_long), width = 9, height = 6, dpi = 200)
} else {
  message("Heatmap skipped (make_heatmap <- FALSE)")
}

cat("Figures saved in", fig_dir, "\n")

