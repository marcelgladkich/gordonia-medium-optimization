# plot_growth_curves.R
# Part 3 of the analysis: first overview figures from the long table
# (output of prepare_plate_long.R).
#   1. growth curves per condition (mean +/- SD), linear y axis
#   2. same curves with logarithmic y axis (display only, data are not changed)
#   3. blank wells over time (contamination check)
# Figures are saved as PNG for the presentation.

library(tidyverse)

# ---- functions --------------------------------------------------------------
# mean +/- SD per condition; optionally the single wells as thin grey lines
plot_curves <- function(plate_long, show_wells = FALSE, log_y = FALSE) {
  
  samples <- plate_long |> filter(type == "sample")
  
  summ <- samples |>
    group_by(condition, time_min) |>
    summarise(mean = mean(od_net), sd = sd(od_net), .groups = "drop")
  
  if (log_y) {
    n_before <- nrow(summ)
    summ <- summ |> filter(mean > 0.003)   # net OD ~ 0 cannot be shown on a log axis
    if (nrow(summ) < n_before) {
      message("Log axis: ", n_before - nrow(summ),
              " points with mean net OD <= 0.003 not shown")
    }
    if (show_wells) {                       # single wells: values <= 0 cannot be logged
      n_wells_before <- nrow(samples)
      samples <- samples |> filter(od_net > 0)
      if (nrow(samples) < n_wells_before) {
        message("Log axis: ", n_wells_before - nrow(samples),
                " single-well points with net OD <= 0 not shown")
      }
    }
  }
  if (show_wells && n_distinct(samples$condition) > 1) {
    message("Note: show_wells = TRUE on a plate with several conditions - the grey ",
            "single-well lines are not assigned to a condition (meant for plate 1)")
  }
  
  p <- ggplot()
  if (show_wells) {
    p <- p + geom_line(data = samples,
                       aes(time_min, od_net, group = well),
                       colour = "grey75", linewidth = 0.3)
  }
  p <- p +
    geom_ribbon(data = summ,
                aes(time_min, ymin = pmax(mean - sd, 0.001), ymax = mean + sd,
                    fill = condition),
                alpha = 0.25) +
    geom_line(data = summ, aes(time_min, mean, colour = condition),
              linewidth = 0.9) +
    labs(x = "Zeit [min]", y = "OD (netto, Blank abgezogen)",
         colour = NULL, fill = NULL) +
    theme_bw(base_size = 14)
  
  if (log_y) p <- p + scale_y_log10()
  if (n_distinct(samples$condition) > 1) {
    p <- p + scale_colour_viridis_d(option = "turbo") +
      scale_fill_viridis_d(option = "turbo")
  } else {
    p <- p + guides(colour = "none", fill = "none")
  }
  p
}

# raw OD of all blank wells; excluded wells are dashed
plot_blanks <- function(plate_long, exclude_blank_wells = character(0)) {
  blanks <- plate_long |>
    filter(type == "blank") |>
    mutate(used = factor(!(well %in% exclude_blank_wells),
                         levels = c(TRUE, FALSE),
                         labels = c("im Blank-Mittel", "ausgeschlossen")))
  ggplot(blanks, aes(time_min, od_raw, colour = well, linetype = used)) +
    geom_line(linewidth = 0.8) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "Zeit [min]", y = "OD (roh)", colour = "Blank-Well",
         linetype = NULL) +
    theme_bw(base_size = 14)
}

# ---- settings for this measurement (edit here) ---------------------------------
long_csv   <- "data/processed/Marcel-Ecoli_Plate2_595nm_20261021_152420_long.csv"
show_wells <- FALSE                              # TRUE: single wells as grey lines
exclude_blank_wells <- character(0)   # only used for the dashed lines
# plate 1: long_csv <- "data/processed/14102026-Ecoli_Platte1-2_long.csv"
#          show_wells <- TRUE;  exclude_blank_wells <- character(0)

# ---- run --------------------------------------------------------------------------
plate_long <- read_csv(long_csv, show_col_types = FALSE)

plate_id <- sub("_long\\.csv$", "", basename(long_csv))
fig_dir  <- "results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

p1 <- plot_curves(plate_long, show_wells = show_wells, log_y = FALSE)
p2 <- plot_curves(plate_long, show_wells = show_wells, log_y = TRUE)
p3 <- plot_blanks(plate_long, exclude_blank_wells)

ggsave(file.path(fig_dir, paste0(plate_id, "_growth_linear.png")), p1,
       width = 9, height = 5.5, dpi = 200)
ggsave(file.path(fig_dir, paste0(plate_id, "_growth_log.png")), p2,
       width = 9, height = 5.5, dpi = 200)
ggsave(file.path(fig_dir, paste0(plate_id, "_blanks.png")), p3,
       width = 9, height = 5.5, dpi = 200)

cat("Figures saved in", fig_dir, "\n")