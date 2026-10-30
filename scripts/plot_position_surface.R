# plot_position_surface.R
# 3D surface plot of the position effect: column and row of the well on the plate
# (x and y axis) against the deviation of the well from the plate mean in % (z axis).
# Same deviation as in the heatmap (plot_linearity_heatmap.R). Meant for plates where
# all sample wells should be identical (plate 1). Base R only, no extra packages
# beyond tidyverse for reading the data.
# Input: long table from prepare_plate_long.R. Figure is saved as PNG.

library(tidyverse)

# ---- function: 12 x 8 matrix (columns x rows) of deviations in % ----------------
deviation_matrix <- function(plate_long, cycle_use = NULL) {
  if (is.null(cycle_use)) cycle_use <- max(plate_long$cycle)
  d <- plate_long |> filter(cycle == cycle_use)
  
  n_cond <- d |> filter(type == "sample") |> distinct(condition) |> nrow()
  if (n_cond > 1) {
    message("Note: ", n_cond, " different conditions on this plate - deviations are ",
            "relative to the mean of each condition")
  }
  
  cond <- d |> filter(type == "sample") |>
    group_by(condition) |>
    summarise(cond_mean = mean(od_net), .groups = "drop")
  d <- d |>
    left_join(cond, by = "condition") |>
    mutate(rel = if_else(type == "sample", 100 * (od_net / cond_mean - 1), NA_real_),
           row_no = match(row, LETTERS[1:8]))
  
  z <- matrix(NA_real_, nrow = 12, ncol = 8)           # z[column, row]
  z[cbind(d$col, d$row_no)] <- d$rel
  list(z = z, time_min = round(unique(d$time_min)))
}

# ---- function: fill gaps (blank / empty wells) for the DISPLAY only --------------
# A surface cannot have holes. Each gap gets the mean of its neighbouring wells.
fill_gaps <- function(z) {
  gaps <- which(is.na(z), arr.ind = TRUE)
  z_new <- z
  for (k in seq_len(nrow(gaps))) {
    i <- gaps[k, 1]; j <- gaps[k, 2]
    nb <- z[max(1, i - 1):min(12, i + 1), max(1, j - 1):min(8, j + 1)]
    z_new[i, j] <- mean(nb, na.rm = TRUE)
  }
  if (anyNA(z_new)) stop("Gaps are too large to be filled from neighbours")
  list(z = z_new, filled = paste0(LETTERS[gaps[, 2]], gaps[, 1]))
}

# ---- function: draw the surface ----------------------------------------------------
draw_surface <- function(z, time_min, theta = -35, phi = 30) {
  lim <- ceiling(max(abs(z)) / 5) * 5                  # symmetric colour range
  zlim <- c(-lim, lim)
  
  # colour of each little square = mean height of its four corners
  zf <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)]) / 4
  pal <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
  cols <- pal[cut(zf, breaks = seq(-lim, lim, length.out = 101), include.lowest = TRUE)]
  
  par(mar = c(1, 1, 3, 1))
  persp(x = 1:12, y = 1:8, z = z, zlim = zlim,
        theta = theta, phi = phi, expand = 0.6, col = cols,
        border = "grey60", lwd = 0.4, ticktype = "detailed", nticks = 6,
        xlab = "Spalte", ylab = "Reihe (1 = A ... 8 = H)",
        zlab = "Abweichung [%]",
        main = paste0("Positionseffekt nach ", time_min, " min"))
}

# ---- settings for this measurement (edit here) ---------------------------------
long_csv  <- "data/processed/14102026-Ecoli_Platte1-2_long.csv"
cycle_use <- NULL        # NULL = last cycle, otherwise a cycle number
theta     <- -35         # turns the plot around the vertical axis (degrees)
phi       <- 30          # viewing height (degrees)

# ---- run --------------------------------------------------------------------------
plate_long <- read_csv(long_csv, show_col_types = FALSE)
plate_id   <- sub("_long\\.csv$", "", basename(long_csv))
fig_dir    <- "results/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

dev_mat <- deviation_matrix(plate_long, cycle_use)
filled  <- fill_gaps(dev_mat$z)
message("Display only: ", length(filled$filled), " wells without sample (blank/empty) ",
        "filled with the mean of their neighbours: ",
        paste(filled$filled, collapse = ", "))

draw_surface(filled$z, dev_mat$time_min, theta, phi)      # shows it in RStudio

out_file <- file.path(fig_dir, paste0(plate_id, "_position_surface.png"))
png(out_file, width = 2000, height = 1600, res = 250)
draw_surface(filled$z, dev_mat$time_min, theta, phi)
invisible(dev.off())
cat("Saved:", out_file, "\n")