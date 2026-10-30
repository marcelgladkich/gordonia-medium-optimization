#Import_lt4500_asc.R

# Reads a kinetic plate reader export (LT4500, .asc, 96-well, one block per cycle)
# and returns a WIDE table: one row per cycle, one column per well
# (Cycle No., Time /s, Temp. /°C, A1, B1, ..., H1, A2, ..., H12).
# RAW VALUES ONLY: no blank correction, no ln transform, no layout information.
# All of that belongs into the analysis script. THe raw file is never modified.

library(tidyverse)

# ---- function ----------------------------------------------------------------
import_lt4500_asc <- function(path, interval_s, decimal_mark = ",", temp_c = NA_real_) {
  
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines, which = "right")       # remove trailing tabs / \r
  
  # every cycle starts with a header line "<>", rows A-H follow
  cycle <- cumsum(startsWith(lines, "<>"))
  is_row <- grepl("^[A-H](\t|$)", lines)        # also rowswith a lable only
  n_rows <- sum(is_row)
  
  # lines that are neither a block header nor a plate row (e.g. date, method name)
  other <- lines[!startsWith(lines, "<>") & !is_row & nzchar(lines)]
  if (length(other) > 0) {
    message("Lines not imported (no plate data): ",
            paste(head(other, 5), collapse = " | "),
            if (length(other) > 5) paste0(" | ... (", length(other), " lines in total)"))
  }
  if (n_rows == 0) stop ("No plate rows (A-H) found in", path)
  
  # split each row at the tabs; keep exactly 12 values (pad with NA if shorter)
  parts <- strsplit(lines[is_row], "\t", fixed = TRUE)
  mat <- t(vapply(parts, function(x) { v <- x[-1]; length(v) <- 12; v },
                  character(12)))
  
  # text -> numbers (decimal comma or point); empty / text cells become NA
  values <- suppressWarnings(
    as.numeric(gsub(decimal_mark, ".", as.vector(mat), fixed = TRUE))
  )
  
  long <- tibble(
    cycle = rep(cycle[is_row], times = 12),
    row = rep(substr(lines[is_row], 1, 1), times = 12),
    col = rep(1:12, each = n_rows),
    od = values
  )
  
  # keep only the complete cycles at the start of the file
  complete <- long |>
    group_by(cycle) |>
    summarise(ok = (n() == 96) & !anyNA(od))
  first_bad <- which(!complete$ok)[1]
  n_valid <- if (is.na(first_bad)) nrow(complete) else first_bad - 1
  if (n_valid == 0) {
    stop("No complete cycle found in ", path,
         " - check decimal_mark (currently'", decimal_mark, "') and the file format")
  }
  
  message("Cycles in file: ", nrow(complete),
          " | complete cycles kept: ", n_valid)
  if (n_valid < nrow(complete)) {
    message("NOT imported: cycles ", n_valid + 1, " to ", nrow(complete),
            " (incomplete block or no plate data) - please check the raw file")
  }
  
  # wells in the usual export order: A1, B1, ..., H1, A2, ...
  well_order <- paste0(rep(LETTERS[1:8], times = 12), rep(1:12, each = 8))
  
  wells <- long |>
    filter(cycle <= n_valid) |>
    mutate(well = paste0(row, col)) |>
    select(cycle, well, od) |>
    pivot_wider(names_from = well, values_from = od) |>
    select(cycle, all_of(well_order))
  
  out <- tibble(
    'Cycle Nr.' = wells$cycle,
    'Time [s]' = (wells$cycle - 1) * interval_s,    # from the set interval
    temp = temp_c     # NOMINAL set temperature, not measured
  ) |>
    bind_cols(select(wells, -cycle))
  
  names(out)[3] <- "Temp. [\u00b0C]"   # escape sequence
  out
}
  
  # ---- settings for this measurement (edit here) -----------------------------
  raw_file <- "data/raw/20102026-Ecoli_Plate2-2.asc"
  interval_s <- 300       # kinetic interval in seconds
  decimal <- ","          # decimal mar used in the .asc file
  temp_c <- 30            # assumed constant temperature in deg D (not measured!)
  
  # ---- run -------------------------------------------------------------------
  plate_wide <- import_lt4500_asc(raw_file, interval_s, decimal, temp_c)
  
  cat("Rows (cycles); ", nrow(plate_wide), "| columns:", ncol(plate_wide), "\n")
  cat("Missing OD values;", sum(is.na(select(plate_wide, A1:H12))), "\n")
  print(plate_wide[1:5, 1:8])
  
  dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path("data/processed",
                        paste0(tools::file_path_sans_ext(basename(raw_file)), "_wide.csv"))
  write_csv(plate_wide, out_file)
  cat("Saved:", out_file, "\n")
  