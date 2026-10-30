# prepare_plate_long.R

# Part 2 of the analysis: combine the raw wide table (output of import_lt4500_asc.R)
# with the plate layout, convert to a LONG table and subtract the blank.
# Output: one row per well and cycle, raw OD and blank-corrected OD side by side.
# No log transform ad no growth rates here.

library(tidyverse)

# ---- function ----------------------------------------------------------------
prepare_plate_long <- function(
    wide_csv,
    layout_csv,
    exclude_blank_wells = character(0),
    drift_threshold = 0.01) {
  
  # ---- 1. raw data (wide) ----------------------------------------------------
  wide <- read_csv(wide_csv, show_col_types = FALSE)
  names(wide)[1:3] <- c("cycle", "time_s", "temp_c")  # simpler names for analysis
  
  # ---- 2. layout: 8 * 12 grid, one label per well ----------------------------
  # works for comma- and semicolon-separated files (German Excel uses ';')
  first_line <- readLines(layout_csv, n = 1, warn = FALSE)
  layout_grid <- if (grepl(";", first_line)) {
    read_csv2(layout_csv, col_types = cols(.default = "c"))
  } else {
    read_csv(layout_csv, col_types = cols(.default = "c"))
  }
  layout <- layout_grid |>
    pivot_longer(-row, names_to = "col", values_to = "condition") |>
    mutate(well = paste0(row, col),
           row = factor(row, levels = LETTERS[1:8]),
           col = as.integer(col),
           type = case_when(
             is.na(condition)                            ~ "empty",
             str_detect(tolower(condition), "blank")     ~ "blank",
             tolower(condition) %in% c("empty", "leer")  ~ "empty",
             TRUE                                        ~ "sample"
           ))
  if (nrow(layout) != 96) stop("layout must contain 96 wells, found", nrow(layout))
  
  # ---- 3. wide -> long, add layout -------------------------------------------
  long <- wide |>
    pivot_longer(A1:H12, names_to = "well", values_to = "od_raw") |>
    left_join(layout, by = "well") |>
    mutate(time_min = time_s / 60)
  
  # ---- 4. blank: mean of the blank wells per cycle --------------------------
  blank_wells <- layout |>
    filter(type == "blank") |>
    pull(well)
  if (length(blank_wells) == 0) stop ("No blank wells found in the layout")
  message("Blank wells in layout: ", length(blank_wells),
          " | excluded from the blank mean: ", length(exclude_blank_wells),
          if (length(exclude_blank_wells) > 0)
            paste0(" (", paste(exclude_blank_wells, collapse = ", "), ")"))
  
  # drift check: blanks should stay constant
  drift <- long |>
    filter(well %in% blank_wells) |>
    group_by(well) |>
    summarise(change = last(od_raw) - first(od_raw), .groups = "drop") |>
    filter(change > drift_threshold)
  if (nrow(drift) > 0) {
    flagged <- drift |>
      arrange(desc(change)) |>
      mutate(txt = paste0(well, "(+", round(change, 3), ")")) |>
      pull(txt)
    message("WARNING: blank wells with no OD increase > ", drift_threshold,
            " (contamination?): ", paste(flagged, collapse = ", "))
    still_used <- setdiff(drift$well, exclude_blank_wells)
    if (length(still_used) > 0) {
      message(" still used for the blank mean: ", paste(still_used, collapse = ", "))
    }
  }
  
  blank_per_cycle <- long |>
    filter(well %in% setdiff(blank_wells, exclude_blank_wells)) |>
    group_by(cycle) |>
    summarise(blank_od = mean(od_raw), .groups = "drop")
  
  long |>
    left_join(blank_per_cycle , by = "cycle") |>
    mutate(od_net = od_raw - blank_od) |>
    select(cycle, time_s, time_min, temp_c, well, row, col, condition, type,
           od_raw, blank_od, od_net)
}

# ---- settings for this measurement (edit here) -------------------------------
wide_csv <- "data/processed/Marcel-Ecoli_Plate2_595nm_20261021_152420_wide.csv"
layout_csv <- "data/metadata/20102026-Ecoli_Plate2-2_layout.csv"
exclude_blank_wells <- character(0)  # contaminated blanks

# ---- run ---------------------------------------------------------------------
plate_long <- prepare_plate_long(wide_csv, layout_csv, exclude_blank_wells)

cat("Rows:", nrow(plate_long), "| cycles:", n_distinct(plate_long$cycle),
    "| wells:", n_distinct(plate_long$well), "\n")
print(count(distinct(plate_long, well, type), type))
print(plate_long |> filter(cycle == 1, well %in% c("A1", "B12", "E1")))

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
out_file <- sub("_wide\\.csv$", "_long.csv", wide_csv)
write_csv(plate_long, out_file)
cat("Saved:", out_file, "\n")