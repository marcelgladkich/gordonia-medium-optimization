# import_spark_xlsx.R
# Reads a kinetic plate reader export from the Tecan Spark (.xlsx, sheet "Result sheet",
# instrument settings in the first rows, data table starting at the row "Cycle Nr.")
# and returns the SAME wide table as import_lt4500_asc.R:
# (Cycle Nr., Time [s], Temp. [deg C], A1, B1, ..., H1, A2, ..., H12).
# Unlike the LT4500 file, the Spark file contains the REAL time and the measured
# temperature of every cycle, these are taken over unchanged.
# RAW VALUES ONLY: no blank correction, no log transform, no layout information.
# The raw file is never modified.

library(tidyverse)
library(readxl)

# ---- function ----------------------------------------------------------------
import_spark_xlsx <- function(path, sheet = 1) {
  
  # everything as text first, so that we can find the header row ourselves
  raw <- read_excel(path, sheet = sheet, col_names = FALSE, col_types = "text",
                    .name_repair = "minimal")
  first_col <- raw[[1]]
  
  # --- instrument settings from the top of the sheet (for the log only) ---
  setting <- function(label) {
    i <- which(startsWith(first_col, label))[1]
    if (is.na(i)) return(NA_character_)
    vals <- as.character(unlist(raw[i, -1]))
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) NA_character_ else vals[1]
  }
  message("Settings: wavelength ", setting("Measurement wavelength"), " nm | target ",
          setting("Target temperature"), " deg C | interval ", setting("Interval time"),
          " | duration ", setting("Kinetic duration"))
  message("Start: ", setting("Start Time"), " | end: ", setting("End Time"))
  warn_txt <- setting("Warning")
  if (!is.na(warn_txt)) {
    # the instrument lists every affected cycle - keep the message short
    message("INSTRUMENT WARNING: ", sub("for cycles.*$", "(for many cycles)", warn_txt))
  }
  
  # --- header row and data block ---
  hdr <- which(first_col == "Cycle Nr.")[1]
  if (is.na(hdr)) stop("Row 'Cycle Nr.' not found in ", path)
  header <- as.character(unlist(raw[hdr, ]))
  header <- header[!is.na(header)]
  
  well_order <- paste0(rep(LETTERS[1:8], times = 12), rep(1:12, each = 8))
  if (!all(well_order %in% header)) {
    stop("Wells missing in the header row: ",
         paste(setdiff(well_order, header), collapse = ", "))
  }
  
  block <- raw[(hdr + 1):nrow(raw), seq_along(header)]
  names(block) <- header
  block <- block |> select(`Cycle Nr.`, `Time [s]`, starts_with("Temp."), all_of(well_order))
  names(block)[1:3] <- c("cycle", "time_s", "temp_c")    # plain names while working
  
  # text -> numbers (also accepts a decimal comma)
  to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
  block <- block |> mutate(across(everything(), to_num))
  
  # keep only the complete cycles at the start of the data block
  complete  <- complete.cases(block)
  first_bad <- which(!complete)[1]
  n_valid   <- if (is.na(first_bad)) nrow(block) else first_bad - 1
  if (n_valid == 0) stop("No complete cycle found in ", path)
  message("Rows in data block: ", nrow(block), " | complete cycles kept: ", n_valid)
  if (n_valid < nrow(block)) {
    others <- block[(n_valid + 1):nrow(block), ]
    others <- others[rowSums(!is.na(others)) > 0, ]   # ignore completely empty rows
    if (nrow(others) > 0) {
      message("NOT imported: ", nrow(others), " incomplete rows after cycle ", n_valid,
              " - please check the raw file")
    }
  }
  out <- block[seq_len(n_valid), ]
  
  # real interval compared with the set interval
  real_s <- median(diff(out$time_s))
  set_hms <- as.numeric(strsplit(setting("Interval time"), ":")[[1]])
  set_s <- if (length(set_hms) == 3 && !anyNA(set_hms)) sum(set_hms * c(3600, 60, 1)) else NA_real_
  message("Real interval (median of Time [s]): ", round(real_s, 1), " s | set: ", set_s, " s")
  if (!is.na(set_s) && abs(real_s / set_s - 1) > 0.02) {
    message("WARNING: real interval differs from the set interval by more than 2 %")
  }
  names(out)[1:3] <- c("Cycle Nr.", "Time [s]", "Temp. [\u00b0C]")   # same names as LT4500 import
  out
}

# ---- settings for this measurement (edit here) ---------------------------------
raw_file <- "data/raw/Marcel-Ecoli_Plate2_595nm_20261021_152420.xlsx"

# ---- run --------------------------------------------------------------------------
plate_wide <- import_spark_xlsx(raw_file)

cat("Rows (cycles):", nrow(plate_wide), "| columns:", ncol(plate_wide), "\n")
cat("Missing OD values:", sum(is.na(select(plate_wide, A1:H12))), "\n")
print(plate_wide[1:5, 1:8])

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
out_file <- file.path("data/processed",
                      paste0(tools::file_path_sans_ext(basename(raw_file)), "_wide.csv"))
write_csv(plate_wide, out_file)
cat("Saved:", out_file, "\n")
