# fetch_sfts.R
# Downloads SFTS weekly case data from JIHS via jpinfect and exports clean CSVs.

suppressPackageStartupMessages({
  library(jpinfect)
  library(dplyr)
})

DISEASE_PATTERN <- "Severe fever with thrombocytopenia"
DATA_DIR        <- "raw_data"
OUTPUT_DIR      <- "data"
CURRENT_YEAR    <- as.integer(format(Sys.Date(), "%Y"))

dir.create(DATA_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== SFTS automated download:", format(Sys.time()), "===\n")

# ── 1. Confirmed historical data (built-in + optional fresh download) ─────────
# Use the built-in dataset first — always available, covers 1999-2024
data("sex_prefecture", package = "jpinfect")

# Optionally refresh from JIHS (comment out if server is slow / unavailable)
tryCatch(
  jpinfect_get_confirmed(type = "sex", dest_dir = DATA_DIR),
  error = function(e) cat("  NOTE: confirmed download skipped:", conditionMessage(e), "\n")
)

confirmed_raw <- tryCatch(
  jpinfect_read_confirmed(path = DATA_DIR, type = "sex"),
  error = function(e) {
    cat("  Falling back to built-in sex_prefecture dataset\n")
    sex_prefecture
  }
)

# ── 2. Provisional data (current year) ────────────────────────────────────────
tryCatch(
  jpinfect_get_bullet(year = CURRENT_YEAR, dest_dir = DATA_DIR),
  error = function(e) cat("  NOTE: bullet download skipped:", conditionMessage(e), "\n")
)

bullet_raw <- tryCatch(
  jpinfect_read_bullet(year = CURRENT_YEAR, directory = DATA_DIR),
  error = function(e) { cat("  WARNING: bullet import failed:", conditionMessage(e), "\n"); NULL }
)

# ── 3. Filter SFTS from confirmed data ────────────────────────────────────────
sfts_confirmed <- NULL
if (!is.null(confirmed_raw)) {
  confirmed_long <- jpinfect_pivot(confirmed_raw)
  
  cat("  Confirmed data columns:", paste(head(colnames(confirmed_long), 8), collapse = ", "), "\n")
  cat("  Sample prefecture values:", paste(head(unique(confirmed_long$prefecture), 5), collapse = ", "), "\n")
  
  sfts_confirmed <- confirmed_long |>
    filter(grepl(DISEASE_PATTERN, disease, ignore.case = TRUE)) |>
    mutate(source = "confirmed")
  
  cat("  Confirmed SFTS rows (all prefectures):", nrow(sfts_confirmed), "\n")
}

# ── 4. Filter SFTS from provisional data ─────────────────────────────────────
sfts_bullet <- NULL
if (!is.null(bullet_raw)) {
  sfts_col <- grep(DISEASE_PATTERN, colnames(bullet_raw), ignore.case = TRUE, value = TRUE)
  
  if (length(sfts_col) == 0) {
    cat("  WARNING: SFTS column not found in bullet data.\n")
    cat("  Available columns:\n", paste(" ", colnames(bullet_raw), collapse = "\n"), "\n")
  } else {
    cat("  SFTS column found:", sfts_col[1], "\n")
    sfts_bullet <- bullet_raw |>
      select(date, prefecture, cases = all_of(sfts_col[1])) |>
      mutate(disease = DISEASE_PATTERN, source = "provisional")
    cat("  Provisional SFTS rows:", nrow(sfts_bullet), "\n")
  }
}

# ── 5. Build national totals ──────────────────────────────────────────────────
# FIX: sum across ALL prefectures per week rather than filtering for "Total"
# (confirmed data has individual prefecture rows; bullet has a "Total" row
#  but we harmonise by always summing to avoid missing either source)

make_national <- function(df, source_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  df |>
    # Drop any pre-computed total rows to avoid double-counting
    filter(!grepl("^total$|^national$|^all$", prefecture, ignore.case = TRUE)) |>
    group_by(date, source) |>
    summarise(cases = sum(cases, na.rm = TRUE), .groups = "drop")
}

sfts_national_confirmed    <- make_national(sfts_confirmed, "confirmed")
sfts_national_provisional  <- make_national(sfts_bullet,   "provisional")

sfts_national <- bind_rows(sfts_national_confirmed, sfts_national_provisional) |>
  arrange(date)

# ── 6. All-prefecture data (keep individual rows, drop pre-computed totals) ───
sfts_all_pref <- bind_rows(
  sfts_confirmed |>
    filter(!grepl("^total$|^national$|^all$", prefecture, ignore.case = TRUE)) |>
    select(date, prefecture, cases, source),
  sfts_bullet |>
    filter(!grepl("^total$|^national$|^all$", prefecture, ignore.case = TRUE)) |>
    select(date, prefecture, cases, source)
) |>
  arrange(date, prefecture)

# ── 7. Export ─────────────────────────────────────────────────────────────────
if (nrow(sfts_national) > 0) {
  write.csv(sfts_national,  file.path(OUTPUT_DIR, "sfts_weekly_national.csv"),        row.names = FALSE)
  write.csv(sfts_all_pref,  file.path(OUTPUT_DIR, "sfts_weekly_all_prefectures.csv"), row.names = FALSE)
  cat("=== Done. National rows:", nrow(sfts_national),
      "| Year range:", format(min(sfts_national$date), "%Y"), "–",
      format(max(sfts_national$date), "%Y"), "===\n")
} else {
  stop("No SFTS data retrieved. Check download errors above.")
}

