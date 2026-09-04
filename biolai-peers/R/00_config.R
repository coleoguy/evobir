## 00_config.R -- shared configuration and helpers for the BiolAI peer pipeline
##
## Every other script sources this file first. It defines the project paths,
## the study window, the seed pool, and three helpers that every pull uses:
##   get_json()   polite GET returning parsed JSON, with retry and cursor paging
##   save_raw()   write a raw pull to data/raw/ with a date stamp
##   log_pull()   append a manifest row (file, url, source, date, notes)
## The manifest is the audit trail: every number in the memo must trace to a
## row here, so no script writes to data/raw/ except through save_raw().
##
## Base R only. jsonlite and curl are the only add-on packages used anywhere
## in the pipeline; both are installed on the TAMU HPRC login nodes.

suppressPackageStartupMessages({
  library(jsonlite)
  library(curl)
})

## ---- paths -----------------------------------------------------------------
## Scripts are run from the biolai-peers/ directory (Rscript R/01_ipeds.R).
proj_dir    <- normalizePath(".")
raw_dir     <- file.path(proj_dir, "data", "raw")
derived_dir <- file.path(proj_dir, "data", "derived")
int_dir     <- file.path(proj_dir, "data", "internal")   # not committed
for (d in c(raw_dir, derived_dir, int_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
manifest_file <- file.path(raw_dir, "manifest.csv")

## ---- study design constants ------------------------------------------------
## Matching covariates are averaged over match_years. Pre-period for the DiD
## is 2019-2025 and BiolAI is treated as starting in treat_year; the event
## study uses the last pre-treatment year as the omitted reference.
match_years <- 2019:2024
pre_years   <- 2019:2025
treat_year  <- 2026
## OpenAlex research-active window: authors with >= 3 works in these years.
oa_years    <- 2021:2025
oa_min_works <- 3

## Contact address for OpenAlex polite pool and NSF/NIH user agents.
## Set BIOLAI_MAILTO in ~/.Renviron; falls back to a placeholder that still
## works but is rate-limited.
mailto <- Sys.getenv("BIOLAI_MAILTO", unset = "biolai-peers@tamu.edu")

## ---- seed pool -------------------------------------------------------------
## One row per candidate institution. Structural fields (colleges, unit lists,
## contamination verdict) were filled by hand from public web sources on
## 2026-09-04; see data/seed_pool_sources.md. Quantitative fields are filled
## by scripts 01-04 and consumed by 05_match.R.
seed <- read.csv(file.path(proj_dir, "data", "seed_pool.csv"),
                 stringsAsFactors = FALSE, na.strings = c("", "NA"))
stopifnot(all(c("inst_key", "institution", "unitid", "openalex_query",
                "dept_regex", "phd_program_regex", "ai_bio_verdict") %in% names(seed)))
tamu_key <- "tamu"

## ---- helpers ---------------------------------------------------------------

## GET a URL and parse JSON. Retries on transient failure with exponential
## backoff because OpenAlex and NSF both throttle bursts. Returns NULL after
## the last failed attempt so callers can decide whether to abort.
## BIOLAI_TRIES=1 in the environment makes offline test runs fail fast.
default_tries <- as.integer(Sys.getenv("BIOLAI_TRIES", unset = "5"))

get_json <- function(url, tries = default_tries, pause = 1, ...) {
  h <- new_handle(useragent = paste0("biolai-peers (mailto:", mailto, ")"))
  for (i in seq_len(tries)) {
    res <- tryCatch(curl_fetch_memory(url, handle = h), error = function(e) NULL)
    if (!is.null(res) && res$status_code == 200) {
      return(fromJSON(rawToChar(res$content), simplifyVector = TRUE, ...))
    }
    Sys.sleep(pause * 2^(i - 1))
  }
  warning("get_json failed after ", tries, " tries: ", url)
  NULL
}

## POST a JSON body (NIH RePORTER uses POST for search).
post_json <- function(url, body, tries = default_tries, pause = 1) {
  h <- new_handle(useragent = paste0("biolai-peers (mailto:", mailto, ")"))
  handle_setopt(h, customrequest = "POST",
                postfields = toJSON(body, auto_unbox = TRUE))
  handle_setheaders(h, "Content-Type" = "application/json",
                    "Accept" = "application/json")
  for (i in seq_len(tries)) {
    res <- tryCatch(curl_fetch_memory(url, handle = h), error = function(e) NULL)
    if (!is.null(res) && res$status_code == 200) {
      return(fromJSON(rawToChar(res$content), simplifyVector = TRUE))
    }
    Sys.sleep(pause * 2^(i - 1))
  }
  warning("post_json failed after ", tries, " tries: ", url)
  NULL
}

## Walk an OpenAlex cursor. `base` is the URL without the cursor parameter.
## Returns a list of page results (each a data.frame from $results) which the
## caller rbinds; we keep pages separate so a single bad page can be inspected.
oa_pages <- function(base, per_page = 200, max_pages = Inf, verbose = TRUE) {
  cursor <- "*"; pages <- list(); n <- 0
  while (!is.null(cursor) && n < max_pages) {
    u <- paste0(base, "&per-page=", per_page, "&cursor=", URLencode(cursor, reserved = TRUE),
                "&mailto=", mailto)
    p <- get_json(u)
    if (is.null(p)) break
    n <- n + 1
    pages[[n]] <- p$results
    cursor <- p$meta$next_cursor
    if (verbose && n %% 10 == 0) message("  page ", n, " (", p$meta$count, " total)")
    if (length(p$results) == 0) break
  }
  pages
}

## Write a raw pull. `x` may be a data.frame (csv) or anything else (json).
## The date stamp in the filename is what lets the annual re-pull coexist with
## the locked baseline instead of overwriting it.
save_raw <- function(x, name, url, source, notes = "", subdir = NULL,
                     date = format(Sys.Date(), "%Y%m%d")) {
  d <- if (is.null(subdir)) raw_dir else file.path(raw_dir, subdir)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  ext <- if (is.data.frame(x)) "csv" else "json"
  f <- file.path(d, sprintf("%s_%s.%s", name, date, ext))
  if (ext == "csv") write.csv(x, f, row.names = FALSE) else write_json(x, f, auto_unbox = TRUE, pretty = TRUE)
  log_pull(f, url, source, notes)
  invisible(f)
}

## Append a manifest row. Creates the manifest on first use.
log_pull <- function(file, url, source, notes = "") {
  row <- data.frame(file = sub(paste0("^", proj_dir, "/"), "", file),
                    url = url, source = source,
                    pulled = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                    n_bytes = file.info(file)$size, notes = notes,
                    stringsAsFactors = FALSE)
  write.table(row, manifest_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(manifest_file), append = file.exists(manifest_file),
              qmethod = "double")
  invisible(row)
}

## Download a binary file (IPEDS zips) with the same retry/backoff.
download_raw <- function(url, dest, tries = default_tries, pause = 2) {
  if (file.exists(dest) && file.info(dest)$size > 0) return(invisible(dest))
  for (i in seq_len(tries)) {
    ok <- tryCatch({ curl_download(url, dest, quiet = TRUE); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(dest))
    Sys.sleep(pause * 2^(i - 1))
  }
  stop("download failed: ", url)
}

## Most recent dated raw file matching a stem, e.g. latest_raw("ipeds_completions").
## Used by downstream scripts so they read the newest pull by default while an
## explicit date can be passed to reproduce the locked baseline.
latest_raw <- function(stem, subdir = NULL, date = NULL) {
  d <- if (is.null(subdir)) raw_dir else file.path(raw_dir, subdir)
  f <- list.files(d, pattern = paste0("^", stem, "_\\d{8}\\.(csv|json)$"), full.names = TRUE)
  if (!is.null(date)) f <- f[grepl(date, f)]
  if (length(f) == 0) stop("no raw file for ", stem)
  sort(f)[length(f)]
}

## Normalise a person name to "lastname f" for PI-to-author matching.
## Diacritics, punctuation and middle names are stripped because NSF, NIH and
## OpenAlex each format names differently.
norm_name <- function(last, first) {
  l <- tolower(iconv(last, to = "ASCII//TRANSLIT")); l <- gsub("[^a-z]", "", l)
  f <- tolower(iconv(first, to = "ASCII//TRANSLIT")); f <- substr(gsub("[^a-z]", "", f), 1, 1)
  paste(l, f)
}
