## 06_outcomes.R -- department-year outcome panel (run annually)
##
## Outcomes, one row per institution-year, for TAMU and every peer:
##   phd        PhDs granted (coordinator reply > ProQuest > OpenAlex)
##   grants3    3-year rolling sum of department-attributed NSF+NIH dollars
##   pubs       department-attributed works (OpenAlex, from 03 raw pulls)
##   new_trainees  authors whose FIRST work anywhere is in year t and carries
##              the department affiliation: the trainee-output proxy
##   placement  TAMU only: share of BS graduates entering PhD programs
##              (data/internal/placement.csv, columns year, n_grads, n_phd)
## Rolling sums and first-publication years need the full work-level history,
## so this script re-reads the raw OpenAlex pulls rather than the counts.

source("R/00_config.R")
peers <- read.csv(file.path(proj_dir, "data", "peers.csv"), stringsAsFactors = FALSE)
keys  <- c(tamu_key, peers$inst_key)
yrs   <- min(pre_years):as.integer(format(Sys.Date(), "%Y"))

## ---- PhDs ----
phd_src <- c(coord = file.path(int_dir, "coordinator_replies.csv"),
             proquest = file.path(derived_dir, "phd_proquest_panel.csv"),
             openalex = file.path(derived_dir, "phd_openalex_panel.csv"))
phd_src <- phd_src[file.exists(phd_src)]
phd <- do.call(rbind, lapply(names(phd_src), function(n) {
  d <- read.csv(phd_src[[n]], stringsAsFactors = FALSE)
  v <- c(coord = "phd_granted", proquest = "phd_proquest", openalex = "diss_dept")[n]
  data.frame(inst_key = d$inst_key, year = d$year, phd = d[[v]], phd_source = n)
}))
phd <- phd[!duplicated(phd[, c("inst_key", "year")]), ]

## ---- grants: 3-year rolling sum (t-2, t-1, t) ----
g <- read.csv(file.path(derived_dir, "grants_panel.csv"), stringsAsFactors = FALSE)
g <- g[g$inst_key %in% keys, ]
g <- g[order(g$inst_key, g$year), ]
roll3 <- function(x) stats::filter(x, rep(1, 3), sides = 1)   # NA for the first two years
g$grants3 <- ave(g$grants_dept, g$inst_key, FUN = roll3)

## ---- publications and new trainees from raw OpenAlex work-author rows ----
## The raw pull covers oa_years only; the annual re-pull of 03 extends it.
## For "first work anywhere" we query each candidate author's earliest year
## from OpenAlex once and cache it, because the pull is per-institution and
## an author's first paper may be elsewhere.
oa <- do.call(rbind, lapply(keys, function(k) {
  f <- tryCatch(latest_raw(paste0("openalex_works_", k), subdir = "openalex"), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  read.csv(f, stringsAsFactors = FALSE)
}))
oa <- oa[oa$in_dept, ]
pubs <- aggregate(work ~ inst_key + year, data = unique(oa[, c("inst_key", "year", "work")]), FUN = length)
names(pubs)[3] <- "pubs"

first_file <- file.path(derived_dir, "openalex_author_first_year.csv")
first <- if (file.exists(first_file)) read.csv(first_file, stringsAsFactors = FALSE) else data.frame(author_id = character(), first_year = integer())
todo <- setdiff(unique(oa$author_id), first$author_id)
if (length(todo) > 0) {
  message("looking up first publication year for ", length(todo), " authors")
  ## OpenAlex authors endpoint accepts up to 50 ids per filter call.
  chunks <- split(todo, ceiling(seq_along(todo) / 50))
  new <- list()
  for (ch in chunks) {
    r <- get_json(sprintf("https://api.openalex.org/authors?filter=ids.openalex:%s&select=id,counts_by_year&per-page=50&mailto=%s",
                          paste(ch, collapse = "|"), mailto))
    if (is.null(r)) { warning("OpenAlex author lookup failed; first_year left NA for remaining authors"); break }
    new[[length(new) + 1]] <- data.frame(
      author_id = sub("https://openalex.org/", "", r$results$id),
      first_year = vapply(r$results$counts_by_year, function(c) if (length(c) == 0) NA_integer_ else as.integer(min(c$year[c$works_count > 0])), integer(1)))
  }
  new <- do.call(rbind, new)
  first <- rbind(first, new)
  write.csv(first, first_file, row.names = FALSE)
}
oa$first_year <- first$first_year[match(oa$author_id, first$author_id)]
## counts_by_year only reaches back ~10 years; an author whose first work is
## older than that shows the window floor, which is fine: they are not new.
nt <- unique(oa[oa$year == oa$first_year, c("inst_key", "year", "author_id")])
nt <- aggregate(author_id ~ inst_key + year, data = nt, FUN = length); names(nt)[3] <- "new_trainees"

## ---- assemble ----
grid <- expand.grid(inst_key = keys, year = yrs, stringsAsFactors = FALSE)
panel <- Reduce(function(a, b) merge(a, b, all.x = TRUE), list(grid, phd, g[, c("inst_key", "year", "grants_dept", "grants3")], pubs, nt))
## An institution-year covered by the OpenAlex pull but with no matching
## works or first-year authors is a zero, not a missing value; leave NA only
## for years the raw pull does not cover.
covered <- paste(panel$inst_key, panel$year) %in% paste(oa$inst_key, oa$year) |
  panel$year %in% unique(oa$year)
panel$pubs[is.na(panel$pubs) & covered] <- 0
panel$new_trainees[is.na(panel$new_trainees) & covered] <- 0
panel$treated <- as.integer(panel$inst_key == tamu_key)
panel$post <- as.integer(panel$year >= treat_year)
panel$weight <- ifelse(panel$treated == 1, 1, peers$weight[match(panel$inst_key, peers$inst_key)])

## TAMU placement, internal rosters only.
pl_file <- file.path(int_dir, "placement.csv")
if (file.exists(pl_file)) {
  pl <- read.csv(pl_file, stringsAsFactors = FALSE)
  pl$inst_key <- tamu_key; pl$placement <- pl$n_phd / pl$n_grads
  panel <- merge(panel, pl[, c("inst_key", "year", "placement")], all.x = TRUE)
}
panel <- panel[order(panel$inst_key, panel$year), ]
write.csv(panel, file.path(derived_dir, sprintf("outcome_panel_%s.csv", format(Sys.Date(), "%Y%m%d"))), row.names = FALSE)
write.csv(panel, file.path(derived_dir, "outcome_panel.csv"), row.names = FALSE)
message("06_outcomes done: ", nrow(panel), " institution-years")
