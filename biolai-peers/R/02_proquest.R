## 02_proquest.R -- PhDs granted per year, by department
##
## ProQuest Dissertations & Theses (PQDT) has no public API. The workflow is:
##   1. This script prints one PQDT advanced-search string per institution.
##   2. Run each in PQDT (TAMU Libraries access), restrict to Doctoral +
##      years 2019-2025, and export the result list as "CSV" (Export/Save ->
##      CSV, all fields). Save each file as
##        data/raw/proquest/pqdt_<inst_key>_<YYYYMMDD>.csv
##   3. Re-run this script: it parses every export, applies the department
##      regex from seed_pool.csv to the "Department" column, and writes a
##      department-year count.
## The exports are raw pulls and are logged to the manifest like any other.
##
## Independent check: OpenAlex indexes dissertations (type = dissertation)
## with raw affiliation strings that usually name the department. Section 4
## pulls that count for the same window so the two can be compared, and so
## the annual re-pull has a fully automated fallback when a PQDT export is
## late. Coordinator-reported counts (data/internal/coordinator_replies.csv)
## are merged in 06_outcomes.R and are the number of record where available.

source("R/00_config.R")
pq_dir <- file.path(raw_dir, "proquest"); dir.create(pq_dir, showWarnings = FALSE)
yrs <- pre_years

## ---- 1. search strings -------------------------------------------------------
## PQDT field codes: SCH() = school/institution, DEP() = department,
## DG() = degree. The department clause is deliberately broad (every unit
## in the aggregation list); the regex filter below does the precise cut.
pq_terms <- vapply(strsplit(seed$proquest_dept_terms, "|", fixed = TRUE),
                   function(v) paste0('"', v, '"', collapse = " OR "), character(1))
pq_strings <- sprintf('SCH("%s") AND DEP(%s) AND DG("Ph.D.")', seed$proquest_school, pq_terms)
writeLines(paste(seed$inst_key, pq_strings, sep = "\t"),
           file.path(derived_dir, "proquest_search_strings.tsv"))
cat("PQDT search strings written to data/derived/proquest_search_strings.tsv\n")

## ---- 2. parse exports --------------------------------------------------------
files <- list.files(pq_dir, pattern = "^pqdt_.*_\\d{8}\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  message("No PQDT exports found in data/raw/proquest/. Run the searches, then re-run.")
} else {
  ## Keep only the newest export per institution so a re-pull supersedes.
  key <- sub("^pqdt_(.*)_\\d{8}\\.csv$", "\\1", basename(files))
  files <- tapply(files, key, function(f) sort(f)[length(f)])
  pq <- do.call(rbind, lapply(names(files), function(k) {
    x <- read.csv(files[[k]], stringsAsFactors = FALSE, check.names = FALSE)
    names(x) <- tolower(gsub("[^A-Za-z]", "", names(x)))
    ## PQDT CSV column names vary by export version; take the first match.
    dep <- x[[grep("^department", names(x))[1]]]
    yr  <- as.integer(substr(x[[grep("^(year|pubyear|publicationyear|degreeyear)", names(x))[1]]], 1, 4))
    deg <- x[[grep("^degree", names(x))[1]]]
    data.frame(inst_key = k, year = yr, department = dep, degree = deg, stringsAsFactors = FALSE)
  }))
  pq <- pq[grepl("ph\\.?d", pq$degree, ignore.case = TRUE) & pq$year %in% yrs, ]
  ## Department match: regex per institution from seed_pool.csv, applied to
  ## the PQDT department string. Aggregated peers list several units joined
  ## by "|"; exclusions (e.g. "Biomedical") are handled by the negative regex.
  ## PQDT "Department" is a bare program name ("Biology", "Botany"), so the
  ## match uses phd_program_regex, not the affiliation-string dept_regex.
  rx  <- setNames(seed$phd_program_regex, seed$inst_key)
  rxn <- setNames(seed$dept_exclude_regex, seed$inst_key)
  keep <- mapply(function(d, k) grepl(rx[k], trimws(d), ignore.case = TRUE) &&
                   (is.na(rxn[k]) || !grepl(rxn[k], d, ignore.case = TRUE)),
                 pq$department, pq$inst_key)
  pq$in_dept <- unname(keep)
  save_raw(pq, "proquest_parsed", "PQDT export (TAMU Libraries)", "ProQuest PQDT",
           "parsed exports; in_dept = matched department regex")
  phd <- aggregate(in_dept ~ inst_key + year, data = pq, FUN = sum)
  names(phd)[3] <- "phd_proquest"
  grid <- expand.grid(inst_key = names(files), year = yrs, stringsAsFactors = FALSE)
  phd <- merge(grid, phd, all.x = TRUE); phd$phd_proquest[is.na(phd$phd_proquest)] <- 0
  write.csv(phd[order(phd$inst_key, phd$year), ],
            file.path(derived_dir, "phd_proquest_panel.csv"), row.names = FALSE)
  ## Unmatched department strings, for regex tuning.
  um <- sort(table(pq$department[!pq$in_dept]), decreasing = TRUE)
  write.csv(data.frame(department = names(um), n = as.integer(um)),
            file.path(derived_dir, "proquest_unmatched_departments.csv"), row.names = FALSE)
}

## ---- 3. OpenAlex dissertations (automated cross-check) -------------------------
## Filter: type dissertation, institution, year. We then apply the same
## department regex to raw_affiliation_strings. OpenAlex coverage of
## dissertations is incomplete and uneven across institutions, so this is a
## check on trend, not a substitute for the PQDT count.
oa_diss <- lapply(seq_len(nrow(seed)), function(i) {
  s <- seed[i, ]
  if (is.na(s$openalex_id)) return(NULL)
  base <- sprintf("https://api.openalex.org/works?filter=type:dissertation,authorships.institutions.id:%s,publication_year:%d-%d&select=id,publication_year,authorships",
                  s$openalex_id, min(yrs), max(yrs))
  message("OpenAlex dissertations: ", s$inst_key)
  pages <- oa_pages(base, verbose = FALSE)
  if (length(pages) == 0) return(NULL)
  do.call(rbind, lapply(pages, function(p) {
    if (length(p) == 0) return(NULL)
    aff <- vapply(p$authorships, function(a) {
      if (is.null(a) || length(a) == 0) return("")
      paste(unlist(a$raw_affiliation_strings), collapse = " || ")
    }, character(1))
    data.frame(inst_key = s$inst_key, id = p$id, year = p$publication_year,
               affiliation = aff, stringsAsFactors = FALSE)
  }))
})
oa_diss <- do.call(rbind, oa_diss)
if (!is.null(oa_diss) && nrow(oa_diss) > 0) {
  rx  <- setNames(seed$dept_regex, seed$inst_key)
  rxn <- setNames(seed$dept_exclude_regex, seed$inst_key)
  oa_diss$in_dept <- mapply(function(d, k) grepl(rx[k], d, ignore.case = TRUE) &&
                              (is.na(rxn[k]) || !grepl(rxn[k], d, ignore.case = TRUE)),
                            oa_diss$affiliation, oa_diss$inst_key)
  save_raw(oa_diss, "openalex_dissertations", "https://api.openalex.org/works?filter=type:dissertation,...",
           "OpenAlex", "dissertations with raw affiliation strings")
  d <- aggregate(cbind(diss_all = 1, diss_dept = in_dept) ~ inst_key + year, data = oa_diss, FUN = sum)
  write.csv(d[order(d$inst_key, d$year), ], file.path(derived_dir, "phd_openalex_panel.csv"), row.names = FALSE)
}
message("02_proquest done")
