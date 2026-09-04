## 04_grants.R -- NSF and NIH award dollars attributed to the department
##
## Sources:
##   NSF Award Search API  https://api.nsf.gov/services/v1/awards.json
##   NIH RePORTER API      https://api.reporter.nih.gov/v2/projects/search
## Both return PI names and the awardee institution but not (for NSF) the
## department. Attribution to the department is by PI: an award counts if its
## PI (or any co-PI, NSF only) matches, by normalised "lastname f", an author
## flagged research-active in the department by 03_openalex_faculty.R.
## NIH also carries organization.dept_type (e.g. "BIOLOGY"); we keep it as a
## second attribution route and report both.
##
## Dollars are recorded by fiscal year of obligation (NSF fundsObligatedAmt
## by award start year; NIH award_amount by fiscal_year). 06_outcomes.R
## builds the 3-year rolling sum, which is the outcome, because single-year
## obligations are lumpy (a 5-year R01 obligates in one year at NSF but
## annually at NIH).
##
## Year window: pre_years plus every later year so the same script serves the
## annual re-pull; the pre-period cut is applied downstream.

source("R/00_config.R")
gr_dir <- file.path(raw_dir, "grants"); dir.create(gr_dir, showWarnings = FALSE)
yrs <- min(pre_years):as.integer(format(Sys.Date(), "%Y"))

## Department-active authors from 03, for PI matching.
au_files <- list.files(derived_dir, pattern = "^openalex_authors_.*\\.csv$", full.names = TRUE)
if (length(au_files) == 0) stop("run 03_openalex_faculty.R first")
au <- do.call(rbind, lapply(au_files, read.csv, stringsAsFactors = FALSE))
au <- au[au$active_dept, ]
## OpenAlex display names are "First M. Last"; take last token as surname.
au$key <- norm_name(sub(".* ", "", au$author), sub(" .*", "", au$author))
dept_keys <- split(au$key, au$inst_key)

## ---- NSF ------------------------------------------------------------------------
nsf_fields <- "id,title,piFirstName,piLastName,coPDPI,awardeeName,awardee,fundsObligatedAmt,estimatedTotalAmt,startDate,expDate,date,primaryProgram,fundProgramDir,cfdaNumber"
nsf_pull <- function(s) {
  out <- list(); off <- 1
  repeat {
    u <- sprintf("https://api.nsf.gov/services/v1/awards.json?awardeeName=%s&dateStart=01/01/%d&dateEnd=12/31/%d&printFields=%s&rpp=25&offset=%d",
                 URLencode(s$nsf_awardee, reserved = TRUE), min(yrs), max(yrs), nsf_fields, off)
    r <- get_json(u)
    a <- r$response$award
    if (is.null(a) || length(a) == 0) break
    out[[length(out) + 1]] <- a
    if (nrow(a) < 25) break
    off <- off + 25
  }
  if (length(out) == 0) return(NULL)
  x <- do.call(rbind, lapply(out, function(a) { a[setdiff(strsplit(nsf_fields, ",")[[1]], names(a))] <- NA; a[, strsplit(nsf_fields, ",")[[1]]] }))
  x$inst_key <- s$inst_key
  x
}
nsf <- do.call(rbind, lapply(seq_len(nrow(seed)), function(i) { message("NSF: ", seed$inst_key[i]); nsf_pull(seed[i, ]) }))
nsf$year <- as.integer(substr(nsf$startDate, 7, 10))          # startDate is MM/DD/YYYY
nsf$amount <- as.numeric(nsf$fundsObligatedAmt)
nsf$pi_key <- norm_name(nsf$piLastName, nsf$piFirstName)
## coPDPI is a list column of "First Last" strings.
nsf$copi_keys <- vapply(nsf$coPDPI, function(v) {
  if (is.null(v) || all(is.na(v))) return("")
  paste(norm_name(sub(".* ", "", v), sub(" .*", "", v)), collapse = ";")
}, character(1))
nsf$pi_in_dept <- mapply(function(k, i) k %in% dept_keys[[i]], nsf$pi_key, nsf$inst_key)
nsf$copi_in_dept <- mapply(function(k, i) any(strsplit(k, ";")[[1]] %in% dept_keys[[i]]), nsf$copi_keys, nsf$inst_key)
nsf$bio_directorate <- grepl("^BIO", nsf$fundProgramDir) | grepl("Biological Sciences", nsf$fundProgramDir)
nsf$coPDPI <- vapply(nsf$coPDPI, function(v) paste(v, collapse = "; "), character(1))
save_raw(nsf, "nsf_awards", "https://api.nsf.gov/services/v1/awards.json", "NSF Award Search API",
         "all awards to seed institutions; pi_in_dept by OpenAlex match", subdir = "grants")

## ---- NIH -------------------------------------------------------------------------
nih_pull <- function(s) {
  out <- list(); off <- 0
  repeat {
    body <- list(criteria = list(org_names = list(s$nih_org), fiscal_years = as.list(yrs)),
                 offset = off, limit = 500,
                 include_fields = list("ProjectNum", "FiscalYear", "AwardAmount", "ContactPiName",
                                       "PrincipalInvestigators", "Organization", "AgencyIcAdmin",
                                       "ProjectStartDate", "ProjectEndDate", "ActivityCode"))
    r <- post_json("https://api.reporter.nih.gov/v2/projects/search", body)
    x <- r$results
    if (is.null(x) || length(x) == 0) break
    out[[length(out) + 1]] <- data.frame(
      project = x$project_num, year = x$fiscal_year, amount = x$award_amount,
      pi = x$contact_pi_name, dept_type = x$organization$org_dept,
      org = x$organization$org_name, ic = x$agency_ic_admin$abbreviation,
      activity = x$activity_code, stringsAsFactors = FALSE)
    off <- off + 500
    if (nrow(x) < 500 || off >= 14999) break   # RePORTER caps offset+limit at 15000
    Sys.sleep(1)
  }
  if (length(out) == 0) return(NULL)
  x <- do.call(rbind, out); x$inst_key <- s$inst_key; x
}
nih <- do.call(rbind, lapply(seq_len(nrow(seed)), function(i) { message("NIH: ", seed$inst_key[i]); nih_pull(seed[i, ]) }))
## contact_pi_name is "LAST, FIRST M".
nih$pi_key <- norm_name(sub(",.*", "", nih$pi), trimws(sub("^[^,]*,", "", nih$pi)))
nih$pi_in_dept <- mapply(function(k, i) k %in% dept_keys[[i]], nih$pi_key, nih$inst_key)
nih$dept_is_bio <- grepl("^BIOLOG", toupper(nih$dept_type)) & !grepl("BIOCHEM|BIOMED", toupper(nih$dept_type))
save_raw(nih, "nih_projects", "https://api.reporter.nih.gov/v2/projects/search", "NIH RePORTER API",
         "all projects at seed institutions; pi_in_dept by OpenAlex match", subdir = "grants")

## ---- department-year totals ----------------------------------------------------------
tot <- function(d, sel, nm) { a <- aggregate(amount ~ inst_key + year, data = d[sel, ], FUN = sum); names(a)[3] <- nm; a }
grants <- Reduce(function(a, b) merge(a, b, all = TRUE), list(
  tot(nsf, nsf$pi_in_dept, "nsf_pi_dept"),
  tot(nsf, nsf$pi_in_dept | nsf$copi_in_dept, "nsf_anypi_dept"),
  tot(nsf, nsf$bio_directorate, "nsf_bio_directorate"),
  tot(nih, nih$pi_in_dept, "nih_pi_dept"),
  tot(nih, nih$dept_is_bio, "nih_dept_bio")))
grid <- expand.grid(inst_key = seed$inst_key, year = yrs, stringsAsFactors = FALSE)
grants <- merge(grid, grants, all.x = TRUE); grants[is.na(grants)] <- 0
grants$grants_dept <- grants$nsf_pi_dept + grants$nih_pi_dept      # headline measure
write.csv(grants[order(grants$inst_key, grants$year), ], file.path(derived_dir, "grants_panel.csv"), row.names = FALSE)

## PI match rate is the diagnostic that decides whether OpenAlex attribution is
## trustworthy: at TAMU it should be well above 50% of NSF BIO awards.
diag <- aggregate(cbind(n = 1, matched = pi_in_dept) ~ inst_key, data = nsf[nsf$bio_directorate, ], FUN = sum)
diag$match_rate <- diag$matched / diag$n
print(diag)
write.csv(diag, file.path(derived_dir, "grants_pi_match_diagnostic.csv"), row.names = FALSE)
message("04_grants done")
