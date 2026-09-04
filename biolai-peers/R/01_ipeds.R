## 01_ipeds.R -- IPEDS completions (CIP 26.01xx bachelor's) and Pell share
##
## Sources (all public, no login):
##   Completions   https://nces.ed.gov/ipeds/datacenter/data/C{YYYY}_A.zip
##   Directory     https://nces.ed.gov/ipeds/datacenter/data/HD{YYYY}.zip
##   Financial aid https://nces.ed.gov/ipeds/datacenter/data/SFA{YY}{YY+1}.zip
## IPEDS labels a completions file by the survey year: C2023_A covers degrees
## conferred 2022-07-01 to 2023-06-30. We label each row by that survey year.
##
## Why CIP 26.01xx only: the peer unit is the general biology department, and
## 26.0101 (Biology, General) is the code every such department reports its
## main BS under. 26.02-26.13 (biochemistry, microbiology, zoology, ecology,
## ...) are excluded per the study design because at ag/vet-heavy campuses
## those codes are dominated by ag- and vet-college majors.
##
## Caution flagged for the PI: 26.0102 "Biomedical Sciences, General" is inside
## the 26.01 group but at TAMU it is the vet college's BIMS major (~1,000
## degrees/yr), which the design explicitly excludes. We therefore report it
## as its own column and use 26.0101 + 26.0199 as the default covariate.
## Split-biology peers (UC Davis, Minnesota, ...) report some general-biology
## degrees under 26.05/26.07/26.13, so the narrow count undercounts them; the
## all-26 column is kept for a sensitivity check in 05_match.R.

source("R/00_config.R")

years   <- match_years                      # 2019:2024
ipeds   <- "https://nces.ed.gov/ipeds/datacenter/data/"
zip_dir <- file.path(raw_dir, "ipeds_zip"); dir.create(zip_dir, showWarnings = FALSE)

## Read the data csv out of an IPEDS zip. NCES ships a revised file
## (suffix _rv) alongside the original when one exists; prefer it.
read_ipeds_zip <- function(zip) {
  files <- unzip(zip, list = TRUE)$Name
  csv <- files[grepl("\\.csv$", files, ignore.case = TRUE)]
  csv <- if (any(grepl("_rv", csv))) csv[grepl("_rv", csv)][1] else csv[1]
  x <- read.csv(unz(zip, csv), stringsAsFactors = FALSE, check.names = FALSE)
  names(x) <- toupper(trimws(names(x)))
  x
}

## ---- 1. institution names (sanity check on UNITIDs) -----------------------
hd_zip <- file.path(zip_dir, "HD2023.zip")
download_raw(paste0(ipeds, "HD2023.zip"), hd_zip)
log_pull(hd_zip, paste0(ipeds, "HD2023.zip"), "IPEDS HD2023", "institution directory")
hd <- read_ipeds_zip(hd_zip)
hd <- hd[hd$UNITID %in% seed$unitid, c("UNITID", "INSTNM", "STABBR", "CONTROL")]
chk <- merge(seed[, c("inst_key", "institution", "unitid")], hd,
             by.x = "unitid", by.y = "UNITID", all.x = TRUE)
print(chk[, c("inst_key", "institution", "INSTNM")])
if (any(is.na(chk$INSTNM))) stop("UNITID not found in HD2023 for: ",
                                 paste(chk$inst_key[is.na(chk$INSTNM)], collapse = ", "))
save_raw(chk, "ipeds_unitid_check", paste0(ipeds, "HD2023.zip"), "IPEDS HD2023")

## ---- 2. completions ----------------------------------------------------------
comp <- lapply(years, function(y) {
  z <- file.path(zip_dir, sprintf("C%d_A.zip", y))
  u <- paste0(ipeds, sprintf("C%d_A.zip", y))
  ok <- tryCatch({ download_raw(u, z); TRUE }, error = function(e) FALSE)
  if (!ok) { warning("completions not available for ", y); return(NULL) }
  log_pull(z, u, sprintf("IPEDS C%d_A", y), "completions by CIP")
  x <- read_ipeds_zip(z)
  x <- x[x$UNITID %in% seed$unitid & x$AWLEVEL == 5, ]   # 5 = bachelor's
  x$CIPCODE <- sprintf("%.4f", as.numeric(x$CIPCODE))    # "26.0101" (some files drop the trailing zero)
  x$year <- y
  x[, c("UNITID", "year", "CIPCODE", "MAJORNUM", "CTOTALT")]
})
comp <- do.call(rbind, comp)
save_raw(comp, "ipeds_completions_cip26_raw", paste0(ipeds, "C{YYYY}_A.zip"),
         "IPEDS Completions", "bachelor's, seed pool, all CIPs, both majors")

## Aggregate. MAJORNUM 1 = first major; IPEDS counts second majors separately
## and we follow the convention of counting first majors only so that a
## double major in Biology + Chemistry is not two biology degrees.
comp1 <- comp[comp$MAJORNUM == 1, ]
agg <- function(sel, nm) {
  a <- aggregate(CTOTALT ~ UNITID + year, data = comp1[sel, ], FUN = sum)
  names(a)[3] <- nm; a
}
is26     <- substr(comp1$CIPCODE, 1, 3) == "26."
narrow   <- comp1$CIPCODE %in% c("26.0101", "26.0199")
bims     <- comp1$CIPCODE == "26.0102"
out <- Reduce(function(a, b) merge(a, b, all = TRUE),
              list(agg(narrow, "bs_2601_narrow"), agg(bims, "bs_260102_biomed"),
                   agg(is26, "bs_26_all")))
out[is.na(out)] <- 0
out <- merge(seed[, c("inst_key", "unitid")], out, by.x = "unitid", by.y = "UNITID")
## Make sure every institution-year is present even if it had zero degrees.
grid <- expand.grid(inst_key = seed$inst_key, year = unique(comp$year), stringsAsFactors = FALSE)
out <- merge(grid, out, all.x = TRUE); out[is.na(out)] <- 0
out <- out[order(out$inst_key, out$year), ]
write.csv(out, file.path(derived_dir, "ipeds_completions_panel.csv"), row.names = FALSE)

## ---- 3. Pell share -----------------------------------------------------------
## SFA{yy}{yy+1}: e.g. SFA2223 = award year 2022-23. UPGRNTP is "percent of
## undergraduate students awarded Pell grants". Older files used PGRNT_P for
## first-time full-time only; we want the all-undergraduate measure.
pell <- lapply(years, function(y) {
  tag <- sprintf("SFA%02d%02d", (y - 1) %% 100, y %% 100)
  z <- file.path(zip_dir, paste0(tag, ".zip")); u <- paste0(ipeds, tag, ".zip")
  ok <- tryCatch({ download_raw(u, z); TRUE }, error = function(e) FALSE)
  if (!ok) { warning("SFA not available for ", y); return(NULL) }
  log_pull(z, u, paste0("IPEDS ", tag), "student financial aid")
  x <- read_ipeds_zip(z)
  v <- intersect(c("UPGRNTP", "PGRNT_P"), names(x))
  if (length(v) == 0) stop("Pell variable not found in ", tag, "; columns: ", paste(head(names(x), 40), collapse = " "))
  x <- x[x$UNITID %in% seed$unitid, c("UNITID", v[1])]
  data.frame(UNITID = x$UNITID, year = y, pell_pct = as.numeric(x[[v[1]]]), pell_var = v[1])
})
pell <- do.call(rbind, pell)
pell <- merge(seed[, c("inst_key", "unitid")], pell, by.x = "unitid", by.y = "UNITID")
save_raw(pell, "ipeds_pell", paste0(ipeds, "SFA{YY}{YY+1}.zip"), "IPEDS SFA", "UPGRNTP, all undergraduates")
write.csv(pell[order(pell$inst_key, pell$year), ], file.path(derived_dir, "ipeds_pell_panel.csv"), row.names = FALSE)

message("01_ipeds done: ", nrow(out), " institution-years of completions, ", nrow(pell), " of Pell")
