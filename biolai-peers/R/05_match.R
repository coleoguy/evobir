## 05_match.R -- Mahalanobis matching on logged, 2019-2024-averaged covariates
##
## Covariates (one value per institution, log of the 2019-2024 mean):
##   bs_2601      IPEDS bachelor's, CIP 26.0101 + 26.0199, first major
##   phd          PhDs granted per year (coordinator reply > ProQuest > OpenAlex)
##   faculty      OpenAlex research-active authors matched to department
##   herd_bio     NSF HERD biological & biomedical sciences R&D expenditures (institution level)
##   pell         IPEDS Pell share of undergraduates (logged like the rest for
##                consistency; it is a proportion, so log is monotone and harmless)
## Distance: D_i = sqrt((x_i - x_T)' S^-1 (x_i - x_T)), with S the covariance of
## the logged covariates across the whole seed pool (TAMU included; it is one
## draw from the same population). Logging first is what makes a pooled S
## sensible when institutions differ tenfold in size.
##
## Hard filters (R1, land-grant/flagship, ag college, vet college) and the
## contamination cull are applied from seed_pool.csv AFTER the distance is
## computed so the full table shows every candidate's distance and the reason
## it was dropped; the Dean's office can see what a filter cost.
##
## HERD is not pulled by script: NCSES publishes it as interactive tables
## (https://ncses.nsf.gov/surveys/higher-education-research-development).
## Download "R&D expenditures by field: biological and biomedical sciences"
## for each fiscal year as CSV to data/raw/herd/herd_bio_<FY>.csv (columns:
## institution, unitid or name, amount_k). The manifest row is written here.

source("R/00_config.R")

## ---- assemble covariates -------------------------------------------------------
avg <- function(d, var, yrs = match_years) {
  d <- d[d$year %in% yrs, ]
  a <- aggregate(d[[var]], list(inst_key = d$inst_key), mean, na.rm = TRUE)
  names(a)[2] <- var; a
}
ipeds <- read.csv(file.path(derived_dir, "ipeds_completions_panel.csv"), stringsAsFactors = FALSE)
pell  <- read.csv(file.path(derived_dir, "ipeds_pell_panel.csv"), stringsAsFactors = FALSE)
fac   <- read.csv(file.path(derived_dir, "openalex_faculty_counts.csv"), stringsAsFactors = FALSE)

## PhDs: prefer coordinator-reported, then ProQuest, then OpenAlex dissertations.
phd_src <- c(coord = file.path(int_dir, "coordinator_replies.csv"),
             proquest = file.path(derived_dir, "phd_proquest_panel.csv"),
             openalex = file.path(derived_dir, "phd_openalex_panel.csv"))
phd_src <- phd_src[file.exists(phd_src)]
phd_list <- lapply(names(phd_src), function(n) {
  d <- read.csv(phd_src[[n]], stringsAsFactors = FALSE)
  v <- c(coord = "phd_granted", proquest = "phd_proquest", openalex = "diss_dept")[n]
  a <- avg(d, v); names(a)[2] <- "phd"; a$phd_source <- n; a
})
phd <- do.call(rbind, phd_list)
phd <- phd[!duplicated(phd$inst_key), ]         # first in priority order wins

## HERD: read whatever fiscal years are present.
herd_files <- list.files(file.path(raw_dir, "herd"), pattern = "^herd_bio_\\d{4}\\.csv$", full.names = TRUE)
if (length(herd_files) == 0) stop("no HERD files in data/raw/herd/; see header of this script")
herd <- do.call(rbind, lapply(herd_files, function(f) {
  x <- read.csv(f, stringsAsFactors = FALSE)
  x$year <- as.integer(sub(".*herd_bio_(\\d{4})\\.csv", "\\1", f))
  log_pull(f, "https://ncses.nsf.gov/surveys/higher-education-research-development", "NSF HERD", "biological and biomedical sciences R&D, $K")
  x
}))
herd <- merge(seed[, c("inst_key", "unitid")], herd, by = "unitid")
herd_a <- avg(herd, "amount_k"); names(herd_a)[2] <- "herd_bio"

cov <- Reduce(function(a, b) merge(a, b, all.x = TRUE), list(
  seed[, c("inst_key", "institution")],
  avg(ipeds, "bs_2601_narrow"), avg(ipeds, "bs_26_all"),
  phd[, c("inst_key", "phd", "phd_source")],
  fac[, c("inst_key", "n_active_dept", "faculty_calibrated")],
  herd_a, avg(pell, "pell_pct")))
names(cov)[names(cov) == "bs_2601_narrow"] <- "bs_2601"
cov$faculty <- if (all(is.na(cov$faculty_calibrated))) cov$n_active_dept else cov$faculty_calibrated
cov$pell <- cov$pell_pct
write.csv(cov, file.path(derived_dir, "covariates_raw.csv"), row.names = FALSE)

## ---- distance ---------------------------------------------------------------------
vars <- c("bs_2601", "phd", "faculty", "herd_bio", "pell")
X <- log(as.matrix(cov[, vars]) + 1)          # +1 guards a zero PhD year at a small unit
rownames(X) <- cov$inst_key
if (any(is.na(X))) {
  message("missing covariates for: ", paste(rownames(X)[!complete.cases(X)], collapse = ", "))
  X <- X[complete.cases(X), ]
}
S <- cov(X)
d2 <- mahalanobis(X, center = X[tamu_key, ], cov = S)
res <- data.frame(inst_key = rownames(X), D = sqrt(d2), stringsAsFactors = FALSE)

## Per-covariate standardised gaps, so the memo can say WHY a peer is near or far.
z <- sweep(X, 2, X[tamu_key, ]) / rep(sqrt(diag(S)), each = nrow(X))
colnames(z) <- paste0("z_", vars)
res <- cbind(res, round(z, 2))

## Sensitivity: replace narrow BS count with all-26 count.
Xb <- X; Xb[, "bs_2601"] <- log(cov$bs_26_all[match(rownames(X), cov$inst_key)] + 1)
res$D_bs26all <- sqrt(mahalanobis(Xb, Xb[tamu_key, ], cov(Xb)))

## ---- filters and cull ------------------------------------------------------------------
res <- merge(res, seed[, c("inst_key", "institution", "public", "carnegie_r1", "land_grant_or_flagship",
                           "has_ag_college", "has_vet_college", "ai_bio_verdict")], by = "inst_key")
res$pass_hard <- with(res, public == "Y" & carnegie_r1 == "Y" & land_grant_or_flagship == "Y" & has_ag_college == "Y")
res$pass_vet  <- res$pass_hard & res$has_vet_college == "Y"
res$pass_ai   <- res$ai_bio_verdict %in% c("NONE", "UNIVERSITY-LEVEL ONLY")
res$eligible  <- res$pass_hard & res$pass_ai & res$inst_key != tamu_key
res$drop_reason <- ifelse(res$inst_key == tamu_key, "treated unit",
                   ifelse(!res$pass_hard, "hard filter",
                   ifelse(res$ai_bio_verdict == "UNKNOWN", "AI program check not yet run",
                   ifelse(!res$pass_ai, "AI/data-science biology program", ""))))
res <- res[order(res$D), ]
res$rank_all <- seq_len(nrow(res)) - 1                          # TAMU is rank 0
res$rank_eligible <- NA; res$rank_eligible[res$eligible] <- seq_len(sum(res$eligible))
write.csv(res, file.path(derived_dir, "mahalanobis_table.csv"), row.names = FALSE)
print(res[, c("rank_all", "inst_key", "D", "pass_vet", "ai_bio_verdict", "drop_reason")], digits = 3)

## Nearest 8-10 eligible. We take 10 when the 9th and 10th are within 10% of
## the 8th's distance, otherwise 8, so the cut falls at a natural gap rather
## than an arbitrary count.
el <- res[res$eligible, ]
n_take <- if (nrow(el) >= 10 && el$D[10] <= 1.10 * el$D[8]) 10 else min(8, nrow(el))
peers <- el[seq_len(n_take), ]
peers$weight <- (1 / peers$D) / sum(1 / peers$D)               # inverse-distance, normalised
peers$vet_pass <- peers$pass_vet
peers$locked_date <- format(Sys.Date())
out <- peers[, c("inst_key", "institution", "D", "weight", "rank_eligible", "vet_pass", "ai_bio_verdict", "locked_date")]
write.csv(out, file.path(proj_dir, "data", "peers.csv"), row.names = FALSE)
message("05_match done: ", n_take, " peers written to data/peers.csv; ",
        sum(peers$vet_pass), " pass the vet-college filter. Tag with: git tag peers-v1")
