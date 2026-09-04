## 03_openalex_faculty.R -- research-active faculty proxy from OpenAlex
##
## Definition (study design): an author is "research-active in the department"
## in the window oa_years if they have >= oa_min_works works in that window
## whose authorship lists the institution AND whose raw affiliation string
## for that author matches the department regex. Works are restricted to the
## life-sciences concept cluster so a chemist with a courtesy biology
## affiliation is not counted.
##
## Why works -> authors rather than the /authors endpoint: OpenAlex author
## records carry only the institution, not the department. Department lives in
## the per-work raw_affiliation_strings, so we must pull works and roll up.
##
## The count cannot separate faculty from postdocs and senior students who
## also clear 3 works; that is why 05_match.R uses the calibration ratio
## from data/internal/faculty_handcount.csv (five departments hand-counted
## from directory pages) to rescale, and reports both raw and rescaled.
##
## Output: data/derived/openalex_authors_<inst_key>.csv (one row per author
## with works in window, matched flag) and openalex_faculty_counts.csv.

source("R/00_config.R")
oa_dir <- file.path(raw_dir, "openalex"); dir.create(oa_dir, showWarnings = FALSE)

## Life-sciences scope: OpenAlex topics carry a field id. 11 = Agricultural
## and Biological Sciences, 13 = Biochemistry, Genetics and Molecular Biology,
## 28 = Neuroscience, 24 = Immunology and Microbiology. Domain 1 = Life
## Sciences covers all four plus 30 (Pharmacology) and 12 (Arts & Humanities
## is domain 2, not included). We use the domain filter for recall and keep
## the field id on each work so the cut can be narrowed later.
oa_filter_domain <- "primary_topic.domain.id:1"

## Resolve OpenAlex institution ids once, from the ROR in seed_pool.csv.
## ROR ids are stable; OpenAlex ids are looked up rather than hard-coded so a
## re-pull survives an OpenAlex id merge.
if (any(is.na(seed$openalex_id))) {
  for (i in which(is.na(seed$openalex_id))) {
    r <- get_json(sprintf("https://api.openalex.org/institutions/%s?mailto=%s", seed$ror[i], mailto))
    if (is.null(r)) { warning("could not resolve ROR for ", seed$inst_key[i]); next }
    ## The ROR ids were entered by hand; refuse to proceed if the resolved
    ## name does not contain the first word of the institution name.
    first_word <- strsplit(seed$institution[i], "[ -]")[[1]][1]
    if (!grepl(first_word, r$display_name, ignore.case = TRUE))
      stop("ROR ", seed$ror[i], " resolves to '", r$display_name, "', not ", seed$institution[i])
    message(seed$inst_key[i], " -> ", r$display_name, " (", sub("https://openalex.org/", "", r$id), ")")
    seed$openalex_id[i] <- sub("https://openalex.org/", "", r$id)
  }
  ## Write ids back into the full roster file, not the subset.
  full <- read.csv(file.path(proj_dir, "data", "seed_pool.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
  full$openalex_id[match(seed$inst_key, full$inst_key)] <- seed$openalex_id
  write.csv(full, file.path(proj_dir, "data", "seed_pool.csv"), row.names = FALSE)
}

pull_inst <- function(s) {
  base <- sprintf("https://api.openalex.org/works?filter=authorships.institutions.id:%s,publication_year:%d-%d,%s&select=id,publication_year,type,primary_topic,authorships",
                  s$openalex_id, min(oa_years), max(oa_years), oa_filter_domain)
  message("OpenAlex works: ", s$inst_key)
  pages <- oa_pages(base)
  ## Flatten to one row per (work, author) keeping only authors whose
  ## institution list includes this institution.
  rows <- lapply(pages, function(p) {
    if (length(p) == 0) return(NULL)
    do.call(rbind, lapply(seq_len(nrow(p)), function(j) {
      a <- p$authorships[[j]]
      if (is.null(a) || length(a) == 0) return(NULL)
      inst_ids <- vapply(a$institutions, function(x) paste(sub("https://openalex.org/", "", x$id), collapse = "|"), character(1))
      here <- grepl(s$openalex_id, inst_ids, fixed = TRUE)
      if (!any(here)) return(NULL)
      data.frame(work = p$id[j], year = p$publication_year[j], type = p$type[j],
                 field = if (is.null(p$primary_topic$field$id[j])) NA else sub("https://openalex.org/fields/", "", p$primary_topic$field$id[j]),
                 author_id = sub("https://openalex.org/", "", a$author$id[here]),
                 author = a$author$display_name[here],
                 affiliation = vapply(a$raw_affiliation_strings[here], function(x) paste(x, collapse = " || "), character(1)),
                 stringsAsFactors = FALSE)
    }))
  })
  x <- do.call(rbind, rows)
  x$inst_key <- s$inst_key
  x
}

counts <- list()
for (i in seq_len(nrow(seed))) {
  s <- seed[i, ]
  x <- pull_inst(s)
  if (is.null(x) || nrow(x) == 0) { warning("no works for ", s$inst_key); next }
  x <- x[x$type %in% c("article", "review", "preprint", "book-chapter"), ]   # drop editorials, errata, datasets
  x$in_dept <- grepl(s$dept_regex, x$affiliation, ignore.case = TRUE) &
    (is.na(s$dept_exclude_regex) | !grepl(s$dept_exclude_regex, x$affiliation, ignore.case = TRUE))
  save_raw(x, paste0("openalex_works_", s$inst_key), sub("&cursor.*", "", sprintf("works?filter=authorships.institutions.id:%s", s$openalex_id)),
           "OpenAlex", "work-author rows, life-sciences domain", subdir = "openalex")
  ## Roll up to author: works in window at this institution, and works with
  ## a department-matching affiliation string. An author is department
  ## research-active if the department-matched count clears the threshold.
  au <- aggregate(cbind(n_works = 1, n_dept = in_dept) ~ author_id + author, data = x, FUN = sum)
  au$active_dept <- au$n_dept >= oa_min_works
  au$active_inst <- au$n_works >= oa_min_works
  au$inst_key <- s$inst_key
  write.csv(au, file.path(derived_dir, paste0("openalex_authors_", s$inst_key, ".csv")), row.names = FALSE)
  counts[[s$inst_key]] <- data.frame(inst_key = s$inst_key,
                                     n_active_dept = sum(au$active_dept),
                                     n_active_inst_lifesci = sum(au$active_inst),
                                     n_works_dept = sum(x$in_dept), n_works = nrow(x),
                                     share_affil_nonempty = mean(nzchar(x$affiliation)))
}
counts <- do.call(rbind, counts)

## ---- calibration ---------------------------------------------------------------
## data/internal/faculty_handcount.csv: inst_key, n_tt_faculty, directory_url,
## date_counted. Five departments (tamu, purdue, florida, colostate, kstate).
## The ratio hand / OpenAlex is applied to every institution; its spread
## across the five is the honest measure of how noisy this covariate is.
hc_file <- file.path(int_dir, "faculty_handcount.csv")
if (file.exists(hc_file)) {
  hc <- read.csv(hc_file, stringsAsFactors = FALSE)
  cal <- merge(counts, hc[, c("inst_key", "n_tt_faculty")], by = "inst_key")
  cal$ratio <- cal$n_tt_faculty / cal$n_active_dept
  print(cal[, c("inst_key", "n_active_dept", "n_tt_faculty", "ratio")])
  k <- median(cal$ratio)
  message("calibration ratio (median hand/OpenAlex) = ", round(k, 3),
          "; range ", round(min(cal$ratio), 3), " to ", round(max(cal$ratio), 3))
  counts$faculty_calibrated <- counts$n_active_dept * k
  write.csv(cal, file.path(derived_dir, "openalex_faculty_calibration.csv"), row.names = FALSE)
} else {
  message("No hand-count file at data/internal/faculty_handcount.csv; calibrated column left NA")
  counts$faculty_calibrated <- NA
}
write.csv(counts, file.path(derived_dir, "openalex_faculty_counts.csv"), row.names = FALSE)
message("03_openalex_faculty done")
