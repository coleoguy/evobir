## 08_directory_counts.R -- count faculty by rank and graduate students from
## department directory pages. Run on a machine with open web access:
##   Rscript R/08_directory_counts.R
## Every page is saved to data/raw/directories/<key>_<page>_<date>.html and its
## visible text to a .txt, so the count can be checked by eye. Counts are
## keyword tallies of rank labels in the page text; directories that load
## names by JavaScript will show near-zero counts, and for those the .txt
## file will make it obvious and a browser "select all, copy" into the same
## file is the fallback.

source("R/00_config.R")
dir_dir <- file.path(raw_dir, "directories"); dir.create(dir_dir, showWarnings = FALSE)

pages <- read.csv(text = '
key,page,url
tamu,faculty,https://www.bio.tamu.edu/full-directory/faculty-directory/
tamu,instructional,https://artsci.tamu.edu/biology/contact/instructional-faculty.html
tamu,grad,https://artsci.tamu.edu/biology/contact/graduate-students.html
purdue,faculty,https://www.bio.purdue.edu/People/faculty/index.html
purdue,grad,https://www.bio.purdue.edu/People/graduate_students.html
ncstate,faculty,https://bio.sciences.ncsu.edu/group/faculty/
ncstate,grad,https://bio.sciences.ncsu.edu/group/graduate-students/
lsu,faculty,https://www.lsu.edu/science/biosci/people/faculty.php
lsu,grad,https://www.lsu.edu/science/biosci/people/graduate-students.php
vt,faculty,https://www.biol.vt.edu/People/Faculty.html
vt,grad,https://www.biol.vt.edu/People/Graduate_Students.html
auburn,faculty,https://www.auburn.edu/cosam/departments/biology/biology-faculty/index.htm
auburn,grad,https://www.auburn.edu/cosam/departments/biology/graduate-student-current.htm
kstate,faculty,https://www.k-state.edu/biology/about/people/faculty/
kstate,grad,https://www.k-state.edu/biology/about/people/graduate/
texastech,faculty,https://www.depts.ttu.edu/biology/people/faculty/
texastech,grad,https://www.depts.ttu.edu/biology/people/graduate-students/
colostate,faculty,https://www.biology.colostate.edu/faculty/
florida,faculty,https://biology.ufl.edu/people/faculty/
', stringsAsFactors = FALSE, strip.white = TRUE)

## Strip tags crudely; base R only. Scripts and styles removed first.
html_text <- function(h) {
  h <- gsub("(?s)<script.*?</script>", " ", h, perl = TRUE)
  h <- gsub("(?s)<style.*?</style>", " ", h, perl = TRUE)
  h <- gsub("<[^>]+>", "\n", h)
  h <- gsub("&amp;", "&", h); h <- gsub("&nbsp;", " ", h)
  h <- gsub("[ \t]+", " ", h); h <- gsub("\n\\s*\n+", "\n", h)
  trimws(h)
}

## Classify each text line by the rank label it carries. Directory pages put
## the rank in its own element, so after tag stripping it is its own line.
## Order matters: emeritus, teaching, research, and adjunct labels are
## checked before the plain tenure-line pattern so "Instructional Assistant
## Professor" is not counted as tenure-line.
classify <- function(lines) {
  l <- trimws(lines)
  is <- function(rx) grepl(rx, l, ignore.case = TRUE, perl = TRUE)
  out <- rep(NA_character_, length(l))
  out[is("Emerit")] <- "emeritus"
  out[is.na(out) & is("Teaching|Instructional|Professor of Practice|Lecturer|Instructor|Professor of the Practice")] <- "teaching_track"
  out[is.na(out) & is("^Research (Assistant |Associate )?Professor|Research Scientist|Research Associate$|Postdoc")] <- "research_track"
  out[is.na(out) & is("Adjunct|Affiliate|Courtesy|Joint|Visiting|Clinical")] <- "adjunct_affiliate"
  out[is.na(out) & is("^(University |Regents |Distinguished |Boyd |Presidential |[A-Z][A-Za-z.'\\- ]{2,40} )?(Assistant |Associate )?Professor\\b(,| and| of |$)")] <- "tenure_line"
  out
}
rank_levels <- c("tenure_line", "teaching_track", "research_track", "adjunct_affiliate", "emeritus")

stamp <- format(Sys.Date(), "%Y%m%d")
out <- lapply(seq_len(nrow(pages)), function(i) {
  p <- pages[i, ]
  f_html <- file.path(dir_dir, sprintf("%s_%s_%s.html", p$key, p$page, stamp))
  ok <- tryCatch({ download_raw(p$url, f_html); TRUE }, error = function(e) FALSE)
  if (!ok) return(data.frame(key = p$key, page = p$page, url = p$url, fetched = FALSE))
  log_pull(f_html, p$url, "department directory", p$page)
  txt <- html_text(paste(readLines(f_html, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  writeLines(txt, sub("\\.html$", ".txt", f_html))
  cls <- classify(strsplit(txt, "\n")[[1]])
  n <- vapply(rank_levels, function(k) sum(cls == k, na.rm = TRUE), integer(1))
  ## Graduate pages: count e-mail addresses or "Ph.D." / "M.S." tags as a
  ## rough headcount; hand-check against the .txt.
  n_email <- length(gregexpr("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.edu", txt)[[1]][gregexpr("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.edu", txt)[[1]] > 0])
  data.frame(key = p$key, page = p$page, url = p$url, fetched = TRUE, n_chars = nchar(txt),
             as.list(n), n_email = n_email)
})
out <- do.call(rbind, lapply(out, function(d) { d[setdiff(c("fetched","n_chars",names(ranks),"n_email"), names(d))] <- NA; d }))
out$tenure_line <- with(out, professor + associate_professor + assistant_professor)
print(out[, c("key", "page", "fetched", "n_chars", "tenure_line", "teaching_track", "research_track", "emeritus", "n_email")])
write.csv(out, file.path(derived_dir, sprintf("directory_counts_%s.csv", stamp)), row.names = FALSE)
message("Pages with n_chars < 2000 or tenure_line == 0 are JavaScript-rendered: open the URL, select all, paste into the .txt, re-run the tally by hand.")
