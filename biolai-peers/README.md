# BiolAI peer departments and evaluation pipeline

Peer selection and difference-in-differences evaluation of BiolAI for TAMU Biology
(College of Arts and Sciences). See `context.md` for status and decisions.

Run from this directory, in order, with base R plus `jsonlite` and `curl`:

```
Rscript R/01_ipeds.R              # IPEDS completions (CIP 26.01xx) and Pell share
Rscript R/02_proquest.R           # PQDT search strings; parses exports in data/raw/proquest/
Rscript R/03_openalex_faculty.R   # research-active faculty proxy per department
Rscript R/04_grants.R             # NSF + NIH dollars attributed to department by PI
Rscript R/05_match.R              # Mahalanobis distance, filters, writes data/peers.csv
Rscript R/06_outcomes.R           # department-year outcome panel (annual)
Rscript R/07_did.R                # DiD, event study, synthetic control (annual)
```

Every raw pull is written to `data/raw/` with a date stamp and a row in
`data/raw/manifest.csv` (file, source URL, timestamp, bytes). Derived tables go to
`data/derived/`. Internal rosters and email replies go to `data/internal/` (not
committed; templates provided). Set `BIOLAI_MAILTO` in `~/.Renviron` for the
OpenAlex polite pool.

Two inputs are manual: ProQuest exports (no API; strings are printed by 02) and NSF
HERD biological-sciences expenditures (download the NCSES table per fiscal year to
`data/raw/herd/herd_bio_<FY>.csv`).

`data/seed_pool.csv` holds the candidate roster and every structural fact used by
the filters, with a `verification` column saying how each row was checked.
