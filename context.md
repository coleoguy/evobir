# evobiR Package — Working Context

## Overview
R package for comparative analyses and teaching evolutionary biology. Maintained by Heath Blackmon.

## Recent Changes

### 2026-09-01: AncCond missing-data support
- Student hit the "All tips must have non-missing values" error with x and y measured on different species. No prior version allowed missing data.
- `AncCond()` now accepts x and y on any tip subsets. `.pruning2` marginalises NA tips; `.map.pruned` replaced by `.anc.subset` (pruned-tree fastAnc, mapped to full tree with BM-correct Brownian-bridge interpolation; verified against a GLS conditional-expectation oracle in tests).
- Null sims masked to the observed y pattern. Type-I error checked at nominal with disjoint subsets (120 reps). Power with disjoint subsets is much lower than with complete data (signal only via phylogenetic signal in x); warn students.
- New file `tests/testthat/test-AncCond.R`; man page updated by hand (AncCond.Rd is not roxygen-generated).
- Not yet done: version bump / NEWS heading for a release; vignette still references `AncCondFast`, which no longer exists in `R/`.

### 2026-03-20: Replaced `Prop_SA` with `pSAF`
- **Old function**: `Prop_SA(Da, scs, mud, type, wA, wX, wY)` — returned a `PropSA` S3 object with `print.PropSA` method.
- **New function**: `pSAF(n_auto, n_x, n_y, n_z, n_w, mu, weighted)` — returns a plain numeric (or numeric vector when `n_auto` is a vector).
- Key API changes:
  - `Da` (diploid autosome count) → `n_auto` (haploid autosome count, i.e., number of pairs)
  - `scs` (string like "XY", "XXO") → explicit `n_x`, `n_y`, `n_z`, `n_w` integer counts
  - `mud` → `mu` (now proportion from **females**, not homogametic sex)
  - `type = "fixation"/"occurrence"` → `weighted = TRUE/FALSE`
  - Removed user-facing `wA`, `wX`, `wY` parameters (hardcoded internally)
  - Removed `PropSA` S3 class and `print.PropSA` method
- Files updated: NAMESPACE, vignette, tests, man pages (old .Rd files deleted, new pSAF.Rd generated via roxygen2)

### 2026-03-20: Added missing functions to vignette
- Added vignette sections for: AncCondFast, GetTipRates, countTrees, WinCalcD, ReOrderAlignment
- Vignette now covers all 18 exported functions

## Package Structure
- Source: `R/` directory
- Tests: `tests/testthat/`
- Vignette: `vignettes/evolutionary.biology.in.R.Rmd`
- Documentation generated via roxygen2

## Key Exported Functions
AncCond, AncCondFast, CalcD, CalcPopD, FuzzyMatch, GetTipRates, pSAF, ReOrderAlignment, ResSel, SampleTrees, SlidingWindow, SuperMatrix, WinCalcD, countTrees, fix.simmap, getNe, make.simmap2, scaleTreeRates
