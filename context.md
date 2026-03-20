# evobiR Package — Working Context

## Overview
R package for comparative analyses and teaching evolutionary biology. Maintained by Heath Blackmon.

## Recent Changes

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

## Package Structure
- Source: `R/` directory
- Tests: `tests/testthat/`
- Vignette: `vignettes/evolutionary.biology.in.R.Rmd`
- Documentation generated via roxygen2

## Key Exported Functions
AncCond, AncCondFast, CalcD, CalcPopD, FuzzyMatch, GetTipRates, pSAF, ReOrderAlignment, ResSel, SampleTrees, SlidingWindow, SuperMatrix, WinCalcD, countTrees, fix.simmap, getNe, make.simmap2, scaleTreeRates
