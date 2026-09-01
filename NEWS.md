# evobiR 2.2.1

## New features

* `AncCond()` accepts missing data: `x` and `y` may each be observed on
  any subset of tips, including disjoint subsets (previously any missing
  value was an error). Tips without `y` are marginalised in the Mk
  likelihood and the stochastic maps, and null simulations are masked to
  the observed missingness pattern. Nodes without a direct continuous
  estimate receive the Brownian-motion conditional expectation given the
  `x`-bearing tips (branch-length interpolation along the pruned-tree
  edge; attachment-point value for dataless subtrees). The result gains an
  `n.obs` element. Note for `prune = TRUE`: nodes in state-1 subtrees now
  take this conditional expectation instead of the previous fallback,
  which used a `fastAnc` fit to all tips including the state-1 lineages
  that `prune` is meant to exclude.

## Bug fixes

* `make.simmap2()` now indexes the transition matrix by character state
  name everywhere instead of treating state labels as numeric row
  indices. Traits coded 0/1 previously failed with "no positive
  probabilities"; states with labels like chromosome numbers ("5",
  "12", ...) could silently index the wrong rows. User-supplied Q
  matrices without dimnames are now given state names automatically.
* `make.simmap2()` now honors `rejmax`, `rejint`, and `monitor` for
  `multiPhylo` input (it previously fell back to
  `phytools::make.simmap`, silently dropping these arguments).
* `countTrees()` now works with single-tree collection or reference
  files (previously indexed tree *components* on `phylo` inputs).
* `SampleTrees()` burn-in no longer retains the last burn-in tree,
  no longer produces fractional/zero indices, validates `burnin` and
  `final.number`, and writes correctly named output files (previously
  "prefix .nwk" with a space).
* `scaleTreeRates()` returns class `c("phyloscaled", "phylo")` so that
  `plot()` dispatches to `plot.phyloscaled()` (scalars were previously
  silently ignored in plots); tie-breaking in the final-likelihood
  selection now uses `which.max()`.
* `GetTipRates()` stops with an informative error naming missing tips
  when tip labels cannot be matched (previously crashed with
  `if (NA)` or proceeded on misordered data after a printed note).
* `ResSel()` now selects the same number of individuals in both tails
  (the low tail previously selected one fewer); plotting no longer
  assumes the two trait columns are adjacent and ascending; an
  unsupported `model` now errors informatively.
* `WinCalcD()` returns properly typed (numeric) columns instead of
  all-character ones and no longer drops the final window.
* `CalcD()`/`CalcPopD()` jackknife now drops exactly `block.size`
  contiguous sites (previously `block.size + 1`) at whole-number
  positions (previously fractional).
* `fix.simmap()` correctly resolves start states for failed tip
  branches whose parent is the root, no longer leaks the working
  table across simulations, and the shared-leading-edge test no
  longer always fires (self-match).
* `SuperMatrix()` ignores its own previous output files on re-runs,
  treats `input` as a literal string rather than a regular expression,
  and errors when no files match.
* `getNe()` errors informatively on an unrecognized `locus` value
  (previously "object 'ne' not found").

## Performance improvements

* `CalcD()`/`WinCalcD()`: the per-site biallelic test is now fully
  vectorized (~20x faster; the gain multiplies under bootstrap
  replication).
* `scaleTreeRates()`: matrix exponentials are memoized by scaled edge
  length once Q is fixed, removing >95% of `expm()` calls.
* `AncCond()`: per-edge transition matrices are precomputed once and
  reused across all null simulations.
* `CalcPopD()`: population row indices are computed once per dataset
  rather than roughly ten times per site.

## Documentation and metadata

* `AncCond()` documentation now matches the code defaults
  (`n.sims = 1000`, `Q = NULL`).
* Non-ASCII characters in R sources replaced with ASCII equivalents.
* `rmarkdown` added to Suggests (vignette builder requirement);
  unused imports removed from NAMESPACE; `stats::rbinom` imported;
  `print.anccond` registered as an S3 method; `Mode()` exported.
* Package `Version` bumped to 2.2.1.
