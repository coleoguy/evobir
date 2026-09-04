## 07_did.R -- difference-in-differences, event study, synthetic control
##
## Model, for each outcome y (log(y+1) for counts and dollars):
##   y_it = a_i + g_t + b * treated_i * post_t + e_it
## estimated by weighted least squares with peers weighted 1/D_i (normalised)
## and TAMU weighted 1, so the comparison group is a distance-weighted
## average of peers rather than a flat mean. Standard errors are clustered by
## institution (CR1 sandwich, computed by hand: no add-on package). With one
## treated unit the cluster count is small and the t-stat should be read
## against the permutation distribution reported alongside, not a normal.
##
## Event study: replace b * treated * post with treated * 1[year = k] for every
## k except the last pre-year (treat_year - 1), the omitted reference.
## Pre-period coefficients near zero are the parallel-trends check.
##
## Synthetic control: weights w >= 0, sum(w) = 1 chosen to minimise pre-period
## squared error of the outcome path (softmax parameterisation solved with
## optim so it stays base R). Inference by in-space placebo: fit every peer as
## if treated and compare TAMU's post/pre RMSPE ratio to the placebo ratios.
##
## This script runs on whatever years exist; before treat_year it only
## produces the pre-trend diagnostics.

source("R/00_config.R")
panel <- read.csv(file.path(derived_dir, "outcome_panel.csv"), stringsAsFactors = FALSE)
outcomes <- intersect(c("phd", "grants3", "pubs", "new_trainees"), names(panel))
fig_dir <- file.path(proj_dir, "figures"); dir.create(fig_dir, showWarnings = FALSE)

## CR1 cluster-robust covariance for an lm fit, clusters = institution.
cluster_vcov <- function(fit, cl) {
  X <- model.matrix(fit)[, !is.na(coef(fit)), drop = FALSE]   # drop aliased columns
  w <- if (is.null(fit$weights)) rep(1, nrow(X)) else fit$weights
  u <- residuals(fit) * w
  Xw <- X * w
  bread <- solve(crossprod(Xw, X))
  meat <- Reduce(`+`, lapply(split(seq_len(nrow(X)), cl), function(i) {
    s <- crossprod(X[i, , drop = FALSE], u[i]); s %*% t(s) }))
  G <- length(unique(cl)); n <- nrow(X); k <- ncol(X)
  adj <- G / (G - 1) * (n - 1) / (n - k)
  adj * bread %*% meat %*% bread
}

## Skip an outcome with too little data to identify anything (no treated
## rows, or fewer than two institutions): the annual re-pull can hit this
## for phd before coordinator replies arrive.
enough <- function(d) nrow(d) > 0 && any(d$treated == 1) && length(unique(d$inst_key)) >= 2

did_one <- function(y) {
  d <- panel[!is.na(panel[[y]]), ]
  if (!enough(d)) return(list(outcome = y, stage = "skipped: insufficient data"))
  d$ly <- log(d[[y]] + 1)
  d$inst <- factor(d$inst_key); d$yr <- factor(d$year)
  if (!any(d$post == 1)) {
    ## Pre-period only: report pre-trend slope difference.
    fit <- lm(ly ~ inst + yr + treated:I(year - min(year)), data = d, weights = weight)
    return(list(outcome = y, stage = "pre-trend only", fit = fit))
  }
  d$did <- d$treated * d$post
  fit <- lm(ly ~ inst + yr + did, data = d, weights = weight)
  V <- cluster_vcov(fit, d$inst)
  b <- coef(fit)["did"]; se <- sqrt(V["did", "did"])
  V_ev <- NULL
  ## Permutation: reassign treatment to each peer in turn.
  perm <- vapply(setdiff(levels(d$inst), tamu_key), function(k) {
    dd <- d; dd$did <- as.integer(dd$inst_key == k) * dd$post
    coef(lm(ly ~ inst + yr + did, data = dd, weights = weight))["did"]
  }, numeric(1))
  p_perm <- mean(abs(perm) >= abs(b))
  ## Event study.
  ref <- treat_year - 1
  ks <- setdiff(sort(unique(d$year)), ref)
  for (k in ks) d[[paste0("ev_", k)]] <- d$treated * (d$year == k)
  ks <- ks[vapply(ks, function(k) any(d[[paste0("ev_", k)]] > 0), logical(1))]  # TAMU must have data in year k
  f_ev <- as.formula(paste("ly ~ inst + yr +", paste(paste0("ev_", ks), collapse = " + ")))
  fit_ev <- lm(f_ev, data = d, weights = weight)
  V_ev <- cluster_vcov(fit_ev, d$inst)
  ev <- data.frame(year = ks, est = coef(fit_ev)[paste0("ev_", ks)],
                   se = sqrt(diag(V_ev)[paste0("ev_", ks)]))
  ev <- rbind(ev, data.frame(year = ref, est = 0, se = 0)); ev <- ev[order(ev$year), ]
  png(file.path(fig_dir, paste0("event_study_", y, ".png")), 900, 600, res = 120)
  plot(ev$year, ev$est, type = "b", pch = 19, ylim = range(c(ev$est - 1.96 * ev$se, ev$est + 1.96 * ev$se)),
       xlab = "Year", ylab = paste0("Effect on log(", y, " + 1)"), main = paste("Event study:", y))
  ok <- ev$se > 0
  arrows(ev$year[ok], (ev$est - 1.96 * ev$se)[ok], ev$year[ok], (ev$est + 1.96 * ev$se)[ok], angle = 90, code = 3, length = 0.03)
  abline(h = 0, lty = 2); abline(v = treat_year - 0.5, col = "grey50")
  dev.off()
  list(outcome = y, stage = "did", est = b, se = se, p_perm = p_perm, event = ev, fit = fit)
}

## Synthetic control for one outcome. Returns weights and the placebo test.
synth_one <- function(y) {
  d <- panel[!is.na(panel[[y]]), ]
  if (!enough(d)) return(list(outcome = y, weights = NA, rmspe_ratio = NA, p_placebo = NA))
  d$ly <- log(d[[y]] + 1)
  Y <- xtabs(ly ~ year + inst_key, data = d)          # years x institutions
  pre <- rownames(Y) %in% as.character(pre_years)
  fit_w <- function(target, donors) {
    y1 <- Y[pre, target]; Y0 <- Y[pre, donors, drop = FALSE]
    obj <- function(theta) { w <- exp(theta) / sum(exp(theta)); sum((y1 - Y0 %*% w)^2) }
    o <- optim(rep(0, length(donors)), obj, method = "BFGS", control = list(maxit = 2000))
    w <- exp(o$par) / sum(exp(o$par)); names(w) <- donors
    gap <- Y[, target] - Y[, donors, drop = FALSE] %*% w
    list(w = w, gap = as.numeric(gap), years = as.integer(rownames(Y)))
  }
  donors <- setdiff(colnames(Y), tamu_key)
  main <- fit_w(tamu_key, donors)
  rmspe_ratio <- function(g, yrs) {
    post <- yrs >= treat_year; if (!any(post)) return(NA)
    sqrt(mean(g[post]^2)) / sqrt(mean(g[!post]^2))
  }
  placebo <- vapply(donors, function(k) rmspe_ratio(fit_w(k, setdiff(donors, k))$gap, main$years), numeric(1))
  r_tamu <- rmspe_ratio(main$gap, main$years)
  png(file.path(fig_dir, paste0("synth_", y, ".png")), 900, 600, res = 120)
  plot(main$years, Y[, tamu_key], type = "b", pch = 19, xlab = "Year", ylab = paste0("log(", y, " + 1)"),
       main = paste("Synthetic control:", y))
  lines(main$years, Y[, tamu_key] - main$gap, lty = 2)
  abline(v = treat_year - 0.5, col = "grey50"); legend("topleft", c("TAMU", "synthetic"), lty = 1:2, bty = "n")
  dev.off()
  list(outcome = y, weights = round(main$w, 3), rmspe_ratio = r_tamu,
       p_placebo = if (is.na(r_tamu)) NA else mean(placebo >= r_tamu, na.rm = TRUE))
}

did <- lapply(outcomes, did_one); names(did) <- outcomes
sc  <- lapply(outcomes, synth_one); names(sc) <- outcomes

summ <- do.call(rbind, lapply(outcomes, function(y) {
  r <- did[[y]]; s <- sc[[y]]
  data.frame(outcome = y, stage = r$stage,
             did_est = if (r$stage == "did") unname(r$est) else NA, did_se = if (r$stage == "did") r$se else NA,
             p_perm = if (r$stage == "did") r$p_perm else NA,
             synth_rmspe_ratio = if (is.null(s$rmspe_ratio)) NA else s$rmspe_ratio,
             synth_p_placebo = if (is.null(s$p_placebo)) NA else s$p_placebo)
}))
print(summ, digits = 3)
write.csv(summ, file.path(derived_dir, sprintf("did_summary_%s.csv", format(Sys.Date(), "%Y%m%d"))), row.names = FALSE)
saveRDS(list(did = did, synth = sc), file.path(derived_dir, "did_fits.rds"))
message("07_did done; figures in figures/")
