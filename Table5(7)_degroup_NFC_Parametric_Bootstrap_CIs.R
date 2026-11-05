# DEGROUPING: deterministic interior spacing of rounded claim amounts
# ============================================================================
# Apply this preprocessing once, before fitting either model. Repeated amounts
# are treated as rounded observations within a bin of width rounding_unit NOK.
# For m observations tied at v, use L = max(d, v - rounding_unit/2) and
# U = v + rounding_unit/2, then replace the group by the m interior points
# L + j*(U-L)/(m+1), j = 1,...,m. Endpoints are excluded deliberately.
# Example: three ties at 601000 with rounding_unit = 1000 become
# 600750, 601000, 601250. Clipping L at d keeps all amounts above truncation.
# Singleton amounts are unchanged; the sample size is preserved. This rule
# assigns plausible distinct values; it does not recover unobserved exact claims.
# The sorted original amounts are retained in k_grouped for tie diagnostics.
# Run-length encoding locates each tied group; k receives the replacement points.
# Assertions check spacing within each bin and the final sample's validity.
# The processed sample is exported to NFC_Degrouped_Data.csv for inspection.
# No seed is needed for this deterministic step. The bootstrap seed controls
# simulated samples only. Bootstrap intervals condition on this processed data;
# they do not account for uncertainty about the original rounding or this rule.
# Continuous bootstrap samples are fitted directly without further degrouping.
#
# Norwegian Fire Claims: parametric bootstrap percentile confidence intervals
# Pareto I and shifted lognormal, both with FIXED ground-up x0 = 100000 NOK.
# Observations are LEFT TRUNCATED at d = 500000 NOK.
# Run via RStudio Source or Rscript. Outputs are written beside this R file
# when its path is detectable; otherwise use the working directory.
# Input: NFC_Data.xlsx in that same folder, with numeric claims and no header.
# Install once: install.packages('readxl'). rstudioapi is optional.
# Monetary calculations use NOK; risk tables convert results to million NOK.

# Configuration: fixed support threshold x0 and observation threshold d.
# Neither threshold is estimated during fitting or bootstrap refitting.
x0 <- 100000
d <- 500000
# beta is upper-tail probability: VaR here is the (1-beta) quantile.
# delta controls GS weighting, r the proportional-hazards transform (PHT),
# and lambda the Wang transform (WT). B sets replications per model.
beta <- 0.10
delta <- 0.25
r <- 0.95
lambda <- 0.25
B <- 5000L
# Fix the random-number sequence for reproducible bootstrap simulation.
seed <- 20261021
# Central confidence levels; each uses equal lower and upper tail probabilities.
levels <- c(0.90, 0.95, 0.99)

rounding_unit <- 1000

# Check core configuration assumptions before reading or fitting any data.
stopifnot(
  x0 > 0, x0 < d, B >= 99,
  beta > 0, beta < 1,
  r > 0, r <= 1,
  rounding_unit > 0
)

# Resolve output folder robustly across Rscript and RStudio Source.
script_path <- NULL
# Rscript supplies --file=; take its first occurrence when present.
args <- grep('^--file=', commandArgs(FALSE), value=TRUE)
if (length(args)) script_path <- sub('^--file=', '', args[1])

# If running in RStudio, try the active source editor path. Catch API errors
# so path detection can continue in other execution environments.
if (is.null(script_path) &&
    requireNamespace('rstudioapi', quietly=TRUE)) {
  script_path <- tryCatch(
    rstudioapi::getSourceEditorContext()$path,
    error=function(e) NULL
  )
}

# As a further fallback, inspect source() frames from newest to oldest
# for an existing source filename stored in ofile.
if (is.null(script_path) ||
    length(script_path) != 1L ||
    is.na(script_path) ||
    !nzchar(script_path) ||
    !file.exists(script_path)) {
  for (f in rev(sys.frames())) {
    if (!is.null(f$ofile) && file.exists(f$ofile)) {
      script_path <- f$ofile
      break
    }
  }
}

# Use the containing folder of a valid script path, or the current working
# directory when no usable path was found.
folder <- if (
  !is.null(script_path) &&
  length(script_path) == 1L &&
  !is.na(script_path) &&
  nzchar(script_path) &&
  file.exists(script_path)
) {
  dirname(normalizePath(script_path))
} else {
  getwd()
}

# Stop early if the input workbook or required reader is unavailable.
input <- file.path(folder, 'NFC_Data.xlsx')
if (!file.exists(input)) {
  stop('NFC_Data.xlsx not found in: ', folder)
}
if (!requireNamespace('readxl', quietly=TRUE)) {
  stop('Install readxl: install.packages("readxl")')
}

# Flatten worksheet cells into one numeric vector, discard non-finite values,
# and sort. Non-numeric cells may produce coercion warnings and be discarded.
k <- as.numeric(unlist(readxl::read_excel(input, col_names=FALSE)))
k <- sort(k[is.finite(k)])
n <- length(k)

# Require at least ten observations and no claims below d; flag deviations
# from the expected sample size without preventing execution.
if (n < 10 || any(k < d)) stop('Check claim data and truncation d.')
if (n != 647) warning('Expected 647 claims, found ', n)

k_grouped <- k
tied <- duplicated(k_grouped) |
  duplicated(k_grouped, fromLast=TRUE)

groups <- rle(k_grouped)
group_end <- cumsum(groups$lengths)
group_start <- group_end - groups$lengths + 1L
tied_groups <- which(groups$lengths > 1L)

for (g in tied_groups) {
  m <- groups$lengths[g]
  value <- groups$values[g]
  
  lower <- max(d, value - rounding_unit/2)
  upper <- value + rounding_unit/2
  
  idx <- seq.int(group_start[g], group_end[g])
  interior_points <- lower +
    (seq_len(m)/(m + 1))*(upper - lower)
  
  stopifnot(
    all(interior_points > lower),
    all(interior_points < upper),
    all(diff(interior_points) > 0)
  )
  
  k[idx] <- interior_points
}

k <- sort(k)

degrouped_out <- file.path(folder, 'NFC_Degrouped_Data.csv')

write.csv(
  data.frame(
    Observation=seq_along(k),
    Claim_NOK=k
  ),
  file=degrouped_out,
  row.names=FALSE
)

stopifnot(length(k) == n, all(is.finite(k)), all(k >= d))
cat(
  'Interior-spacing de-grouping applied to', sum(tied),
  'observations in', length(tied_groups), 'tied groups.\n'
)

# Model 1: conditional Pareto likelihood. With d fixed, its shape MLE is
# n / sum(log(x/d)); x0 cancels from the conditional density above d.
# Truncated Pareto I: f(x|X>=d) = alpha*d^alpha/x^(alpha+1).
fit_pa <- function(x) length(x)/sum(log(x/d))

# Shifted LN: X=x0+exp(Y), Y~N(mu,sigma^2), conditioned on X>=d.
# Model 2: maximize the shifted-lognormal likelihood conditional on X >= d.
# Optimize (mu, log(sigma)) so sigma remains positive automatically.
fit_ln <- function(x, start=NULL) {
  # Transform claims onto the normal scale; yd is its truncation boundary.
  y <- log(x-x0)
  yd <- log(d-x0)
  nn <- length(x)
  
  # Negative conditional log-likelihood: subtract normal log densities and
  # add n*log(P(Y >= yd)). The x-to-y Jacobian is parameter-independent
  # for fixed x0, so it can be omitted for this optimization.
  nll <- function(par) {
    mu <- par[1]
    s <- exp(par[2])
    # Use log survival probabilities to avoid loss of precision in the tail.
    z <- (yd-mu)/s
    val <- -sum(dnorm(y, mu, s, log=TRUE)) +
      nn*pnorm(z, lower.tail=FALSE, log.p=TRUE)
    # Return a large penalty for numerical failures encountered by optim().
    if (!is.finite(val)) 1e100 else val
  }
  
  # Use transformed sample moments only as initial guesses, not final fits.
  if (is.null(start)) {
    start <- c(mean(y), log(max(sd(y), 0.01)))
  }
  
  # Run two optimization methods from the same start and retain the smaller
  # objective. Selection checks numerical validity below, not convergence codes.
  candidates <- list(
    optim(start, nll, method='BFGS', control=list(maxit=1500)),
    optim(start, nll, method='Nelder-Mead', control=list(maxit=2000))
  )
  
  # Extract scalar objective values and choose the best candidate.
  vals <- vapply(candidates, function(z) z$value, numeric(1))
  best <- candidates[[which.min(vals)]]
  
  # Reject unusable estimates before downstream risk and simulation calculations.
  if (!is.finite(best$value) ||
      any(!is.finite(best$par)) ||
      !is.finite(exp(best$par[2])) ||
      exp(best$par[2]) <= 0) {
    stop('Shifted LN fitting failed or produced invalid parameters.')
  }
  
  # Return consistently named parameters on their natural scales.
  c(
    mu=unname(best$par[1]),
    sigma=unname(exp(best$par[2]))
  )
}

# Ground-up risk functionals, not conditional risk measures above d.
# Quantile-based integrals evaluated with midpoint quadrature on v=-log(t).
# t denotes survival probability. WT and GS are integrated numerically.
# For survival probability t = exp(-v), dt contributes exp(-v) dv.
# Midpoint sums approximate the integrals up to finite cutoff vmax with step h;
# these are numerical approximations, with no adaptive error assessment.
vmax <- 180
h <- 0.01
v_all <- seq(h/2, vmax-h/2, by=h)
# Tail integrals start at v = -log(beta), corresponding to t = beta.
v_tail <- seq(-log(beta)+h/2, vmax-h/2, by=h)

# Lognormal survival quantile Q(t)=x0+exp(mu+sigma*z_t).
q_ln <- function(t, mu, sigma) {
  x0+exp(mu+sigma*qnorm(t, lower.tail=FALSE))
}
# Ground-up Pareto survival quantile (defined here for reference).
q_pa <- function(t, a) x0*exp(-log(t)/a)

# Return ground-up VaR, conditional tail expectation (CTE), GS,
# proportional-hazards transform (PHT), and Wang transform (WT).
# GS uses the explicitly coded tail weight controlled by delta.
risk_ln <- function(mu, sigma) {
  # Remove inherited names so output names remain exactly as specified.
  mu <- unname(as.numeric(mu))
  sigma <- unname(as.numeric(sigma))
  
  # Require finite scalar location and positive finite scalar scale.
  if (length(mu) != 1L ||
      length(sigma) != 1L ||
      !is.finite(mu) ||
      !is.finite(sigma) ||
      sigma <= 0) {
    stop('Invalid shifted LN parameters')
  }
  
  # Closed-form VaR and CTE at upper-tail probability beta. The additive
  # shift x0 contributes once to each monetary risk measure.
  zb <- qnorm(beta, lower.tail=FALSE)
  VaR <- x0+exp(mu+sigma*zb)
  CTE <- x0+exp(mu+sigma^2/2)*pnorm(sigma-zb)/beta
  
  # GS: integrate the survival quantile over t in (0,beta), with linear
  # weight 1 + 2*delta - 4*delta*t/beta, then divide by beta.
  tt <- exp(-v_tail)
  qq <- q_ln(tt, mu, sigma)
  GS <- sum(qq*(1+2*delta-4*delta*tt/beta)*tt)*h/beta
  
  # Evaluate the survival quantile on the full integration grid.
  t <- exp(-v_all)
  q <- q_ln(t, mu, sigma)
  # PHT uses density weight r*t^(r-1); after transformation its weight is
  # r*exp(-r*v), including the change-of-variable factor.
  PHT <- sum(q*r*exp(-r*v_all))*h
  
  z <- qnorm(t, lower.tail=FALSE)
  # WT multiplies the transformed quantile by exp(lambda*z-lambda^2/2),
  # where z is the standard-normal survival quantile.
  WT <- sum(q*exp(-v_all+lambda*z-lambda^2/2))*h
  
  # Strip intermediate names and assemble the agreed risk-measure columns.
  ans <- c(
    VaR=unname(VaR),
    CTE=unname(CTE),
    GS=unname(GS),
    PHT=unname(PHT),
    WT=unname(WT)
  )
  
  # Treat overflow or other non-finite LN results as calculation failures.
  if (any(!is.finite(ans))) {
    stop('Non-finite shifted LN risk estimate; check fit')
  }
  ans
}

# Pareto risks use closed forms where available. Infinite values represent
# theoretical divergence when the required tail moment does not exist.
risk_pa <- function(a) {
  # CTE and GS require a > 1; PHT requires r*a > 1.
  VaR <- x0*beta^(-1/a)
  CTE <- if (a > 1) VaR*a/(a-1) else Inf
  GS <- if (a > 1) CTE*(1+2*delta/(2*a-1)) else Inf
  PHT <- if (r*a > 1) x0*r*a/(r*a-1) else Inf
  
  # WT finite only for a>1 when lambda>0.
  if (a <= 1) {
    WT <- Inf
  } else {
    # Numerically integrate Pareto WT on the same transformed grid.
    t <- exp(-v_all)
    z <- qnorm(t, lower.tail=FALSE)
    WT <- sum(
      x0*exp(-v_all*(1-1/a)+lambda*z-lambda^2/2)
    )*h
  }
  
  c(VaR=VaR, CTE=CTE, GS=GS, PHT=PHT, WT=WT)
}

# Simulate from fitted conditional distributions, respecting truncation.
# Inverse-survival sampling for the conditional Pareto distribution: its
# threshold is d because only claims exceeding d are observed.
sim_pa <- function(nn, a) d*exp(-log(runif(nn))/a)

# Sample shifted LN conditional on X >= d by drawing survival probabilities
# uniformly between zero and P(X >= d), then applying the survival quantile.
sim_ln <- function(nn, mu, sigma) {
  z_d <- (log(d-x0)-mu)/sigma
  tail_d <- pnorm(z_d, lower.tail=FALSE)
  
  # Fail explicitly if the conditional probability is lost to underflow.
  if (tail_d <= 0) {
    stop('LN conditional survival probability underflow')
  }
  
  z <- qnorm(runif(nn)*tail_d, lower.tail=FALSE)
  x0+exp(mu+sigma*z)
}

# Fit both models to k and calculate their ground-up point risk estimates.
alpha_hat <- fit_pa(k)
ln_hat <- fit_ln(k)
pa_est <- risk_pa(alpha_hat)
ln_est <- risk_ln(unname(ln_hat['mu']), unname(ln_hat['sigma']))

# Verify the risk vector has the exact names used for bootstrap table lookup.
stopifnot(
  identical(names(ln_est), c('VaR', 'CTE', 'GS', 'PHT', 'WT'))
)

# Initialize bootstrap randomness once, so replications form a reproducible
# sequence rather than restarting the generator inside the loop.
set.seed(seed)
measures <- names(pa_est)

# Preallocate B-by-measure result matrices. Missing LN results remain NA;
# Inf Pareto risks are retained as meaningful outputs.
boot_pa <- matrix(
  NA_real_, B, length(measures),
  dimnames=list(NULL, measures)
)
boot_ln <- matrix(
  NA_real_, B, length(measures),
  dimnames=list(NULL, measures)
)

# Also retain refitted parameters for diagnostics; failed_ln counts failed
# LN fits or risk evaluations, at most one increment per replication.
par_pa <- numeric(B)
par_ln <- matrix(
  NA_real_, B, 2,
  dimnames=list(NULL, c('mu', 'sigma'))
)
failed_ln <- 0L

# Each model gets its own simulated sample of size n from its fitted
# conditional distribution, followed by refitting and risk recalculation.
for (b in seq_len(B)) {
  # Pareto branch: simulate, re-estimate shape, then recompute all risks.
  xp <- sim_pa(n, alpha_hat)
  par_pa[b] <- fit_pa(xp)
  boot_pa[b, ] <- risk_pa(par_pa[b])
  
  # LN branch: warm-start optimization at the original fitted parameters;
  # the second optimizer coordinate must be log(sigma).
  xl <- sim_ln(n, ln_hat['mu'], ln_hat['sigma'])
  fitted <- tryCatch(
    fit_ln(xl, c(ln_hat['mu'], log(ln_hat['sigma']))),
    error=function(e) NULL
  )
  
  # A failed fit leaves this replication missing. Otherwise store parameters
  # and attempt risk evaluation independently of optimization success.
  if (is.null(fitted)) {
    failed_ln <- failed_ln+1L
  } else {
    stopifnot(
      identical(names(fitted), c('mu', 'sigma')),
      all(is.finite(fitted))
    )
    par_ln[b, ] <- fitted
    boot_ln[b, ] <- tryCatch(
      risk_ln(unname(fitted['mu']), unname(fitted['sigma'])),
      error=function(e) {
        # Update the enclosing failure counter from inside the error handler.
        failed_ln <<- failed_ln+1L
        rep(NA_real_, length(measures))
      }
    )
  }
  
  # Print progress every 100 replications for this potentially long computation.
  if (b %% 100 == 0) {
    cat('Completed', b, '/', B, 'replications\n')
  }
}

# Percentile intervals: retain Inf for divergent Pareto risks;
# exclude only NA/NaN numerical failures. Group tables by CI level.
# Build one row per confidence level and risk measure, preserving that order.
# These are percentile intervals from refitted risks, not parameter intervals.
make_table <- function(model, point, boots) {
  do.call(rbind, lapply(levels, function(level) {
    do.call(rbind, lapply(measures, function(m) {
      vals <- boots[, m]
      # Drop NA and NaN only; retain both finite and infinite numerical results.
      vals <- vals[!is.na(vals)]
      
      # Warn if fewer than 90% of requested replications are usable for this row.
      if (length(vals) < ceiling(0.90*B)) {
        warning('Few valid bootstrap results: ', model, ' ', m)
      }
      
      # Central endpoints are empirical quantiles at (1-level)/2 and (1+level)/2.
      # R type 7 uses interpolation; infinite values can affect interpolated endpoints.
      q <- if (length(vals)) {
        quantile(
          vals,
          c((1-level)/2, (1+level)/2),
          na.rm=TRUE,
          names=FALSE,
          type=7
        )
      } else {
        c(NA_real_, NA_real_)
      }
      
      # Convert point estimates and interval endpoints from NOK to million NOK;
      # include usable and infinite replication counts so diagnostics stay visible.
      data.frame(
        Model=model,
        CI_Level=paste0(round(level*100), '%'),
        Risk_Measure=m,
        Estimate_million_NOK=unname(point[m])/1e6,
        CI_Lower_million_NOK=q[1]/1e6,
        CI_Upper_million_NOK=q[2]/1e6,
        Valid_Bootstrap=length(vals),
        Infinite_Bootstrap=sum(is.infinite(vals)),
        stringsAsFactors=FALSE
      )
    }))
  }))
}

# Create model-specific tables using the same confidence levels and measures.
pa_tab <- make_table('Pareto I', pa_est, boot_pa)
ln_tab <- make_table('Shifted Lognormal', ln_est, boot_ln)

# One CSV with separate labeled sections ordered by 90%, 95%, 99%.
out <- file.path(folder, 'NFC_Parametric_Bootstrap_Risk_CIs.csv')
# Open one text connection for both labeled CSV sections. This is a report
# with titles and a blank separator, rather than one rectangular CSV table.
con <- file(out, 'wt')
writeLines(
  'Table 1. Pareto I - parametric bootstrap percentile CIs', con
)
write.table(
  pa_tab, con, sep=',', row.names=FALSE,
  col.names=TRUE, quote=TRUE, na=''
)
writeLines(
  c('', 'Table 2. Shifted Lognormal - parametric bootstrap percentile CIs'),
  con
)
write.table(
  ln_tab, con, sep=',', row.names=FALSE,
  col.names=TRUE, quote=TRUE, na=''
)
# Flush and close the combined confidence-interval output.
close(con)

# Record sample size, fixed thresholds, and point-fit parameters. NA marks
# parameters that are not part of the corresponding model.
params <- data.frame(
  Model=c('Pareto I', 'Shifted Lognormal'),
  n=n,
  x0=x0,
  d=d,
  Alpha=c(alpha_hat, NA),
  Mu=c(NA, unname(ln_hat['mu'])),
  Sigma=c(NA, unname(ln_hat['sigma']))
)

# Write a separate, rectangular CSV of fitted parameters.
params_out <- file.path(
  folder, 'NFC_Parametric_Bootstrap_Fitted_Parameters.csv'
)
write.csv(params, params_out, row.names=FALSE)

# Redirect console-style output into a text report with settings, diagnostic
# counts, methodological notes, and confidence-level-specific tables.
report <- file.path(folder, 'NFC_Parametric_Bootstrap_CI_Report.txt')
sink(report)
cat('NFC parametric bootstrap percentile confidence intervals\n')
cat(
  'Ground-up x0=', x0,
  '; left truncation d=', d,
  '; n=', n,
  '; B=', B,
  '; bootstrap seed=', seed, '\n',
  sep=''
)
cat(
  'Risk settings: beta=', beta,
  ', delta=', delta,
  ', r=', r,
  ', lambda=', lambda, '\n',
  sep=''
)
cat(
  'Pareto alpha=', alpha_hat,
  '; LN mu=', ln_hat['mu'],
  '; LN sigma=', ln_hat['sigma'], '\n'
)
cat('Failed LN bootstrap fits:', failed_ln, '\n')
cat(
  'Interior-spacing de-grouping: rounding unit=', rounding_unit,
  ' NOK; observations in tied groups=', sum(tied),
  '; tied groups=', length(tied_groups), '\n',
  sep=''
)
cat('Rule: L + j*(U-L)/(m+1), j=1,...,m, for each group of m ties.\n')
cat('Only tied groups were processed; singleton values were unchanged.\n')
cat('Rounding intervals were clipped at d before interior spacing.\n')
cat('Degrouping was deterministic; no degrouping seed was required.\n')
cat('Intervals condition on this interior-spaced dataset.\n')
cat('Risk measures are ground-up (untruncated); bootstrap sampling is conditional on X>=d.\n')
cat('Inf denotes theoretical divergence; NA denotes missing/failed calculations.\n')
cat('Intervals are percentile bootstrap intervals, not bias-corrected or studentized.\n')

# Print both model tables for each requested confidence level.
for (level in paste0(round(levels*100), '%')) {
  cat('\n============ ', level, ' ============\n', sep='')
  cat('\nPareto I\n')
  print(pa_tab[pa_tab$CI_Level == level, ], row.names=FALSE)
  cat('\nShifted Lognormal\n')
  print(ln_tab[ln_tab$CI_Level == level, ], row.names=FALSE)
}
# Restore ordinary console output after completing the report.
sink()

# List all output paths in the console for easy retrieval.
cat(
  'Saved:\n',
  degrouped_out, '\n',
  out, '\n',
  params_out, '\n',
  report, '\n'
)