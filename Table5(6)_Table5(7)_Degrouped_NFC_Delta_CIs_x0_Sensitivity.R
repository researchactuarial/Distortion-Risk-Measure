# ---- How degrouping is applied ----
# Repeated claim amounts (ties) may arise because claims were rounded/grouped.
# Before fitting either model, this script deterministically spreads each tied
# group within a rounding interval; it does not add random jitter.
# The default interval width is 1000 NOK (see degrouping_width below).
# For m identical claims at v, define the interval endpoints as
#   L = max(500000, v - 500), U = v + 500  (with the default width).
# More generally, the half-width is degrouping_width/2 and the lower bound
# is the truncation setting, which defaults to 500000 NOK.
# Replace the m ties by the m equally spaced INTERIOR points
#   L + j*(U-L)/(m+1), for j = 1, ..., m.
# Thus neither endpoint is used. For three ties at 601000 NOK, the interval
# is [600500, 601500] and the replacements are 600750, 601000, 601250 NOK.
# At the observation threshold, the lower interval endpoint is clipped:
# three ties at 500000 use [500000, 500500] and become
# 500125, 500250, 500375 NOK. Clipping can shift a tied group's mean upward.
# Singleton amounts are left unchanged, including a singleton at 500000.
# The number of claims is preserved, and the result is sorted again.
# All model fits, confidence intervals, x0 sensitivity calculations and QQ
# plots use these same degrouped observations. Values depend on the chosen
# rounding width; degrouping does not recover the original unrounded claims.

# Norwegian Fire Claims (1986): combined Pareto I and shifted-lognormal analysis
# Reproduces the STRUCTURE of Table 5.3 and Figure 5.1 in the manuscript.
# Run with RStudio Source or: Rscript Norwegian_Fire_Claims_Three_CIs_X0_Sensitivity.R
# Install once if needed: install.packages('readxl')
# Six labeled tables in one CSV and optional QQ PDF saved beside this script.
# Paper: x0=100000 NOK, d=500000 NOK, n=647, beta=.10, delta=.25,
# PHT r=.95, WT lambda=.25. No right censoring (u=Inf).

# ---- Editable settings ----
# x0 is the assumed ground-up Pareto scale / lognormal location, in NOK.
# d (truncation) is the reporting threshold, not a fitted location parameter.
x0 <- 100000
truncation <- 500000
# beta is the upper-tail probability: VaR therefore uses percentile 1-beta.
# delta controls GS tail weighting; r and lambda control PHT and Wang distortion.
beta <- 0.10
delta <- 0.25
r <- 0.95  # Paper specification: use identically for LN and Pareto I
stopifnot(isTRUE(all.equal(r, 0.95)))
lambda <- 0.25
# Two-sided normal delta-method intervals will be calculated at these levels.
confidence_levels <- c(0.90, 0.95, 0.99)
# Sensitivity to the assumed ground-up location x0 (NOK), not truncation d.
x0_sensitivity <- c(10000, 50000, 100000, 150000, 200000)
# Set FALSE to suppress PDF creation; numerical tables are still produced.
make_qq_plot <- TRUE
degrouping_width <- 1000  # NOK: rounding interval centered on each tied value

# ---- Find this script's folder ----
# Resolve output/input paths independently of the current working directory.
# Try the active RStudio editor first, source() call frames next, then Rscript.
script_path <- NULL
if (requireNamespace('rstudioapi', quietly=TRUE)) {
  script_path <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                          error=function(e) NULL)
  if (length(script_path)!=1L || !nzchar(script_path)) script_path <- NULL
}
if (is.null(script_path)) {
  # Search newest source frames first for the saved script filename.
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      script_path <- frames[[i]]$ofile; break
    }
  }
}
if (is.null(script_path)) {
  # Rscript exposes its script path through the --file= command-line argument.
  arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  if (length(arg)) script_path <- sub('^--file=', '', arg[1])
}
if (is.null(script_path) || !nzchar(script_path))
  stop('Cannot locate script. Use RStudio Source or Rscript with a saved .R file.')
# All generated files are written beside this script and the input workbook.
folder <- dirname(normalizePath(script_path, mustWork=TRUE))
input <- file.path(folder, 'NFC_Data.xlsx')
if (!file.exists(input)) stop('Put NFC_Data.xlsx in the same folder as this R script.')
if (!requireNamespace('readxl', quietly=TRUE)) stop('Install readxl: install.packages("readxl")')
# Read every workbook cell as a claim candidate; remove missing/nonfinite values.
# The workbook must contain claim amounts in NOK, without unrelated numeric cells.
k <- as.numeric(unlist(readxl::read_excel(input, col_names=FALSE)))
k <- sort(k[is.finite(k)])
# n is the number of retained claims; sorting also prepares the empirical QQ axis.
n <- length(k)
if (n < 10L || any(k < truncation)) stop('Check data units and truncation threshold.')
if (n != 647L) warning('Paper uses 647 claims; loaded ', n)

# ---- Deterministic degrouping: interior-spacing convention ----
degroup_interior <- function(values, width=1000, lower=500000) {
  # Validate the input vector and scalar interval settings before replacement.
  stopifnot(is.numeric(values), all(is.finite(values)),
            length(width)==1L, is.finite(width), width>0,
            length(lower)==1L, is.finite(lower), all(values>=lower))
  if (!length(values)) return(values)
  values <- sort(values)
  # Runs in the sorted vector identify tied groups; cumulative lengths locate them.
  groups <- rle(values)
  ends <- cumsum(groups$lengths)
  # Visit only repeated values; singleton observations bypass the loop.
  for (g in which(groups$lengths > 1L)) {
    m <- groups$lengths[g]
    v <- groups$values[g]
    L <- max(lower, v-width/2)
    U <- v+width/2
    # Locate this group in the original run structure and assign its interior grid.
    indices <- seq.int(ends[g]-m+1L, ends[g])
    values[indices] <- L + seq_len(m)*(U-L)/(m+1)
  }
  # Re-sort because intervals from nearby groups can overlap.
  sort(values)
}
# Use the same degrouped observations for all fits, CIs, sensitivity and QQ plots.
k <- degroup_interior(k, width=degrouping_width, lower=truncation)

# ---- Fit truncated Pareto and shifted lognormal ----
# Conditioning Pareto I on X>=d removes x0 from the shape likelihood.
# This closed-form MLE uses the observed excess log-ratios above d.
alpha_hat <- n / sum(log(k/truncation))
# Shifted lognormal: X=x0+exp(Y), Y~N(mu,sigma^2), conditioned on X>=d.
# Transform to a normal variate; yd is its left-truncation boundary.
y <- log(k-x0)
yd <- log(truncation-x0)
# Optimize in (mu, log(sigma)) so sigma remains strictly positive.
# The normal log-survival term corrects for conditioning on X>=d.
# The omitted log(k-x0) Jacobian is constant in mu and sigma for fixed x0.
ln_nll <- function(par) {
  mu <- par[1]; sigma <- exp(par[2])
  z <- (yd-mu)/sigma
  -sum(dnorm(y, mean=mu, sd=sigma, log=TRUE)) + n*pnorm(z, lower.tail=FALSE, log.p=TRUE)
}
# Start from transformed sample moments; request a tightly converged BFGS fit.
fit <- optim(c(mean(y),log(sd(y))), ln_nll, method='BFGS', hessian=TRUE,
             control=list(reltol=1e-13, maxit=2000))
if (fit$convergence!=0) warning('Lognormal optimizer convergence code: ',fit$convergence)
mu_hat <- fit$par[1]; sigma_hat <- exp(fit$par[2])
# Invert observed information in (mu, log sigma), then transform to (mu, sigma).
# Recompute curvature at the optimum and require positive observed information.
info <- optimHess(fit$par, ln_nll)
if (min(eigen(info, symmetric=TRUE, only.values=TRUE)$values) <= 0)
  stop('Lognormal observed information is not positive definite.')
cov_log <- solve(info)
# Delta transformation: d sigma / d log(sigma) = sigma; mu is unchanged.
J <- diag(c(1,sigma_hat))
cov_ln <- J %*% cov_log %*% J
# Approximate variance of the uncensored conditional Pareto shape MLE.
cov_pareto <- matrix(alpha_hat^2/n, 1, 1) # uncensored, left-truncated Pareto

# ---- Ground-up risk measures: robust, deterministic finite quadrature ----
# No calls to integrate(): eliminates adaptive-integration divergence errors.
# Write t=exp(-v), v=-log(t), so dt=exp(-v)dv and v is finite.
# A midpoint grid avoids t=0 and t=1 (and hence infinite qnorm values).
# Increase quadrature_max_v or reduce quadrature_step for accuracy checks.
quadrature_max_v <- 500
quadrature_step <- 0.02
# Build equal-width midpoint bins in v=-log(t); h is their actual spacing.
# The upper limit truncates the infinite v domain, so integration is approximate.
logtail_grid <- function(v_start=0) {
  count <- ceiling((quadrature_max_v-v_start)/quadrature_step)
  if (count < 1L) stop('Invalid quadrature grid settings')
  h <- (quadrature_max_v-v_start)/count
  list(v=v_start + (seq_len(count)-0.5)*h, h=h)
}
# Full-tail distortion measures start at t=1 (v=0); GS uses only t<=beta.
grid_all <- logtail_grid()
grid_cte <- logtail_grid(-log(beta))
# Multiply midpoint sums by bin width and stop on numerical overflow/invalid terms.
weighted_sum <- function(values, h, label) {
  if (any(!is.finite(values))) stop('Nonfinite quadrature terms in ', label)
  ans <- sum(values)*h
  if (!is.finite(ans)) stop('Nonfinite quadrature result in ', label)
  ans
}
# Wang transform used here: g(t)=Phi(Phi^{-1}(t)+lambda).
# g'(t)=exp(-lambda*Phi^{-1}(t)-lambda^2/2).
# Since t is upper-tail probability, Phi^{-1}(t)=-qnorm(t,lower.tail=FALSE),
# the weight is exp(lambda*z-lambda^2/2), z=qnorm(t,lower.tail=FALSE).
# Return five unconditional (ground-up) lognormal risk measures in NOK.
# Fitting is conditional on reporting, but these functionals use the full model.
risk_ln <- function(mu,sigma) {
  # Upper-tail normal quantile defines VaR; CTE has an analytic lognormal formula.
  zbeta <- qnorm(beta,lower.tail=FALSE)
  var <- x0 + exp(mu+sigma*zbeta)
  cte <- x0 + exp(mu+sigma^2/2)*pnorm(sigma-zbeta)/beta
  # GS integrates upper-tail quantiles against its linear tail-probability weight.
  vc <- grid_cte$v
  tc <- exp(-vc)
  zc <- qnorm(tc,lower.tail=FALSE)
  qc <- x0+exp(mu+sigma*zc)
  gs <- weighted_sum(qc*(1+2*delta-4*delta*tc/beta)*tc,
                     grid_cte$h,'lognormal GS')/beta
  # Use the full probability range for PHT and Wang distortion integrals.
  v <- grid_all$v
  t <- exp(-v)
  z <- qnorm(t,lower.tail=FALSE)
  q <- x0+exp(mu+sigma*z)
  # PHT derivative r*t^(r-1), multiplied by dt, becomes r*exp(-r*v).
  pht <- weighted_sum(q*r*exp(-r*v),grid_all$h,'lognormal PHT')
  # Wang derivative times dt produces exp(-v + lambda*z - lambda^2/2).
  wt <- weighted_sum(q*exp(-v+lambda*z-lambda^2/2),
                     grid_all$h,'lognormal WT')
  c(VaR=var,CTE=cte,GS=gs,PHT=pht,WT=wt)
}
# Pareto ground-up risk measures use x0 as the scale and a as the shape.
# Require a>1 for finite CTE and r*a>1 for finite PHT.
risk_pa <- function(a) {
  if (!is.finite(a) || a<=1 || r*a<=1)
    stop('Pareto shape outside risk-measure finiteness range.')
  # VaR, CTE, GS and PHT have closed forms; WT uses the common midpoint grid.
  var <- x0*beta^(-1/a)
  cte <- var*a/(a-1)
  gs <- cte*(1+2*delta/(2*a-1))
  pht <- x0*r*a/(r*a-1)
  # Use the full probability range for PHT and Wang distortion integrals.
  v <- grid_all$v
  t <- exp(-v)
  z <- qnorm(t,lower.tail=FALSE)
  # q(t)=x0*t^(-1/a); combine exponentials before evaluating.
  wt <- weighted_sum(x0*exp(-v*(1-1/a)+lambda*z-lambda^2/2),
                     grid_all$h,'Pareto WT')
  c(VaR=var,CTE=cte,GS=gs,PHT=pht,WT=wt)
}

# ---- Delta-method confidence intervals at 90%, 95%, 99% ----
# Central finite differences, with both mu and sigma uncertainty included.
# Numerical Jacobian: rows are risk measures and columns are model parameters.
num_jac <- function(fn,theta) {
  base <- fn(theta)
  out <- matrix(NA_real_,length(base),length(theta),
                dimnames=list(names(base),names(theta)))
  for (j in seq_along(theta)) {
    # Scale the central-difference step to the parameter magnitude.
    h <- 1e-4*max(1,abs(theta[j]))
    plus <- minus <- theta
    plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
    out[,j] <- (fn(plus)-fn(minus))/(2*h)
  }
  out
}
# Refit at one assumed x0 and return estimates, SEs, and convergence details.
# x0 is treated as fixed in each scenario; its uncertainty is not in the CIs.
compute_at_x0 <- function(location, initial_par=NULL) {
  if (!is.finite(location) || location >= truncation || location >= min(k))
    stop('Each x0 must be below the truncation point and all observed claims.')
  # The Pareto shape MLE is unchanged by x0 when d and the data are fixed.
  # The shifted-lognormal likelihood and covariance MUST be refitted for each x0.
  # Update the enclosing x0 used by risk_ln() and risk_pa(); reset after scenarios.
  x0 <<- location
  yloc <- log(k-location)
  ydloc <- log(truncation-location)
  # Rebuild the transformed sample and truncation correction for this location.
  local_nll <- function(par) {
    mu <- par[1]; sigma <- exp(par[2])
    -sum(dnorm(yloc,mean=mu,sd=sigma,log=TRUE)) +
      n*pnorm((ydloc-mu)/sigma,lower.tail=FALSE,log.p=TRUE)
  }
  start <- if (is.null(initial_par)) c(mean(yloc),log(sd(yloc))) else initial_par
  opt <- optim(start,local_nll,method='BFGS',
               control=list(reltol=1e-11,maxit=2500))
  if (opt$convergence != 0) warning('Lognormal optimization at x0=',location,
                                    ' returned code ',opt$convergence)
  # Invert local likelihood curvature and transform covariance to (mu, sigma).
  hess <- optimHess(opt$par,local_nll)
  if (any(!is.finite(hess)) || min(eigen(hess,symmetric=TRUE,
                                       only.values=TRUE)$values) <= 0)
    stop('Non-positive observed information at x0=',location)
  mu <- opt$par[1]; sigma <- exp(opt$par[2])
  jac <- diag(c(1,sigma))
  cov <- jac %*% solve(hess) %*% jac
  # Evaluate both model functionals at this scenario and differentiate numerically.
  le <- risk_ln(mu,sigma)
  pe <- risk_pa(alpha_hat)
  lg <- num_jac(function(par) risk_ln(par[1],par[2]),c(mu=mu,sigma=sigma))
  pg <- num_jac(function(par) risk_pa(par[1]),c(alpha=alpha_hat))
  # Var(R) = gradient(R)^T Cov(theta) gradient(R), one row per measure.
  # pmax removes small negative roundoff before taking square roots.
  ls <- sqrt(pmax(0,rowSums((lg%*%cov)*lg)))
  ps <- sqrt(pmax(0,rowSums((pg%*%cov_pareto)*pg)))
  list(x0=location,mu=mu,sigma=sigma,ln=le,pa=pe,ln_se=ls,pa_se=ps,
       convergence=opt$convergence,fit_par=opt$par)
}

cat(sprintf('PHT distortion parameter r = %.2f for BOTH LN and Pareto I, all x0 scenarios.\n', r))
cat('Fitting baseline and sensitivity models...\n')
# Preserve the requested baseline before scenario fits mutate the enclosing x0.
baseline_x0 <- x0
# Each scenario is refitted from its own transformed sample mean and SD.
# initial_par is available in compute_at_x0(), but this call does not supply it.
scenarios <- sort(unique(c(x0_sensitivity,baseline_x0)))
fits <- lapply(scenarios,function(xx) compute_at_x0(xx))
# Named lookup avoids depending on the order of the sorted scenario vector.
names(fits) <- as.character(scenarios)
base <- fits[[as.character(baseline_x0)]]
x0 <- baseline_x0
mu_hat <- base$mu
sigma_hat <- base$sigma

# Baseline parameter table: NA marks parameters that do not belong to a model.
params <- data.frame(
  Model=c('Pareto_I','Shifted_Lognormal'),
  Number_of_Claims=rep(n,2L),
  Ground_Up_Threshold_NOK=rep(baseline_x0,2L),
  Left_Truncation_NOK=rep(truncation,2L),
  Alpha=c(alpha_hat,NA_real_),Mu=c(NA_real_,mu_hat),
  Sigma=c(NA_real_,sigma_hat),check.names=FALSE
)

# Form symmetric Wald intervals estimate +/- z*SE, then convert NOK to millions.
# Endpoints are not clipped to a risk measure's admissible range.
make_ci_table <- function(level, fit=base) {
  z <- qnorm((1+level)/2)
  data.frame(
    Risk_Measure=names(fit$ln),
    LN_Estimate_million_NOK=unname(fit$ln)/1e6,
    LN_CI_Lower_million_NOK=unname(fit$ln-z*fit$ln_se)/1e6,
    LN_CI_Upper_million_NOK=unname(fit$ln+z*fit$ln_se)/1e6,
    Pareto_Estimate_million_NOK=unname(fit$pa)/1e6,
    Pareto_CI_Lower_million_NOK=unname(fit$pa-z*fit$pa_se)/1e6,
    Pareto_CI_Upper_million_NOK=unname(fit$pa+z*fit$pa_se)/1e6,
    check.names=FALSE
  )
}
# One table per confidence level; estimates repeat while interval widths vary.
ci_tables <- lapply(confidence_levels,make_ci_table)
names(ci_tables) <- paste0(round(100*confidence_levels),'pct')

# The sensitivity table reports both models' point estimates and changes
# relative to x0=100000. The observed claims and d are held fixed.
# Stack one row per scenario, model, and risk measure; percentage changes compare
# each model with its own fitted baseline, not with the competing model.
sensitivity <- do.call(rbind,lapply(fits,function(fit) {
  measures <- names(fit$ln)
  rbind(data.frame(
    X0_NOK=rep(fit$x0,length(measures)),
    Model=rep('Shifted_Lognormal',length(measures)),
    Risk_Measure=measures,
    Estimate_million_NOK=unname(fit$ln)/1e6,
    Change_from_baseline_percent=100*(unname(fit$ln)/unname(base$ln)-1),
    Alpha=NA_real_,Mu=fit$mu,Sigma=fit$sigma,
    stringsAsFactors=FALSE
  ),data.frame(
    X0_NOK=rep(fit$x0,length(measures)),
    Model=rep('Pareto_I',length(measures)),
    Risk_Measure=measures,
    Estimate_million_NOK=unname(fit$pa)/1e6,
    Change_from_baseline_percent=100*(unname(fit$pa)/unname(base$pa)-1),
    Alpha=alpha_hat,Mu=NA_real_,Sigma=NA_real_,
    stringsAsFactors=FALSE
  ))
}))
rownames(sensitivity) <- NULL

# ---- Ground-up probability P(X < 500000 NOK) ----
# These are unconditional fitted probabilities, extrapolated below the
# observation threshold; the left-truncated sample cannot estimate them empirically.
# Shifted LN: Phi((log(c-x0)-mu)/sigma), for c>x0.
# Pareto I: 1-(x0/c)^alpha, for c>x0.
# Evaluate the full fitted CDF at a common monetary cutoff in every scenario.
probability_cutoff <- 500000
probability_table <- do.call(rbind,lapply(fits,function(fit) {
  ln_probability <- if (probability_cutoff <= fit$x0) 0 else
    plnorm(probability_cutoff-fit$x0,meanlog=fit$mu,sdlog=fit$sigma)
  # -expm1(log survival) computes 1-survival accurately near probability zero.
  pa_probability <- if (probability_cutoff <= fit$x0) 0 else
    -expm1(alpha_hat*log(fit$x0/probability_cutoff))
  data.frame(
    X0_NOK=fit$x0,
    Cutoff_NOK=probability_cutoff,
    LN_Probability_X_Less_Than_Cutoff=ln_probability,
    Pareto_Probability_X_Less_Than_Cutoff=pa_probability,
    check.names=FALSE
  )
}))
rownames(probability_table) <- NULL

# This CSV is a report of labeled sections, not one rectangular data frame.
output_file <- file.path(folder,'NFC_Combined_Parameters_CIs_and_X0_Sensitivity.csv')
# First section overwrites an earlier report; later sections append a title,
# header, and rows. NA cells are written blank and row numbers are suppressed.
write_section <- function(title,dat,first=FALSE) {
  if (first) writeLines(title,output_file)
  else cat('\n',title,'\n',sep='',file=output_file,append=TRUE)
  write.table(dat,file=output_file,sep=',',row.names=FALSE,col.names=TRUE,
              quote=TRUE,na='',append=TRUE)
}
write_section('Table 1. Norwegian Fire Claims - Fitted Model Parameters (baseline x0)',params,TRUE)
# Number the interval sections after the parameter table, preserving level order.
for (i in seq_along(ci_tables)) {
  lev <- round(100*confidence_levels[i])
  write_section(paste0('Table ',i+1L,'. Norwegian Fire Claims - ',lev,
                       '% Confidence Intervals (million NOK; baseline x0)'),
                ci_tables[[i]])
}
write_section('Table 5. Sensitivity to Assumed Ground-Up Threshold x0 (NOK; fixed d and claims)',
              sensitivity)
write_section('Table 6. Ground-Up P(X < 500000 NOK) by Assumed x0',
              probability_table)

# ---- LaTeX sensitivity table generated from the SAME calculated results ----
# This avoids transcription/rounding differences between CSV and manuscript.
tex_file <- file.path(folder, 'NFC_X0_Sensitivity_From_R_Results.tex')
# Pull full-precision point estimates from fitted scenarios and convert to millions.
get_est <- function(xx, model, measure) {
  f <- fits[[as.character(xx)]]
  unname(if (model == 'LN') f$ln[measure] else f$pa[measure]) / 1e6
}
# Apply three-decimal display rounding only when constructing LaTeX cells.
fmt <- function(x) sprintf('%.3f', x)
# The manuscript layout uses these three x0 columns for each model.
# Additional scenarios remain in the CSV; these three must exist in scenarios.
scenario_cols <- c(50000, 100000, 200000)
if (!all(scenario_cols %in% scenarios))
  stop('LaTeX output requires x0_sensitivity = c(50000,100000,200000).')
# These manuscript labels are literal strings matching the default settings.
# If beta/delta/r/lambda change, update these labels to keep the captioned values accurate.
labels <- c(
  VaR = '$\\mathrm{VaR}[F,\\beta=0.10]$',
  CTE = '$\\mathrm{cte}[F,\\beta=0.10]$',
  GS  = '$\\mathrm{gs}[F,\\beta=0.10,\\delta=0.25]$',
  PHT = '$\\mathrm{pht}[F,r=0.95]$',
  WT  = '$\\mathrm{wt}[F,\\lambda=0.25]$'
)
# Double-backslashes in R source string literals are converted to single
# LaTeX backslashes when written to disk.
lines <- c(
  '% Add \\usepackage{xcolor} to the manuscript preamble if needed.',
  '\\begingroup', '\\color{blue}',
  '\\begin{table}[htbp]', '\\centering',
  '\\caption{Sensitivity analysis of point estimates of selected risk measures for Norwegian Fire Claims data based on LN and Pa I models for different choices of $x_0$ (millions nok).}',
  '\\label{tab:nfc-x0-sensitivity}',
  '\\small', '\\setlength{\\tabcolsep}{4pt}',
  '\\begin{tabular}{lrrrrrr}', '\\hline',
  ' & \\multicolumn{3}{c}{LN distribution} & \\multicolumn{3}{c}{Pa I distribution} \\\\',
  '\\cline{2-4}\\cline{5-7}',
  'Risk Measure $R[F]$ & $5\\times10^4$ & $10^5$ & $2\\times10^5$ & $5\\times10^4$ & $10^5$ & $2\\times10^5$ \\\\',
  '\\hline'
)
# Assemble each row with all LN scenarios first, followed by all Pareto scenarios.
for (m in names(labels)) {
  values <- c(vapply(scenario_cols, get_est, numeric(1), model='LN', measure=m),
              vapply(scenario_cols, get_est, numeric(1), model='Pa', measure=m))
  formatted <- fmt(values)
  # Emphasize the x0=100000 baseline column in each model block.
  formatted[c(2,5)] <- paste0('\\textbf{', formatted[c(2,5)], '}')
  lines <- c(lines, paste0(labels[[m]], ' & ', paste(formatted, collapse=' & '), ' \\\\'))
}
lines <- c(lines, '\\hline', '\\end{tabular}', '\\end{table}', '\\endgroup')
# Write a complete table fragment for inclusion in a LaTeX manuscript.
writeLines(lines, tex_file)
cat('LaTeX sensitivity table:', tex_file, '\\n')

# Restore baseline for Q-Q plots and printed diagnostics.
x0 <- baseline_x0
# Select the default 90% interval table for the final console summary.
# This lookup assumes confidence_levels includes 0.90.
result <- ci_tables[['90pct']]
# ---- Figure 5.1: fitted versus observed truncated quantiles ----
if (make_qq_plot) {
  # Midpoint plotting positions avoid zero/one probabilities and infinite quantiles.
  p <- (seq_len(n)-0.5)/n
  # Compare with conditional quantiles above d, matching the observed sample.
  qpa <- truncation*(1-p)^(-1/alpha_hat)
  # conditional CDF of X at truncation under fitted lognormal
  Fd <- pnorm((log(truncation-x0)-mu_hat)/sigma_hat)
  # Map conditional p to unconditional F(d)+(1-F(d))*p before inversion.
  qln <- x0+exp(mu_hat+sigma_hat*qnorm(Fd+(1-Fd)*p))
  # Open one PDF page with two panels; both axes use natural logarithms.
  pdf(file.path(folder,'NFC_Figure_5_1_QQ.pdf'),width=10,height=5)
  par(mfrow=c(1,2),mar=c(4.3,4.3,2.3,1.0))
  plot(log(qln),log(k),xlab='LOG ( Fitted Lognormal Quantiles )',
       ylab='LOG ( Empirical Quantiles )',main='Lognormal QQ-plot',pch=16,cex=.4)
  # The identity line marks equality of fitted and empirical log-quantiles.
  abline(0,1,lty=2,col='gray40')
  plot(log(qpa),log(k),xlab='LOG ( Fitted Pareto Quantiles )',
       ylab='LOG ( Empirical Quantiles )',main='Pareto QQ-plot',pch=16,cex=.4)
  # The identity line marks equality of fitted and empirical log-quantiles.
  abline(0,1,lty=2,col='gray40')
  # Close the graphics device to finish writing the PDF.
  dev.off()
}
# Print baseline diagnostics, the selected CI table, output paths and probabilities.
cat('\nNorwegian Fire Claims n=',n,'\n',sep='')
cat('Pareto alpha=',format(alpha_hat,digits=9),'\n')
cat('Lognormal mu=',format(mu_hat,digits=9),' sigma=',format(sigma_hat,digits=9),'\n')
print(round(result[,-1],3))
cat('\nCombined CSV:',output_file,'\n')
cat('QQ plot (if enabled):',file.path(folder,'NFC_Figure_5_1_QQ.pdf'),'\n')

cat('\nGround-up P(X < 500000 NOK), baseline and x0 sensitivity:\n')
print(probability_table, row.names=FALSE, digits=8)
