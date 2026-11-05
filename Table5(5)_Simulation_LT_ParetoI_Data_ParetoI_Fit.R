# CORRECT SPECIFICATION: Pareto I data fitted with Pareto I.
# True survival: P(X>x)=(x0/x)^alpha for x>=x0.
# Fit S_P(x)=(x0/x)^alpha with fixed x0 set in the configuration below.
# Outer data: conditional on X>=d, fully observed without right censoring.
# Coverage targets the true GROUND-UP Pareto I risk measures.
# Model-based Delta, Pareto parametric bootstrap, ordinary chi-square LR.
# No sandwich correction or misspecification-adjusted LR calibration.
# Base R only. R execution unavailable during authoring.
# Run with Rscript using the filename of this script.

integration_tolerance <- 1e-8

# STUDY DESIGN AND USER SETTINGS
# Pareto I generates the outer samples and supplies fitted risks and intervals.
# The target is the unconditional (ground-up) risk, even though fitting uses
# only observations above d. Coverage assesses interval calibration and
# ground-up extrapolation under a correctly specified model.
# MODE selects the replication budget; use validation for a small local pilot.
MODE <- 'production'  # Set to 'production' after a local pilot
# N_MC is the number of independent outer datasets for EACH sample size.
# B is the number of parametric bootstrap refits WITHIN each usable outer fit.
# Increasing B improves bootstrap quantile resolution; increasing N_MC reduces
# Monte Carlo uncertainty in reported performance measures.
N_MC <- if(MODE=='production') 100000L else 10L
B <- if(MODE=='production') 500L else 30L
sample_sizes <- c(647L)
# x0 is the true and fitted fixed Pareto minimum; alpha is the true tail index.
# d is the left-truncation threshold: all n observed claims satisfy X>=d.
# No censoring, upper loss limit, or cap is imposed.
x0 <- 100000
alpha <- 1.13
d <- 500000
# confidence is the nominal two-sided confidence level for all three methods.
# beta is the UPPER-tail probability: VaR corresponds to quantile 1-beta.
# delta controls the tail Gini shortfall loading; r is the PHT distortion power;
# lambda is the positive Wang distortion shift. These are fixed, not estimated.
confidence <- .90
beta <- .10; delta <- .25; r <- .95; lambda <- .25
# Require this fraction of B bootstrap alpha estimates to be finite and positive.
# A value of 1 requires every refit; failed refits are never replaced by redraws.
# bootstrap_block limits the temporary n-by-block matrix, not the total B.
bootstrap_required_fraction <- 1  # .979: at least 979 valid of B=1000
bootstrap_block <- 250L
# Outer and bootstrap random streams use separate seed bases. Within each,
# cell index and replication index determine a reproducible seed.
# checkpoint_every controls how much completed work could be lost on interruption.
# resume loads a compatible saved state; RUN_CHECKS enables deterministic audits.
seed <- 20261022L
checkpoint_every <- 100L
resume <- TRUE
RUN_CHECKS <- TRUE
# The version enters the output folder name. It separates publication-summary
# checkpoints from older accumulator layouts. outer_rng_kind records the active
# R RNG algorithms; reproducibility also depends on the R numerical environment.
version <- 'ParetoI_to_ParetoI_left_truncated_v1_publication'
bootstrap_seed <- 20261023L
outer_rng_kind <- RNGkind()
delta_information <- 'expected_Fisher_per_observation'
# All risk vectors use this exact order; endpoint vectors and truth must align.
# VaR: value at risk; CTE: conditional tail expectation; GS: Gini shortfall;
# PHT: proportional hazards transform; WT: Wang transform.
measures <- c('VaR','CTE','GS','PHT','WT')
methods <- c('Delta','Bootstrap_percentile','Profile_likelihood')

# LOCATE THE SCRIPT AND PREPARE OUTPUTS
# Try Rscript --file first, then source() frames, then an optional RStudio API.
# Interactive sessions without a discoverable script fall back to getwd().
# The simulation itself uses base R; rstudioapi is optional path discovery only.
script <- NULL
args <- grep('^--file=',commandArgs(FALSE),value=TRUE)
if(length(args)) script <- sub('^--file=','',args[1])
if(is.null(script)) for(fr in rev(sys.frames())) {
  if(!is.null(fr$ofile) && file.exists(fr$ofile)) {script <- fr$ofile; break}
}
if(is.null(script) && requireNamespace('rstudioapi',quietly=TRUE))
  script <- tryCatch(if(rstudioapi::isAvailable())
    rstudioapi::getSourceEditorContext()$path else NULL,error=function(e) NULL)
folder <- if(length(script)==1 && !is.na(script) && nzchar(script) && file.exists(script))
  dirname(normalizePath(script)) else getwd()
# Each design/version gets a descriptive output directory beside the script.
# Checkpoint configuration comparison checks additional settings not in its name;
# changing settings may therefore require a new version or output folder.
outdir <- file.path(folder,sprintf('MC_%s_%s_x0%g_alpha%g_d%g_r%g_lambda%g',
  version,MODE,x0,alpha,d,r,lambda))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
# Reject settings outside the implemented parameter domains before simulation.
# The normal critical value is for symmetric Delta intervals; the chi-square
# critical value is for a one-parameter ordinary likelihood-ratio interval.
stopifnot(MODE %in% c('production','validation'),N_MC>0,B>=20,
          all(sample_sizes>=2),is.finite(alpha),alpha>1,alpha*r>1,x0>0,d>x0,
          beta>0,beta<1,delta>=0,delta<.5,r>0,r<=1,lambda>0,
          confidence>0,confidence<1,bootstrap_block>=1,
          bootstrap_required_fraction>0,bootstrap_required_fraction<=1)
zcrit <- qnorm((1+confidence)/2)
cutoff <- qchisq(confidence,1)

# Numerical diagnostics, checkpointed with the streaming summaries.
# Diagnostic counters record evaluation events, not necessarily distinct datasets:
# checks, endpoints and derivatives can evaluate a functional more than once.
# An environment lets nested numerical routines increment shared counters.
# run_cell restores/saves these counters with its own checkpoint state.
diag_env <- new.env(parent=emptyenv())
diag_env$counts <- c(WT_Integration_Failures=0,WT_Finite_Overflow=0,
                     Risk_Finite_Overflow=0,CI_Boundary_Exclusions=0,CI_Infinite_Exclusions=0)
count_diagnostic <- function(name) {
  diag_env$counts[name] <- diag_env$counts[name]+1
}

# Stable log upper inverse Mills ratio log(phi(z)/S(z)).
log_mills <- function(z) {
  # For moderate z, subtract log survival from log density directly.
  # For large positive z, use an asymptotic expansion of the inverse Mills ratio;
  # this avoids subtracting nearly equal, very large negative log probabilities.
  ans <- numeric(length(z)); high <- z>20
  if(any(high)) {
    t <- 1/z[high]
    excess <- t*(1+t^2*(-2+t^2*(10+t^2*(-74+t^2*(706-8162*t^2)))))
    ans[high] <- log(z[high])+log1p(excess/z[high])
  }
  if(any(!high)) ans[!high] <- dnorm(z[!high],log=TRUE)-
    pnorm(z[!high],lower.tail=FALSE,log.p=TRUE)
  ans
}
# WT and optional derivative with respect to log(alpha).
# Normal-coordinate quadrature with square completion avoids a long tail
# near alpha=1. A finite-domain overflow is NA, not theoretical infinity.
wang <- function(a,gradient=FALSE) {
  if(a<=1) return(list(value=Inf,derivative=NA_real_))
  tryCatch({
    # a>1 is the finite WT domain. Completing the square gives a Gaussian
    # kernel in t centered at lambda/q with scale 1/sqrt(q).
    # As a approaches 1 from above, the center and WT can become extremely large.
    q <- 1-1/a; center <- lambda/q; scale <- 1/sqrt(q)
    # Subtract a constant from the log integrand to keep quadrature values
    # numerically manageable, then restore it in logvalue.
    shift <- log_mills(center)/a
    integrand <- function(t) exp(-t^2/2+log_mills(center+scale*t)/a-shift)
    I <- integrate(integrand,-Inf,Inf,rel.tol=1e-8,subdivisions=1000)$value
    logvalue <- log(x0)+lambda^2/(2*q)-lambda^2/2-
      q*log(2*pi)/2-log(q)/2+shift+log(I)
    # A finite mathematical WT that exceeds double precision is a numerical
    # failure (NA), distinct from theoretical divergence (Inf at a<=1).
    if(!is.finite(logvalue) || logvalue>log(.Machine$double.xmax)) {
      count_diagnostic('WT_Finite_Overflow')
      return(list(value=NA_real_,derivative=NA_real_))
    }
    value <- exp(logvalue)
    derivative <- NA_real_
    # Differentiate with respect to eta=log(a), not directly with respect to a.
    # J/I is a weighted integral ratio; the negative sign reflects decreasing risk
    # as the Pareto tail becomes lighter. A failed derivative leaves the value usable.
    if(gradient) {
      J <- integrate(function(t) {
        integrand(t)*(-pnorm(center+scale*t,lower.tail=FALSE,log.p=TRUE))/a
      },-Inf,Inf,rel.tol=1e-8,subdivisions=1000)$value
      derivative <- -value*J/I
      if(!is.finite(derivative)) derivative <- NA_real_
    }
    list(value=value,derivative=derivative)
  # An integration error is recorded and converted to missing output so the
  # outer study can count the failure and continue rather than silently omit it.
  },error=function(e) {
    count_diagnostic('WT_Integration_Failures')
    list(value=NA_real_,derivative=NA_real_)
  })
}
# Risk vector and analytic derivatives with respect to eta=log(alpha).
risk_all <- function(a,gradient=FALSE) {
  if(!is.finite(a) || a<=0)
    return(list(value=setNames(rep(NA_real_,5),measures),derivative=rep(NA_real_,5)))
  # Initialize risks at theoretical infinity; overwrite only inside each
  # functional's finite domain. Invalid alpha estimates instead return all NA.
  # Derivatives are with respect to log(alpha), so d/deta = alpha*d/dalpha.
  value <- setNames(rep(Inf,5),measures)
  derivative <- setNames(rep(NA_real_,5),measures)
  # Pareto survival inversion gives VaR_beta=x0*beta^(-1/a).
  # CTE is its tail mean; GS multiplies CTE by its tail Gini loading.
  value['VaR'] <- x0*beta^(-1/a)
  derivative['VaR'] <- value['VaR']*log(beta)/a
  if(a>1) {
    value['CTE'] <- value['VaR']*a/(a-1)
    factor <- 1+2*delta/(2*a-1)
    value['GS'] <- value['CTE']*factor
    derivative['CTE'] <- value['CTE']*(log(beta)/a-1/(a-1))
    derivative['GS'] <- value['GS']*(log(beta)/a-1/(a-1)-
                                      4*delta*a/((2*a-1)^2*factor))
  }
  # Under the power distortion S(x)^r, the effective Pareto exponent is r*a.
  # Its integral is finite exactly when r*a>1.
  if(r*a>1) {
    value['PHT'] <- x0*r*a/(r*a-1)
    derivative['PHT'] <- -x0*r*a/(r*a-1)^2
  }
  wt <- wang(a,gradient)
  value['WT'] <- wt$value; derivative['WT'] <- wt$derivative
  # VaR is theoretically finite for a>0. Other finite domains are explicit.
  finite_domain <- c(TRUE,a>1,a>1,r*a>1,a>1)
  overflow <- finite_domain & is.infinite(value)
  if(any(overflow)) {
    count_diagnostic('Risk_Finite_Overflow'); value[overflow] <- NA_real_
  }
  derivative[!is.finite(derivative)] <- NA_real_
  list(value=value,derivative=derivative)
}
# Delta-method components matching the supplied table.
# All conditional observations are exact: p(a)=1.
# Per-observation I_alpha=1/a^2 and I_log(alpha)=1; variance divides by n.
# D terms are chosen with positive signs; the table uses their squares.
# WT: C=a*(WT/x0-1), D=C-a*dC/da=-a*(dWT/dlog(a))/x0.
delta_table_components <- function(a,point=NULL) {
  if(is.null(point)) point <- risk_all(a,gradient=TRUE)
  # probability refers to observing an exact claim in the conditional sampling
  # scheme; it is NOT the population retention probability P(X>=d).
  # With n fully observed retained claims it equals one in these table formulas.
  probability <- 1
  coefficient <- setNames(rep(NA_real_,length(measures)),measures)
  D_CTE <- D_GS <- C_WT <- D_WT <- NA_real_
  if(!is.finite(probability) || probability<=0)
    return(list(probability=probability,coefficient=coefficient,
                D_CTE=D_CTE,D_GS=D_GS,C_WT=C_WT,D_WT=D_WT))
  V <- point$value
  if(is.finite(V['VaR']))
    coefficient['VaR'] <- log(beta)^2*V['VaR']^2/(a^2*probability)
  if(a>1) {
    D_CTE <- -log(beta)+a/(a-1)
    D_GS <- D_CTE+4*delta*a^2/((2*a-1)*(2*a-1+2*delta))
    if(is.finite(V['CTE'])) coefficient['CTE'] <- D_CTE^2*V['CTE']^2/(a^2*probability)
    if(is.finite(V['GS'])) coefficient['GS'] <- D_GS^2*V['GS']^2/(a^2*probability)
    if(is.finite(V['WT']) && is.finite(point$derivative['WT'])) {
      C_WT <- a*(V['WT']/x0-1)
      D_WT <- -a*point$derivative['WT']/x0
      coefficient['WT'] <- x0^2*D_WT^2/(a^2*probability)
    }
  }
  if(r*a>1 && is.finite(V['PHT']))
    coefficient['PHT'] <- V['PHT']^2/((r*a-1)^2*probability)
  # The coefficient is the asymptotic variance of sqrt(n)*(estimated risk-risk)
  # under the working Pareto model. The actual squared standard error is coefficient/n.
  # Overflow in a squared term makes this Delta interval unavailable.
  coefficient[!is.finite(coefficient)] <- NA_real_
  list(probability=probability,coefficient=coefficient,
       D_CTE=unname(D_CTE),D_GS=unname(D_GS),C_WT=unname(C_WT),D_WT=unname(D_WT))
}
# Independent survival-integral definitions for checking the WT table constants.
# Let t=exp(-v). C=integral exp(v/a)*g(exp(-v)) dv; D adds factor 1+v/a.
# This audit is performed at an arbitrary working-model shape, not each fit.
wang_table_constants_check <- function(a) {
  stopifnot(a>1)
  # This separate integral representation audits WT constants independently
  # of the transformed quadrature used in wang(); agreement checks algebra and
  # numerics under Pareto, not the adequacy of Pareto for the true distribution.
  logbase <- function(v) v/a+pnorm(qnorm(-v,log.p=TRUE)+lambda,log.p=TRUE)
  C <- integrate(function(v) exp(logbase(v)),0,Inf,
                 rel.tol=1e-8,subdivisions=1000)$value
  D <- integrate(function(v) exp(logbase(v)+log1p(v/a)),0,Inf,
                 rel.tol=1e-8,subdivisions=1000)$value
  c(C_WT=C,D_WT=D)
}
# True ground-up Pareto I targets, as in the original Pareto-data script.
pareto_true_risks <- function() {
  tail_quantile <- x0*beta^(-1/alpha)
  # Use p=Phi(z) in the distorted-quantile representation. The factor
  # exp(dnorm(z+lambda,log=TRUE)) is phi(z+lambda), and the Pareto quantile
  # factor is Phi(z)^(-1/alpha). Log probabilities protect the lower tail.
  wang_integrand <- function(z) {
    ans <- exp(-pnorm(z,log.p=TRUE)/alpha+dnorm(z+lambda,log=TRUE))
    ans[is.infinite(z)] <- 0
    ans
  }
  wang <- x0*integrate(wang_integrand,-Inf,Inf,
    rel.tol=integration_tolerance,subdivisions=1000L)$value
  c(VaR=tail_quantile,
    CTE=tail_quantile*alpha/(alpha-1),
    GS=tail_quantile*((1+2*delta)/(1-1/alpha)-4*delta/(2-1/alpha)),
    PHT=x0*alpha*r/(alpha*r-1),WT=wang)
}

true_risks <- pareto_true_risks()
stopifnot(all(is.finite(true_risks)))
write.csv(data.frame(Risk_Measure=measures,True_Risk=unname(true_risks)),
          file.path(outdir,"True_Ground_Up_Risks.csv"),row.names=FALSE)
true_retention_probability <- (x0/d)^alpha
write.csv(data.frame(Truncation_Probability=1-true_retention_probability,
 Retention_Probability=true_retention_probability),
 file.path(outdir,"True_Truncation_Probability.csv"),row.names=FALSE)

# Original Pareto sampler: identical rexp draws and loss construction.
# The Pareto likelihood uses log(X/d)=log_ratio rather than the LN coordinate.
yd <- log(d-x0)
simulate_pareto_statistics <- function(n) {
  log_ratio <- rexp(n,rate=alpha)
  log_x <- log(d)+log_ratio
  # log(X-x0)=log(X)+log1p(-x0/X), stable even for extreme X.
  exposure <- log_x+log1p(-exp(log(x0)-log_x))-yd
  if(any(!is.finite(exposure)) || any(exposure<0))
    stop('Invalid uncensored Pareto sample')
  T <- sum(log_ratio)
  c(n=n,T=T,ahat=if(is.finite(T) && T>0) n/T else NA_real_)
}
# Outer draws have a separate deterministic seed for every replication.
# Validation generates only N_MC datasets, without a production-sized warmup.
outer_configuration <- list(distribution='Pareto I',seed=seed,
  RNGkind=outer_rng_kind,x0=x0,alpha=alpha,
  n=sample_sizes,N_MC=N_MC,d=d,observation_scheme='left truncation only')
outer_path <- file.path(outdir,'ParetoI_Outer_Statistics.rds')
# Reuse sufficient statistics only when the whole outer design matches exactly.
# Otherwise create N_MC rows for each n. Storing n, T and alpha_hat avoids keeping
# all raw claims and lets resumed interval calculations use identical datasets.
if(file.exists(outer_path)) {
  outer_cache <- readRDS(outer_path)
  if(!identical(outer_cache$configuration,outer_configuration))
    stop('Outer sample cache configuration differs; use a new output folder')
} else {
  outer_cells <- vector('list',length(sample_sizes))
  for(ii in seq_along(sample_sizes)) {
    nn <- sample_sizes[ii]
    stats <- matrix(NA_real_,N_MC,3L,dimnames=list(NULL,c('n','T','ahat')))
    for(rr in seq_len(N_MC)) {
      # Seed offsets identify each dataset separately; restarting or changing
      # checkpoint frequency therefore does not shift later outer draws.
      set.seed(seed+ii*1000000L+rr)
      stats[rr,] <- simulate_pareto_statistics(nn)
      if(rr %% 10000L==0L) {
        cat('Preparing Pareto I samples: n=',nn,', ',rr,'/',N_MC,'\n',sep='')
        flush.console()
      }
    }
    outer_cells[[ii]] <- stats
  }
  outer_cache <- list(configuration=outer_configuration,cells=outer_cells)
  saveRDS(outer_cache,outer_path)
}
# Boundary counts include invalid/nonpositive working-model shapes. Keep them as an
# explicit audit; failed fits remain separate from undefined delta intervals.
outer_boundary_audit <- function(n,cell,completed,acc,fit_failures) {
  aa <- cell[seq_len(completed),'ahat']
  invalid <- !is.finite(aa) | aa<=0
  # Audit CTE/GS at alpha<=1 and PHT at r*alpha<=1, including invalid fits.
  # WT shares the alpha=1 boundary but is not a separate row in this legacy audit.
  # These counts concern fitted-risk boundaries, not all interval endpoint failures.
  masks <- list(CTE=invalid | aa<=1,GS=invalid | aa<=1,
                PHT=invalid | r*aa<=1)
  ans <- data.frame(Risk_Measure=names(masks),Sample_Size=n,
    Replications=completed,Violation_Count=vapply(masks,sum,integer(1)))
  ans$Violation_Rate <- ans$Violation_Count/completed
  ans$Violation_Percent <- 100*ans$Violation_Rate
  ans$Delta_Boundary_Excluded <- vapply(ans$Risk_Measure,function(risk)
    acc$Delta_Boundary_Excluded[acc$Method=='Delta' & acc$Risk_Measure==risk],integer(1))
  ans$Invalid_Shape_Count <- sum(invalid)
  ans$Outer_Fit_Failures <- fit_failures
  # With valid positive shapes, the definitions agree replication by replication.
  stopifnot(all(ans$Violation_Count==ans$Delta_Boundary_Excluded+sum(invalid)))
  ans
}
# PROFILE LIKELIHOOD ENDPOINTS
# Only n and the confidence cutoff determine the alpha endpoint ratios, so cache
# them across replications. alpha_hat supplies the dataset-specific scale.
lr_cache <- new.env(parent=emptyenv())
lr_ratios <- function(n) {
  if(n<=0) return(c(NA_real_,NA_real_))
  key <- as.character(n)
  if(exists(key,envir=lr_cache,inherits=FALSE)) return(get(key,envir=lr_cache))
  # Write alpha=alpha_hat*exp(t). Twice the log-likelihood loss is
  # 2*n*(exp(t)-1-t). expm1(t) is accurate near zero.
  # Find one negative and one positive root; exponentiate to get alpha ratios.
  f <- function(t) 2*n*(expm1(t)-t)-cutoff
  left <- -1; right <- 1
  while(f(left)<0) left <- left*2
  while(f(right)<0) right <- right*2
  ans <- exp(c(uniroot(f,c(left,0),tol=1e-10)$root,
               uniroot(f,c(0,right),tol=1e-10)$root))
  assign(key,ans,envir=lr_cache)
  ans
}
# Vectorized Pareto bootstrap alpha estimates, conditional on X>=d.
# No censoring is applied. At least ceiling(B*fraction) valid alpha refits are required.
bootstrap_alpha <- function(a,n) {
  # Each column represents one bootstrap dataset of n exponential log-exposures
  # under the fitted Pareto model. Column sums yield B MLE refits without generating
  # raw Pareto claims. Blocking bounds memory while retaining all B alpha estimates.
  estimates <- rep(NA_real_,B)
  for(first in seq.int(1L,B,by=bootstrap_block)) {
    last <- min(B,first+bootstrap_block-1L); nb <- last-first+1L
    e <- matrix(rexp(n*nb,rate=a),nrow=n,ncol=nb)
    T <- colSums(e)
    good <- is.finite(T) & T>0
    values <- rep(NA_real_,nb); values[good] <- n/T[good]
    estimates[first:last] <- values
  }
  good <- is.finite(estimates) & estimates>0
  # Leave both risk-bound alpha values missing if too few refits succeed.
  # The success threshold applies to alpha fits, not to finite risk functionals:
  # a valid positive alpha can still imply an infinite CTE/GS/PHT/WT.
  selected <- c(NA_real_,NA_real_)
  if(sum(good)>=ceiling(B*bootstrap_required_fraction)) {
    aa <- sort(estimates[good]); K <- length(aa)
    i_low <- max(1L,ceiling(K*(1-confidence)/2))
    i_high <- min(K,ceiling(K*(1+confidence)/2))
    # Risk decreases in alpha: reverse the order-statistic indices.
    # This exactly reproduces quantile(risk(alpha_boot),type=1), including
    # integer-index cases where mapping complementary quantiles is off by one.
    selected <- c(aa[K-i_low+1L],aa[K-i_high+1L])
  }
  list(alpha_for_risk_bounds=selected,valid=sum(good),failed=sum(!good))
}
# DETERMINISTIC FORMULA CHECKS
# Verify the LR cutoff, log-alpha derivatives, percentile ordering, and agreement
# between the table and gradient expressions. These are implementation audits;
# they do not establish nominal coverage under misspecification.
run_checks <- function() {
  n_check <- 20; ahat <- 2; T <- n_check/ahat
  ll <- function(a) n_check*log(a)-a*T
  ratios <- lr_ratios(n_check)
  stopifnot(all(vapply(ratios,function(q)
    abs(2*(ll(ahat)-ll(ahat*q))-cutoff)<1e-6,logical(1))))
  # Central differences perturb log(alpha), matching the analytic derivative's
  # parameterization. Relative tolerances accommodate different risk scales.
  eps <- 1e-5; point <- risk_all(ahat,TRUE)
  numerical <- (risk_all(ahat*exp(eps))$value-risk_all(ahat*exp(-eps))$value)/(2*eps)
  stopifnot(all(abs(numerical-point$derivative)/(1+abs(numerical))<1e-5))
  # Include a PHT value outside its finite domain to check percentile mapping
  # when the transformed sample contains theoretical infinities.
  aa <- c(.8,1.1,1.4,1.8,2,2.5,3,4,5,6)
  risks <- vapply(aa,function(a) risk_all(a)$value['PHT'],numeric(1))
  K <- length(aa); ps <- c(.1,.9); idx <- ceiling(K*ps)
  transformed <- vapply(aa[K-idx+1L],function(a) risk_all(a)$value['PHT'],numeric(1))
  stopifnot(isTRUE(all.equal(unname(quantile(risks,ps,type=1)),unname(transformed))))
  # Compare every table coefficient with the derivative/Fisher-information form.
  for(aa in c(1.5,2,3)) {
    pp <- risk_all(aa,gradient=TRUE)
    tt <- delta_table_components(aa,pp)
    via_gradient <- pp$derivative^2/tt$probability
    stopifnot(all(abs(tt$coefficient-via_gradient)/(1+abs(via_gradient))<1e-7))
  }
  ww <- wang_table_constants_check(2)
  tt <- delta_table_components(2)
  stopifnot(abs(tt$C_WT-ww['C_WT'])/(1+abs(ww['C_WT']))<1e-6,
            abs(tt$D_WT-ww['D_WT'])/(1+abs(ww['D_WT']))<1e-6)
  cat('LR, gradients, bootstrap ordering and expected-information table checks passed.\n')
}
if(RUN_CHECKS) run_checks()

# Only finite, ordered intervals enter Valid, Covered and length summaries.
# Theoretical risk-boundary exclusions are tracked separately from numerical
# failures and outer fits that did not exist.
# ACCUMULATOR LAYOUT
# One row per risk/method combination, reused for a single sample size.
# Valid and Covered count only ordered finite intervals with representable length.
# Boundary_Excluded combines infinite non-Delta intervals and undefined Delta
# boundary cases; Numerical_CI_Failures excludes failures of the outer fit itself.
new_accumulator <- function() {
  a <- expand.grid(Risk_Measure=measures,Method=methods,stringsAsFactors=FALSE)
  a$Valid <- a$Covered <- a$Finite_Length_Count <- 0L
  a$Boundary_Excluded <- a$Infinite_CI_Excluded <- a$Delta_Boundary_Excluded <- 0L
  a$Numerical_CI_Failures <- 0L
  a$Sum_Finite_Length <- 0
  # Exact quantiles: retain finite lengths in checkpointed list columns.
  # I() preserves each list as one data-frame column. Exact lengths cost memory
  # proportional to N_MC times the number of risk/method rows; they support exact
  # sample quantiles rather than an approximate streaming quantile algorithm.
  # Largest_WT separately retains at most ten complete finite interval records per row.
  a$Finite_Lengths <- I(lapply(seq_len(nrow(a)),function(i) numeric(0)))
  a$Largest_WT <- I(lapply(seq_len(nrow(a)),function(i) NULL))
  a
}
# Classify each interval once, using the priority below:
# (1) explicit theoretical boundary or usable ordered infinite endpoints;
# (2) missing, nonfinite, reversed, or overflowing-length numerical failures;
# (3) finite valid intervals, eligible for coverage and length summaries.
# Input lower/upper vectors must follow measures; boundary_violation is explicit
# because an undefined Delta functional is different from an infinite endpoint.
update_accumulator <- function(a,method,lower,upper,
                               boundary_violation=rep(FALSE,length(measures))) {
  for(j in seq_along(measures)) {
    i <- which(a$Method==method & a$Risk_Measure==measures[j])
    ci <- c(lower[j],upper[j])
    # For delta CIs, an infinite computed endpoint can be arithmetic overflow.
    # A theoretical delta boundary failure is passed explicitly below.
    # Mixed NA/infinite or reversed bounds are numerical failures.
    infinite_endpoint <- method!='Delta' && !anyNA(ci) &&
      ci[1]<=ci[2] && any(is.infinite(ci))
    if(isTRUE(boundary_violation[j]) || infinite_endpoint) {
      a$Boundary_Excluded[i] <- a$Boundary_Excluded[i]+1L
      count_diagnostic('CI_Boundary_Exclusions')
      if(infinite_endpoint) {
        a$Infinite_CI_Excluded[i] <- a$Infinite_CI_Excluded[i]+1L
        count_diagnostic('CI_Infinite_Exclusions')
      } else if(method=='Delta') {
        # A delta CI is undefined when the fitted functional is infinite;
        # it is not itself an interval with an infinite endpoint.
        a$Delta_Boundary_Excluded[i] <- a$Delta_Boundary_Excluded[i]+1L
      }
      next
    }
    if(anyNA(ci) || !all(is.finite(ci)) || ci[1]>ci[2] || !is.finite(ci[2]-ci[1])) {
      a$Numerical_CI_Failures[i] <- a$Numerical_CI_Failures[i]+1L
      next
    }
    # Coverage includes equality at either endpoint. Length summaries and coverage
    # share the same finite-valid subset, so every accepted interval contributes once.
    len <- ci[2]-ci[1]
    a$Valid[i] <- a$Valid[i]+1L
    a$Covered[i] <- a$Covered[i]+as.integer(ci[1]<=true_risks[j] && true_risks[j]<=ci[2])
    a$Finite_Length_Count[i] <- a$Finite_Length_Count[i]+1L
    a$Sum_Finite_Length[i] <- a$Sum_Finite_Length[i]+len
    a$Finite_Lengths[[i]] <- c(a$Finite_Lengths[[i]],len)
    # Keep extreme finite WT intervals with replication IDs so the same outer
    # dataset can be reconstructed from its deterministic seed. Alpha distance from
    # one helps diagnose near-boundary instability but is not a significance measure.
    if(measures[j]=='WT') {
      context <- attr(a,'replication_context')
      if(is.null(context)) context <- c(replication=NA_real_,ahat=NA_real_)
      record <- data.frame(Replication=unname(context['replication']),
        Alpha_Hat=unname(context['ahat']),Alpha_Distance_From_One=unname(context['ahat'])-1,
        Lower=ci[1],Upper=ci[2],Finite_Length=len,
        Covered=ci[1]<=true_risks[j] && true_risks[j]<=ci[2])
      # Sort within this risk/method row and retain the ten largest observed lengths.
      # If fewer than ten finite WT intervals exist, retain all of them.
      largest <- rbind(a$Largest_WT[[i]],record)
      a$Largest_WT[[i]] <- head(largest[order(largest$Finite_Length,decreasing=TRUE),],10L)
    }
  }
  a
}
# Check boundary/NA classification without running the Monte Carlo study.
if(RUN_CHECKS) {
  # Synthetic cases exercise the classifier independently of the Monte Carlo loop.
  # Four non-VaR Pareto functionals diverge at alpha=.9; Delta marks those risks
  # undefined, while non-Delta transformations may have infinite endpoints.
  test <- new_accumulator()
  lo <- risk_all(2)$value; hi <- risk_all(.9)$value
  test <- update_accumulator(test,'Profile_likelihood',lo,hi)
  idx <- test$Method=='Profile_likelihood'
  stopifnot(sum(test$Valid[idx])==1L,sum(test$Infinite_CI_Excluded[idx])==4L,
            sum(test$Boundary_Excluded[idx])==4L)
  test <- update_accumulator(test,'Delta',rep(NA_real_,5),rep(NA_real_,5),
                             boundary_violation=is.infinite(hi))
  idx <- test$Method=='Delta'
  stopifnot(sum(test$Delta_Boundary_Excluded[idx])==4L,
            sum(test$Numerical_CI_Failures[idx])==1L)
  # A mixed missing/infinite interval is a numerical failure, not a usable infinite CI.
  mixed <- update_accumulator(new_accumulator(),'Profile_likelihood',
    rep(NA_real_,5),rep(Inf,5))
  stopifnot(sum(mixed$Numerical_CI_Failures)==5L,
            sum(mixed$Infinite_CI_Excluded)==0L)
  lengths_test <- new_accumulator()
  for(width in c(1,2,1000)) lengths_test <- update_accumulator(lengths_test,
    'Delta',rep(0,5),rep(width,5))
  idx <- which(lengths_test$Method=='Delta' & lengths_test$Risk_Measure=='WT')
  stopifnot(identical(lengths_test$Finite_Lengths[[idx]],c(1,2,1000)),
    median(lengths_test$Finite_Lengths[[idx]])==2,
    lengths_test$Largest_WT[[idx]]$Finite_Length[1]==1000)
  cat('Interval classification, retained lengths and largest-WT checks passed.\n')
}
# CONFIGURATION AND CHECKPOINT PROVENANCE
# Store settings and the discoverable script source with each state. Exact source
# comparison means even a comment-only edit makes an existing checkpoint
# incompatible; use a new output folder/version when continuing from edited source.
source_text <- if(length(script)==1 && file.exists(script)) readLines(script,warn=FALSE) else NULL
configuration <- list(version=version,MODE=MODE,N_MC=N_MC,B=B,n=sample_sizes,
  data_distribution='Pareto I',fitted_distribution='Pareto I',
  observation_scheme='left truncation only',
  x0=x0,alpha=alpha,d=d,confidence=confidence,
  beta=beta,delta=delta,r=r,lambda=lambda,bootstrap_block=bootstrap_block,
  bootstrap_required_fraction=bootstrap_required_fraction,seed=seed,source=source_text,
  delta_information=delta_information,outer_configuration=outer_configuration,
  bootstrap_seed=bootstrap_seed,RNGkind=outer_rng_kind)
saveRDS(configuration,file.path(outdir,'Simulation_Configuration.rds'))
# Write the new RDS state to a temporary file before moving the old state to
# backup. If promotion fails, try to restore the prior checkpoint. This reduces
# risk of losing the last saved state; it is not a guarantee against every crash.
save_checkpoint <- function(state,path) {
  tmp <- paste0(path,'.tmp'); backup <- paste0(path,'.bak'); saveRDS(state,tmp)
  if(file.exists(backup)) unlink(backup)
  if(file.exists(path) && !file.rename(path,backup)) stop('Cannot back up checkpoint')
  if(!file.rename(tmp,path)) {
    if(file.exists(backup)) file.rename(backup,path)
    stop('Cannot replace checkpoint')
  }
  if(file.exists(backup)) unlink(backup)
}
# RUN ONE SAMPLE-SIZE CELL
# The returned data frame contains scalar publication and audit columns only.
# The checkpoint retains the richer accumulator, lengths, extremes and counts.
run_cell <- function(n,cell_id) {
  path <- file.path(outdir,sprintf('checkpoint_ParetoI_n%d.rds',n)); backup <- paste0(path,'.bak')
  # done is the last processed outer replication, including failed fits.
  # Bootstrap counters count refits, whereas fit_failures counts outer datasets.
  # Numerical counters belong to this cell rather than accumulating across n.
  state <- list(done=0L,fit_failures=0L,
    bootstrap_attempted=0,bootstrap_valid=0,bootstrap_failed=0,
    numerical_counts=setNames(rep(0,length(diag_env$counts)),names(diag_env$counts)),
    acc=new_accumulator(),configuration=configuration)
  # Prefer the current checkpoint, using backup only if the current file is absent.
  # Configuration equality protects against mixing different designs or source.
  # An unreadable existing checkpoint is not automatically replaced by the backup.
  if(resume && (file.exists(path) || file.exists(backup))) {
    state <- readRDS(if(file.exists(path)) path else backup)
    if(!identical(state$configuration,configuration)) stop('Checkpoint configuration differs')
  }
  diag_env$counts <- state$numerical_counts
  # Skip the loop if the cell is already complete. Otherwise restart at done+1;
  # per-replication seeds make results independent of the resumed RNG state.
  if(state$done<N_MC) for(replication in seq.int(state$done+1L,N_MC)) {
    row <- outer_cache$cells[[cell_id]][replication,]
    s <- list(n=n,T=unname(row['T']))
    if(!is.finite(s$T) || s$T<=0 || !is.finite(row['ahat']) || row['ahat']<=0) {
      state$fit_failures <- state$fit_failures+1L
    } else {
      ahat <- unname(row['ahat'])
      # Attach metadata for WT extreme records; all three methods use this same fit.
      # The attribute is checkpointed along with the accumulator.
      attr(state$acc,'replication_context') <- c(replication=replication,ahat=ahat)
      point <- risk_all(ahat,gradient=TRUE)
      # Table coefficient is the asymptotic variance for sqrt(n)*(Rhat-R).
      # Divide by n; all observations are exact, so probability=1.
      tab <- delta_table_components(ahat,point)
      se <- sqrt(tab$coefficient/n)
      # The Delta interval is symmetric on the risk scale, with no positivity clipping.
      # A negative finite lower endpoint is accepted. Infinite fitted risks are flagged
      # as theoretical Delta boundary failures; missing values/SEs are numerical failures.
      lower <- point$value-zcrit*se; upper <- point$value+zcrit*se
      bad <- !is.finite(point$value) | !is.finite(se)
      lower[bad] <- upper[bad] <- NA_real_
      state$acc <- update_accumulator(state$acc,'Delta',lower,upper,
                                      boundary_violation=is.infinite(point$value))
      # Transform the alpha LR set through each decreasing risk functional.
      # If an endpoint crosses a finite-risk boundary, its transformed risk may be Inf;
      # that interval is excluded from finite-valid coverage and counted separately.
      ratios <- lr_ratios(n)
      # Decreasing functionals: high alpha gives the lower risk endpoint.
      state$acc <- update_accumulator(state$acc,'Profile_likelihood',
        risk_all(ahat*ratios[2])$value,risk_all(ahat*ratios[1])$value)
      # A deterministic bootstrap-only seed makes resume independent of RNG state.
      # Outer statistics were already generated and cannot be changed by this.
      set.seed(bootstrap_seed+cell_id*1000000L+replication)
      boot <- bootstrap_alpha(ahat,n)
      # Count bootstrap alpha successes separately from interval availability:
      # even successful alpha refits may produce infinite risk percentile endpoints.
      state$bootstrap_attempted <- state$bootstrap_attempted+B
      state$bootstrap_valid <- state$bootstrap_valid+boot$valid
      state$bootstrap_failed <- state$bootstrap_failed+boot$failed
      state$acc <- update_accumulator(state$acc,'Bootstrap_percentile',
        risk_all(boot$alpha_for_risk_bounds[1])$value,
        risk_all(boot$alpha_for_risk_bounds[2])$value)
    }
    # Mark the replication complete only after its fit/interval updates.
    # Checkpoint at the configured cadence and at the final replication so all
    # finished cells can be reconstructed without repeating bootstrap calculations.
    state$done <- replication
    if(replication %% checkpoint_every==0 || replication==N_MC) {
      state$numerical_counts <- diag_env$counts; save_checkpoint(state,path)
      cat('Pareto I n=',n,': ',replication,'/',N_MC,'\n',sep=''); flush.console()
    }
  }
  # FINALIZE THE CELL
  # The audit checks fitted-risk boundary exclusions against the stored outer MLEs.
  # All percentages below use completed outer replications, unless explicitly
  # labelled as conditioning on finite valid intervals.
  audit <- outer_boundary_audit(n,outer_cache$cells[[cell_id]],state$done,
                                state$acc,state$fit_failures)
  write.csv(audit,file.path(outdir,sprintf('Working_Model_Boundary_Audit_n%d.csv',n)),
            row.names=FALSE)
  a <- state$acc; a$Data_Model <- 'Pareto_I'; a$Model <- 'Pareto_I'; a$n <- n; a$Completed <- state$done
  # Failed is the complement of finite validity, so it includes infinite intervals.
  # Other_Failed groups outer fit failures and numerical interval failures only.
  # The first partition identity checks that each replication was classified once.
  a$Failed <- state$done-a$Valid; a$Fit_Failures <- state$fit_failures
  a$Other_Failed <- a$Fit_Failures+a$Numerical_CI_Failures
  a$Boundary_Exclusion_Rate <- a$Boundary_Excluded/state$done
  a$Infinite_CI_Exclusion_Rate <- a$Infinite_CI_Excluded/state$done
  stopifnot(all(a$Completed==a$Valid+a$Boundary_Excluded+a$Other_Failed))
  # Coverage_Success_Over_All is a finite-interval success rate over all attempts:
  # excluded/infinite/failed intervals count as noncoverage by this convention.
  # It does not evaluate containment by infinite intervals.
  # Coverage_Among_Finite_Valid is conditional empirical coverage, Covered/Valid.
  # Its binomial plug-in MCSE uses Valid, the observed subset size, not Completed;
  # when Valid=0 coverage and MCSE are NA. A boundary estimate 0 or 1 has MCSE 0,
  # which does not imply absence of uncertainty.
  a$Coverage_Success_Over_All <- a$Covered/state$done
  a$Coverage_Among_Finite_Valid <- ifelse(a$Valid>0,a$Covered/a$Valid,NA_real_)
  a$MCSE_Coverage_Among_Finite_Valid <- sqrt(a$Coverage_Among_Finite_Valid*(1-a$Coverage_Among_Finite_Valid)/a$Valid)
  # Keep legacy mean columns for compatibility. Recompute final means from
  # retained values below using mean(), avoiding avoidable overflow in a raw sum.
  # Sum_Finite_Length remains an audit total and can overflow for extreme values.
  a$Mean_Length_Valid <- ifelse(a$Valid>0,a$Sum_Finite_Length/a$Valid,NA_real_)
  a$Mean_Finite_Length <- ifelse(a$Finite_Length_Count>0,a$Sum_Finite_Length/a$Finite_Length_Count,NA_real_)
  # Apply each statistic to exactly the retained finite-valid lengths.
  # Return NA for empty subsets. The median measures typical precision; mean
  # measures average precision and can be driven by rare enormous WT intervals.
  # P95/P99 use R type-7 interpolation, unlike type-1 bootstrap endpoint percentiles.
  length_stat <- function(fun) vapply(a$Finite_Lengths,function(x)
    if(length(x)) as.numeric(fun(x)) else NA_real_,numeric(1))
  a$Mean_Finite_Length <- length_stat(mean)
  a$Mean_Length_Valid <- a$Mean_Finite_Length
  a$Median_Finite_Length <- length_stat(median)
  a$P95_Finite_Length <- length_stat(function(x) quantile(x,.95,type=7,names=FALSE))
  a$P99_Finite_Length <- length_stat(function(x) quantile(x,.99,type=7,names=FALSE))
  a$Maximum_Finite_Length <- length_stat(max)
  # These two failure percentages have denominator Completed and are disjoint.
  # Undefined/numerical combines undefined Delta boundaries, numerical CI failures
  # and outer fit failures. Together with Valid they partition all replications.
  a$Percent_Infinite_Intervals <- 100*a$Infinite_CI_Excluded/a$Completed
  a$Undefined_Or_Numerical_Failures <- a$Delta_Boundary_Excluded+a$Other_Failed
  a$Percent_Undefined_Or_Numerical_Failures <- 100*a$Undefined_Or_Numerical_Failures/a$Completed
  stopifnot(all(a$Completed==a$Valid+a$Infinite_CI_Excluded+a$Undefined_Or_Numerical_Failures))
  # Export the retained WT records with sample size, method and rank.
  # An explicit empty schema still produces a readable CSV if no finite WT CI exists.
  wt_rows <- which(a$Risk_Measure=='WT')
  largest_wt <- do.call(rbind,lapply(wt_rows,function(i) {
    records <- a$Largest_WT[[i]]
    if(is.null(records)) return(NULL)
    cbind(data.frame(n=n,Method=a$Method[i],Rank=seq_len(nrow(records))),records)
  }))
  if(is.null(largest_wt)) largest_wt <- data.frame(n=integer(),Method=character(),
    Rank=integer(),Replication=integer(),Alpha_Hat=double(),Alpha_Distance_From_One=double(),
    Lower=double(),Upper=double(),Finite_Length=double(),Covered=logical())
  write.csv(largest_wt,file.path(outdir,sprintf('Largest_Finite_WT_Intervals_n%d.csv',n)),row.names=FALSE)
  # Remove nested columns only from this output copy: CSV requires scalar cells.
  # The checkpoint state still contains exact lengths and largest-WT records.
  # Mean_Length_All is reported only when every attempted interval was finite valid;
  # it is otherwise NA, not an average silently computed after exclusions.
  a$Finite_Lengths <- NULL; a$Largest_WT <- NULL
  a$Mean_Length_All <- ifelse(a$Failed==0,a$Mean_Length_Valid,NA_real_)
  a$Mean_Length_Valid_Thousand <- a$Mean_Length_Valid/1000
  write.csv(a,file.path(outdir,sprintf('Summary_ParetoI_n%d.csv',n)),row.names=FALSE,na='NA')
  # Export outer-fit and bootstrap-refit counts with numerical event counters.
  # These event totals need not equal the per-interval failure counts in the summary.
  diagnostics <- data.frame(n=n,Completed=state$done,Outer_Fit_Failures=state$fit_failures,
    Bootstrap_Attempted=state$bootstrap_attempted,Bootstrap_Valid=state$bootstrap_valid,
    Bootstrap_Failed=state$bootstrap_failed)
  diagnostics <- cbind(diagnostics,as.data.frame(as.list(diag_env$counts)))
  write.csv(diagnostics,file.path(outdir,sprintf('Diagnostics_ParetoI_n%d.csv',n)),row.names=FALSE)
  a
}

# COMBINED OUTPUTS
# Run sample sizes sequentially and refresh the combined summary after each cell.
# Thus completed cells have CSV output even if a later sample size is interrupted.
summaries <- vector('list',length(sample_sizes))
for(i in seq_along(sample_sizes)) {
  summaries[[i]] <- run_cell(sample_sizes[i],i)
  write.csv(do.call(rbind,summaries[seq_len(i)]),
    file.path(outdir,'Monte_Carlo_Summary.csv'),row.names=FALSE,na='NA')
}
summary <- do.call(rbind,summaries)
# The main publication table pairs conditional finite-valid coverage and its MCSE
# with median finite precision and the percentages of infinite/undefined intervals.
# Completed and Valid expose the denominators needed to interpret these results.
main_columns <- c('Data_Model','Model','n','Risk_Measure','Method','Completed','Valid',
  'Coverage_Among_Finite_Valid','MCSE_Coverage_Among_Finite_Valid','Median_Finite_Length',
  'Percent_Infinite_Intervals','Percent_Undefined_Or_Numerical_Failures')
# The supplement describes tail sensitivity of finite length distributions and
# retains detailed exclusion counts. Inspect largest-WT CSVs alongside these tails;
# a single maximum alone cannot explain why an interval became enormous.
supplementary_columns <- c('Data_Model','Model','n','Risk_Measure','Method','Completed',
  'Finite_Length_Count','Mean_Finite_Length','P95_Finite_Length','P99_Finite_Length',
  'Maximum_Finite_Length','Covered','Valid','Failed','Infinite_CI_Excluded',
  'Delta_Boundary_Excluded','Boundary_Excluded','Fit_Failures','Numerical_CI_Failures',
  'Undefined_Or_Numerical_Failures','Coverage_Success_Over_All')
write.csv(summary[,main_columns],file.path(outdir,'Main_Results.csv'),row.names=FALSE,na='NA')
write.csv(summary[,supplementary_columns],file.path(outdir,'Supplementary_Results.csv'),row.names=FALSE,na='NA')

# Write a self-contained results README with design, interpretation and denominator
# conventions. It accompanies the exported CSVs so coverage cannot be mistaken
# for unconditional coverage of every generated interval.
writeLines(c(
  'CORRECT SPECIFICATION: Pareto I truth and fitted model; fixed x0.',
  sprintf('Truth: x0=%g, alpha=%g; d=%g; no right censoring.',x0,alpha,d),
  'Risk coverage targets true ground-up Pareto I risks.',
  'n counts fully observed claims conditional on X>=d; no upper cap.',
  'Pareto conditional likelihood: alpha=n/T, T=sum(log(X/d)).',
  'The fixed Pareto x0 cancels from conditional likelihood; it enters ground-up extrapolation.',
  'Delta uses working-model expected Fisher information, without a sandwich correction.',
  'Without censoring, I_alpha=1/alpha^2 and I_log(alpha)=1 per observation.',
  'Delta SE=abs(dR/dlog(alpha))/sqrt(n); LR uses n exact observations.',
  'Parametric bootstrap generates fitted Pareto I left-truncated-only samples.',
  'Profile likelihood uses ordinary chi-square(1) calibration without misspecification correction.',
  'Finite-sample calibration error can affect coverage of true Pareto I risks.',
  'Outer samples and bootstrap draws use separate deterministic per-replication seeds.',
  'Finite positive alpha bootstrap refits are required at the configured fraction; failures are not redrawn.',
  'Bootstrap percentiles use type=1 decreasing-order-statistic mapping.',
  'CTE, GS and WT are infinite for alpha<=1; PHT is infinite for r*alpha<=1.',
  'Final CIs with infinite endpoints are excluded and counted separately.',
  'Delta intervals are unavailable when the fitted risk is infinite.',
  'Primary: Coverage_Among_Finite_Valid = coverage among finite valid intervals; denominator is Valid.',
  'MCSE_Coverage_Among_Finite_Valid = sqrt(p*(1-p)/Valid); NA if Valid=0.',
  'Median finite interval length describes typical precision; mean finite length describes average precision.',
  'Extreme WT lengths can dominate the mean; neither measure replaces the other.',
  'Main_Results.csv reports finite-valid coverage with MCSE, median finite length and failure percentages.',
  'Infinite and undefined/numerical percentages use Completed, including outer fit failures.',
  'Supplementary_Results.csv reports mean/P95/P99/maximum finite lengths and detailed exclusion counts.',
  'P95/P99 use R quantile type=7; all length statistics condition on finite valid intervals.',
  'Largest_Finite_WT_Intervals_n*.csv retains the ten largest finite WT intervals per method and n.',
  'Those diagnostics include replication, fitted alpha, distance from one, endpoints, length and coverage.',
  'Exact finite lengths and largest-WT diagnostics are checkpointed; v5 checkpoints cannot be reused.',
  'Coverage_Success_Over_All treats excluded or failed intervals as noncoverage.',
  'Completed=Valid+Boundary_Excluded+Other_Failed.',
  'Other_Failed=Fit_Failures+Numerical_CI_Failures.',
  'Mean_Length_Valid_Thousand=Mean_Length_Valid/1000.',
  'Working_Model_Boundary_Audit checks exclusions, including invalid fitted shapes.',
  'The derivative/Fisher-information checks validate working-model formulas only.',
  'No production simulation was executed during authoring; run validation locally first.',
  sprintf('N_MC=%d per n; B=%d; confidence=%.2f; r=%.2f; lambda=%.2f.',
          N_MC,B,confidence,r,lambda)
),file.path(outdir,'README_Results.txt'))
# Display the compact main table in the console; detailed outputs are on disk.
print(summary[,main_columns],row.names=FALSE)
cat('Results saved in:',outdir,'\n')
