# PARETO I ONLY: left-truncated/right-censored Monte Carlo study.
# Ground-up S(x)=(x0/x)^alpha, x>=x0, fixed x0.
# n counts observed conditional claims above d, including censored claims.
# Closed-form MLE, vectorized bootstrap, cached LR ratios.
# Bootstrap risk endpoints use decreasing-order statistics exactly (type=1).
# Delta CIs use the table variance coefficient divided by n:
# expected I_alpha=(1-(d/censor_at)^alpha)/alpha^2 per observation.
# Outer samples are generated in thousand-unit scale and cached before fitting.
# Bootstrap draws are isolated from the outer RNG; checkpoints can be resumed.
# Run in a fresh R session; the RNG kind is recorded for reproducibility.
# Run MODE='validation' to check the pipeline before a production run.

# -----------------------------------------------------------------------------
# User settings: edit parameter values only in this section.
# -----------------------------------------------------------------------------
# Execution guide: choose validation for a small pipeline smoke check, then
# production for the planned Monte Carlo size. Validation is too small for
# publication estimates. N_MC is the number of independent outer datasets;
# B is the number of parametric bootstrap datasets within each successful fit.
MODE <- "production"  # Use "validation" for an initial pipeline check.
mode_settings <- list(
  production = list(N_MC = 10000L, B = 500L),
  validation = list(N_MC = 10L, B = 30L)
)

# Sampling settings. Generate and evaluate N_MC datasets per sample size.
# Each sample size defines a separate simulation cell. The observed sample
# is conditional on X >= d: losses below d are absent, not zero-valued records.
# Losses at or above censor_at contribute a right-censored observation.
# unit_scale changes generation arithmetic only; reported risks use x0's units.
sample_sizes <- c(50L, 100L, 500L)
d <- 4000
censor_at <- 14000
unit_scale <- 1000

# Ground-up model and risk parameters (x0 is in the reporting scale).
# x0 is the known lower support of the ground-up Pareto distribution;
# alpha_true controls tail heaviness (smaller alpha means a heavier tail).
# confidence is the nominal two-sided CI level. beta is the upper-tail
# probability for VaR/CTE/GS; delta modifies GS; r is the proportional-hazards
# power; lambda shifts normal quantiles in the Wang distortion (WT).
# These risk targets describe the ground-up distribution, although fitting
# uses only the truncated/censored observations.
x0 <- 1000
alpha_true <- 2
confidence <- 0.90
beta <- 0.10
delta <- 0.25
r <- 0.75
lambda <- 0.50

# Bootstrap controls. A fraction of 1 requires all B refits to be valid;
# for example, 0.979 would require at least 979 valid refits when B = 1000.
# bootstrap_block limits the temporary n-by-block matrix size. It affects
# memory use, while B determines the requested number of bootstrap refits.
# The validity threshold concerns alpha estimates, not finite risk values:
# a positive bootstrap alpha can still imply a theoretically infinite risk.
bootstrap_required_fraction <- 1
bootstrap_block <- 250L

# Reproducibility and checkpoint controls.
# seed controls only outer datasets; bootstrap_seed starts a separate,
# deterministic seed schedule indexed by sample-size cell and replication.
# checkpoint_every controls saved progress. resume loads compatible progress;
# RUN_CHECKS enables mathematical and interval-classification sanity checks.
seed <- 200L
bootstrap_seed <- 20261023L
checkpoint_every <- 100L
resume <- TRUE
RUN_CHECKS <- TRUE

# -----------------------------------------------------------------------------
# Derived settings and study metadata: no repeated parameter values to edit.
# -----------------------------------------------------------------------------
# Validate the mode before indexing its settings. Record RNGkind so that
# reproducibility includes the generator, normal generator and sampling kind.
# The version label separates this summary/checkpoint schema from older runs.
stopifnot(length(MODE) == 1L, MODE %in% names(mode_settings))
N_MC <- mode_settings[[MODE]]$N_MC
B <- mode_settings[[MODE]]$B
scaled_d <- d / unit_scale
scaled_u <- censor_at / unit_scale
outer_rng_kind <- RNGkind()
version <- "ParetoI_expected_information_v6_publication"
delta_information <- "expected_Fisher_per_observation"
measures <- c("VaR", "CTE", "GS", "PHT", "WT")
methods <- c("Delta", "Bootstrap_percentile", "Profile_likelihood")

# Locate the current script under Rscript, source(), or RStudio, in that
# order. If none supplies a usable path, write results beneath getwd().
# Output folders include mode and selected parameters; full configuration
# checks below prevent silently reusing caches with different settings.
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
outdir <- file.path(folder,sprintf('MC_ParetoI_%s_r%.2f_lambda%.2f_%s',MODE,r,lambda,version))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
# Reject settings outside the mathematical or computational domains before
# starting a potentially expensive simulation. The true risk vector is also
# checked below because publication coverage requires finite target values.
stopifnot(MODE %in% c('production','validation'),N_MC>0,B>=20,
          all(sample_sizes>=2),alpha_true>0,x0>0,d>x0,censor_at>d,
          beta>0,beta<1,delta>=0,delta<.5,r>0,r<=1,lambda>0,
          confidence>0,confidence<1,bootstrap_block>=1,
          unit_scale>0,checkpoint_every>=1,
          bootstrap_required_fraction>0,bootstrap_required_fraction<=1)
# For Y=log(X/d) conditional on X>=d, Y is exponential with rate alpha.
# H is its censoring time. zcrit is the central normal critical value;
# cutoff is the asymptotic one-parameter likelihood-ratio critical value.
H <- log(censor_at/d)
zcrit <- qnorm((1+confidence)/2)
cutoff <- qchisq(confidence,1)

# Numerical diagnostics, checkpointed with the streaming summaries.
# An environment lets nested numerical helpers increment counters without
# passing them through every return value. Counts are numerical events, not
# necessarily distinct replications: one fit can trigger several evaluations.
# Per-risk/per-method exclusion counts are stored separately in the accumulator.
diag_env <- new.env(parent=emptyenv())
diag_env$counts <- c(WT_Integration_Failures=0,WT_Finite_Overflow=0,
                     Risk_Finite_Overflow=0,CI_Boundary_Exclusions=0,CI_Infinite_Exclusions=0)
count_diagnostic <- function(name) {
  diag_env$counts[name] <- diag_env$counts[name]+1
}

# Stable log upper inverse Mills ratio log(phi(z)/S(z)).
log_mills <- function(z) {
  # For large positive z, use an asymptotic expansion instead of forming a
  # ratio of tiny tail probabilities. For other z, subtract log probabilities.
  # Both branches compute the same log(phi(z) / P(Z>z)) quantity.
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
    # q is positive only on the finite WT domain. Completing the square gives
    # an integral centered at lambda/q with scale 1/sqrt(q). Subtracting shift
    # inside the exponent and restoring it in logvalue limits intermediate size.
    q <- 1-1/a; center <- lambda/q; scale <- 1/sqrt(q)
    shift <- log_mills(center)/a
    integrand <- function(t) exp(-t^2/2+log_mills(center+scale*t)/a-shift)
    # Adaptive quadrature evaluates the transformed integral on the real line.
    # As alpha approaches 1 from above, finite WT values can become enormous;
    # assemble their logarithms before checking the floating-point range.
    I <- integrate(integrand,-Inf,Inf,rel.tol=1e-8,subdivisions=1000)$value
    logvalue <- log(x0)+lambda^2/(2*q)-lambda^2/2-
      q*log(2*pi)/2-log(q)/2+shift+log(I)
    if(!is.finite(logvalue) || logvalue>log(.Machine$double.xmax)) {
      count_diagnostic('WT_Finite_Overflow')
      return(list(value=NA_real_,derivative=NA_real_))
    }
    value <- exp(logvalue)
    derivative <- NA_real_
    # Only delta intervals need the derivative. The additional integral gives
    # dWT/deta, eta=log(alpha); its negative sign reflects decreasing risk with
    # increasing shape. Nonfinite derivatives make the delta SE unavailable.
    if(gradient) {
      J <- integrate(function(t) {
        integrand(t)*(-pnorm(center+scale*t,lower.tail=FALSE,log.p=TRUE))/a
      },-Inf,Inf,rel.tol=1e-8,subdivisions=1000)$value
      derivative <- -value*J/I
      if(!is.finite(derivative)) derivative <- NA_real_
    }
    list(value=value,derivative=derivative)
  # A quadrature error is a numerical failure, returned as NA. It must not
  # be reclassified as theoretical infinity, which has a distinct scientific
  # meaning and is counted separately in the results.
  },error=function(e) {
    count_diagnostic('WT_Integration_Failures')
    list(value=NA_real_,derivative=NA_real_)
  })
}
# Risk vector and analytic derivatives with respect to eta=log(alpha).
risk_all <- function(a,gradient=FALSE) {
  if(!is.finite(a) || a<=0)
    return(list(value=setNames(rep(NA_real_,5),measures),derivative=rep(NA_real_,5)))
  # Initialize risks at theoretical infinity and overwrite only where their
  # finite-domain formulas apply. Derivatives are with respect to log(alpha),
  # so dR/deta = alpha*dR/dalpha. Named vectors keep measure order consistent.
  value <- setNames(rep(Inf,5),measures)
  derivative <- setNames(rep(NA_real_,5),measures)
  # VaR is the (1-beta) quantile. For a>1, CTE is the tail conditional mean
  # above VaR and GS multiplies CTE by its specified loading factor.
  # PHT is finite when r*a>1; WT is computed by numerical quadrature.
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
# Delta-method variance components.
# p(a)=P(exact | X>=d)=1-(d/censor_at)^a.
# Per-observation I_alpha=p(a)/a^2; the full variance divides by n.
# D terms are chosen with positive signs; the table uses their squares.
# WT: C=a*(WT/x0-1), D=C-a*dC/da=-a*(dWT/dlog(a))/x0.
delta_table_components <- function(a,point=NULL) {
  if(is.null(point)) point <- risk_all(a,gradient=TRUE)
  # Compute 1-exp(a*log(d/censor_at)) with expm1 to reduce cancellation.
  # The returned coefficient is the variance of sqrt(n)*(Rhat-R), not the
  # variance of Rhat itself. Callers divide this coefficient by n.
  probability <- -expm1(a*log(d/censor_at))
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
  # Overflow or an unavailable derivative leaves the coefficient undefined.
  # Returning NA avoids reporting a numerical overflow as an infinite-risk CI.
  coefficient[!is.finite(coefficient)] <- NA_real_
  list(probability=probability,coefficient=coefficient,
       D_CTE=unname(D_CTE),D_GS=unname(D_GS),C_WT=unname(C_WT),D_WT=unname(D_WT))
}
# Independent survival-integral definitions for checking the WT table constants.
# Let t=exp(-v). C=integral exp(v/a)*g(exp(-v)) dv; D adds factor 1+v/a.
# This audit is performed at the true shape in startup checks, not each fit.
wang_table_constants_check <- function(a) {
  stopifnot(a>1)
  # Use log probabilities for the distorted survival integrand. This is an
  # independent integral representation, so agreement with the main WT helper
  # checks more than repeated use of the same transformation.
  logbase <- function(v) v/a+pnorm(qnorm(-v,log.p=TRUE)+lambda,log.p=TRUE)
  C <- integrate(function(v) exp(logbase(v)),0,Inf,
                 rel.tol=1e-8,subdivisions=1000)$value
  D <- integrate(function(v) exp(logbase(v)+log1p(v/a)),0,Inf,
                 rel.tol=1e-8,subdivisions=1000)$value
  c(C_WT=C,D_WT=D)
}
# Evaluate the fixed population targets once. Coverage comparisons always
# use these ground-up targets, never a conditional-on-d risk or a fitted risk.
# The audit CSVs expose true risks, censoring probability, variance coefficients,
# and sample-size-specific asymptotic SEs for independent inspection.
true_point <- risk_all(alpha_true,gradient=TRUE)
true_risks <- true_point$value
if(any(!is.finite(true_risks))) stop('This evaluation requires finite true risks')
write.csv(data.frame(Risk_Measure=measures,True_Risk=unname(true_risks)),
          file.path(outdir,'True_Ground_Up_Risks.csv'),row.names=FALSE)
write.csv(data.frame(Conditional_Censoring_Probability=(d/censor_at)^alpha_true),
          file.path(outdir,'True_Censoring_Probability.csv'),row.names=FALSE)
true_table <- delta_table_components(alpha_true,true_point)
table_audit <- data.frame(Risk_Measure=measures,True_Risk=unname(true_risks),
  Uncensored_Probability=true_table$probability,
  Variance_Coefficient=unname(true_table$coefficient))
for(nn in sample_sizes) table_audit[[paste0('Asymptotic_SE_n',nn)]] <-
  sqrt(unname(true_table$coefficient)/nn)
write.csv(table_audit,file.path(outdir,'True_Risks_And_Delta_Table.csv'),row.names=FALSE)
write.csv(data.frame(Alpha=alpha_true,D_CTE=true_table$D_CTE,D_GS=true_table$D_GS,
  C_WT=true_table$C_WT,D_WT=true_table$D_WT),
  file.path(outdir,'Table_Auxiliary_Constants.csv'),row.names=FALSE)


# Precompute the sufficient statistics BEFORE any bootstrap draws.
# Generate N_MC datasets per sample size using a single outer seed.
# Sampling uses thousand-unit arithmetic. Each mode has its own output folder.
# The cache is independent of bootstrap settings and CI computations.
# Only quantities determining the outer sufficient statistics belong in
# this cache key. CI settings do not change the generated outer datasets.
# The stricter checkpoint key later also includes bootstrap and source details.
outer_configuration <- list(seed=seed,RNGkind=outer_rng_kind,
  alpha=alpha_true,n=sample_sizes,N_MC=N_MC,
  d=scaled_d,u=scaled_u)
outer_path <- file.path(outdir,'Outer_Statistics.rds')
if(file.exists(outer_path)) {
  outer_cache <- readRDS(outer_path)
  if(!identical(outer_cache$configuration,outer_configuration))
    stop('Outer sample cache configuration differs; use a new output folder')
} else {
  # Generate every outer dataset before any bootstrap RNG calls. Caching the
  # sufficient statistics rather than all observations saves memory and ensures
  # all three CI methods are evaluated on exactly the same outer datasets.
  set.seed(seed)
  outer_cells <- vector('list',length(sample_sizes))
  for(ii in seq_along(sample_sizes)) {
    nn <- sample_sizes[ii]
    stats <- matrix(NA_real_,N_MC,3L,
                    dimnames=list(NULL,c('m','T','ahat')))
    for(rr in seq_len(N_MC)) {
      # Inverse survival sampling produces conditional Pareto losses above d.
      # pmin implements right censoring. m counts exact observations; T includes
      # both exact and censored log-exposures. Ignoring constants, the conditional
      # log likelihood is m*log(alpha)-alpha*T, giving alpha_hat=m/T when m>0.
      # Scaling d and censor_at together leaves log ratios and alpha_hat unchanged.
      b <- sort(pmin(scaled_d*runif(nn)^(-1/alpha_true),scaled_u))
      m <- sum(b<scaled_u)
      T <- sum(log(b/scaled_d))
      stats[rr,] <- c(m,T,m/T)
      if(rr %% 10000L==0L) {
        cat('Preparing outer samples: n=',nn,', ',rr,'/',
            N_MC,'\n',sep=''); flush.console()
      }
    }
    outer_cells[[ii]] <- stats
  }
  outer_cache <- list(configuration=outer_configuration,cells=outer_cells)
  saveRDS(outer_cache,outer_path)
}
# Audit invalid/nonpositive shapes and risk-domain violations.
# Failed fits remain separate from undefined delta intervals.
outer_boundary_audit <- function(n,cell,completed,acc,fit_failures) {
  # Use only completed replications when resuming. This audit reports CTE,
  # GS and PHT domains; WT shares the a>1 domain but is not a separate row here.
  # Invalid shapes are counted as failed outer fits, not delta boundary exclusions.
  aa <- cell[seq_len(completed),'ahat']
  invalid <- !is.finite(aa) | aa<=0
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
# For q=alpha/alpha_hat, twice the log-likelihood drop is
# 2*m*(q-1-log(q)). It depends on m alone, allowing roots to be reused across
# replications. Solve in t=log(q) to keep both alpha endpoints positive.
lr_cache <- new.env(parent=emptyenv())
lr_ratios <- function(m) {
  if(m<=0) return(c(NA_real_,NA_real_))
  key <- as.character(m)
  if(exists(key,envir=lr_cache,inherits=FALSE)) return(get(key,envir=lr_cache))
  f <- function(t) 2*m*(expm1(t)-t)-cutoff
  # Expand brackets until each contains one root on its side of zero.
  # Exponentiating the roots returns multiplicative lower/upper alpha limits;
  # these are subsequently mapped through decreasing risk functionals.
  left <- -1; right <- 1
  while(f(left)<0) left <- left*2
  while(f(right)<0) right <- right*2
  ans <- exp(c(uniroot(f,c(left,0),tol=1e-10)$root,
               uniroot(f,c(0,right),tol=1e-10)$root))
  assign(key,ans,envir=lr_cache)
  ans
}
# Vectorized bootstrap alpha estimates. Every sample repeats truncation and
# censoring. At least ceiling(B*fraction) valid alpha refits are required.
bootstrap_alpha <- function(a,n) {
  estimates <- rep(NA_real_,B)
  for(first in seq.int(1L,B,by=bootstrap_block)) {
    last <- min(B,first+bootstrap_block-1L); nb <- last-first+1L
    # Columns are independent bootstrap datasets for log(X/d). Censor at H,
    # then refit using the same m/T likelihood as in the outer study. Blocking
    # avoids allocating all n*B observations simultaneously.
    e <- matrix(rexp(n*nb,rate=a),nrow=n,ncol=nb)
    m <- colSums(e<H); T <- colSums(pmin(e,H))
    good <- m>0 & is.finite(T) & T>0
    values <- rep(NA_real_,nb); values[good] <- m[good]/T[good]
    estimates[first:last] <- values
  }
  good <- is.finite(estimates) & estimates>0
  # Leave endpoint shapes as NA if too few refits pass the threshold.
  # Otherwise choose empirical order statistics without interpolation (type=1).
  # Only the two selected shapes need expensive risk/WT evaluations.
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
# Startup checks cover LR cutoff equations, derivatives on the log-shape
# scale, reversed bootstrap ordering, and two independent delta/WT formulas.
# They are numerical sanity checks, not evidence of nominal finite-sample coverage.
run_checks <- function() {
  m <- 20; ahat <- 2; T <- m/ahat
  ll <- function(a) m*log(a)-a*T
  ratios <- lr_ratios(m)
  stopifnot(all(vapply(ratios,function(q)
    abs(2*(ll(ahat)-ll(ahat*q))-cutoff)<1e-6,logical(1))))
  # Multiplicative perturbations alpha*exp(+/-eps) approximate derivatives
  # with respect to eta=log(alpha), matching the analytical derivative scale.
  eps <- 1e-5; point <- risk_all(ahat,TRUE)
  numerical <- (risk_all(ahat*exp(eps))$value-risk_all(ahat*exp(-eps))$value)/(2*eps)
  stopifnot(all(abs(numerical-point$derivative)/(1+abs(numerical))<1e-5))
  # Include a shape outside the finite PHT domain to verify that infinity
  # is preserved correctly by the decreasing-order-statistic mapping.
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
  ww <- wang_table_constants_check(alpha_true)
  tt <- delta_table_components(alpha_true)
  stopifnot(abs(tt$C_WT-ww['C_WT'])/(1+abs(ww['C_WT']))<1e-6,
            abs(tt$D_WT-ww['D_WT'])/(1+abs(ww['D_WT']))<1e-6)
  cat('LR, gradients, bootstrap ordering and expected-information table checks passed.\n')
}
if(RUN_CHECKS) run_checks()

# Only finite, ordered intervals enter Valid, Covered and length summaries.
# Theoretical risk-boundary exclusions are tracked separately from numerical
# failures and outer fits that did not exist.
# One row represents one risk/method combination. Valid and Covered refer
# only to finite ordered intervals with representable finite lengths.
# Boundary_Excluded combines infinite endpoint CIs and undefined delta risks;
# Numerical_CI_Failures excludes those theoretical boundary cases.
new_accumulator <- function() {
  a <- expand.grid(Risk_Measure=measures,Method=methods,stringsAsFactors=FALSE)
  a$Valid <- a$Covered <- a$Finite_Length_Count <- 0L
  a$Boundary_Excluded <- a$Infinite_CI_Excluded <- a$Delta_Boundary_Excluded <- 0L
  a$Numerical_CI_Failures <- 0L
  a$Sum_Finite_Length <- 0
  # Retain exact finite lengths and endpoints for quantiles and WT diagnostics.
  # A list-column stores each row's finite interval records. I() preserves
  # the list structure in the data frame. Exact records support exact medians,
  # tail quantiles and ranked WT diagnostics, at a memory/checkpoint-size cost.
  a$Finite_Intervals <- I(lapply(seq_len(nrow(a)),function(i)
    data.frame(Replication=integer(),Alpha_Hat=numeric(),Uncensored=numeric(),
               T=numeric(),Lower=numeric(),Upper=numeric(),Length=numeric())))
  a
}
# Inputs lower/upper follow the measures vector. boundary_violation flags
# undefined delta risks explicitly; replication, ahat, m and T are retained
# for tracing extreme intervals back to their outer dataset. Default NA
# metadata permits the synthetic startup classification checks.
update_accumulator <- function(a,method,lower,upper,
                               boundary_violation=rep(FALSE,length(measures)),
                               replication=NA_integer_,ahat=NA_real_,
                               m=NA_real_,T=NA_real_) {
  for(j in seq_along(measures)) {
    i <- which(a$Method==method & a$Risk_Measure==measures[j])
    ci <- c(lower[j],upper[j])
    # For delta CIs, an infinite computed endpoint can be arithmetic overflow.
    # A theoretical delta boundary failure is passed explicitly below.
    infinite_endpoint <- method!='Delta' && !anyNA(ci) &&
      ci[1]<=ci[2] && any(is.infinite(ci))
    # Classification is mutually exclusive: theoretical boundary first,
    # numerical failure second, finite valid interval last. An endpoint containing
    # NA or an unordered pair is not accepted as a theoretical infinite interval.
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
    # Even finite endpoints can have a difference that overflows. Such a pair
    # cannot enter finite-length summaries and is recorded as a numerical failure.
    # Raw-scale delta lower bounds are not clipped to zero in this study.
    if(anyNA(ci) || !all(is.finite(ci)) || ci[1]>ci[2] || !is.finite(ci[2]-ci[1])) {
      a$Numerical_CI_Failures[i] <- a$Numerical_CI_Failures[i]+1L
      next
    }
    # Coverage is inclusive at both endpoints. Update the count, sum and exact
    # record together so all finite-valid length and coverage denominators agree.
    len <- ci[2]-ci[1]
    a$Valid[i] <- a$Valid[i]+1L
    a$Covered[i] <- a$Covered[i]+as.integer(ci[1]<=true_risks[j] && true_risks[j]<=ci[2])
    a$Finite_Length_Count[i] <- a$Finite_Length_Count[i]+1L
    a$Sum_Finite_Length[i] <- a$Sum_Finite_Length[i]+len
    a$Finite_Intervals[[i]] <- rbind(a$Finite_Intervals[[i]],
      data.frame(Replication=replication,Alpha_Hat=ahat,Uncensored=unname(m),
                 T=unname(T),Lower=ci[1],Upper=ci[2],Length=len))
  }
  a
}
# Check boundary/NA classification without running the Monte Carlo study.
if(RUN_CHECKS) {
  # At alpha=.9 all non-VaR risks are infinite: the profile interval therefore
  # has four theoretical infinite exclusions and one finite valid interval.
  # The delta test distinguishes four undefined risks from a numerical NA case.
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
  cat('Infinite-boundary exclusion and separate numerical-failure checks passed.\n')
}
# Store the source alongside all settings as a reproducibility record.
# Because source text is part of configuration, even comment edits cause an
# existing checkpoint compatibility check to fail; use a fresh output folder
# or retain the exact original script when resuming a previous run.
source_text <- if(length(script)==1 && file.exists(script)) readLines(script,warn=FALSE) else NULL
configuration <- list(version=version,MODE=MODE,N_MC=N_MC,B=B,n=sample_sizes,
  x0=x0,alpha_true=alpha_true,d=d,censor_at=censor_at,confidence=confidence,
  beta=beta,delta=delta,r=r,lambda=lambda,bootstrap_block=bootstrap_block,
  bootstrap_required_fraction=bootstrap_required_fraction,seed=seed,source=source_text,
  delta_information=delta_information,outer_configuration=outer_configuration,
  bootstrap_seed=bootstrap_seed,unit_scale=unit_scale,RNGkind=outer_rng_kind)
saveRDS(configuration,file.path(outdir,'Simulation_Configuration.rds'))
# Write a temporary RDS before rotating the current checkpoint to .bak.
# If the final rename fails, attempt to restore the backup and stop. The
# backup is also considered at startup when the primary file is absent.
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
# Run or resume one sample-size cell and return its flat summary table.
# State contains progress, fitting/bootstrap counts, numerical event counters,
# finite interval records and configuration; all are saved together.
run_cell <- function(n,cell_id) {
  path <- file.path(outdir,sprintf('checkpoint_ParetoI_n%d.rds',n)); backup <- paste0(path,'.bak')
  state <- list(done=0L,fit_failures=0L,censored_total=0,
    bootstrap_attempted=0,bootstrap_valid=0,bootstrap_failed=0,
    numerical_counts=setNames(rep(0,length(diag_env$counts)),names(diag_env$counts)),
    acc=new_accumulator(),configuration=configuration)
  if(resume && (file.exists(path) || file.exists(backup))) {
    state <- readRDS(if(file.exists(path)) path else backup)
    if(!identical(state$configuration,configuration)) stop('Checkpoint configuration differs')
  }
  # Restore this cell's counters before continuing, preventing startup-check
  # or previous-cell events from leaking into the cell's diagnostics.
  diag_env$counts <- state$numerical_counts
  if(state$done<N_MC) for(replication in seq.int(state$done+1L,N_MC)) {
    # Fetch cached sufficient statistics; no new outer observations are drawn.
    # A dataset with no exact claims has no positive interior MLE and contributes
    # one fit failure to every risk/method combination through the final accounting.
    row <- outer_cache$cells[[cell_id]][replication,]
    s <- list(m=row['m'],T=row['T'],k=n-row['m'])
    state$censored_total <- state$censored_total+s$k
    if(s$m==0L || !is.finite(s$T) || s$T<=0) {
      state$fit_failures <- state$fit_failures+1L
    } else {
      ahat <- unname(row['ahat'])
      point <- risk_all(ahat,gradient=TRUE)
      # Table coefficient is the asymptotic variance for sqrt(n)*(Rhat-R).
      # Divide by n, not m. The fitted probability supplies the censoring term.
      tab <- delta_table_components(ahat,point)
      se <- sqrt(tab$coefficient/n)
      # Construct symmetric Wald limits on the raw risk scale. If risk or SE is
      # unavailable, pass NA endpoints while separately flagging a theoretical
      # infinite fitted risk. This separates domain failure from numerical failure.
      lower <- point$value-zcrit*se; upper <- point$value+zcrit*se
      bad <- !is.finite(point$value) | !is.finite(se)
      lower[bad] <- upper[bad] <- NA_real_
      state$acc <- update_accumulator(state$acc,'Delta',lower,upper,
                                      boundary_violation=is.infinite(point$value),
        replication=replication,ahat=ahat,m=s$m,T=s$T)
      ratios <- lr_ratios(s$m)
      # Decreasing functionals: high alpha gives the lower risk endpoint.
      state$acc <- update_accumulator(state$acc,'Profile_likelihood',
        risk_all(ahat*ratios[2])$value,risk_all(ahat*ratios[1])$value,
        replication=replication,ahat=ahat,m=s$m,T=s$T)
      # A deterministic bootstrap-only seed makes resume independent of RNG state.
      # Outer statistics were already generated and cannot be changed by this.
      set.seed(bootstrap_seed+cell_id*1000000L+replication)
      # Bootstrap diagnostics count alpha refits, so their totals are in units
      # of bootstrap datasets rather than outer intervals. The selected endpoint
      # risks may still fail even after the required fraction of refits is valid.
      boot <- bootstrap_alpha(ahat,n)
      state$bootstrap_attempted <- state$bootstrap_attempted+B
      state$bootstrap_valid <- state$bootstrap_valid+boot$valid
      state$bootstrap_failed <- state$bootstrap_failed+boot$failed
      state$acc <- update_accumulator(state$acc,'Bootstrap_percentile',
        risk_all(boot$alpha_for_risk_bounds[1])$value,
        risk_all(boot$alpha_for_risk_bounds[2])$value,
        replication=replication,ahat=ahat,m=s$m,T=s$T)
    }
    # Advance progress even for failed fits: N_MC denotes generated datasets,
    # not a target number of successful intervals. Save at the chosen frequency
    # and at the final replication so a completed cell can be resumed directly.
    state$done <- replication
    if(replication %% checkpoint_every==0 || replication==N_MC) {
      state$numerical_counts <- diag_env$counts; save_checkpoint(state,path)
      cat('Pareto I n=',n,': ',replication,'/',N_MC,'\n',sep=''); flush.console()
    }
  }
  audit <- outer_boundary_audit(n,outer_cache$cells[[cell_id]],state$done,
                                state$acc,state$fit_failures)
  write.csv(audit,file.path(outdir,sprintf('Boundary_Audit_n%d.csv',n)),
            row.names=FALSE)
  # Flatten the completed state into reported counts. Failed means every
  # interval excluded from finite-valid summaries, including theoretical infinity.
  # Fit failures are common to all rows; other exclusions depend on risk/method.
  a <- state$acc; a$Model <- 'Pareto_I'; a$n <- n; a$Completed <- state$done
  a$Failed <- state$done-a$Valid; a$Fit_Failures <- state$fit_failures
  a$Other_Failed <- a$Fit_Failures+a$Numerical_CI_Failures
  a$Boundary_Exclusion_Rate <- a$Boundary_Excluded/state$done
  a$Infinite_CI_Exclusion_Rate <- a$Infinite_CI_Excluded/state$done
  stopifnot(all(a$Completed==a$Valid+a$Boundary_Excluded+a$Other_Failed))
  # Two denominators answer different questions. Covered/Completed is the
  # success rate when every excluded interval counts as noncoverage. Covered/Valid
  # is empirical coverage conditional on obtaining a finite valid interval.
  # The latter must be labeled explicitly in publication. Its plug-in binomial
  # MCSE is sqrt(p_hat*(1-p_hat)/Valid); it is undefined if Valid is zero.
  a$Coverage_Success_Over_All <- a$Covered/state$done
  a$Coverage_Among_Finite_Valid <- ifelse(a$Valid>0,a$Covered/a$Valid,NA_real_)
  a$MCSE_Coverage_Among_Finite_Valid <- sqrt(a$Coverage_Among_Finite_Valid*(1-a$Coverage_Among_Finite_Valid)/a$Valid)
  # Both finite-mean columns use the same finite valid intervals. The legacy
  # Mean_Length_All is provided only when no interval was excluded; it does
  # not impute lengths for missing intervals. Thousand is a display rescaling.
  a$Mean_Length_Valid <- ifelse(a$Valid>0,a$Sum_Finite_Length/a$Valid,NA_real_)
  a$Mean_Finite_Length <- ifelse(a$Finite_Length_Count>0,a$Sum_Finite_Length/a$Finite_Length_Count,NA_real_)
  a$Mean_Length_All <- ifelse(a$Failed==0,a$Mean_Length_Valid,NA_real_)
  a$Mean_Length_Valid_Thousand <- a$Mean_Length_Valid/unit_scale
  # Exact conditional length summaries; type=7 is R's default quantile rule.
  length_stat <- function(fun) vapply(a$Finite_Intervals,function(x)
    if(nrow(x)>0L) unname(fun(x$Length)) else NA_real_,numeric(1))
  # Median describes typical finite precision. P95/P99 and maximum reveal
  # extreme finite lengths that can dominate the mean, particularly for WT.
  # These conditional summaries exclude theoretical infinite lengths and failures.
  a$Median_Finite_Length <- length_stat(median)
  a$P95_Finite_Length <- length_stat(function(x) quantile(x,.95,type=7))
  a$P99_Finite_Length <- length_stat(function(x) quantile(x,.99,type=7))
  a$Maximum_Finite_Length <- length_stat(max)
  # Percentages use all Completed datasets. The undefined/numerical category
  # combines delta domain failures, numerical CI failures and outer fit failures.
  # The identity below verifies the three-way partition into finite valid,
  # infinite, and undefined/numerical outcomes for every risk/method row.
  a$Percent_Infinite_Intervals <- 100*a$Infinite_CI_Excluded/a$Completed
  a$Undefined_Or_Numerical_Failures <- a$Delta_Boundary_Excluded+
    a$Numerical_CI_Failures+a$Fit_Failures
  a$Percent_Undefined_Or_Numerical_Failures <-
    100*a$Undefined_Or_Numerical_Failures/a$Completed
  stopifnot(all(a$Completed==a$Valid+a$Infinite_CI_Excluded+
                  a$Undefined_Or_Numerical_Failures))
  # Ten largest finite WT intervals per method, with outer-fit diagnostics.
  # Rank each WT method separately, retaining at most ten records. Alpha_Hat
  # and m/T help investigate tail instability; these records contain finite
  # intervals only, so numerical overflow and infinite intervals remain in counts.
  wt_diagnostics <- do.call(rbind,lapply(which(a$Risk_Measure=='WT'),function(i) {
    x <- a$Finite_Intervals[[i]]
    x <- x[head(order(x$Length,decreasing=TRUE),10L),,drop=FALSE]
    data.frame(n=rep(n,nrow(x)),Method=rep(a$Method[i],nrow(x)),x)
  }))
  write.csv(wt_diagnostics,file.path(outdir,sprintf('Largest_Finite_WT_Intervals_n%d.csv',n)),
            row.names=FALSE,na='NA')
  # Remove the nested records from the CSV-facing table after diagnostics
  # are written. The saved state still contains them for resume and exact summaries.
  # Observed censoring uses all n*Completed observed conditional claims.
  a$Finite_Intervals <- NULL
  a$Observed_Censoring_Fraction <- state$censored_total/(state$done*n)
  write.csv(a,file.path(outdir,sprintf('Summary_ParetoI_n%d.csv',n)),row.names=FALSE,na='NA')
  # This file combines cell-level fit/bootstrap counts and numerical helper
  # events. Consult summary rows for distinct interval exclusions; helper event
  # counts can exceed interval counts because endpoints may be evaluated repeatedly.
  diagnostics <- data.frame(n=n,Completed=state$done,Outer_Fit_Failures=state$fit_failures,
    Bootstrap_Attempted=state$bootstrap_attempted,Bootstrap_Valid=state$bootstrap_valid,
    Bootstrap_Failed=state$bootstrap_failed)
  diagnostics <- cbind(diagnostics,as.data.frame(as.list(diag_env$counts)))
  write.csv(diagnostics,file.path(outdir,sprintf('Diagnostics_ParetoI_n%d.csv',n)),row.names=FALSE)
  a
}

# Process sample sizes in their configured order. Rewrite the combined CSV
# after each completed cell so results already finished remain accessible if
# a later cell stops. Each cell also has its own summary and checkpoint.
summaries <- vector('list',length(sample_sizes))
for(i in seq_along(sample_sizes)) {
  summaries[[i]] <- run_cell(sample_sizes[i],i)
  write.csv(do.call(rbind,summaries[seq_len(i)]),
    file.path(outdir,'Monte_Carlo_Summary.csv'),row.names=FALSE,na='NA')
}
summary <- do.call(rbind,summaries)
# Primary publication table: conditional finite-valid coverage with MCSE,
# median finite length, and percentages of infinite and undefined outcomes.
# The full supplementary table retains mean/tail lengths and exclusion counts.
# Coverage and MCSE are proportions; columns prefixed Percent use a 0-100 scale.
main_columns <- c('Model','n','Risk_Measure','Method','Completed','Valid',
  'Coverage_Among_Finite_Valid','MCSE_Coverage_Among_Finite_Valid',
  'Median_Finite_Length','Percent_Infinite_Intervals',
  'Percent_Undefined_Or_Numerical_Failures')
write.csv(summary[,main_columns],file.path(outdir,'Main_Results.csv'),row.names=FALSE,na='NA')
write.csv(summary,file.path(outdir,'Supplementary_Results.csv'),row.names=FALSE,na='NA')
# Write a plain-text data dictionary and interpretation guide next to the
# CSV outputs. These notes make the denominator and numerical conventions
# available to readers who receive results without the simulation source.
writeLines(c(
  'Pareto I only; fixed x0; ground-up risk targets.',
  'Primary measures: coverage among finite valid intervals (with MCSE), median finite length,',
  'percentage infinite intervals, and percentage undefined or numerical failures.',
  'Median finite length describes typical precision; mean finite length describes average length.',
  'Rare enormous WT lengths can dominate the mean; neither measure replaces the other.',
  'Main_Results.csv contains the primary measures; Supplementary_Results.csv contains full counts and tails.',
  'Median/P95/P99/max use finite valid lengths only; quantiles use R type=7.',
  'Length summaries use reporting units; only columns ending Thousand divide by unit_scale.',
  'All percentages use Completed as denominator; undefined/numerical includes failed outer fits.',
  'Largest_Finite_WT_Intervals files retain the ten largest finite intervals per method.',
  'Coverage_Among_Finite_Valid is coverage among finite valid intervals, not unconditional coverage.',
  'Its binomial plug-in MCSE uses Valid as denominator and is NA when Valid=0.',
  'Outer draws use runif, a fixed seed, sample-size order and thousand-unit arithmetic.',
  'Exactly N_MC datasets are generated and evaluated per sample size.',
  'Bootstrap draws use separate per-replication seeds after the outer samples are cached.',
  'Reproducibility requires the recorded R RNG kind and unchanged simulation settings.',
  'Boundary_Audit files check Violation_Count = Delta_Boundary_Excluded + Invalid_Shape_Count.',
  'Cache and checkpoint configurations are checked before resuming.',
  sprintf('N_MC=%d per n; B=%d; confidence=%.2f; r=%.2f; lambda=%.2f.',N_MC,B,confidence,r,lambda),
  'n includes exact and right-censored records conditional on X>=d.',
  'MLE alpha=m/T, T=sum(log(w/d)); observed information for log alpha equals m.',
  'Delta uses EXPECTED Fisher information: p=1-(d/censor_at)^alpha, I_alpha=p/alpha^2 per observation.',
  'Delta variance equals the supplied table coefficient divided by n.',
  'Delta SE=abs(dR/dlog(alpha))/sqrt(n*p) at alpha_hat, not divided by sqrt(m).',
  'True_Risks_And_Delta_Table.csv audits table coefficients and SEs at the true parameters.',
  'Delta: raw-scale Wald with expected information; unavailable when fitted risk is infinite.',
  'Bootstrap: percentile type=1 quantiles; exact decreasing-order-statistic mapping.',
  'Only selected bootstrap risk endpoints are evaluated; risks are not integrated B times.',
  'Bootstrap validity fraction applies to finite interior alpha estimates.',
  'A risk interval is still invalid if a selected numerical risk endpoint fails.',
  'CTE, GS and WT require alpha>1; PHT requires r*alpha>1.',
  'CIs with theoretical infinite endpoints, including [Inf,Inf], are excluded.',
  'Infinite_CI_Excluded counts intervals with an infinite endpoint by risk and method.',
  'Delta_Boundary_Excluded counts undefined delta CIs due to an infinite fitted risk.',
  'Boundary_Excluded is the combined theoretical risk-boundary exclusion count.',
  'Numerical_CI_Failures excludes theoretical boundary exclusions and failed outer fits.',
  'Completed=Valid+Boundary_Excluded+Other_Failed; Other_Failed=Fit_Failures+Numerical_CI_Failures.',
  'Mean_Length_Valid and Coverage_Among_Finite_Valid use only finite valid CIs.',
  'Finite-domain overflow remains a numerical failure (NA), not theoretical infinity.',
  'Mean_Length_All is NA when any interval failed.',
  sprintf('Mean_Length_Valid_Thousand=Mean_Length_Valid/%.0f.',unit_scale),
  'Coverage_Among_Finite_Valid conditions on finite valid intervals; report exclusions with coverage.',
  'Coverage_Success_Over_All counts failed intervals as noncoverage.',
  'Profile likelihood uses asymptotic chi-square(1) cutoffs.',
  'Finite-sample 90% coverage is not guaranteed; no uncertainty in x0 included.'
),file.path(outdir,'README_Results.txt'))
# Display the same primary columns written to Main_Results.csv and report
# the output directory. The script creates files when executed; editing this
# source alone does not regenerate any Monte Carlo results.
print(summary[,main_columns],row.names=FALSE)
cat('Results saved in:',outdir,'\n')
