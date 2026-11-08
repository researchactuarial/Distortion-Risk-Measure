# MISSPECIFICATION STUDY: Pareto I truth, shifted lognormal working model.
# Truth: P(X>x)=(x0/x)^alpha, x>=x0. x0 is the Pareto minimum.
# Working model: X=x0+exp(Y); its shift is fixed at the same numeric x0.
# Outer samples are Pareto I conditional on X>=d, with no right censoring.
# Bootstrap samples come from the fitted conditional shifted LN model.
# Coverage targets TRUE GROUND-UP PARETO risks, not fitted-model risks.
# Delta information and LR calibration remain working-model based;
# no sandwich or misspecification correction is applied.
# X=x0+exp(Y), Y~N(mu,sigma^2), observed conditional on X>=d,
# with no right censoring. n counts fully observed conditional claims.
# Fits the finite conditional model and its exponential log-excess limit.
# Retains verified limiting bootstrap fits; never redraws failed datasets.
# Delta uses expected information, profile checks limiting likelihoods.
# This changes the estimator to an EXTENDED MLE when a finite MLE is absent.
# Boundary-aware estimation does NOT guarantee nominal bootstrap/LR coverage.
# Base R only. R execution was unavailable during authoring.
# Run: Rscript Simulation_ParetoI_Fit_ShiftedLN_LeftTruncated.R

# SECTION 1. Study settings and interpretation.
# Run this file with Rscript or source it in R. All calculations use base R;
# rstudioapi is optional and is used only to locate an open script.
# The outer experiment assesses a deliberately misspecified working model.
# Consequently, increasing sample size need not restore nominal coverage.
# A pilot uses 10 outer replications and 30 bootstrap draws; production uses
# 1000 and 500 respectively. Select the mode before starting a fresh run.
MODE <- 'production'  # Use 'production' after local pilot validation.
N_MC <- if(MODE=='production') 5000L else 10L
B <- if(MODE=='production') 500L else 30L
# Each n is the number of observed claims above the truncation threshold,
# not the number of claims generated before truncation. There is no censoring.
sample_sizes <- c(647L)
# x0 is both the true Pareto minimum and the fixed lognormal shift.
# d is the observation threshold; extrapolation below d is model dependent.
x0 <- 100000; d <- 500000
# The Pareto tail index controls heaviness. alpha>1 ensures a finite mean;
# alpha*r>1 additionally ensures that the true PHT risk is finite.
alpha <- 1.13
# Risk settings: beta is an upper-tail probability (VaR level 1-beta);
# delta controls the GS weighting; r is the proportional-hazards exponent;
# lambda is the positive Wang-transform shift. CTE averages the upper beta tail.
beta <- .10; delta <- .25; r <- .95; lambda <- .25
# Two-sided confidence level shared by Delta, percentile bootstrap and LR CIs.
confidence <- .90
# Required share of nonmissing bootstrap risk values, assessed separately per
# measure. A value of 1 requires all B draws; genuine limiting Inf values count
# as available. Failed bootstrap draws are never replaced by fresh draws.
bootstrap_required_fraction <- 1
EXCLUDE_INFINITE_CIS <- TRUE # Retain boundary bootstrap values; exclude only final infinite CIs.
# Save progress after this many outer replications. The seed is deterministic
# per sample-size cell and replication, so restarting preserves that sequence.
checkpoint_every <- 10L
seed <- 20261022L
# Resume only when the complete saved configuration matches, including source
# text. Even comment changes make earlier checkpoints incompatible; use a fresh
# output directory or archive old checkpoints before running this edited file.
resume <- TRUE
# Execute numerical consistency checks before the expensive simulation begins.
RUN_CHECKS <- TRUE
# Algorithm/output schema tag used in the output directory name.
version <- 'ParetoI_to_ShiftedLN_left_truncated_v5'
# Engineering guards: reaching them returns UNRESOLVED, not a fitted bound.
# Search limits are numerical stopping rules. They are not statistical parameter
# constraints and must not be interpreted as evidence of a finite endpoint.
max_t_bound <- 80
max_log_sigma <- 16
min_log_sigma <- -12
# Accept finite fits only with a sufficiently small per-observation score.
# Likelihood tolerance permits small floating-point differences when comparing
# finite candidates with the analytic boundary maximum. Integration tolerance
# is passed to adaptive quadrature; it is not an accuracy guarantee for all CIs.
score_tolerance <- 1e-6
likelihood_tolerance <- 1e-6
integration_tolerance <- 1e-8
measures <- c('VaR','CTE','GS','PHT','WT')
methods <- c('Delta','Bootstrap_percentile','Profile_likelihood')
# Precompute the log-excess threshold and common normal/chi-square cutoffs.
# cutoff is the ordinary one-degree-of-freedom LR threshold; no boundary or
# misspecification calibration is substituted in this implementation.
yd <- log(d-x0)
zb <- qnorm(beta,lower.tail=FALSE)
zcrit <- qnorm((1+confidence)/2); cutoff <- qchisq(confidence,1)
stopifnot(x0>0,d>x0,alpha>1,alpha*r>1,beta>0,beta<1,
          delta>=0,delta<.5,r>0,r<=1,lambda>0,B>=20,N_MC>0,
          bootstrap_required_fraction>0,bootstrap_required_fraction<=1)

# SECTION 2. Locate the script and prepare the output directory.
# Try Rscript's --file argument, then source() frames, then an optional RStudio
# editor path. If none succeeds, write beneath the current working directory.
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
# Directory names encode selected settings; the saved Configuration.rds captures
# the remaining settings and exact source. Checkpoint matching prevents mixing
# incompatible runs even when their directory names happen to coincide.
outdir <- file.path(folder,sprintf('MC_%s_%s_x0%g_alpha%g_d%g_r%g_lambda%g',
  version,MODE,x0,alpha,d,r,lambda))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# Aggregated diagnostics. Counts of optimizer attempts are not fit counts.
D <- new.env(parent=emptyenv()); D$counts <- numeric(0)
# Increment a named diagnostic counter in a shared environment. scope separates
# outer fits, bootstrap fits, profile searches and CI processing. Several
# optimizer events can occur for one dataset, so these totals are not rates.
record <- function(scope,reason,count=1) {
  key <- paste(scope,reason,sep='|'); v <- D$counts
  if(!key %in% names(v)) v[key] <- 0
  v[key] <- v[key]+count; D$counts <- v
}
# Export the currently accumulated counters for a single sample-size cell.
# An empty counter collection needs no diagnostic file.
write_diagnostics <- function(n) {
  if(!length(D$counts)) return(invisible(NULL))
  keys <- do.call(rbind,strsplit(names(D$counts),'|',fixed=TRUE))
  write.csv(data.frame(n=n,Scope=keys[,1],Reason=keys[,2],Count=unname(D$counts)),
            file.path(outdir,sprintf('Diagnostics_LN_n%d.csv',n)),row.names=FALSE)
}

# Natural conditional density:
# f(U)=exp(-a*U-b*U^2/2)/Z(a,b), U>=0,
# b=1/sigma^2, a=(yd-mu)/sigma^2.
# b=0, a>0 is its exponential limiting family.
# Return log Z and the first two moments of U.
# SECTION 3. Conditional model normalization, likelihood and sampling.
# Inputs a and b are natural parameters. For b>0, sigma=1/sqrt(b) and
# mu=yd-a/b. Return logZ, E[U] and E[U^2] for the conditional density.
# For b=0 only a>0 yields an integrable exponential density.
normalizer <- function(a,b) {
  if(b==0) {
    if(a<=0) stop('Improper limiting density')
    return(list(logZ=-log(a),M1=1/a,M2=2/a^2))
  }
  if(!is.finite(a) || !is.finite(b) || b<0) stop('Invalid natural parameters')
  z <- a/sqrt(b)
  # For a large positive standardized threshold, direct inverse-Mills
  # subtraction loses precision. Use a truncated asymptotic expansion instead.
  # The threshold 20 selects a numerical evaluation formula, not a data cap.
  if(z>20) {
    w <- b/a^2
    e <- 1+w*(-2+w*(10+w*(-74+w*(706-8162*w))))
    return(list(logZ=-log(a)-log1p(w*e),M1=e/a,
      M2=(2+w*(-10+w*(74+w*(-706+8162*w))))/a^2))
  }
  # Compute the inverse Mills ratio on the log scale to avoid tail underflow.
  # Its excess over z determines the conditional mean and second moment.
  lm <- dnorm(z,log=TRUE)-pnorm(z,lower.tail=FALSE,log.p=TRUE)
  excess <- exp(lm)-z
  list(logZ=-log(b)/2-lm,M1=excess/sqrt(b),M2=(1-z*excess)/b)
}
# Fully observed conditional log excesses U=log((X-x0)/(d-x0)).
# T1=sum(U), T2=sum(U^2). No censoring terms enter the likelihood.
# The parameter-independent Jacobian is omitted.
# s contains only n, T1=sum(U), T2=sum(U^2), which are sufficient here.
# Return log likelihood and analytic derivatives with respect to (a,b).
# The omitted Jacobian does not depend on these parameters and therefore
# cancels in maximization, score checks and likelihood-ratio comparisons.
likelihood <- function(a,b,s) {
  z0 <- normalizer(a,b)
  ll <- -a*s$T1-b*s$T2/2-s$n*z0$logZ
  ga <- -s$T1+s$n*z0$M1
  gb <- -s$T2/2+s$n*z0$M2/2
  list(ll=ll,score=c(ga,gb))
}
# Generate a parametric bootstrap dataset from a fitted conditional model.
# A verified boundary fit samples exponential U directly. A finite fit samples
# a normal distribution truncated below zd, then converts to log excess U.
# Only sufficient statistics are retained; raw losses need not be exponentiated.
simulate_statistics <- function(fit,n) {
  if(fit$type=='limit') exposure <- rexp(n,rate=fit$a) else {
    sigma <- 1/sqrt(fit$b); zd <- fit$a/sqrt(fit$b)
    exponential <- rexp(n)
    if(zd>40) {
      # Invert the conditional cumulative hazard directly. This avoids
      # subtracting almost equal enormous normal quantiles near b=0.
      logZ0 <- normalizer(fit$a,fit$b)$logZ
      # Initialize the hazard inversion using its local exponential approximation.
      # Newton steps use the reciprocal hazard; a negative proposal is halved back
      # into the positive support. Nonconvergence is a failure, not a clipped draw.
      exposure <- exponential*exp(logZ0)
      for(iter in 1:50) {
        logZh <- vapply(fit$a+fit$b*exposure,function(a)
          normalizer(a,fit$b)$logZ,numeric(1))
        residual <- fit$a*exposure+fit$b*exposure^2/2+logZ0-logZh-exponential
        if(max(abs(residual)/(1+exponential))<1e-10) break
        proposed <- exposure-residual*exp(logZh)
        exposure <- ifelse(proposed>0,proposed,exposure/2)
      }
      if(max(abs(residual)/(1+exponential))>=1e-8) stop('Conditional sampler unresolved')
    } else {
      # If E~Exp(1), multiplying a survival probability by exp(-E) samples
      # uniformly below that survival probability. Log-tail arguments keep this
      # inverse-normal construction stable for small tail probabilities.
      z <- qnorm(pnorm(zd,lower.tail=FALSE,log.p=TRUE)-exponential,
                 lower.tail=FALSE,log.p=TRUE)
      exposure <- sigma*(z-zd)
    }
  }
  if(any(!is.finite(exposure)) || any(exposure<0))
    stop('Invalid uncensored conditional sample')
  list(n=n,T1=sum(exposure),T2=sum(exposure^2))
}

# Pareto I is scale invariant: X | X>=d = d*exp(E/alpha), E~Exp(1).
# Map to U=log((X-x0)/(d-x0)), the working likelihood coordinate.
# Work entirely on the log scale: no cap and no exponentiation of large losses.
# Outer data come from the Pareto truth, unlike the bootstrap working model.
# Conditional Pareto scale invariance is exact because d exceeds x0. The
# working-model U transformation is applied after drawing log(X/d).
simulate_pareto_statistics <- function(n) {
  log_ratio <- rexp(n,rate=alpha)
  log_x <- log(d)+log_ratio
  # log(X-x0)=log(X)+log1p(-x0/X), stable even for extreme X.
  exposure <- log_x+log1p(-exp(log(x0)-log_x))-yd
  if(any(!is.finite(exposure)) || any(exposure<0))
    stop('Invalid uncensored Pareto sample')
  list(n=n,T1=sum(exposure),T2=sum(exposure^2))
}

# Ground-up Pareto I risks use the SAME distortion conventions as the LN
# risk functions. beta is the upper-tail probability, not the CDF level.
# GS = integral_0^beta Q(1-p)*(1+2*delta-4*delta*p/beta) dp/beta.
# PHT = integral_0^Inf S(x)^r dx.
# WT = integral_0^Inf Phi(Phi^-1(S(x))+lambda) dx.
# SECTION 4. True ground-up targets for coverage evaluation.
# These targets describe the full Pareto loss distribution starting at x0,
# not the observed conditional distribution above d. VaR, CTE, GS and PHT
# have closed forms below; WT is evaluated by numerical quadrature.
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

# Optimization coordinates (t,b): a=exp(t)-b*exp(-t).
# For b>0 every real a is reachable; at b=0 a=exp(t)>0.
# SECTION 5. Extended maximum likelihood estimation.
# Invert the optimization coordinate map by solving for y=exp(t)>0.
# The alternative expression for negative a avoids subtractive cancellation.
a_to_t <- function(a,b) {
  if(b==0) return(log(a))
  root <- sqrt(a*a+4*b)
  log(if(a>=0) (a+root)/2 else 2*b/(root-a))
}
# Fit one dataset; return type, a, b and maximized log likelihood, or NULL.
# Compare a multistart finite search with the analytic exponential boundary
# MLE a0=n/T1. warm optionally reuses an outer finite fit for bootstrap starts.
# Acceptance requires likelihood and score checks, not just an optimizer code.
fit_model <- function(s,warm=NULL,scope='outer') {
  if(s$n==0 || !is.finite(s$T1) || s$T1<=0) {
    record(scope,'unresolved_invalid_statistics'); return(NULL)
  }
  a0 <- s$n/s$T1
  boundary <- likelihood(a0,0,s)
  # Different starting b values explore different curvature regimes. Including
  # the exact boundary lets the constrained optimizer approach either family.
  starts <- list(c(log(a0),0),c(a_to_t(a0,.5),.5),c(a_to_t(a0,2),2))
  if(!is.null(warm) && warm$type=='finite')
    starts <- c(list(c(a_to_t(warm$a,warm$b),warm$b)),starts)
  candidates <- list()
  bound <- 20
  repeat {
    last <- NULL; cached <- NULL
    # Cache the last objective and gradient because nlminb may request both at
    # the same parameter vector. Chain-rule derivatives convert (a,b) scores to
    # (t,b); minimizing negative average log likelihood stabilizes scaling across n.
    evaluate <- function(par) {
      if(is.null(last) || !identical(last,par)) {
        y <- exp(par[1]); b <- par[2]; a <- y-b/y
        tmp <- likelihood(a,b,s)
        cached <<- list(a=a,b=b,ll=tmp$ll,
          gradient=-c(tmp$score[1]*(y+b/y),tmp$score[2]-tmp$score[1]/y)/s$n)
        last <<- par
      }
      cached
    }
    # Run constrained optimization with b>=0 and the current expandable t range.
    # Record convergence codes for diagnosis, but validate candidates independently.
    for(start in starts) {
      start[1] <- min(bound-1e-5,max(-bound+1e-5,start[1]))
      opt <- tryCatch(nlminb(start,function(p) -evaluate(p)$ll/s$n,
        gradient=function(p) evaluate(p)$gradient,
        lower=c(-bound,0),upper=c(bound,Inf),
        control=list(iter.max=1000,eval.max=2000,rel.tol=1e-11)),error=function(e) NULL)
      if(is.null(opt)) {record(scope,'optimizer_exception'); next}
      record(scope,paste0('optimizer_code_',opt$convergence))
      v <- tryCatch(evaluate(opt$par),error=function(e) NULL)
      if(!is.null(v) && is.finite(v$ll))
        candidates[[length(candidates)+1L]] <- c(v,list(par=opt$par))
    }
    # Discard finite candidates pressed against the current artificial t boundary.
    # Among interior finite candidates, retain the highest log likelihood.
    finite <- Filter(function(v) v$b>0 && abs(v$par[1])<bound-1e-4,candidates)
    best <- if(length(finite)) finite[[which.max(vapply(finite,function(v) v$ll,numeric(1)))]] else NULL
    bestll <- if(is.null(best)) -Inf else best$ll
    # An exact boundary maximum must satisfy its one-sided score condition.
    if(boundary$score[2]<=0 && boundary$ll>=bestll-likelihood_tolerance) {
      record(scope,'accepted_limit')
      return(list(type='limit',a=a0,b=0,ll=boundary$ll))
    }
    if(!is.null(best) && bestll>=boundary$ll-likelihood_tolerance) {
      score <- likelihood(best$a,best$b,s)$score/s$n
      if(all(is.finite(score)) && max(abs(score))<score_tolerance) {
        record(scope,'accepted_finite')
        return(list(type='finite',a=best$a,b=best$b,ll=best$ll))
      }
    }
    # If acceptance checks fail, expand the t search and warm-start from the best
    # candidate. Exhausting the engineering guard produces an unresolved fit.
    if(bound>=max_t_bound) break
    bound <- min(max_t_bound,bound*2)
    if(!is.null(best)) starts <- c(list(best$par),starts)
  }
  record(scope,'unresolved_optimizer_or_score'); NULL
}

# SECTION 6. Stable risk evaluation and derivatives.
# Return the inverse Mills ratio phi(z)/Pr(Z>z), its logarithm, and its excess
# over z. The excess has its own tail expansion to avoid cancellation.
mills <- function(z) {
  if(z>20) {
    t <- 1/z; e <- t*(1+t^2*(-2+t^2*(10+t^2*(-74+t^2*(706-8162*t^2)))))
    return(list(log=log(z)+log1p(e/z),excess=e,value=z+e))
  }
  lm <- dnorm(z,log=TRUE)-pnorm(z,lower.tail=FALSE,log.p=TRUE)
  list(log=lm,excess=exp(lm)-z,value=exp(lm))
}
# log A(sigma)=coefficient*sigma^2+remainder.
# Separate the leading term to avoid cancellation in ground-up extrapolation.
# Represent log risk excess as mu+log A(sigma). Return the quadratic
# coefficient, remaining term, and derivative of that remainder with respect
# to sigma. The derivative of coefficient*sigma^2 is handled separately by
# the risk gradient. The last branch is PHT, given the fixed measures list.
factor_parts <- function(sigma,m,derivative=FALSE) {
  if(m=='VaR') return(list(coefficient=0,remainder=sigma*zb,prime=zb))
  if(m=='CTE') return(list(coefficient=.5,
    remainder=pnorm(sigma-zb,log.p=TRUE)-log(beta),prime=mills(zb-sigma)$value))
  if(m=='WT') return(list(coefficient=.5,remainder=lambda*sigma,prime=lambda))
  # The GS branch integrates its weighted normal-tail expression. Splitting
  # quadrature near zero helps resolve normal mass when the lower limit is far
  # negative. prime is the analytic derivative of the logarithmic remainder.
  if(m=='GS') {
    L <- zb-sigma
    f <- function(t) dnorm(t)*(1+2*delta-4*delta*pnorm(t+sigma,lower.tail=FALSE)/beta)
    # A normal tail below -12 is negligible at the requested accuracy.
    lower <- if(L< -12) -Inf else L
    I <- if(lower<0) integrate(f,lower,0,rel.tol=integration_tolerance)$value+
      integrate(f,0,Inf,rel.tol=integration_tolerance)$value else
      integrate(f,lower,Inf,rel.tol=integration_tolerance)$value
    prime <- (1-2*delta)*dnorm(sigma-zb)/I+
      2*delta/sqrt(pi)*exp(-sigma^2/4+pnorm(sqrt(2)*(sigma/2-zb),log.p=TRUE))/(beta*I)
    return(list(coefficient=.5,remainder=log(I)-log(beta),prime=prime))
  }
  # For r=1, PHT reduces to the mean of the shifted lognormal. Otherwise
  # recenter/rescale the PHT integral and subtract a log shift before numerical
  # integration, restoring that shift in the final logarithmic expression.
  if(r==1) return(list(coefficient=.5,remainder=0,prime=0))
  center <- sigma/r; scale <- 1/sqrt(r); shift <- (1-r)*mills(center)$log
  f <- function(t) vapply(t,function(tt)
    exp(-tt^2/2+(1-r)*mills(center+scale*tt)$log-shift),numeric(1))
  I <- integrate(f,-Inf,Inf,rel.tol=integration_tolerance)$value
  prime <- 0
  if(derivative) {
    J <- integrate(function(t) f(t)*vapply(center+scale*t,
      function(z) mills(z)$excess,numeric(1)),-Inf,Inf,rel.tol=integration_tolerance)$value
    prime <- (1-r)*J/(r*I)
  }
  list(coefficient=1/(2*r),remainder=log(r)/2-r*log(2*pi)/2+shift+log(I),prime=prime)
}
# As b tends to zero, the sign of coefficient-a determines whether the
# extrapolated risk tends to x0 or infinity. VaR has coefficient zero; PHT
# has 1/(2*r); the remaining measures have 1/2.
critical_slope <- function(m) if(m=='VaR') 0 else if(m=='PHT') 1/(2*r) else .5
# Return the risk, its optional (a,b) gradient, and log(R-x0). Boundary risks
# are limits of ground-up extrapolations; their gradients are unavailable.
# Numerical overflow for a finite fit is NA, distinct from a verified Inf limit.
risk_value <- function(fit,m,gradient=FALSE) {
  if(fit$type=='limit') {
    critical <- critical_slope(m)
    # At exact critical slopes the ground-up limit can depend on the path.
    if(m!='VaR' && abs(fit$a-critical)<1e-10*(1+critical))
      return(list(value=NA_real_,gradient=c(NA_real_,NA_real_),logexcess=NA_real_))
    return(list(value=if(fit$a>critical) x0 else Inf,
                gradient=c(NA_real_,NA_real_),logexcess=NA_real_))
  }
  sigma <- 1/sqrt(fit$b)
  f <- factor_parts(sigma,m,gradient)
  # Substitute mu=yd-a/b into mu+log A(sigma). Separating the quadratic
  # term avoids subtracting enormous mu and sigma-squared contributions later.
  q <- yd+(f$coefficient-fit$a)/fit$b+f$remainder
  if(!is.finite(q) || q>log(.Machine$double.xmax))
    return(list(value=NA_real_,gradient=c(NA_real_,NA_real_),logexcess=q))
  excess <- exp(q)
  # Differentiate R=x0+exp(q) by the chain rule. Underflow to zero excess
  # returns x0 but leaves the gradient unavailable under the current convention.
  g <- if(gradient && excess>0) excess*c(-1/fit$b,
    (fit$a-f$coefficient)/fit$b^2-f$prime/(2*fit$b^(3/2))) else c(NA_real_,NA_real_)
  list(value=x0+excess,gradient=g,logexcess=q)
}
# Evaluate every risk independently; one failed integral becomes a missing
# value for that measure rather than erasing other bootstrap risk values.
all_risks <- function(fit) vapply(measures,function(m)
  tryCatch(risk_value(fit,m)$value,error=function(e) NA_real_),numeric(1))

# Expected information per fully observed left-truncated record in (a,b).
# Integrate centered score products over the entire conditional support.
# A location/scale change keeps the integration near its probability mass.
# SECTION 7. Model-based Delta confidence intervals.
# Compute per-observation expected Fisher information as the covariance of
# conditional scores under the fitted working model. This is not the sandwich
# covariance required for a general misspecification-robust Delta procedure.
expected_information <- function(fit) {
  a <- fit$a; b <- fit$b
  z0 <- normalizer(a,b)
  # Center at the mode when it is positive. Scale by an exponential or normal
  # width, then transform U=center+scale*t. log(scale) supplies the Jacobian.
  center <- if(b>0 && a<0) -a/b else 0
  scale <- if(b==0 || a>sqrt(b)) 1/a else 1/sqrt(b)
  lower <- -center/scale
  I <- matrix(0,2,2)
  for(i in 1:2) for(j in i:2) {
    integrand <- function(t) {
      exposure <- center+scale*t
      logdensity <- -a*exposure-b*exposure^2/2-z0$logZ+log(scale)
      # Scores are centered using the fitted conditional moments. Integrating
      # pairwise products gives a symmetric 2-by-2 information matrix.
      scores <- cbind(z0$M1-exposure,(z0$M2-exposure^2)/2)
      ans <- exp(logdensity)*scores[,i]*scores[,j]
      ans[!is.finite(t)] <- 0
      ans
    }
    # Split at the mode to resolve both sides of a narrow normal peak.
    value <- if(lower<0)
      integrate(integrand,lower,0,rel.tol=integration_tolerance,subdivisions=1000L)$value+
      integrate(integrand,0,Inf,rel.tol=integration_tolerance,subdivisions=1000L)$value else
      integrate(integrand,lower,Inf,rel.tol=integration_tolerance,subdivisions=1000L)$value
    I[i,j] <- I[j,i] <- value
  }
  if(any(!is.finite(I))) stop('Nonfinite expected information')
  I
}
# Use the first-order variance grad(R)' I^(-1) grad(R)/n and a symmetric
# normal interval on the original monetary scale. Endpoints are not clipped
# to the risk support. A limiting original fit is explicitly unavailable.
delta_ci <- function(fit,n,m,I=NULL) {
  if(fit$type=='limit') return(list(ci=c(NA_real_,NA_real_),reason='delta_limit_unavailable'))
  rv <- risk_value(fit,m,TRUE)
  if(any(!is.finite(c(rv$value,rv$gradient)))) stop('Nonfinite delta derivative')
  if(is.null(I)) I <- expected_information(fit)
  # Cholesky solves avoid explicitly inverting a nearly singular matrix.
  R <- chol(I); v <- forwardsolve(t(R),rv$gradient)
  se <- sqrt(sum(v^2)/n)
  if(!is.finite(se)) stop('Nonfinite delta standard error')
  list(ci=rv$value+c(-1,1)*zcrit*se,reason='ok')
}
# SECTION 8. Parametric percentile bootstrap confidence intervals.
# Rows are bootstrap datasets; columns are risk measures. A failed sampler,
# fit or risk evaluation stays missing. Successful boundary fits are retained.
bootstrap_cis <- function(fit,n) {
  values <- matrix(NA_real_,B,length(measures))
  for(b in seq_len(B)) {
    sb <- tryCatch(simulate_statistics(fit,n),error=function(e) {
      record('bootstrap','conditional_sampler_unresolved'); NULL
    })
    if(is.null(sb)) next
    fb <- fit_model(sb,warm=fit,scope='bootstrap')
    if(!is.null(fb)) values[b,] <- all_risks(fb)
  }
  # Each output row contains lower/upper endpoints for one measure. Type 1
  # quantiles select empirical order statistics, avoiding interpolation between
  # finite and infinite draws. Final infinite-CI policy is applied only in update().
  result <- matrix(NA_real_,length(measures),2)
  for(j in seq_along(measures)) {
    good <- !is.na(values[,j]) # verified limiting Inf values are retained
    record('bootstrap',paste0(measures[j],'_missing_values'),sum(!good))
    record('bootstrap',paste0(measures[j],'_infinite_values'),sum(is.infinite(values[,j])))
    if(sum(good)>=ceiling(B*bootstrap_required_fraction))
      result[j,] <- quantile(values[good,j],c((1-confidence)/2,(1+confidence)/2),type=1,names=FALSE)
  }
  result
}

# Risk profiles include finite parameters and the relevant likelihood limits.
# At a fixed finite q=log(R-x0), sigma->Inf leads to a=critical_slope(m).
# SECTION 9. Profile likelihood and confidence-set inversion.
# For fixed q=log(R-x0), eliminate a algebraically and maximize over log sigma.
# Include the analytic sigma->infinity limit where the conditional density is
# proper. Return the largest resolved likelihood value across these candidates.
profile_ll <- function(q,m,s) {
  critical <- critical_slope(m)
  limiting <- if(critical>0) likelihood(critical,0,s)$ll else -Inf
  # Use log sigma to preserve positivity. A nonfinite likelihood receives a
  # large minimizing penalty; exceptions are caught by the surrounding grid scan.
  objective <- function(ls) {
    sigma <- exp(ls); b <- exp(-2*ls); f <- factor_parts(sigma,m)
    a <- f$coefficient+b*(yd-q+f$remainder)
    val <- likelihood(a,b,s)$ll
    if(is.finite(val)) -val else 1e100
  }
  low <- -3; high <- 3
  repeat {
    # Scan a grid and refine each apparent local minimum of negative likelihood.
    # This multistart strategy reduces dependence on a single nuisance optimum;
    # it is a numerical search, not a proof of global optimality.
    grid <- seq(low,high,length.out=max(31,ceiling((high-low)/.25)+1))
    vals <- vapply(grid,function(z) tryCatch(objective(z),error=function(e) 1e100),numeric(1))
    idx <- unique(c(which(vals<=c(Inf,head(vals,-1)) & vals<=c(tail(vals,-1),Inf)),which.min(vals)))
    opts <- lapply(idx,function(i) tryCatch(optimize(objective,
      grid[c(max(1,i-1),min(length(grid),i+1))],tol=1e-7),error=function(e) NULL))
    opts <- Filter(Negate(is.null),opts)
    best <- max(c(limiting,-vals,if(length(opts))
      -vapply(opts,function(o) o$objective,numeric(1)) else -Inf))
    # Require each search edge to be clearly below the best log likelihood or,
    # for the upper edge, sufficiently close to the analytic limiting likelihood.
    # Otherwise expand the log-sigma range until resolved or a guard is reached.
    lower_ok <- -vals[1]<best-10
    upper_ok <- (is.finite(limiting) && abs(-tail(vals,1)-limiting)<1e-5) ||
      (-tail(vals,1)<best-10)
    if(lower_ok && upper_ok) return(best)
    if((!lower_ok && low<=min_log_sigma) || (!upper_ok && high>=max_log_sigma))
      stop('Unresolved profile search guard')
    if(!lower_ok) low <- max(min_log_sigma,low-2)
    if(!upper_ok) high <- min(max_log_sigma,high+2)
  }
}
# Invert 2*(fitted maximum-profile maximum)<=chi-square cutoff. Before
# finite root searches, evaluate whether risk=x0 and risk=infinity are admitted
# through the relevant boundary families. Boundary inclusion is likelihood based.
profile_ci <- function(fit,m,s) {
  a0 <- s$n/s$T1; critical <- critical_slope(m)
  # Optimize the exponential boundary likelihood over a on the side yielding
  # risk x0, and separately on the side yielding infinite risk. At a critical
  # slope these are limiting envelope values for endpoint inclusion.
  zero_ll <- likelihood(max(a0,critical),0,s)$ll
  infinity_ll <- if(critical>0) likelihood(min(a0,critical),0,s)$ll else -Inf
  admitted <- function(ll) is.finite(ll) && 2*(fit$ll-ll)<=cutoff+1e-7
  lower_limit <- admitted(zero_ll); upper_limit <- admitted(infinity_ll)
  if(lower_limit && upper_limit) return(c(x0,Inf))
  threshold <- fit$ll-cutoff/2
  # A data-based central q; no generating parameters enter profiling.
  ref <- if(fit$type=='finite') risk_value(fit,m)$logexcess else yd
  if(!is.finite(ref)) ref <- yd
  # A nonpositive dev denotes admission to the LR confidence set. A profile
  # above the purported fitted maximum signals inconsistent optimization.
  dev <- function(q) {
    ll <- profile_ll(q,m,s)
    if(ll>fit$ll+1e-4) stop('Profile exceeds fitted maximum')
    threshold-ll
  }
  # Find an admitted finite point, including when the MLE is a limiting fit.
  center <- NULL
  for(k in 0:18) {
    step <- if(k==0) 0 else .25*2^(k-1)
    candidates <- if(step==0) ref else c(ref-step,ref+step)
    for(q in candidates) if(is.finite(q) && abs(q)<=1e4 && dev(q)<=0) {center <- q; break}
    if(!is.null(center)) break
  }
  if(is.null(center)) stop('No admitted finite profile point resolved')
  # Expand outward from an admitted center until dev changes sign, then solve
  # for a finite endpoint. This finds the component reached from that center;
  # additional disconnected components require separate profile inspection.
  endpoint <- function(side) {
    step <- .25
    for(j in 1:60) {
      q <- center+side*step
      if(abs(q)>1e4) stop('Unresolved log-risk search guard')
      if(dev(q)>0) return(uniroot(dev,sort(c(center,q)),tol=1e-7)$root)
      step <- step*1.5
    }
    stop('Unresolved profile endpoint')
  }
  # Convert finite log-excess endpoints back to monetary units. Guard against
  # floating-point overflow rather than declaring an unresolved finite endpoint
  # to be a mathematically verified infinite endpoint.
  lower <- if(lower_limit) x0 else {
    q <- endpoint(-1)
    if(q>log(.Machine$double.xmax)) stop('Finite endpoint exceeds floating-point range')
    x0+exp(q)
  }
  upper <- if(upper_limit) Inf else {
    q <- endpoint(1)
    if(q>log(.Machine$double.xmax)) stop('Finite endpoint exceeds floating-point range')
    x0+exp(q)
  }
  c(lower,upper)
}

# Numerical unit checks; these stop the run on disagreement.
# SECTION 10. Internal numerical consistency checks.
# Use deterministic small examples to compare scores with finite differences,
# verify the exponential boundary, check risk gradients and validate information.
# These are implementation checks, not evidence of nominal statistical coverage.
run_checks <- function() {
  exposure <- c(.1,.3,.6,1,2)
  ss <- list(n=length(exposure),T1=sum(exposure),T2=sum(exposure^2))
  for(par in list(c(1,.5),c(.2,.1),c(-.5,2),c(4,.01))) {
    analytical <- likelihood(par[1],par[2],ss)$score
    numerical <- vapply(1:2,function(j) {
      eps <- 1e-5; plus <- minus <- par; plus[j] <- plus[j]+eps; minus[j] <- minus[j]-eps
      (likelihood(plus[1],plus[2],ss)$ll-likelihood(minus[1],minus[2],ss)$ll)/(2*eps)
    },numeric(1))
    stopifnot(max(abs(analytical-numerical)/(1+abs(numerical)))<1e-5)
  }
  a0 <- ss$n/ss$T1
  stopifnot(abs(likelihood(a0,0,ss)$ll-(ss$n*log(a0)-a0*ss$T1))<1e-10)
  # At b=0 use a forward difference because b<0 is outside the model.
  numerical_b0 <- (likelihood(a0,1e-6,ss)$ll-likelihood(a0,0,ss)$ll)/1e-6
  stopifnot(abs(numerical_b0-likelihood(a0,0,ss)$score[2])/(1+abs(numerical_b0))<1e-4)
  # Arbitrary finite working-model parameters for derivative checks only.
  true <- list(type='finite',a=.7,b=.4)
  for(m in measures) {
    rv <- risk_value(true,m,TRUE); eps <- 1e-5
    num <- vapply(1:2,function(j) {
      plus <- minus <- true
      key <- c('a','b')[j]; plus[[key]] <- plus[[key]]+eps; minus[[key]] <- minus[[key]]-eps
      (risk_value(plus,m)$value-risk_value(minus,m)$value)/(2*eps)
    },numeric(1))
    stopifnot(max(abs(rv$gradient-num)/(1+abs(num)))<1e-4)
  }
  stopifnot(all(all_risks(list(type='limit',a=2,b=0))==x0),
            is.infinite(risk_value(list(type='limit',a=.2,b=0),'WT')$value))
  # Cholesky must succeed for the selected finite example. Then compare the
  # boundary information matrix against exact exponential moment identities.
  chol(expected_information(true))
  check_a <- 1.7
  limit_I <- expected_information(list(type='limit',a=check_a,b=0))
  analytic_I <- matrix(c(1/check_a^2,2/check_a^3,
                         2/check_a^3,5/check_a^4),2,2)
  stopifnot(max(abs(limit_I-analytic_I)/(1+abs(analytic_I)))<1e-7)
  cat('Likelihood scores, risk derivatives, limit values and expected-information checks passed.\n')
}
if(RUN_CHECKS) run_checks()

# SECTION 11. Interval classification and retained replication diagnostics.
# Create one row per risk-method pair. Valid and Covered refer to intervals
# accepted under EXCLUDE_INFINITE_CIS; finite-specific counts define the primary
# publication denominator independently. List columns retain exact finite lengths
# and WT records for quantiles and diagnostics; they are kept in checkpoints.
new_accumulator <- function() {
  g <- expand.grid(Risk_Measure=measures,Method=methods,stringsAsFactors=FALSE)
  for(key in c('Valid','Covered','Finite_Length_Count','Infinite_Length_Count','Undefined_Length_Count')) g[[key]] <- 0
  g$Infinite_CI_Excluded <- g$Delta_Limit_Unavailable <- g$Unresolved_CI <- 0
  g$Sum_Finite_Length <- 0
  g$Finite_Covered <- g$Infinite_CI_Count <- 0
  g$Finite_Lengths <- I(rep(list(numeric(0)),nrow(g)))
  g$WT_Details <- I(rep(list(list()),nrow(g)))
  g
}
# Process one method/measure interval. Missing, malformed and reversed intervals
# are classified first. Unavailable boundary Delta intervals get their own count;
# other invalid intervals increment Unresolved_CI. Returns the modified table.
update <- function(acc,method,m,ci,reason='ok',replication=NA_integer_,fit=NULL,s=NULL) {
  i <- which(acc$Method==method & acc$Risk_Measure==m)
  # Classify malformed/undefined intervals before counting infinite intervals.
  if(length(ci)!=2 || anyNA(ci) || ci[1]>ci[2]) {
    key <- if(reason=='delta_limit_unavailable') 'Delta_Limit_Unavailable' else 'Unresolved_CI'
    acc[[key]][i] <- acc[[key]][i]+1
    record('CI',paste0(method,'_',m,'_',reason)); return(acc)
  }
  # Count an ordered, nonmissing infinite-endpoint interval even if it will be
  # excluded. With exclusion disabled it contributes to accepted coverage; finite
  # truth is not covered by [Inf,Inf], whose subtraction also has undefined length.
  if(any(is.infinite(ci))) {
    acc$Infinite_CI_Count[i] <- acc$Infinite_CI_Count[i]+1
    if(EXCLUDE_INFINITE_CIS) {
      acc$Infinite_CI_Excluded[i] <- acc$Infinite_CI_Excluded[i]+1
      record('CI',paste0(method,'_',m,'_infinite_excluded')); return(acc)
    }
  }
  # Inclusive endpoint comparisons define coverage against the Pareto target.
  # Retained nonfinite lengths remain separate from the finite-length distribution.
  acc$Valid[i] <- acc$Valid[i]+1
  acc$Covered[i] <- acc$Covered[i]+as.integer(ci[1]<=true_risks[m] && true_risks[m]<=ci[2])
  len <- ci[2]-ci[1]
  if(is.na(len)) acc$Undefined_Length_Count[i] <- acc$Undefined_Length_Count[i]+1 else
    if(is.infinite(len)) acc$Infinite_Length_Count[i] <- acc$Infinite_Length_Count[i]+1 else {
      acc$Finite_Length_Count[i] <- acc$Finite_Length_Count[i]+1
      acc$Sum_Finite_Length[i] <- acc$Sum_Finite_Length[i]+len
      # Store individual lengths so median and extreme quantiles can be computed
      # without a distributional assumption. Rare enormous finite WT lengths can
      # strongly influence the mean even when the median remains moderate.
      acc$Finite_Lengths[[i]] <- c(acc$Finite_Lengths[[i]],len)
      acc$Finite_Covered[i] <- acc$Finite_Covered[i]+
        as.integer(ci[1]<=true_risks[m] && true_risks[m]<=ci[2])
      # Retain finite WT interval endpoints, coverage and fitting context. Replication
      # and cell seeds identify the simulated dataset for a later targeted rerun;
      # T1 and T2 are sufficient for refitting the conditional likelihood.
      if(m=='WT') acc$WT_Details[[i]][[length(acc$WT_Details[[i]])+1L]] <-
        data.frame(Replication=replication,Lower=ci[1],Upper=ci[2],Length=len,
          Covered=ci[1]<=true_risks[m] && true_risks[m]<=ci[2],
          Fit_Type=if(is.null(fit)) NA_character_ else fit$type,
          Fit_a=if(is.null(fit)) NA_real_ else fit$a,
          Fit_b=if(is.null(fit)) NA_real_ else fit$b,
          T1=if(is.null(s)) NA_real_ else s$T1,
          T2=if(is.null(s)) NA_real_ else s$T2)
    }
  acc
}
# SECTION 12. Reproducibility metadata and checkpoint handling.
# Evaluate true targets once and require finite values before any coverage run.
true_risks <- pareto_true_risks()
# True risks are computed from the current parameters; no fixed numerical targets.
stopifnot(all(is.finite(true_risks)))
write.csv(data.frame(Risk_Measure=measures,True_Risk=unname(true_risks)),
          file.path(outdir,'True_Risks.csv'),row.names=FALSE)
# Keep the exact source alongside settings. The strict checkpoint comparison
# below includes this text, so edited comments also prevent accidental resumption.
source_text <- if(length(script)==1 && file.exists(script)) readLines(script,warn=FALSE) else NULL
configuration <- list(version=version,MODE=MODE,N_MC=N_MC,B=B,n=sample_sizes,
  data_distribution='Pareto I',fitted_distribution='shifted lognormal',
  observation_scheme='left truncation only',x0=x0,alpha=alpha,d=d,beta=beta,delta=delta,r=r,lambda=lambda,
  confidence=confidence,bootstrap_required_fraction=bootstrap_required_fraction,
  exclude_infinite_CIs=EXCLUDE_INFINITE_CIS,seed=seed,source=source_text,
  max_t_bound=max_t_bound,max_log_sigma=max_log_sigma,min_log_sigma=min_log_sigma,
  score_tolerance=score_tolerance,likelihood_tolerance=likelihood_tolerance,
  integration_tolerance=integration_tolerance)
saveRDS(configuration,file.path(outdir,'Configuration.rds'))
# Write a temporary RDS before replacing the main checkpoint. Preserve a backup
# during renaming and attempt restoration if installing the new file fails.
# This reduces interruption risk; it does not provide a transactional filesystem.
save_checkpoint <- function(state,path) {
  tmp <- paste0(path,'.tmp'); backup <- paste0(path,'.bak'); saveRDS(state,tmp)
  if(file.exists(backup)) unlink(backup)
  if(file.exists(path) && !file.rename(path,backup)) stop('Cannot back up checkpoint')
  if(!file.rename(tmp,path)) {if(file.exists(backup)) file.rename(backup,path); stop('Cannot save checkpoint')}
  if(file.exists(backup)) unlink(backup)
}
# SECTION 13. Outer Monte Carlo experiment for one sample size.
# state holds replication progress, original-fit counts, the accumulator,
# diagnostics and configuration. A missing original fit affects every method and
# measure in that replication and is accounted for at final summarization.
run_cell <- function(n,cell_id) {
  path <- file.path(outdir,sprintf('checkpoint_n%d.rds',n)); backup <- paste0(path,'.bak')
  state <- list(done=0,finite_fits=0,limit_fits=0,unresolved_fits=0,
                acc=new_accumulator(),diagnostics=numeric(0),configuration=configuration)
  if(resume && (file.exists(path) || file.exists(backup))) {
    state <- readRDS(if(file.exists(path)) path else backup)
    if(!identical(state$configuration,configuration)) stop('Checkpoint configuration differs')
  }
  # Restore only this cell's diagnostics. Deterministic per-replication seeding
  # means the RNG state itself does not need to be saved in the checkpoint.
  D$counts <- state$diagnostics
  if(state$done<N_MC) for(replication in seq.int(state$done+1,N_MC)) {
    # cell_id depends on the order of sample_sizes. Preserve that order to
    # reproduce the same datasets; all methods within a replication share the data.
    set.seed(seed+cell_id*1000000+replication)
    s <- simulate_pareto_statistics(n); fit <- fit_model(s)
    if(is.null(fit)) state$unresolved_fits <- state$unresolved_fits+1 else {
      if(fit$type=='finite') state$finite_fits <- state$finite_fits+1 else state$limit_fits <- state$limit_fits+1
      # Compute information once per finite original fit and bootstrap all risks
      # together. Information failure blocks Delta CIs but allows other methods.
      I <- if(fit$type=='finite') tryCatch(expected_information(fit),error=function(e) NULL) else NULL
      boot <- bootstrap_cis(fit,n)
      for(j in seq_along(measures)) {
        m <- measures[j]
        # Convert numerical Delta errors to missing endpoints with a reason code.
        # This preserves the attempted replication instead of silently resampling it.
        dc <- if(fit$type=='finite' && is.null(I))
          list(ci=c(NA_real_,NA_real_),reason='delta_expected_information_failed') else
          tryCatch(delta_ci(fit,n,m,I),error=function(e)
            list(ci=c(NA_real_,NA_real_),reason='delta_numerical_failure'))
        record('outer',dc$reason)
        # Profile failures are recorded with a sanitized message; CI classification
        # then counts the missing interval. Bootstrap failures were recorded internally.
        pc <- tryCatch(profile_ci(fit,m,s),error=function(e) {
          record('profile',gsub('[|\r\n]',' ',conditionMessage(e))); c(NA_real_,NA_real_)
        })
        state$acc <- update(state$acc,'Delta',m,dc$ci,dc$reason,replication,fit,s)
        state$acc <- update(state$acc,'Bootstrap_percentile',m,boot[j,],replication=replication,fit=fit,s=s)
        state$acc <- update(state$acc,'Profile_likelihood',m,pc,replication=replication,fit=fit,s=s)
      }
    }
    # Mark completion only after all available methods have been processed, then
    # save at the requested cadence and at the final replication.
    state$done <- replication
    if(replication %% checkpoint_every==0 || replication==N_MC) {
      state$diagnostics <- D$counts; save_checkpoint(state,path); write_diagnostics(n)
      cat('Pareto I -> shifted LN n=',n,': ',replication,'/',N_MC,'\n',sep=''); flush.console()
    }
  }
  # SECTION 14. Primary and supplementary publication summaries.
  # Summarize a copy of the checkpoint accumulator; keep its raw list columns
  # intact in state for reproducibility. Completed includes unresolved outer fits.
  acc <- state$acc
  acc$n <- n; acc$Completed <- state$done
  acc$Original_Finite_Fits <- state$finite_fits; acc$Original_Limit_Fits <- state$limit_fits
  acc$Original_Unresolved_Fits <- state$unresolved_fits; acc$Failed <- state$done-acc$Valid
  # Audit an exhaustive partition of attempts: accepted intervals, excluded
  # infinite intervals, unavailable boundary Delta CIs, unresolved CIs, and failed
  # original fits. Infinite/undefined length subcounts are not additional partitions.
  stopifnot(all(acc$Completed==acc$Valid+acc$Infinite_CI_Excluded+
    acc$Delta_Limit_Unavailable+acc$Unresolved_CI+acc$Original_Unresolved_Fits))
  # Supplementary accepted coverage uses Valid, so its interpretation changes
  # when infinite intervals are retained. Operational success uses all attempts;
  # it treats excluded or unavailable outcomes as unsuccessful operationally.
  acc$Coverage_Among_Accepted_Intervals <- ifelse(acc$Valid>0,acc$Covered/acc$Valid,NA_real_)
  acc$Coverage_Success_Over_All <- acc$Covered/state$done
  acc$MCSE_Coverage <- sqrt(acc$Coverage_Among_Accepted_Intervals*(1-acc$Coverage_Among_Accepted_Intervals)/acc$Valid)
  # The accepted-interval mean is undefined if any accepted length is undefined,
  # and infinite if any accepted length is infinite. The finite-only mean discards
  # those lengths; it can still be dominated by rare very large finite values.
  # Sum_Finite_Length itself can overflow even when each individual length is finite.
  acc$Mean_Length_Valid <- ifelse(acc$Valid==0 | acc$Undefined_Length_Count>0,NA_real_,
    ifelse(acc$Infinite_Length_Count>0,Inf,acc$Sum_Finite_Length/acc$Valid))
  acc$Mean_Finite_Length <- ifelse(acc$Finite_Length_Count>0,
    acc$Sum_Finite_Length/acc$Finite_Length_Count,NA_real_)
  acc$Mean_Length_Valid_Thousand <- acc$Mean_Length_Valid/1000
  acc$Mean_Finite_Length_Thousand <- acc$Mean_Finite_Length/1000
  # Primary coverage always conditions on finite valid endpoints AND finite length,
  # regardless of EXCLUDE_INFINITE_CIS. Finite-endpoint subtraction overflow
  # is reported separately, not silently included in finite-length summaries.
  acc$Coverage_Among_Finite_Valid_Intervals <- ifelse(acc$Finite_Length_Count>0,
    acc$Finite_Covered/acc$Finite_Length_Count,NA_real_)
  acc$MCSE_Coverage_Among_Finite_Valid_Intervals <- with(acc,
    sqrt(Coverage_Among_Finite_Valid_Intervals*
      (1-Coverage_Among_Finite_Valid_Intervals)/Finite_Length_Count))
  # Primary median describes typical finite precision; P95, P99 and maximum
  # expose the right tail. Type 7 interpolated sample quantiles summarize lengths,
  # whereas bootstrap CI endpoints deliberately use type 1 order statistics.
  for(stat in c('Median','P95','P99','Maximum')) {
    acc[[paste0(stat,'_Finite_Length')]] <- vapply(acc$Finite_Lengths,function(x) {
      if(!length(x)) return(NA_real_)
      switch(stat,Median=median(x),P95=unname(quantile(x,.95,type=7)),
        P99=unname(quantile(x,.99,type=7)),Maximum=max(x))
    },numeric(1))
    acc[[paste0(stat,'_Finite_Length_Thousand')]] <- acc[[paste0(stat,'_Finite_Length')]]/1000
  }
  # This composite includes statistical unavailability of Delta at the boundary
  # as well as numerical/undefined outcomes. All percentages divide by Completed.
  # Infinite and undefined categories can overlap if infinite CIs are retained;
  # finite-endpoint length overflow is recorded in Infinite_Length_Count.
  acc$Undefined_Or_Numerical_Failure_Count <- with(acc,
    Delta_Limit_Unavailable+Unresolved_CI+Original_Unresolved_Fits+Undefined_Length_Count)
  acc$Percent_Infinite_Intervals <- 100*acc$Infinite_CI_Count/acc$Completed
  acc$Percent_Undefined_Or_Numerical_Failures <-
    100*acc$Undefined_Or_Numerical_Failure_Count/acc$Completed
  # Rank finite WT intervals separately within each method and retain up to ten.
  # No WT diagnostic CSV is written when no eligible finite WT records exist.
  wt <- list()
  for(i in which(acc$Risk_Measure=='WT')) {
    rows <- acc$WT_Details[[i]]
    if(length(rows)) {
      tab <- do.call(rbind,rows)
      tab <- head(tab[order(tab$Length,decreasing=TRUE),],10L)
      tab$n <- n; tab$Method <- acc$Method[i]; tab$Rank <- seq_len(nrow(tab))
      wt[[length(wt)+1L]] <- tab
    }
  }
  if(length(wt)) write.csv(do.call(rbind,wt),
    file.path(outdir,sprintf('Largest_WT_Intervals_n%d.csv',n)),row.names=FALSE)
  # Remove nested/list data before CSV export; their full contents remain in the
  # checkpoint. The main table contains primary measures; the full summary retains
  # supplementary measures and detailed denominator/exclusion counts.
  acc$Finite_Lengths <- acc$WT_Details <- NULL
  primary <- c('n','Risk_Measure','Method','Completed','Finite_Length_Count',
    'Coverage_Among_Finite_Valid_Intervals','MCSE_Coverage_Among_Finite_Valid_Intervals',
    'Median_Finite_Length','Median_Finite_Length_Thousand','Percent_Infinite_Intervals',
    'Percent_Undefined_Or_Numerical_Failures')
  write.csv(acc[,primary],file.path(outdir,sprintf('Main_Results_n%d.csv',n)),
    row.names=FALSE,na='NA')
  write.csv(acc,file.path(outdir,sprintf('Summary_LN_n%d.csv',n)),row.names=FALSE,na='NA')
  write_diagnostics(n); acc
}
# SECTION 15. Run sample sizes, combine tables and write interpretation notes.
# Update the combined CSV after each completed cell so finished cells remain
# available even if a later sample size is interrupted.
results <- vector('list',length(sample_sizes))
for(i in seq_along(sample_sizes)) {
  results[[i]] <- run_cell(sample_sizes[i],i)
  write.csv(do.call(rbind,results[seq_len(i)]),file.path(outdir,'Monte_Carlo_Summary.csv'),row.names=FALSE,na='NA')
}
# Write a companion README documenting assumptions, denominators, numerical
# limitations and output interpretation beside the result CSVs.
writeLines(c('MISSPECIFICATION: Pareto I truth; shifted lognormal working model.',
 sprintf('Pareto minimum x0=%g; alpha=%g; truncation d=%g; no right censoring.',x0,alpha,d),
 'n counts fully observed records conditional on X>=d; no upper cap is applied.',
 'Coverage targets true ground-up Pareto risks; it does not target a pseudo-true LN risk.',
 'Parametric bootstrap draws from fitted LN, never from the Pareto truth.',
 'Delta is model based, without a sandwich correction.',
 'Ordinary chi-square profile calibration is not corrected for misspecification.',
 'Misspecification bias and calibration error can both reduce coverage.',
 'Extended conditional LN MLE: finite and verified b=0 limits.',
 'No arbitrary upper sigma bound; numerical guard hits remain unresolved.',
 'Critical-slope risk limits may be path-dependent and are marked NA.',
 'Delta uses expected information and is unavailable for limiting original fits.',
 'Bootstrap retains verified limit values; failed samples are never redrawn.',
 'Profile endpoints at x0 or Inf require an admitted limiting likelihood.',
 'Profile inversion returns the component found from an admitted point; inspect profile shape if nonconnectedness is suspected.',
 'Boundary bootstrap risk values, including Inf, are retained when computing percentiles.',
 'EXCLUDE_INFINITE_CIS controls excluding only FINAL infinite CIs; exclusions are counted separately.',
 'With exclusions enabled, reported valid coverage/lengths condition on finite accepted intervals.',
 'If EXCLUDE_INFINITE_CIS is FALSE, [Inf,Inf] contributes noncoverage of finite truth and has undefined length.',
 'Main_Results_n*.csv reports primary publication measures; Summary_LN_n*.csv and Monte_Carlo_Summary.csv include supplementary measures and exclusion counts.',
 'Median finite interval length describes typical precision; mean finite length describes average length and may be dominated by rare enormous intervals.',
 'Primary coverage is coverage among finite valid intervals, not unconditional coverage probability. Its denominator is Finite_Length_Count.',
 'Primary coverage MCSE is sqrt(p*(1-p)/Finite_Length_Count); no eligible intervals gives NA.',
 'Coverage_Among_Accepted_Intervals and MCSE_Coverage use Valid as denominator and depend on EXCLUDE_INFINITE_CIS.',
 'Coverage_Success_Over_All counts covered accepted intervals over all completed replications; it is operational success, not unconditional coverage.',
 'All percentages use Completed as denominator. Infinite_CI_Count counts ordered nonmissing intervals with infinite endpoints, including undefined-length [Inf,Inf].',
 'Undefined_Or_Numerical_Failure_Count includes unresolved original fits, unavailable boundary Delta intervals, unresolved CIs, and accepted undefined lengths.',
 'Infinite and undefined percentages may overlap when an infinite interval has undefined length and infinite CIs are retained.',
 'Infinite_Length_Count includes finite-endpoint subtraction overflow; such intervals are omitted from finite coverage/length summaries.',
 'Supplementary finite length quantiles use R type=7; lengths are in original currency units and Thousand columns divide by 1000.',
 'Largest_WT_Intervals_n*.csv contains the top 10 finite lengths per method, endpoints, coverage, replication, fit parameters and sufficient statistics; checkpoint seeds permit reproduction.',
 'Detailed exclusion counts are retained in summary tables; optimizer and numerical reasons are in Diagnostics_LN_n*.csv.',
 'This implementation requires local R validation and a paired pilot before production.',
 'Standard percentile bootstrap and chi-square LR calibration may still be inaccurate near boundaries.'),
 file.path(outdir,'README.txt'))
# Display a compact primary-results view in the console. Full-precision CSV
# outputs and saved metadata remain the authoritative detailed run record.
print(do.call(rbind,results)[,c('n','Risk_Measure','Method','Valid','Failed',
 'Coverage_Among_Finite_Valid_Intervals','MCSE_Coverage_Among_Finite_Valid_Intervals',
 'Median_Finite_Length_Thousand','Percent_Infinite_Intervals','Percent_Undefined_Or_Numerical_Failures',
 'Infinite_CI_Excluded','Delta_Limit_Unavailable','Unresolved_CI','Original_Limit_Fits')],row.names=FALSE)
cat('Outputs:',outdir,'\n')
