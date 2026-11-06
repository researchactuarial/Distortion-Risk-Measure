# SHIFTED LOGNORMAL: boundary-aware Monte Carlo study.
# X=x0+exp(Y), Y~N(mu,sigma^2), observed conditional on X>=d,
# right-censored at u. n includes exact AND censored claims.
# Fits the finite conditional model and its exponential log-excess limit.
# Retains verified limiting bootstrap fits; never redraws failed datasets.
# Delta uses expected information, profile checks limiting likelihoods.
# This changes the estimator to an EXTENDED MLE when a finite MLE is absent.
# Boundary-aware estimation does NOT guarantee nominal bootstrap/LR coverage.
# Base R only. R execution was unavailable during authoring.

# 1. STUDY SETTINGS AND DERIVED CONSTANTS
# Run a small pilot by setting MODE to a value other than 'production'. The
# current production settings request 5,000 outer datasets and 500 bootstrap
# refits per dataset, so a complete run can be computationally expensive.
# N_MC controls Monte Carlo precision; B controls bootstrap quantile resolution.
MODE <- 'production'  # Change to 'production' AFTER pilot checks.
N_MC <- if(MODE=='production') 5000L else 10L
B <- if(MODE=='production') 500L else 30L
# Each n is the number of observed conditional claims, including censored ones.
# x0 is the fixed location shift, d the observation threshold, and u the
# right-censoring threshold. The generating mu/sigma describe log(X-x0).
# The risk targets below refer to the ground-up shifted distribution; fitting
# uses only the distribution conditional on exceeding d.
sample_sizes <- c(50L,100L,500L)
x0 <- 1000; d <- 4000; u <- 14000
mu_true <- 6; sigma_true <- sqrt(2)
# Risk settings: beta is the upper-tail probability for VaR/CTE and the GS
# formula; delta shapes GS, r shapes the proportional hazards transform (PHT),
# and lambda is the Wang transform (WT) normal-quantile shift.
# These are fixed study inputs, rather than estimated parameters.
beta <- .10; delta <- .25; r <- .75; lambda <- .50
confidence <- .90
# A value of 1 requires a nonmissing risk estimate from every bootstrap draw
# for that measure. Verified infinite boundary values count as available.
# Reducing this fraction allows percentile intervals based on fewer draws and
# changes the treatment of bootstrap failures; document any such change.
bootstrap_required_fraction <- 1
EXCLUDE_INFINITE_CIS <- TRUE # Retain boundary bootstrap values; exclude only final infinite CIs.
# Checkpoint after this many completed outer replications. resume reloads a
# compatible saved state; RUN_CHECKS runs deterministic numerical checks first.
# A separate seed is assigned to each sample-size cell and replication below.
checkpoint_every <- 10L
seed <- 20261022L
resume <- TRUE
RUN_CHECKS <- TRUE
version <- 'LN_boundary_aware_v5_publication'
# Engineering guards: reaching them returns UNRESOLVED, not a fitted bound.
max_t_bound <- 80
max_log_sigma <- 16
min_log_sigma <- -12
# Acceptance tolerances apply to the per-record likelihood score and the
# absolute log-likelihood comparison, respectively. integration_tolerance is
# the relative tolerance supplied to numerical quadrature.
score_tolerance <- 1e-6
likelihood_tolerance <- 1e-6
integration_tolerance <- 1e-8
measures <- c('VaR','CTE','GS','PHT','WT')
methods <- c('Delta','Bootstrap_percentile','Profile_likelihood')
# Work with U=log(X-x0)-yd conditional on U>=0. Exact U values are observed
# below H; an observation at H represents censoring. zb is the (1-beta) normal
# quantile; zcrit and cutoff calibrate normal and one-df LR intervals.
yd <- log(d-x0); yc <- log(u-x0); H <- yc-yd
zb <- qnorm(beta,lower.tail=FALSE)
zcrit <- qnorm((1+confidence)/2); cutoff <- qchisq(confidence,1)
stopifnot(x0>0,d>x0,u>d,sigma_true>0,beta>0,beta<1,
          delta>=0,delta<.5,r>0,r<=1,lambda>0,B>=20,N_MC>0,
          bootstrap_required_fraction>0,bootstrap_required_fraction<=1)

# 2. OUTPUT LOCATION
# Locate the script under Rscript, source(), or RStudio, falling back to the
# working directory. rstudioapi is optional; all statistical calculations use
# base R. Output files are placed in a version/mode-specific subdirectory.
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
outdir <- file.path(folder,paste0('MC_',version,'_',MODE))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

# Aggregated diagnostics. Counts of optimizer attempts are not fit counts.
D <- new.env(parent=emptyenv()); D$counts <- numeric(0)
# Increment a diagnostic event counter. scope distinguishes outer fitting,
# bootstrap fitting and interval construction. A dataset can generate several
# optimizer events, so these totals must not be interpreted as dataset counts.
# The environment permits diagnostic updates inside nested fitting functions.
record <- function(scope,reason,count=1) {
  key <- paste(scope,reason,sep='|'); v <- D$counts
  if(!key %in% names(v)) v[key] <- 0
  v[key] <- v[key]+count; D$counts <- v
}
# Export the accumulated scope/reason counters for the current sample-size
# cell. They are also stored in checkpoints so a resumed run preserves them.
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
# 3. CONDITIONAL MODEL AND SAMPLING
# Z(a,b) is the integral of exp(-a*U-b*U^2/2) over U>=0.
# M1 and M2 are raw conditional moments E[U] and E[U^2]. These moments
# provide derivatives of log Z for the likelihood score.
normalizer <- function(a,b) {
  if(b==0) {
    if(a<=0) stop('Improper limiting density')
    return(list(logZ=-log(a),M1=1/a,M2=2/a^2))
  }
  if(!is.finite(a) || !is.finite(b) || b<0) stop('Invalid natural parameters')
  z <- a/sqrt(b)
  if(z>20) {
    # Use a large-positive-z expansion to avoid cancellation in the inverse
    # Mills ratio and its excess over z near the exponential boundary.
    w <- b/a^2
    e <- 1+w*(-2+w*(10+w*(-74+w*(706-8162*w))))
    return(list(logZ=-log(a)-log1p(w*e),M1=e/a,
      M2=(2+w*(-10+w*(74+w*(-706+8162*w))))/a^2))
  }
  # Compute the inverse Mills ratio on the log scale to retain normal-tail
  # accuracy. Subtracting z gives the standardized conditional excess mean.
  lm <- dnorm(z,log=TRUE)-pnorm(z,lower.tail=FALSE,log.p=TRUE)
  excess <- exp(lm)-z
  list(logZ=-log(b)/2-lm,M1=excess/sqrt(b),M2=(1-z*excess)/b)
}
# T1 and T2 include censored exposure H and H^2, respectively.
# s contains n records, m exact observations, k censored observations and
# T1/T2, the sums of observed U and U^2. Exact terms use the density; censored
# terms use P(U>=H). Return both total log likelihood and its (a,b) score.
likelihood <- function(a,b,s) {
  z0 <- normalizer(a,b); zh <- normalizer(a+b*H,b)
  ll <- -a*s$T1-b*s$T2/2-s$m*z0$logZ+s$k*(zh$logZ-z0$logZ)
  ga <- -s$T1+s$n*z0$M1-s$k*zh$M1
  gb <- -s$T2/2+s$n*z0$M2/2-s$k*(zh$M2+2*H*zh$M1)/2
  list(ll=ll,score=c(ga,gb))
}
# Generate a conditional sample under either a finite or limiting fitted
# model, then retain only the sufficient statistics needed for estimation.
# This same sampler serves the outer study and the parametric bootstrap.
simulate_statistics <- function(fit,n) {
  if(fit$type=='limit') exposure <- rexp(n,rate=fit$a) else {
    sigma <- 1/sqrt(fit$b); zd <- fit$a/sqrt(fit$b)
    exponential <- rexp(n)
    if(zd>40) {
      # Invert the conditional cumulative hazard directly. This avoids
      # subtracting almost equal enormous normal quantiles near b=0.
      logZ0 <- normalizer(fit$a,fit$b)$logZ
      # For Exp(1) draws E, solve cumulative_hazard(U)=E by Newton updates.
      # Start from the local exponential approximation and halve any proposed
      # nonpositive exposure. The residual criterion verifies the inversion.
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
      # Ordinary inverse-survival sampling: log survival decreases by an
      # Exp(1) draw. Using log probabilities avoids underflow in small tails.
      z <- qnorm(pnorm(zd,lower.tail=FALSE,log.p=TRUE)-exponential,
                 lower.tail=FALSE,log.p=TRUE)
      exposure <- sigma*(z-zd)
    }
  }
  # Apply censoring only after sampling the conditional exposure. An exact
  # observation contributes its U; a censored record contributes H to T1 and
  # H^2 to T2, together with its censoring indicator through k.
  observed <- pmin(exposure,H); exact <- exposure<H
  list(n=n,m=sum(exact),k=sum(!exact),T1=sum(observed),T2=sum(observed^2))
}

# Optimization coordinates (t,b): a=exp(t)-b*exp(-t).
# For b>0 every real a is reachable; at b=0 a=exp(t)>0.
# 4. EXTENDED MAXIMUM LIKELIHOOD ESTIMATION
# Invert the optimization coordinate transformation. The alternative formula
# for negative a avoids cancellation when sqrt(a^2+4b) is close to -a.
a_to_t <- function(a,b) {
  if(b==0) return(log(a))
  root <- sqrt(a*a+4*b)
  log(if(a>=0) (a+root)/2 else 2*b/(root-a))
}
# Return list(type,a,b,ll) for an accepted finite or boundary fit; return NULL
# when the fit cannot be resolved. A boundary fit is an explicit model limit,
# not a finite sigma selected at an artificial optimization cap.
fit_model <- function(s,warm=NULL,scope='outer') {
  if(s$m==0 || !is.finite(s$T1) || s$T1<=0) {
    record(scope,'unresolved_all_censored_or_invalid'); return(NULL)
  }
  # At b=0, the censored exponential likelihood is maximized by a=m/T1.
  # Use its exact likelihood as a candidate alongside multiple finite starts.
  a0 <- s$m/s$T1
  boundary <- likelihood(a0,0,s)
  starts <- list(c(log(a0),0),c(a_to_t(a0,.5),.5),c(a_to_t(a0,2),2))
  # Bootstrap fitting can start from the original finite fit. Keep the other
  # starts as well, because a warm start alone need not locate the maximum.
  if(!is.null(warm) && warm$type=='finite')
    starts <- c(list(c(a_to_t(warm$a,warm$b),warm$b)),starts)
  candidates <- list()
  bound <- 20
  repeat {
    last <- NULL; cached <- NULL
    # Cache the likelihood and transformed gradient because nlminb may request
    # both at the same point. The chain rule converts the (a,b) score to (t,b);
    # negation and division by n match the minimized average negative likelihood.
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
    # Search with b>=0 and a temporary bound on t. Optimizer status is recorded,
    # but acceptance below also checks likelihood and score conditions.
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
    # Discard finite candidates touching the current t guard; retain the best
    # interior candidate across starts and all completed search expansions.
    finite <- Filter(function(v) v$b>0 && abs(v$par[1])<bound-1e-4,candidates)
    best <- if(length(finite)) finite[[which.max(vapply(finite,function(v) v$ll,numeric(1)))]] else NULL
    bestll <- if(is.null(best)) -Inf else best$ll
    # An exact boundary maximum must satisfy its one-sided score condition.
    if(boundary$score[2]<=0 && boundary$ll>=bestll-likelihood_tolerance) {
      record(scope,'accepted_limit')
      return(list(type='limit',a=a0,b=0,ll=boundary$ll))
    }
    # A finite candidate must match or improve the boundary likelihood within
    # tolerance and have a sufficiently small score in natural coordinates.
    if(!is.null(best) && bestll>=boundary$ll-likelihood_tolerance) {
      score <- likelihood(best$a,best$b,s)$score/s$n
      if(all(is.finite(score)) && max(abs(score))<score_tolerance) {
        record(scope,'accepted_finite')
        return(list(type='finite',a=best$a,b=best$b,ll=best$ll))
      }
    }
    # Expand the temporary search region when acceptance remains unresolved.
    # Hitting the final engineering guard records a failure rather than treating
    # the guard value as an estimate.
    if(bound>=max_t_bound) break
    bound <- min(max_t_bound,bound*2)
    if(!is.null(best)) starts <- c(list(best$par),starts)
  }
  record(scope,'unresolved_optimizer_or_score'); NULL
}

# 5. GROUND-UP RISK FUNCTIONALS AND DERIVATIVES
# Return M(z)=phi(z)/P(Z>z), log M(z), and M(z)-z. The expansion preserves
# the small excess when direct subtraction would lose precision.
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
# For R=x0+exp(mu)*A(sigma), return the leading quadratic coefficient of
# log A, its remainder, and the derivative of that remainder with respect
# to sigma. Keeping the leading term separate stabilizes boundary calculations.
# The derivative is used in the delta method through risk_value().
factor_parts <- function(sigma,m,derivative=FALSE) {
  if(m=='VaR') return(list(coefficient=0,remainder=sigma*zb,prime=zb))
  if(m=='CTE') return(list(coefficient=.5,
    remainder=pnorm(sigma-zb,log.p=TRUE)-log(beta),prime=mills(zb-sigma)$value))
  if(m=='WT') return(list(coefficient=.5,remainder=lambda*sigma,prime=lambda))
  # GS is evaluated by normal-density quadrature. Splitting at zero helps the
  # integrator capture the central mass when the lower limit is far negative.
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
  # The remaining branch is PHT (callers supply one of the configured measures).
  # For r=1 it reduces to the ordinary lognormal mean. Otherwise center and
  # rescale the integrand and remove a log shift before integration.
  if(r==1) return(list(coefficient=.5,remainder=0,prime=0))
  center <- sigma/r; scale <- 1/sqrt(r); shift <- (1-r)*mills(center)$log
  f <- function(t) vapply(t,function(tt)
    exp(-tt^2/2+(1-r)*mills(center+scale*tt)$log-shift),numeric(1))
  I <- integrate(f,-Inf,Inf,rel.tol=integration_tolerance)$value
  prime <- 0
  # Compute the extra PHT integral only when derivatives are requested;
  # ordinary risk evaluation and profiling do not need this work.
  if(derivative) {
    J <- integrate(function(t) f(t)*vapply(center+scale*t,
      function(z) mills(z)$excess,numeric(1)),-Inf,Inf,rel.tol=integration_tolerance)$value
    prime <- (1-r)*J/(r*I)
  }
  list(coefficient=1/(2*r),remainder=log(r)/2-r*log(2*pi)/2+shift+log(I),prime=prime)
}
# The leading quadratic coefficient determines whether the ground-up risk
# collapses to x0 or diverges as b approaches zero at fixed positive a.
critical_slope <- function(m) if(m=='VaR') 0 else if(m=='PHT') 1/(2*r) else .5
# Return the risk value, its gradient with respect to (a,b), and log(R-x0).
# Boundary gradients are unavailable. At critical slopes return NA rather
# than assign a potentially path-dependent limiting risk.
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
  # Substitute mu=yd-a/b and sigma^2=1/b into the ground-up risk formula.
  # A finite fit whose exponent exceeds floating-point range is unresolved
  # (NA), distinct from a verified mathematical infinite boundary value.
  q <- yd+(f$coefficient-fit$a)/fit$b+f$remainder
  if(!is.finite(q) || q>log(.Machine$double.xmax))
    return(list(value=NA_real_,gradient=c(NA_real_,NA_real_),logexcess=q))
  excess <- exp(q)
  # Differentiate log(R-x0) with respect to a and b, then multiply by the
  # excess risk. Zero-underflow excesses do not receive a fabricated gradient.
  g <- if(gradient && excess>0) excess*c(-1/fit$b,
    (fit$a-f$coefficient)/fit$b^2-f$prime/(2*fit$b^(3/2))) else c(NA_real_,NA_real_)
  list(value=x0+excess,gradient=g,logexcess=q)
}
# Evaluate all configured measures independently. An integration or numerical
# exception for one measure becomes NA without discarding the other measures.
all_risks <- function(fit) vapply(measures,function(m)
  tryCatch(risk_value(fit,m)$value,error=function(e) NA_real_),numeric(1))

# Expected information per conditional/censored record in (a,b) coordinates.
# Integrate exact-record score products, then add the censored point mass.
# 6. DELTA AND PARAMETRIC BOOTSTRAP INTERVALS
# The returned information matrix is per observed conditional record. Its
# inverse divided by n approximates the covariance of the finite MLE.
# INFORMATION CHOICE: expected Fisher information, evaluated at the fitted
# parameters (a_hat,b_hat), is used here. This function takes the expectation
# of score outer products under the fitted conditional/censored distribution:
# integrate over exact outcomes, then add the censored outcome contribution.
# It does NOT use observed information, i.e., the negative Hessian of the
# realized sample log likelihood. The expectation is over possible observed
# records under the fitted model, rather than over the realized dataset alone.
expected_information <- function(fit) {
  a <- fit$a; b <- fit$b
  z0 <- normalizer(a,b); zh <- normalizer(a+b*H,b)
  logpc <- -a*H-b*H^2/2+zh$logZ-z0$logZ
  # cens is the score for the discrete censored outcome, whose probability is
  # exp(logpc). sc(t) is the score for an exact observation at exposure t.
  cens <- c(z0$M1-H-zh$M1,(z0$M2-H^2-2*H*zh$M1-zh$M2)/2)
  sc <- function(t) cbind(z0$M1-t,(z0$M2-t^2)/2)
  I <- matrix(0,2,2)
  # Integrate only the upper triangle and reflect it to ensure symmetry.
  # Both the continuous exact component and censored point mass are required.
  for(i in 1:2) for(j in i:2) {
    value <- integrate(function(t) {
      s <- sc(t); exp(-a*t-b*t^2/2-z0$logZ)*s[,i]*s[,j]
    },0,H,rel.tol=integration_tolerance)$value+exp(logpc)*cens[i]*cens[j]
    I[i,j] <- I[j,i] <- value
  }
  I
}
# Construct a symmetric Wald interval Rhat +/- zcrit*SE on the risk scale.
# No clipping to x0 is applied. Boundary original fits are marked unavailable;
# nonfinite gradients, singular information or nonfinite SEs are failures.
# With g=gradient R(a_hat,b_hat) and per-record expected information I,
# SE(R_hat) = sqrt(g' %*% solve(I) %*% g / n).
# Thus the total-sample information is n*I, evaluated at the finite fitted
# parameters. This information matrix is used only for delta intervals:
# bootstrap intervals use simulated refits, and profile likelihood intervals
# use likelihood-ratio inversion, without this information-matrix calculation.
delta_ci <- function(fit,n,m,I=NULL) {
  if(fit$type=='limit') return(list(ci=c(NA_real_,NA_real_),reason='delta_limit_unavailable'))
  rv <- risk_value(fit,m,TRUE)
  if(any(!is.finite(c(rv$value,rv$gradient)))) stop('Nonfinite delta derivative')
  if(is.null(I)) I <- expected_information(fit)
  # Cholesky solves avoid explicitly inverting a nearly singular matrix.
  # If I=t(R)%*%R, solving t(R)*v=gradient gives
  # sum(v^2)=gradient' * solve(I) * gradient, the delta variance multiplier.
  R <- chol(I); v <- forwardsolve(t(R),rv$gradient)
  se <- sqrt(sum(v^2)/n)
  if(!is.finite(se)) stop('Nonfinite delta standard error')
  list(ci=rv$value+c(-1,1)*zcrit*se,reason='ok')
}
# For each original dataset, simulate B datasets from its accepted fit,
# refit the same extended model, and evaluate all risks. Failed draws remain
# missing and are never replaced by newly simulated datasets.
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
  # Rows correspond to measures and columns to lower/upper endpoints.
  # Quantile type 1 selects order statistics without interpolating finite and
  # infinite values. This differs intentionally from type 7 length summaries.
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
# 7. PROFILE LIKELIHOOD INTERVALS
# For a fixed log excess risk q, solve the risk constraint for a at each sigma,
# then maximize over log sigma. Include the admissible limiting likelihood
# so that a finite search region cannot hide a boundary supremum.
profile_ll <- function(q,m,s) {
  critical <- critical_slope(m)
  limiting <- if(critical>0) likelihood(critical,0,s)$ll else -Inf
  objective <- function(ls) {
    sigma <- exp(ls); b <- exp(-2*ls); f <- factor_parts(sigma,m)
    # The constraint q=yd+(coefficient-a)/b+remainder determines a.
    # A large objective penalty lets the search skip nonfinite evaluations.
    a <- f$coefficient+b*(yd-q+f$remainder)
    val <- likelihood(a,b,s)$ll
    if(is.finite(val)) -val else 1e100
  }
  low <- -3; high <- 3
  repeat {
    # Scan a grid for candidate local optima, refine their neighborhoods, and
    # compare every result with the limiting likelihood. This reduces dependence
    # on a single unimodal optimization assumption.
    grid <- seq(low,high,length.out=max(31,ceiling((high-low)/.25)+1))
    vals <- vapply(grid,function(z) tryCatch(objective(z),error=function(e) 1e100),numeric(1))
    idx <- unique(c(which(vals<=c(Inf,head(vals,-1)) & vals<=c(tail(vals,-1),Inf)),which.min(vals)))
    opts <- lapply(idx,function(i) tryCatch(optimize(objective,
      grid[c(max(1,i-1),min(length(grid),i+1))],tol=1e-7),error=function(e) NULL))
    opts <- Filter(Negate(is.null),opts)
    best <- max(c(limiting,-vals,if(length(opts))
      -vapply(opts,function(o) o$objective,numeric(1)) else -Inf))
    # Stop extending each side only when its edge is well below the best
    # likelihood, or the upper side has converged to the known limit.
    # Failure to establish this before the guards is reported as unresolved.
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
# Invert the one-df likelihood-ratio cutoff. First test whether risk limits
# x0 and infinity are admitted, then search for any required finite endpoints.
# The returned interval follows the component located from an admitted point;
# this is not an exhaustive search for every disconnected profile component.
profile_ci <- function(fit,m,s) {
  a0 <- s$m/s$T1; critical <- critical_slope(m)
  # Maximize the boundary likelihood over slopes associated with the lower
  # and upper risk limits. VaR has critical slope zero and no proper positive-a
  # exponential boundary yielding infinite VaR in this parameterization.
  zero_ll <- likelihood(max(a0,critical),0,s)$ll
  infinity_ll <- if(critical>0) likelihood(min(a0,critical),0,s)$ll else -Inf
  admitted <- function(ll) is.finite(ll) && 2*(fit$ll-ll)<=cutoff+1e-7
  lower_limit <- admitted(zero_ll); upper_limit <- admitted(infinity_ll)
  if(lower_limit && upper_limit) return(c(x0,Inf))
  # The LR region satisfies profile_ll >= fitted_ll-cutoff/2.
  # A profile exceeding the fitted maximum beyond tolerance signals an
  # inconsistent optimization result and is rejected by dev().
  threshold <- fit$ll-cutoff/2
  # A data-based central q; no generating parameters enter profiling.
  ref <- if(fit$type=='finite') risk_value(fit,m)$logexcess else yd
  if(!is.finite(ref)) ref <- yd
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
  # Expand away from the admitted center until the LR threshold is crossed,
  # then solve the crossing with uniroot. Guards prevent endless expansion;
  # endpoint overflow is unresolved rather than relabeled as a verified limit.
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
# 8. DETERMINISTIC NUMERICAL CHECKS
# Check likelihood derivatives by central differences at several finite
# parameter pairs, the exact boundary likelihood and its one-sided b score,
# risk gradients, representative limits, and positive-definite information.
# These checks assess numerical identities, not Monte Carlo coverage.
run_checks <- function() {
  ss <- list(n=5,m=4,k=1,T1=sum(c(.1,.3,.6,1,H)),T2=sum(c(.1,.3,.6,1,H)^2))
  for(par in list(c(1,.5),c(.2,.1),c(-.5,2),c(4,.01))) {
    analytical <- likelihood(par[1],par[2],ss)$score
    numerical <- vapply(1:2,function(j) {
      eps <- 1e-5; plus <- minus <- par; plus[j] <- plus[j]+eps; minus[j] <- minus[j]-eps
      (likelihood(plus[1],plus[2],ss)$ll-likelihood(minus[1],minus[2],ss)$ll)/(2*eps)
    },numeric(1))
    stopifnot(max(abs(analytical-numerical)/(1+abs(numerical)))<1e-5)
  }
  a0 <- ss$m/ss$T1
  stopifnot(abs(likelihood(a0,0,ss)$ll-(ss$m*log(a0)-a0*ss$T1))<1e-10)
  numerical_b0 <- (likelihood(a0,1e-6,ss)$ll-likelihood(a0,0,ss)$ll)/1e-6
  stopifnot(abs(numerical_b0-likelihood(a0,0,ss)$score[2])/(1+abs(numerical_b0))<1e-4)
  true <- list(type='finite',a=(yd-mu_true)/sigma_true^2,b=1/sigma_true^2)
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
  chol(expected_information(true))
  cat('Likelihood scores, risk derivatives, limit values and expected-information checks passed.\n')
}
if(RUN_CHECKS) run_checks()

# 9. INTERVAL CLASSIFICATION AND CHECKPOINTED SUMMARY STATE
# One row is maintained for each risk/method combination. Valid/Covered follow
# the final infinite-CI exclusion switch; Finite_Valid/Finite_Covered always
# restrict to intervals with finite endpoints. Length counts are separate:
# even finite endpoints can yield an overflowing computed difference.
new_accumulator <- function() {
  g <- expand.grid(Risk_Measure=measures,Method=methods,stringsAsFactors=FALSE)
  for(key in c('Valid','Covered','Finite_Length_Count','Infinite_Length_Count','Undefined_Length_Count')) g[[key]] <- 0
  g$Infinite_CI_Excluded <- g$Delta_Limit_Unavailable <- g$Unresolved_CI <- 0
  g$Finite_Valid <- g$Finite_Covered <- g$Infinite_CI_Count <- 0
  g$Sum_Finite_Length <- 0
  # Store lengths and WT records in checkpoints for exact quantiles/diagnostics.
  # List columns preserve all finite computed lengths for exact empirical
  # quantiles. WT_Records stores only the ten largest finite lengths per row,
  # with generating replication and fitted parameters for later inspection.
  g$Finite_Lengths <- I(rep(list(numeric(0)),nrow(g)))
  g$WT_Records <- I(rep(list(NULL),nrow(g)))
  g
}
# Classify one method/measure interval. Invalid shape, NA/NaN or reversed
# endpoints take priority over infinite-endpoint classification. Each outer
# fit failure is handled separately in run_cell and never enters update().
update <- function(acc,method,m,ci,reason='ok',replication=NA_integer_,fit=NULL) {
  i <- which(acc$Method==method & acc$Risk_Measure==m)
  if(length(ci)!=2 || anyNA(ci) || ci[1]>ci[2]) {
    key <- if(reason=='delta_limit_unavailable') 'Delta_Limit_Unavailable' else 'Unresolved_CI'
    acc[[key]][i] <- acc[[key]][i]+1
    record('CI',paste0(method,'_',m,'_',reason)); return(acc)
  }
  # Count well-formed infinite-endpoint intervals before applying the optional
  # exclusion. This makes the reported infinite percentage independent of the
  # switch that governs the broader accepted-interval coverage calculation.
  infinite <- any(is.infinite(ci))
  if(infinite) acc$Infinite_CI_Count[i] <- acc$Infinite_CI_Count[i]+1
  if(EXCLUDE_INFINITE_CIS && infinite) {
    acc$Infinite_CI_Excluded[i] <- acc$Infinite_CI_Excluded[i]+1
    record('CI',paste0(method,'_',m,'_infinite_excluded')); return(acc)
  }
  # Coverage includes both endpoints. Infinite accepted intervals can cover a
  # finite truth, but never enter the finite-valid coverage denominator.
  acc$Valid[i] <- acc$Valid[i]+1
  acc$Covered[i] <- acc$Covered[i]+as.integer(ci[1]<=true_risks[m] && true_risks[m]<=ci[2])
  if(!infinite) {
    acc$Finite_Valid[i] <- acc$Finite_Valid[i]+1
    acc$Finite_Covered[i] <- acc$Finite_Covered[i]+
      as.integer(ci[1]<=true_risks[m] && true_risks[m]<=ci[2])
  }
  # Record computed lengths separately as undefined (for example Inf-Inf),
  # infinite, or finite. Only finite lengths enter the median, tail quantiles,
  # maximum and finite mean. Infinite exclusions return before length counting.
  len <- ci[2]-ci[1]
  if(is.na(len)) acc$Undefined_Length_Count[i] <- acc$Undefined_Length_Count[i]+1 else
    if(is.infinite(len)) acc$Infinite_Length_Count[i] <- acc$Infinite_Length_Count[i]+1 else {
      acc$Finite_Lengths[[i]] <- c(acc$Finite_Lengths[[i]],len)
      # Retain the ten largest finite WT intervals for each method. Re-sort after
      # each insertion; checkpointing these records makes resume deterministic.
      if(m=='WT') {
        row <- data.frame(Replication=replication,Lower=ci[1],Upper=ci[2],Length=len,
          Covered=ci[1]<=true_risks[m] && true_risks[m]<=ci[2],
          Fit_Type=if(is.null(fit)) NA_character_ else fit$type,
          a=if(is.null(fit)) NA_real_ else fit$a,b=if(is.null(fit)) NA_real_ else fit$b)
        rows <- rbind(acc$WT_Records[[i]],row)
        acc$WT_Records[[i]] <- head(rows[order(rows$Length,decreasing=TRUE),],10L)
      }
      acc$Finite_Length_Count[i] <- acc$Finite_Length_Count[i]+1
      acc$Sum_Finite_Length[i] <- acc$Sum_Finite_Length[i]+len
    }
  acc
}
# 10. GENERATING TRUTH AND REPRODUCIBILITY METADATA
# Compute the ground-up truth once under the finite generating distribution.
# It is used only to score coverage, not to construct fitted intervals.
true_fit <- list(type='finite',a=(yd-mu_true)/sigma_true^2,b=1/sigma_true^2)
true_risks <- all_risks(true_fit)
stopifnot(all(is.finite(true_risks)))
write.csv(data.frame(Risk_Measure=measures,True_Risk=unname(true_risks)),
          file.path(outdir,'True_Risks.csv'),row.names=FALSE)
# Store all study settings plus the script text, when the script is located.
# Checkpoint compatibility uses exact equality, so even comment edits can
# invalidate a prior checkpoint when source text is available.
source_text <- if(length(script)==1 && file.exists(script)) readLines(script,warn=FALSE) else NULL
configuration <- list(version=version,MODE=MODE,N_MC=N_MC,B=B,n=sample_sizes,
  x0=x0,d=d,u=u,mu=mu_true,sigma=sigma_true,beta=beta,delta=delta,r=r,lambda=lambda,
  confidence=confidence,bootstrap_required_fraction=bootstrap_required_fraction,
  exclude_infinite_CIs=EXCLUDE_INFINITE_CIS,seed=seed,source=source_text,
  max_t_bound=max_t_bound,max_log_sigma=max_log_sigma,min_log_sigma=min_log_sigma,
  score_tolerance=score_tolerance,likelihood_tolerance=likelihood_tolerance,
  integration_tolerance=integration_tolerance)
saveRDS(configuration,file.path(outdir,'Configuration.rds'))
# Write the new state to a temporary RDS file, back up the prior checkpoint,
# then rename the temporary file into place. If promotion fails, attempt to
# restore the backup. These operations are not a substitute for disk backups.
save_checkpoint <- function(state,path) {
  tmp <- paste0(path,'.tmp'); backup <- paste0(path,'.bak'); saveRDS(state,tmp)
  if(file.exists(backup)) unlink(backup)
  if(file.exists(path) && !file.rename(path,backup)) stop('Cannot back up checkpoint')
  if(!file.rename(tmp,path)) {if(file.exists(backup)) file.rename(backup,path); stop('Cannot save checkpoint')}
  if(file.exists(backup)) unlink(backup)
}
# 11. OUTER MONTE CARLO LOOP FOR ONE SAMPLE SIZE
# State contains completed replications, original-fit classifications,
# interval accumulators, diagnostics and the exact configuration. Reloading a
# compatible checkpoint resumes at the next uncompleted replication.
run_cell <- function(n,cell_id) {
  path <- file.path(outdir,sprintf('checkpoint_n%d.rds',n)); backup <- paste0(path,'.bak')
  state <- list(done=0,finite_fits=0,limit_fits=0,unresolved_fits=0,
                acc=new_accumulator(),diagnostics=numeric(0),configuration=configuration)
  if(resume && (file.exists(path) || file.exists(backup))) {
    state <- readRDS(if(file.exists(path)) path else backup)
    if(!identical(state$configuration,configuration)) stop('Checkpoint configuration differs')
  }
  D$counts <- state$diagnostics
  if(state$done<N_MC) for(replication in seq.int(state$done+1,N_MC)) {
    # Each cell/replication has its own deterministic seed. Replaying a
    # replication therefore does not require storing the previous RNG state.
    # Bootstrap randomness is consumed within that replication's random stream.
    set.seed(seed+cell_id*1000000+replication)
    s <- simulate_statistics(true_fit,n); fit <- fit_model(s)
    # An unresolved original fit contributes to every method/measure failure
    # percentage. Do not redraw it: replacing difficult datasets would alter the
    # sampling distribution and conceal operational failures.
    if(is.null(fit)) state$unresolved_fits <- state$unresolved_fits+1 else {
      if(fit$type=='finite') state$finite_fits <- state$finite_fits+1 else state$limit_fits <- state$limit_fits+1
      # Reuse one expected-information matrix and one joint bootstrap simulation
      # across all risk measures. Each measure still has its own CI outcome.
      I <- if(fit$type=='finite') tryCatch(expected_information(fit),error=function(e) NULL) else NULL
      boot <- bootstrap_cis(fit,n)
      for(j in seq_along(measures)) {
        m <- measures[j]
        # Convert delta failures to a missing interval with a reason label. Profile
        # exceptions are also retained as missing outcomes and their messages are
        # sanitized for the scope/reason diagnostic key.
        dc <- if(fit$type=='finite' && is.null(I))
          list(ci=c(NA_real_,NA_real_),reason='delta_expected_information_failed') else
          tryCatch(delta_ci(fit,n,m,I),error=function(e)
            list(ci=c(NA_real_,NA_real_),reason='delta_numerical_failure'))
        record('outer',dc$reason)
        pc <- tryCatch(profile_ci(fit,m,s),error=function(e) {
          record('profile',gsub('[|\r\n]',' ',conditionMessage(e))); c(NA_real_,NA_real_)
        })
        state$acc <- update(state$acc,'Delta',m,dc$ci,dc$reason,replication,fit)
        state$acc <- update(state$acc,'Bootstrap_percentile',m,boot[j,],replication=replication,fit=fit)
        state$acc <- update(state$acc,'Profile_likelihood',m,pc,replication=replication,fit=fit)
      }
    }
    # Advance completion even when fitting or interval construction fails.
    # Persist the diagnostics alongside the interval and fit counts.
    state$done <- replication
    if(replication %% checkpoint_every==0 || replication==N_MC) {
      state$diagnostics <- D$counts; save_checkpoint(state,path); write_diagnostics(n)
      cat('LN n=',n,': ',replication,'/',N_MC,'\n',sep=''); flush.console()
    }
  }
  # 12. PUBLICATION SUMMARIES AND DIAGNOSTIC EXPORTS
  # Extract list columns before CSV export; the checkpoint retains them.
  # Length quantiles use R type 7 interpolation and remain in original monetary
  # units unless a column explicitly ends in _Thousand.
  acc <- state$acc
  finite_lengths <- acc$Finite_Lengths
  wt_records <- acc$WT_Records
  acc$Finite_Lengths <- acc$WT_Records <- NULL
  qlength <- function(prob) vapply(finite_lengths,function(v)
    if(length(v)) unname(quantile(v,prob,type=7)) else NA_real_,numeric(1))
  acc$Median_Finite_Length <- qlength(.5)
  acc$P95_Finite_Length <- qlength(.95)
  acc$P99_Finite_Length <- qlength(.99)
  acc$Max_Finite_Length <- vapply(finite_lengths,function(v)
    if(length(v)) max(v) else NA_real_,numeric(1))
  # Primary coverage is conditional on a finite valid interval being produced:
  # p=Finite_Covered/Finite_Valid. Its plug-in binomial Monte Carlo SE uses that
  # same effective denominator. With no finite valid intervals both are NA.
  # This conditional coverage alone does not measure unconditional reliability.
  acc$Coverage_Among_Finite_Valid_Intervals <- ifelse(acc$Finite_Valid>0,
    acc$Finite_Covered/acc$Finite_Valid,NA_real_)
  acc$MCSE_Coverage_Among_Finite_Valid_Intervals <- with(acc,
    sqrt(Coverage_Among_Finite_Valid_Intervals*(1-Coverage_Among_Finite_Valid_Intervals)/Finite_Valid))
  # Use all completed outer replications for both percentages. Failure here
  # means an unresolved original fit, unavailable delta boundary CI, or invalid/
  # unresolved CI. Undefined lengths of accepted infinite intervals are recorded
  # separately and are not added again to this disjoint failure category.
  acc$Percent_Infinite_Intervals <- 100*acc$Infinite_CI_Count/state$done
  acc$Undefined_Or_Numerical_Failure_Count <- acc$Delta_Limit_Unavailable+
    acc$Unresolved_CI+state$unresolved_fits
  acc$Percent_Undefined_Or_Numerical_Failures <-
    100*acc$Undefined_Or_Numerical_Failure_Count/state$done
  # Export up to ten ranked finite WT intervals per method when any exist.
  # These records diagnose large lengths; they do not change exclusions.
  wt <- lapply(which(acc$Risk_Measure=='WT'),function(i) {
    rows <- wt_records[[i]]
    if(is.null(rows)) return(NULL)
    cbind(n=n,Method=acc$Method[i],Rank=seq_len(nrow(rows)),rows)
  })
  wt <- do.call(rbind,wt)
  if(!is.null(wt)) write.csv(wt,file.path(outdir,sprintf('WT_Largest_Finite_Intervals_n%d.csv',n)),row.names=FALSE)
  # Attach original-fit counts and verify that every completed replication is
  # accounted for exactly once in each row: accepted, infinite-excluded,
  # delta-unavailable, unresolved CI, or unresolved original fit.
  # Failed is an operational nonacceptance count and includes exclusions.
  acc$n <- n; acc$Completed <- state$done
  acc$Original_Finite_Fits <- state$finite_fits; acc$Original_Limit_Fits <- state$limit_fits
  acc$Original_Unresolved_Fits <- state$unresolved_fits; acc$Failed <- state$done-acc$Valid
  stopifnot(all(acc$Completed==acc$Valid+acc$Infinite_CI_Excluded+
    acc$Delta_Limit_Unavailable+acc$Unresolved_CI+acc$Original_Unresolved_Fits))
  # Supplementary coverage follows the exclusion switch. Operational success
  # uses all completed replications, treating exclusions and failures as misses.
  # Report these with their full labels because their denominators differ.
  acc$Coverage_Among_Accepted_Intervals <- ifelse(acc$Valid>0,acc$Covered/acc$Valid,NA_real_)
  acc$Coverage_Success_Over_All <- acc$Covered/state$done
  acc$MCSE_Coverage_Among_Accepted_Intervals <- sqrt(acc$Coverage_Among_Accepted_Intervals*(1-acc$Coverage_Among_Accepted_Intervals)/acc$Valid)
  # The accepted-interval mean is undefined if any accepted length is undefined,
  # and infinite if an accepted length is infinite. The finite mean always uses
  # only finite computed lengths; rare enormous lengths can dominate it.
  acc$Mean_Length_Valid <- ifelse(acc$Valid==0 | acc$Undefined_Length_Count>0,NA_real_,
    ifelse(acc$Infinite_Length_Count>0,Inf,acc$Sum_Finite_Length/acc$Valid))
  acc$Mean_Finite_Length <- ifelse(acc$Finite_Length_Count>0,
    acc$Sum_Finite_Length/acc$Finite_Length_Count,NA_real_)
  acc$Mean_Length_Valid_Thousand <- acc$Mean_Length_Valid/1000
  acc$Mean_Finite_Length_Thousand <- acc$Mean_Finite_Length/1000
  write.csv(acc,file.path(outdir,sprintf('Summary_LN_n%d.csv',n)),row.names=FALSE,na='NA')
  # The compact main table emphasizes finite-valid coverage/MCSE, median
  # precision and outcome percentages. Supplementary and legacy Summary files
  # retain the full scalar table, including detailed counts and tail lengths.
  main <- c('n','Risk_Measure','Method','Completed','Finite_Valid',
    'Coverage_Among_Finite_Valid_Intervals','MCSE_Coverage_Among_Finite_Valid_Intervals',
    'Median_Finite_Length','Percent_Infinite_Intervals','Percent_Undefined_Or_Numerical_Failures')
  write.csv(acc[,main],file.path(outdir,sprintf('Main_Results_LN_n%d.csv',n)),row.names=FALSE,na='NA')
  write.csv(acc,file.path(outdir,sprintf('Supplementary_Results_LN_n%d.csv',n)),row.names=FALSE,na='NA')
  write_diagnostics(n); acc
}
# 13. RUN ALL CELLS AND WRITE COMBINED OUTPUT
# Refresh the combined scalar summary after each completed sample-size cell.
# Per-cell main, supplementary, WT and event-diagnostic files remain separate.
results <- vector('list',length(sample_sizes))
for(i in seq_along(sample_sizes)) {
  results[[i]] <- run_cell(sample_sizes[i],i)
  write.csv(do.call(rbind,results[seq_len(i)]),file.path(outdir,'Monte_Carlo_Summary.csv'),row.names=FALSE,na='NA')
}
# Write a plain-text guide beside the results so exported tables retain
# their denominator definitions, numerical caveats and calibration limits.
writeLines(c(
 'Primary precision measure: median finite interval length; mean finite length is supplementary.',
 'Median describes typical precision; the mean can be dominated by rare enormous intervals.',
 'Main results: coverage among finite valid intervals with plug-in binomial Monte Carlo SE, median finite length, infinite and failure percentages.',
 'Finite coverage denominator is Finite_Valid; infinite intervals and failures are excluded regardless of EXCLUDE_INFINITE_CIS.',
 'MCSE is sqrt(p*(1-p)/Finite_Valid); NA when no finite valid intervals exist.',
 'Infinite and failure percentages use Completed, including unresolved original fits, as denominator.',
 'Failures include unresolved original fits, unavailable limiting delta CIs and unresolved/invalid CIs; these categories are disjoint from valid infinite CIs.',
 'Supplementary results contain mean, P95/P99 (R quantile type 7), maximum finite length and detailed exclusion counts.',
 'WT_Largest_Finite_Intervals lists the ten largest finite intervals per method with replication, endpoints, coverage and fitted natural parameters.',
 'Coverage_Among_Accepted_Intervals follows EXCLUDE_INFINITE_CIS; Coverage_Success_Over_All treats every excluded or failed interval as unsuccessful.',
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
 'Conditional coverage and operational success are separate columns.',
 'This implementation requires local R validation and a paired pilot before production.',
 'Standard percentile bootstrap and chi-square LR calibration may still be inaccurate near boundaries.'),
 file.path(outdir,'README.txt'))
# Show the primary estimates and key failure counts in the R console.
# All fuller summaries and diagnostic records are available in outdir.
print(do.call(rbind,results)[,c('n','Risk_Measure','Method','Valid','Failed',
 'Coverage_Among_Finite_Valid_Intervals','MCSE_Coverage_Among_Finite_Valid_Intervals',
 'Median_Finite_Length','Percent_Infinite_Intervals','Percent_Undefined_Or_Numerical_Failures',
 'Infinite_CI_Excluded','Delta_Limit_Unavailable','Unresolved_CI','Original_Limit_Fits')],row.names=FALSE)
cat('Outputs:',outdir,'\n')
