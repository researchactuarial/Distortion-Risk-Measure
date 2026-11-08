# Deterministic uniform de-grouping of TIES ONLY, assuming rounding to 1,000 NOK.
# k tied observations use equally spaced interior points L+i*(U-L)/(k+1).
# Three 750000 claims become 749750, 750000, and 750250.
# Two 601000 claims become approximately 600833.333 and 601166.667.
# Two 500000 claims become approximately 500166.667 and 500333.333.
# The lower bin endpoint is capped at truncation; no value falls below it.
# Singletons remain unchanged. No random seed is required.
# This is a preprocessing assumption, not additional observed precision;
# the following CIs do not include uncertainty about the rounding mechanism.
# The preprocessing function below implements this rule before either model fit.
# The audit export records the original amount, multiplicity, and rounding bin;
# the claim-level export preserves the pairing of recorded and adjusted values.
# An odd-sized tied group can retain its central recorded value, so the number
# of changed claims can be smaller than the number belonging to tied groups.

# Norwegian Fire Claims: combined Pareto I and shifted-lognormal profile CIs
# Both models use the SAME deterministically de-grouped observations.
# Place beside NFC_Data.xlsx and Source in RStudio or run with Rscript.
# X = x0 + exp(Y), Y ~ N(mu, sigma^2), observed conditional on X >= d.
# x0 is FIXED. For each proposed risk R, maximize the truncated likelihood
# over sigma subject to R(mu,sigma)=R, then invert the likelihood-ratio test.
# Run via RStudio Source or Rscript; results saved beside this file.

# SETTINGS AND MODEL CONVENTIONS
# x0 is the fixed ground-up support/location parameter, whereas d is the
# reporting threshold used to condition the likelihood. They play different
# roles: fitting uses claims above d, but reported risks refer to ground-up X.
x0 <- 100000                  # fixed ground-up location (NOK)
d <- 500000                   # left-truncation threshold (NOK)
# beta is the upper-tail probability for VaR and CTE (0.10 means the
# 90th percentile). delta controls the GS (Gini shortfall) adjustment;
# r is the proportional-hazards transform (PHT) exponent, and lambda is
# the Wang transform (WT) normal-shift parameter. Monetary units are NOK.
beta <- 0.10; delta <- 0.25; r <- 0.95; lambda <- 0.25
# Confidence levels are probabilities, converted to chi-square LR cutoffs.
# The rounding unit determines the assumed width of each recorded-value bin.
levels <- c(0.90, 0.95, 0.99)
rounding_unit <- 1000
# Fail early for the explicit parameter restrictions checked here.
# These checks do not exhaust every possible invalid configuration.
stopifnot(x0 > 0, x0 < d, beta > 0, beta < 1,
          r > 0, r <= 1, rounding_unit > 0)
# Numerical profiling domain for sigma. Enlarge if boundary warnings occur.
sigma_limits <- c(0.02, 8)

# Find file location, falling back to working directory if unavailable.
# INPUT DISCOVERY AND VALIDATION
# Try the active RStudio source document, then sourced-file stack frames,
# then Rscript's --file argument. If none identifies an existing script,
# use the current working directory for both input and generated outputs.
script <- NULL
if (requireNamespace('rstudioapi', quietly=TRUE)) {
  script <- tryCatch(if (rstudioapi::isAvailable())
    rstudioapi::getSourceEditorContext()$path else NULL, error=function(e) NULL)
}
if (length(script)!=1L || is.na(script) || !nzchar(script) || !file.exists(script)) script <- NULL
if (is.null(script)) {
  # Search frames from most recent to oldest for a source() file path.
  for (fr in rev(sys.frames())) if (!is.null(fr$ofile) && file.exists(fr$ofile)) {
    script <- fr$ofile; break
  }
}
if (is.null(script)) {
  # Rscript supplies its script filename as a command-line option.
  a <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  if (length(a)) script <- sub('^--file=', '', a[1])
}
folder <- if (!is.null(script) && file.exists(script)) dirname(normalizePath(script)) else getwd()
input <- file.path(folder, 'NFC_Data.xlsx')
if (!file.exists(input)) stop('Place NFC_Data.xlsx beside the script, or set the working directory to its folder: ', folder)
if (!requireNamespace('readxl', quietly=TRUE)) stop('Install readxl: install.packages("readxl")')
# Read every worksheet cell returned by readxl as a claim candidate.
# Flatten to one numeric vector, remove nonfinite entries, and sort;
# this assumes the workbook's first sheet contains the intended claim data.
x <- as.numeric(unlist(readxl::read_excel(input, col_names=FALSE)))
x <- sort(x[is.finite(x)])
n <- length(x)
# Require sufficient observations and claims inside both model supports.
# The expected-count warning is informational; a different n still runs.
if (n < 10 || !is.finite(x0) || x0 <= 0 || x0 >= d || any(x < d) || any(x <= x0)) stop('Invalid data or x0/d settings')
if (n != 647) warning('Expected 647 claims, found ', n)

# PREPROCESSING IMPLEMENTATION (full assumptions are at the file beginning)
# Return a list containing adjusted observations and one audit row per
# unique recorded amount. Input order and observation count are preserved.
de_group_uniform <- function(values, threshold, rounding_unit=1000) {
  stopifnot(length(values)>0L, all(is.finite(values)),
            is.finite(threshold), is.finite(rounding_unit), rounding_unit>0,
            all(values>=threshold))
  # Group indices by exact numeric equality, retaining access to original rows.
  groups <- split(seq_along(values), match(values, unique(values)))
  adjusted <- values
  # Process each amount once and assemble its bin metadata for the audit.
  audit <- lapply(groups, function(ii) {
    amount <- values[ii[1]]; k <- length(ii)
    # Clip only the lower bin boundary to the reporting threshold.
    lo <- max(threshold, amount-rounding_unit/2)
    hi <- amount+rounding_unit/2
    # Interior points avoid both bin edges. The enclosing assignment (<<-)
    # updates adjusted in the function's environment from inside lapply.
    if (k>1L) adjusted[ii] <<- lo+seq_len(k)*(hi-lo)/(k+1)
    data.frame(Recorded_NOK=amount, Count=k, Bin_Lower_NOK=lo,
               Bin_Upper_NOK=hi, Adjusted=(k>1L))
  })
  stopifnot(length(adjusted)==length(values), all(adjusted>=threshold))
  list(values=adjusted, audit=do.call(rbind,audit))
}
# Keep the recorded vector for a row-by-row comparison with adjusted x.
# All subsequent estimators use the replacement x, with the original n.
raw_claims <- x
degrouped <- de_group_uniform(raw_claims,d,rounding_unit)
x <- degrouped$values
# Export paired observations and the separate amount-level audit;
# existing files with these names in folder will be overwritten.
write.csv(data.frame(Claim_ID=seq_along(x),Recorded_NOK=raw_claims,
                     Degrouped_NOK=x),
          file.path(folder,'NFC_Degrouped_Claims.csv'),row.names=FALSE)
write.csv(degrouped$audit,file.path(folder,'NFC_Degrouping_Audit.csv'),row.names=FALSE)
cat('Uniformly de-grouped',sum(raw_claims!=x),'of',n,'claims before BOTH fits.\n')

# Parameter fitting copied from the parametric bootstrap script.
# For a Pareto I tail conditional on X >= d, log(X/d) is exponential
# with rate alpha. Its MLE is n / sum(log(X/d)); fixed x0 cancels
# from this conditional likelihood and cannot be estimated from it.
fit_pa <- function(x) length(x)/sum(log(x/d))

# Shifted LN: X=x0+exp(Y), Y~N(mu,sigma^2), conditioned on X>=d.
# Estimate the shifted-lognormal parameters by minimizing the negative
# conditional log-likelihood. Optimize log(sigma) to enforce sigma > 0.
# The transformation Jacobian depends only on data and fixed x0, so it
# can be omitted when comparing parameter values.
fit_ln <- function(x, start=NULL) {
  y <- log(x-x0)
  yd <- log(d-x0)
  nn <- length(x)
  
  # Subtract the log-density sum and add n times the log survival at d.
  # Use log probabilities for stability when truncation is deep in the tail.
  nll <- function(par) {
    mu <- par[1]
    s <- exp(par[2])
    z <- (yd-mu)/s
    val <- -sum(dnorm(y, mu, s, log=TRUE)) +
      nn*pnorm(z, lower.tail=FALSE, log.p=TRUE)
    if (!is.finite(val)) 1e100 else val
  }
  
  # Initialize from the empirical log-shifted claims. These are starting
  # values only; truncation is handled by the objective, not this initializer.
  if (is.null(start)) {
    start <- c(mean(y), log(max(sd(y), 0.01)))
  }
  
  # Run both gradient-based BFGS and derivative-free Nelder-Mead from
  # the same start, then select the smaller objective. The code validates
  # finite parameters/objective, but does not require convergence code zero.
  candidates <- list(
    optim(start, nll, method='BFGS', control=list(maxit=1500)),
    optim(start, nll, method='Nelder-Mead', control=list(maxit=2000))
  )
  
  vals <- vapply(candidates, function(z) z$value, numeric(1))
  best <- candidates[[which.min(vals)]]
  
  if (!is.finite(best$value) ||
      any(!is.finite(best$par)) ||
      !is.finite(exp(best$par[2])) ||
      exp(best$par[2]) <= 0) {
    stop('Shifted LN fitting failed or produced invalid parameters.')
  }
  
  c(
    mu=unname(best$par[1]),
    sigma=unname(exp(best$par[2]))
  )
}


# ---- Pareto I model: fixed x0, conditional on X >= d ----
# PARETO LIKELIHOOD INTERVALS AND GROUND-UP RISK FUNCTIONALS
# This function uses the script-level n, d, x0, levels, and risk settings;
# its claims argument is expected to be the full validated sample.
run_pareto <- function(claims) {
  alpha_hat <- fit_pa(claims)
  if (!is.finite(alpha_hat) || alpha_hat<=0) stop('Invalid Pareto estimate')
  # Twice the loss of log likelihood relative to the MLE. For this
  # one-parameter model, compare it with chi-square(1) quantiles.
  deviance <- function(a) {
    ifelse(a>0,2*n*(log(alpha_hat/a)+a/alpha_hat-1),Inf)
  }
  # Find one LR crossing on each side of alpha_hat. The upper bracket
  # expands geometrically until it lies outside the confidence region.
  alpha_ci <- function(level) {
    f <- function(a) deviance(a)-qchisq(level,1)
    lo <- uniroot(f,c(alpha_hat*1e-8,alpha_hat),tol=1e-12)$root
    hi <- alpha_hat*2
    while (f(hi)<0) hi <- hi*2
    c(lower=lo,upper=uniroot(f,c(alpha_hat,hi),tol=1e-12)$root)
  }
  # Adaptive integral over the entire normal quantile domain avoids truncating
  # the Wang tail near alpha=1. Completing the square stabilizes quadrature.
  # Wang premium for the ground-up Pareto quantile. With the positive
  # lambda used here, alpha <= 1 gives a divergent integral.
  wt <- function(a) {
    if (a<=1) return(Inf)  # lambda is positive in these settings
    # Split quadrature near the asymptotic tail peak to help resolve
    # large contributions when alpha approaches its integrability boundary.
    c1 <- 1-1/a
    mode_z <- lambda/c1
    # Direct quantile-density identity: Q(t)=x0*t^(-1/a),
    # Wang weight g'(t)=exp(lambda*qnorm(t,lower.tail=FALSE)-lambda^2/2).
    integrand <- function(z) exp(dnorm(z,log=TRUE)+lambda*z-lambda^2/2-
                                pnorm(z,lower.tail=FALSE,log.p=TRUE)/a)
    ans <- integrate(integrand,-Inf,mode_z,rel.tol=1e-8,subdivisions=1000)$value+
           integrate(integrand,mode_z,Inf,rel.tol=1e-8,subdivisions=1000)$value
    x0*ans
  }
  # VaR is the beta upper-tail quantile; CTE averages beyond that quantile.
  # GS adjusts CTE using delta. PHT has a finite mean only if r*alpha > 1;
  # CTE and GS require alpha > 1. Preserve theoretical divergence as Inf.
  risk <- function(a) {
    VaR <- x0*beta^(-1/a)
    CTE <- if(a>1) VaR*a/(a-1) else Inf
    c(VaR=VaR,CTE=CTE,GS=if(a>1) CTE*(1+2*delta/(2*a-1)) else Inf,
      PHT=if(r*a>1) x0*r*a/(r*a-1) else Inf,WT=wt(a))
  }
  measures <- c('VaR','CTE','GS','PHT','WT')
  point <- risk(alpha_hat)
  # Build one alpha confidence interval per confidence level.
  shape <- do.call(rbind,lapply(levels,function(lev) {
    ci <- alpha_ci(lev)
    data.frame(Confidence_Level=lev,Alpha_Estimate=alpha_hat,
               Alpha_Lower=unname(ci[1]),Alpha_Upper=unname(ci[2]))
  }))
  out <- do.call(rbind,lapply(seq_along(levels),function(j) {
    # Each implemented risk decreases with alpha under these settings,
    # so swap the alpha endpoints when mapping them to risk endpoints.
    lower <- risk(shape$Alpha_Upper[j]); upper <- risk(shape$Alpha_Lower[j])
    stopifnot(all(lower<=upper),!anyNA(lower),!anyNA(upper))
    data.frame(Model='Pareto_I',Confidence_Level=levels[j],Risk_Measure=measures,
      Estimate_NOK=unname(point),CI_Lower_NOK=unname(lower),CI_Upper_NOK=unname(upper),
      Estimate_million_NOK=unname(point)/1e6,CI_Lower_million_NOK=unname(lower)/1e6,
      CI_Upper_million_NOK=unname(upper)/1e6,
      Upper_Bound_Infinite=is.infinite(upper),Estimate_Infinite=is.infinite(point))
  }))
  # Save parameter/risk tables, then draw the LR curve and its cutoffs.
  # Risk tables contain both NOK and million NOK plus infinity indicators.
  write.csv(shape,file.path(folder,'NFC_ParetoI_Profile_Alpha_CIs.csv'),row.names=FALSE)
  write.csv(out,file.path(folder,'NFC_ParetoI_Profile_Risk_CIs.csv'),row.names=FALSE)
  pdf(file.path(folder,'NFC_ParetoI_Profile_Likelihood.pdf'),width=8,height=5)
  xx <- seq(min(shape$Alpha_Lower)*.8,max(shape$Alpha_Upper)*1.15,length.out=700)
  plot(xx,deviance(xx),type='l',lwd=2,xlab=expression(alpha),
       ylab='Likelihood-ratio deviance',main='Pareto I: de-grouped fire claims')
  for(j in seq_along(levels)) abline(h=qchisq(levels[j],1),lty=j+1,col=j+1)
  abline(v=alpha_hat,lty=3)
  legend('topright',legend=c('MLE',paste0(100*levels,'% LR cutoff')),
         lty=c(3,seq_along(levels)+1),col=c(1,seq_along(levels)+1),bty='n')
  dev.off()
  # Write a readable record of the fixed settings, confidence intervals,
  # and interpretation limits alongside the machine-readable exports.
  writeLines(c('Pareto I profile-likelihood intervals; uniformly de-grouped claims.',
    sprintf('n=%d; d=%g; fixed x0=%g; alpha_hat=%.10g',n,d,x0,alpha_hat),
    'Ground-up risks. Inf means a theoretically infinite functional.',
    'Asymptotic chi-square(1) cutoffs; rounding uncertainty is excluded.',
    'Fixed x0 is not identifiable from these left-truncated Pareto observations.',
    capture.output(print(shape,row.names=FALSE)),
    capture.output(print(out,row.names=FALSE))),
    file.path(folder,'NFC_ParetoI_Profile_CI_Report.txt'))
  out
}
# Execute the complete Pareto analysis before the shifted-lognormal fit.
pareto_results <- run_pareto(x)

# ---- Shifted-lognormal model, using the same de-grouped x ----
y <- log(x-x0); yd <- log(d-x0)

# Conditional likelihood: normal log-density of log(X-x0), divided by
# Pr[X >= d]. The log-Jacobian -sum(log(X-x0)) is constant in parameters,
# so it cancels from every likelihood-ratio comparison.
loglik <- function(mu, sigma) {
  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) return(-Inf)
  z <- (yd-mu)/sigma
  ans <- sum(dnorm(y,mu,sigma,log=TRUE)) -
    n*pnorm(z, lower.tail=FALSE, log.p=TRUE)
  if (is.finite(ans)) ans else -Inf
}
# Retain the unrestricted MLE and its log likelihood as the common
# reference for every constrained risk profile and parameter projection.
ln_hat <- fit_ln(x)
mu_hat <- unname(ln_hat['mu'])
sigma_hat <- unname(ln_hat['sigma'])
llmax <- loglik(mu_hat,sigma_hat)
if (sigma_hat <= sigma_limits[1] || sigma_hat >= sigma_limits[2])
  stop('MLE outside sigma profiling domain; increase sigma_limits')

# Risk measures have the form x0 + exp(mu)*A(sigma). Numerical quadrature
# is used only for GS and PHT. Other three factors are analytic.
# Upper-tail probability t=exp(-v), with midpoint integration in v.
# NUMERICAL RISK FACTORS
# Use a midpoint grid after t = exp(-v). h controls resolution and vmax
# the finite integration range; GS uses only tail probabilities t <= beta.
# These numerical settings approximate integrals, not sampling uncertainty.
h <- 0.025; vmax <- 140
v <- seq(h/2,vmax-h/2,by=h)
t <- exp(-v); z <- qnorm(t,lower.tail=FALSE)
vc <- seq(-log(beta)+h/2,vmax-h/2,by=h)
tc <- exp(-vc); zc <- qnorm(tc,lower.tail=FALSE)

# Evaluate a log of summed exponentials after subtracting their maximum
# to reduce overflow/underflow for the finite values used in quadrature.
logsumexp <- function(a) { m <- max(a); m + log(sum(exp(a-m))) }
# Return log A(sigma), where each risk equals x0 + exp(mu)*A(sigma).
# VaR, CTE, and WT have analytic factors. PHT and GS use weighted
# midpoint sums, including the change-of-variable factor exp(-v).
log_A <- function(sigma, measure) {
  if (measure=='VaR') return(sigma*qnorm(beta,lower.tail=FALSE))
  if (measure=='CTE') return(sigma^2/2 +
    pnorm(sigma-qnorm(beta,lower.tail=FALSE),log.p=TRUE)-log(beta))
  if (measure=='WT') return(sigma^2/2 + lambda*sigma)
  if (measure=='PHT') return(log(r*h)+
    logsumexp(sigma*z-r*v))
  if (measure=='GS') {
    # GS weights modify the conditional-tail average. Positive weights are
    # required here because their logarithms enter the quadrature sum.
    w <- 1+2*delta-4*delta*tc/beta
    if (any(w<=0)) stop('GS weight must be positive')
    return(log(h/beta)+logsumexp(sigma*zc-vc+log(w)))
  }
  stop('Unknown measure')
}
measures <- c('VaR','CTE','GS','PHT','WT')
# Evaluate all five factors in a fixed order and restore monetary units.
risks <- function(mu,sigma) {
  a <- vapply(measures,function(m) log_A(sigma,m),numeric(1))
  setNames(x0+exp(mu+a),measures)
}
point <- risks(mu_hat,sigma_hat)

# At a specified log(R-x0)=q, the risk constraint uniquely determines
# mu=q-log_A(sigma). Profile out sigma numerically on a log scale.
# CONSTRAINED RISK LIKELIHOOD
# For a proposed q, eliminate mu exactly via the risk equation; only
# log(sigma) remains to be optimized within the declared finite domain.
# With details=TRUE, also return the best sigma for boundary diagnostics.
profile_at <- function(q, measure, details=FALSE) {
  # Negate the likelihood for minimization and penalize invalid evaluations
  # with a large finite objective value accepted by the optimizer.
  obj <- function(ls) {
    s <- exp(ls); mu <- q-log_A(s,measure)
    val <- loglik(mu,s)
    if (is.finite(val)) -val else .Machine$double.xmax/1e100
  }
  lo <- log(sigma_limits[1]); hi <- log(sigma_limits[2])
  # Multistart one-dimensional bounded minimization to avoid local maxima.
  grid <- seq(lo,hi,length.out=31)
  # Identify grid-local minima, refine their neighboring intervals, and
  # retain the best solution. This helps with multiple modes but is not
  # a mathematical guarantee of a global optimum.
  values <- vapply(grid,obj,numeric(1))
  candidates <- which(values <= c(Inf,head(values,-1)) &
                      values <= c(tail(values,-1),Inf))
  candidates <- unique(c(candidates,which.min(values)))
  opts <- lapply(candidates,function(i) {
    a <- grid[max(1,i-1)]; b <- grid[min(length(grid),i+1)]
    if (a==b) return(list(minimum=grid[i],objective=values[i]))
    optimize(obj,c(a,b),tol=1e-7)
  })
  # Include both domain endpoints explicitly; a boundary optimum can
  # signal that sigma_limits needs expansion.
  opts <- c(opts,list(list(minimum=lo,objective=obj(lo)),
                      list(minimum=hi,objective=obj(hi))))
  best <- opts[[which.min(vapply(opts,function(o) o$objective,numeric(1)))]]
  if (details) return(c(loglik=-best$objective,sigma=exp(best$minimum)))
  -best$objective
}

# Invert profile-likelihood ratio, with logarithmic risk coordinate q.
# This is profiling the RISK itself, not merely mapping marginal CIs for
# mu and sigma. The result accounts for their joint likelihood dependence.
# RISK CONFIDENCE INTERVAL INVERSION
# Solve 2*(unrestricted logLik - constrained logLik) = chi-square(1).
# Clamp small negative deviances to zero for numerical roundoff.
profile_ci <- function(measure, level) {
  qhat <- log(point[[measure]]-x0)
  cutoff <- qchisq(level,df=1)
  dev <- function(q) max(0,2*(llmax-profile_at(q,measure)))
  # Search outward separately below and above the fitted log risk.
  # Once the cutoff is crossed, refine that bracket with uniroot.
  root_side <- function(side) {
    step <- 0.20; limit <- 80
    for (j in 1:80) {
      qtest <- qhat + side*step
      # The finite q limit and iteration cap are numerical safeguards.
      # Returning Inf or x0 when the search stops is a fallback endpoint;
      # it does not alone establish a mathematically unbounded LR region.
      if (!is.finite(qtest) || abs(qtest)>limit) return(if (side>0) Inf else x0)
      dtest <- dev(qtest)
      if (is.finite(dtest) && dtest >= cutoff) {
        bounds <- sort(c(qhat,qtest))
        root <- uniroot(function(q) dev(q)-cutoff,bounds,tol=1e-7)$root
        return(x0+exp(root))
      }
      step <- step*1.5
    }
    if (side>0) Inf else x0
  }
  c(lower=unname(root_side(-1)),upper=unname(root_side(1)))
}

cat('Fitted mu:',mu_hat,' sigma:',sigma_hat,' n:',n,'\n')
# Compute one row for each measure/level pair; preserve asymmetric
# endpoints and express estimates in NOK and million NOK.
rows <- list(); idx <- 0L
for (m in measures) {
  for (lev in levels) {
    cat('Profiling ',m,' at ',100*lev,'% ...\n',sep='')
    ci <- profile_ci(m,lev)
    idx <- idx+1L
    rows[[idx]] <- data.frame(Model='Shifted_Lognormal',Risk_Measure=m,
      Confidence_Level=lev,Estimate_NOK=unname(point[m]),
      CI_Lower_NOK=unname(ci['lower']),CI_Upper_NOK=unname(ci['upper']),
      Estimate_million_NOK=unname(point[m])/1e6,
      CI_Lower_million_NOK=unname(ci['lower'])/1e6,
      CI_Upper_million_NOK=unname(ci['upper'])/1e6)
  }
}
out <- do.call(rbind,rows)
# Group output by confidence level, then by risk measure.
out <- out[order(out$Confidence_Level, match(out$Risk_Measure, measures)), ]
rownames(out) <- NULL
write.csv(out,file.path(folder,'NFC_Shifted_LN_Profile_Risk_CIs.csv'),row.names=FALSE,na='NA')
# Export the fitted sample size, fixed support/threshold, and LN MLE.
# The saved logLik omits the parameter-constant transformation Jacobian.
parameters <- data.frame(n=n,x0=x0,d=d,mu_hat=mu_hat,sigma_hat=sigma_hat,logLik=llmax)
write.csv(parameters,file.path(folder,'NFC_Shifted_LN_Fitted_Parameters.csv'),row.names=FALSE)
# Same result schema as the Pareto I shape-parameter CI file, but with
# shifted-lognormal mu and sigma estimates. These are joint 2-parameter
# likelihood-region projections, NOT separate 1-df profile CIs.
# PARAMETER PROJECTIONS OF A JOINT LIKELIHOOD REGION
# Use chi-square(2), since the underlying region concerns (mu, sigma).
# Project that region by profiling the other coordinate at each trial value;
# these bounds differ from individual chi-square(1) profile intervals.
parameter_rows <- lapply(levels, function(lev) {
  cutoff <- qchisq(lev, df=2)
  # Along each coordinate, optimize over the other coordinate.
  # For fixed mu, optimize over log(sigma) within sigma_limits.
  mu_prof <- function(mu) {
    f <- function(ls) -loglik(mu, exp(ls))
    opt <- optimize(f, log(sigma_limits), tol=1e-8)
    2*(llmax+opt$objective)
  }
  # For fixed log(sigma), optimize mu over a broad data-based interval.
  # These scalar optimize calls assume the objective is suitable for bounded
  # one-dimensional optimization; they do not use the multistart grid above.
  sigma_prof <- function(ls) {
    f <- function(mu) -loglik(mu,exp(ls))
    opt <- optimize(f, c(min(y)-20,max(y)+20), tol=1e-8)
    2*(llmax+opt$objective)
  }
  # Find the joint-region cutoff crossing within the specified coordinate
  # limits. Return NA if the limit is still inside the region; this marks
  # an unresolved endpoint within the chosen search domain.
  find_bound <- function(fn, center, side, lower_limit, upper_limit) {
    end <- if (side < 0) lower_limit else upper_limit
    if (fn(end) < cutoff) return(NA_real_)
    uniroot(function(v) fn(v)-cutoff,sort(c(center,end)),tol=1e-7)$root
  }
  mu_bounds <- c(find_bound(mu_prof,mu_hat,-1,min(y)-20,max(y)+20),
                 find_bound(mu_prof,mu_hat,1,min(y)-20,max(y)+20))
  # Locate sigma endpoints on its log scale, then exponentiate for output.
  ls_bounds <- c(find_bound(sigma_prof,log(sigma_hat),-1,log(sigma_limits[1]),log(sigma_limits[2])),
                 find_bound(sigma_prof,log(sigma_hat),1,log(sigma_limits[1]),log(sigma_limits[2])))
  data.frame(Confidence_Level=lev,Mu_Estimate=mu_hat,Mu_Lower=mu_bounds[1],
    Mu_Upper=mu_bounds[2],Sigma_Estimate=sigma_hat,
    Sigma_Lower=exp(ls_bounds[1]),Sigma_Upper=exp(ls_bounds[2]))
})
parameter_ci <- do.call(rbind, parameter_rows)
write.csv(parameter_ci,file.path(folder,'NFC_Shifted_LN_Profile_Parameter_CIs.csv'),row.names=FALSE,na='NA')


# One page per confidence level, grouped consistently with the CSV.
# INTERVAL PLOTS
# Each page has five risk panels plus one title panel. Finite values set
# the vertical range; error bars are drawn only when both endpoints are
# finite. A panel annotation flags nonfinite endpoints.
pdf(file.path(folder,'NFC_Shifted_LN_Profile_Risk_CIs.pdf'),width=9,height=8)
for (lev in levels) {
  par(mfrow=c(3,2),mar=c(4,4,3,1))
  for (m in measures) {
    rr <- subset(out,Risk_Measure==m & abs(Confidence_Level-lev)<1e-9)
    vals <- c(rr$CI_Lower_million_NOK,rr$Estimate_million_NOK,rr$CI_Upper_million_NOK)
    vals <- vals[is.finite(vals)]
    ylim <- if (length(vals)) range(c(0,vals)) else c(0,1)
    # Expand a degenerate range so the plotting device has a valid axis.
    if (diff(ylim)==0) ylim <- ylim+c(-0.1,0.1)
    plot(1,rr$Estimate_million_NOK,xlim=c(.7,1.3),ylim=ylim,
         xaxt='n',xlab='',ylab='Million NOK',
         main=sprintf('%s (%.0f%% profile CI)',m,100*lev),pch=19)
    axis(1,at=1,labels='Shifted LN')
    if (all(is.finite(c(rr$CI_Lower_million_NOK,rr$CI_Upper_million_NOK))))
      arrows(1,rr$CI_Lower_million_NOK,1,rr$CI_Upper_million_NOK,
             angle=90,code=3,length=.07,lwd=2)
    else mtext('Nonfinite confidence endpoint',side=3,cex=.7)
  }
  plot.new();text(.5,.5,sprintf('NFC shifted lognormal\n%.0f%% confidence intervals',100*lev))
}
dev.off()

# Direct likelihood-ratio profile plots, one panel per risk measure.
# Risk coordinate is log(R-x0), allowing wide and asymmetric intervals.
# PROFILE CURVE PLOTS
# Use the widest requested confidence interval to choose the plotting
# range, with finite fallbacks when endpoints are unavailable. Evaluate
# 65 profile points and display the horizontal risk axis in million NOK.
pdf(file.path(folder,'NFC_Shifted_LN_Profile_Likelihood.pdf'),width=10,height=8)
par(mfrow=c(3,2),mar=c(4,4,3,1))
for (m in measures) {
  qhat <- log(unname(point[m])-x0)
  rr <- subset(out,Risk_Measure==m & Confidence_Level==max(levels))
  qlo <- if (is.finite(rr$CI_Lower_NOK) && rr$CI_Lower_NOK>x0)
    log(rr$CI_Lower_NOK-x0) else qhat-1.5
  qhi <- if (is.finite(rr$CI_Upper_NOK)) log(rr$CI_Upper_NOK-x0) else qhat+1.5
  qs <- seq(qlo,qhi,length.out=65)
  ds <- vapply(qs,function(q) max(0,2*(llmax-profile_at(q,m))),numeric(1))
  plot((x0+exp(qs))/1e6,ds,type='l',lwd=2,
       xlab='Risk (million NOK)',ylab='Likelihood-ratio deviance',
       main=paste(m,'risk profile'))
  abline(v=unname(point[m])/1e6,lty=3)
  for (j in seq_along(levels)) abline(h=qchisq(levels[j],1),lty=j+1,col=j+1)
  legend('topright',legend=paste0(100*levels,'% cutoff'),
         lty=seq_along(levels)+1,col=seq_along(levels)+1,bty='n',cex=.75)
}
plot.new();text(.5,.5,'NFC shifted-lognormal\nRisk profile likelihoods')
dev.off()

# DETAILED SHIFTED-LOGNORMAL REPORT
# Capture printed tables together with settings and interpretation notes.
# The report distinguishes one-risk LR intervals from joint parameter
# projections and records the excluded sources of uncertainty.
report <- c('Norwegian Fire Claims: shifted-lognormal profile-likelihood risk intervals',
  sprintf('n=%d, fixed x0=%g NOK, truncation d=%g NOK',n,x0,d),
  sprintf('MLE: mu=%.9f, sigma=%.9f',mu_hat,sigma_hat),
  'Risk profile: maximize left-truncated log-likelihood over sigma, with mu fixed by R.',
  'Risk CIs use chi-square(1) cutoffs. Parameter bounds use chi-square(2) joint-region projections.',
  'All confidence intervals are asymptotic, not bootstrap calibrated.',
  'Ties were uniformly de-grouped before fitting; rounding uncertainty is excluded.',
  'x0 is fixed, so its uncertainty is not included.',
  'All risks are ground-up, not conditional-on-exceedance risks.',
  'Numerical sigma profiling limits should be checked for boundary solutions.',
  '', 'PARAMETER CONFIDENCE BOUNDS:',
  capture.output(print(parameter_ci,row.names=FALSE)), '')
for (lev in levels) {
  report <- c(report,sprintf('===== %.0f%% CONFIDENCE LEVEL =====',100*lev),
    capture.output(print(subset(out,Confidence_Level==lev)[,
      c('Risk_Measure','Estimate_million_NOK','CI_Lower_million_NOK',
        'CI_Upper_million_NOK')],row.names=FALSE)), '')
}
writeLines(report,file.path(folder,'NFC_Shifted_LN_Profile_CI_Report.txt'))
cat('Results saved in:',folder,'\n')

# Comparable rows for both models in one CSV (identical units and data).
# COMBINED MODEL EXPORT
# Select the shared columns so both model results have identical units
# and schema before binding rows; model-specific flags stay in their files.
common_columns <- c('Model','Confidence_Level','Risk_Measure','Estimate_NOK',
 'CI_Lower_NOK','CI_Upper_NOK','Estimate_million_NOK',
 'CI_Lower_million_NOK','CI_Upper_million_NOK')
combined <- rbind(pareto_results[,common_columns],out[,common_columns])
combined <- combined[order(combined$Confidence_Level,
 match(combined$Risk_Measure,measures),combined$Model),]
write.csv(combined,file.path(folder,'NFC_Combined_Profile_Risk_CIs.csv'),row.names=FALSE)
cat('Both models fitted to identical de-grouped claims. All results in:',folder,'\n')


# Bootstrap-style summary tables, retaining PROFILE LIKELIHOOD intervals.
# Bootstrap replicate counts do not apply to these likelihood intervals.
# PRESENTATION TABLES AND FINAL SUMMARY
# Reformat model labels and confidence percentages for readable tables.
# This changes presentation only; estimates and interval endpoints are copied.
make_profile_table <- function(results, model) {
  data.frame(
    Model=model,
    CI_Level=paste0(round(results$Confidence_Level*100), '%'),
    Risk_Measure=results$Risk_Measure,
    Estimate_million_NOK=results$Estimate_million_NOK,
    CI_Lower_million_NOK=results$CI_Lower_million_NOK,
    CI_Upper_million_NOK=results$CI_Upper_million_NOK,
    stringsAsFactors=FALSE
  )
}
pa_tab <- make_profile_table(pareto_results, 'Pareto I')
ln_tab <- make_profile_table(out, 'Shifted Lognormal')
summary_out <- file.path(folder, 'NFC_Profile_Likelihood_Risk_CIs.csv')
# Write two labeled CSV table blocks through one text connection.
# Because title lines and a blank separator are included, this file is a
# presentation export rather than one rectangular CSV data table.
con <- file(summary_out, 'wt')
writeLines('Table 1. Pareto I - profile likelihood CIs', con)
write.table(pa_tab, con, sep=',', row.names=FALSE,
            col.names=TRUE, quote=TRUE, na='')
writeLines(c('', 'Table 2. Shifted Lognormal - profile likelihood CIs'), con)
write.table(ln_tab, con, sep=',', row.names=FALSE,
            col.names=TRUE, quote=TRUE, na='')
close(con)

# Assemble one parameter row per model. NA marks a parameter that is
# not part of that model, rather than an unsuccessful estimate.
alpha_hat <- fit_pa(x)
params <- data.frame(
  Model=c('Pareto I', 'Shifted Lognormal'), n=n, x0=x0, d=d,
  Alpha=c(alpha_hat, NA), Mu=c(NA, mu_hat), Sigma=c(NA, sigma_hat)
)
params_out <- file.path(folder, 'NFC_Profile_Likelihood_Fitted_Parameters.csv')
write.csv(params, params_out, row.names=FALSE)
# Same degrouped-data layout as the bootstrap script.
degrouped_out <- file.path(folder, 'NFC_Degrouped_Data.csv')
write.csv(data.frame(Observation=seq_along(x), Claim_NOK=sort(x)),
          degrouped_out, row.names=FALSE)

# Capture console-style summary text, settings, and both model tables
# into one report. Confidence levels determine the report grouping.
summary_report <- file.path(folder, 'NFC_Profile_Likelihood_CI_Report.txt')
report_lines <- capture.output({
  cat('NFC profile likelihood confidence intervals\n')
  cat('Ground-up x0=', x0, '; left truncation d=', d, '; n=', n, '\n', sep='')
  cat('Risk settings: beta=', beta, ', delta=', delta,
      ', r=', r, ', lambda=', lambda, '\n', sep='')
  cat('Pareto alpha=', alpha_hat, '; LN mu=', mu_hat,
      '; LN sigma=', sigma_hat, '\n')
  cat('Fitting matches bootstrap: closed-form Pareto MLE; best of BFGS and Nelder-Mead for LN.\n')
  cat('Interior-spacing de-grouping: rounding unit=', rounding_unit,
      ' NOK; observations in tied groups=',
      sum(degrouped$audit$Count[degrouped$audit$Adjusted]),
      '; tied groups=', sum(degrouped$audit$Adjusted), '\n', sep='')
  cat('Rule: L + j*(U-L)/(m+1), j=1,...,m, for each group of m ties.\n')
  cat('Only ties were processed; lower rounding endpoints were clipped at d.\n')
  cat('Intervals condition on the deterministic interior-spaced dataset.\n')
  cat('Risk measures are ground-up; fitting is conditional on X>=d.\n')
  cat('x0 is fixed; its uncertainty and rounding uncertainty are excluded.\n')
  cat('Risk intervals use asymptotic chi-square(1) likelihood-ratio cutoffs.\n')
  cat('LN parameter bounds use chi-square(2) joint-region projections.\n')
  cat('Inf denotes theoretical divergence; NA denotes missing calculations.\n')
  cat('Risk evaluation retains the profile script quadrature and analytic formulas;\n')
  cat('small numerical differences from bootstrap risk estimates may remain.\n')
  cat('Numerical sigma profiling limits:', sigma_limits, '\n')
  cat('\nFITTED PARAMETERS\n')
  print(params, row.names=FALSE)
  for (level in paste0(round(levels*100), '%')) {
    cat('\n============ ', level, ' ============\n', sep='')
    cat('\nPareto I\n')
    print(pa_tab[pa_tab$CI_Level == level, ], row.names=FALSE)
    cat('\nShifted Lognormal\n')
    print(ln_tab[ln_tab$CI_Level == level, ], row.names=FALSE)
  }
})
# Persist the final report and announce the principal output paths.
writeLines(report_lines, summary_report)
cat('Saved:\n', degrouped_out, '\n', summary_out, '\n',
    params_out, '\n', summary_report, '\n')
