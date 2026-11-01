# SHIFTED LOGNORMAL SIMULATION -- BASE R
# Units: thousands. Known shift x0.
# ML targets ground-up X. EMP targets Y=min(X,u) conditional on X>d.
# These are DIFFERENT estimands.
# Finite-n ML summaries retain finite risks from accepted finite AND limiting fits.
# Infinite, ambiguous limiting and numerical failures are counted separately.
# Expected-information raw-scale delta CIs are computed for interior fits only.
# B below is the number of Monte Carlo replications, not bootstrap samples.
# nInf: Fisher-information/delta-method limits, evaluated numerically.

x0 <- 1
sigma <- 2
mu_original <- 5  # log(1000) - sigma^2/2  # 5.933607313867482
unit_scale <- 1000
mu <- mu_original - log(unit_scale)  # internal meanlog in thousands
# X = x0 + exp(mu_original + sigma*Z)/unit_scale
# E[X]=2; P(X>14 | X>4)=(4/14)^2.
d <- 4
u <- 14
beta <- 0.10
Delta <- 0.25
r <- 0.75
lambda <- 0.50
sample_sizes <- c(50L,100L,500L)
B <- 100000L  # Try 100 first. Full run is computationally intensive.
random_seed <- 200L
# NULL saves beside this script; alternatively supply an existing folder.
output_dir <- NULL

if (is.null(output_dir)) {
  script_file <- NULL
  if (requireNamespace('rstudioapi',quietly=TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(),error=function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) script_file <- ctx$path
  }
  if (is.null(script_file)) {
    paths <- Filter(function(p) is.character(p)&&length(p)==1L&&nzchar(p),
                    lapply(sys.frames(),function(f) f$ofile))
    if (length(paths)) script_file <- tail(paths,1L)[[1L]]
  }
  if (is.null(script_file)) {
    a <- grep('^--file=',commandArgs(FALSE),value=TRUE)
    if (length(a)) script_file <- sub('^--file=','',a[1L])
  }
  if (is.null(script_file)||!file.exists(script_file))
    stop('Save and Source this script, or set output_dir explicitly.')
  output_dir <- dirname(normalizePath(script_file,mustWork=TRUE))
}
# Noninteger n*beta is supported by exact distortion increments in emp_risk.
stopifnot(dir.exists(output_dir),x0>=0,sigma>0,d>x0,u>d,
          beta>0,beta<1,Delta>=0,Delta<=0.5,r>0,r<=1,lambda>=0,
          B>=2,B==as.integer(B),all(sample_sizes>=2),
          all(sample_sizes==as.integer(sample_sizes)))
output_file <- file.path(output_dir,'Shifted_LN_Combined_Results_v4.csv')
measures <- c('VaR','CTE','GS','PHT','WT')
ld <- log(d-x0)
lu <- log(u-x0)
cap_log <- lu-ld
logS <- function(t,m=mu,s=sigma)
  pnorm((log(t-x0)-m)/s,lower.tail=FALSE,log.p=TRUE)
logSd <- logS(d)
atom <- exp(logS(u)-logSd)
mills <- function(z)
  exp(dnorm(z,log=TRUE)-pnorm(z,lower.tail=FALSE,log.p=TRUE))

# Stable ground-up risks: adaptive quadrature, with exponential factors
# removed before integration. Numerical overflow remains explicitly Inf.
risk <- function(par) {
  m <- par[1L]; s <- exp(par[2L])
  if (!is.finite(m)||!is.finite(s)||s<=0) return(rep(NA_real_,5L))
  zb <- qnorm(beta,lower.tail=FALSE)
  logEX <- m+s*s/2
  VaR <- x0+exp(m+s*zb)
  CTE <- x0+exp(logEX+pnorm(zb-s,lower.tail=FALSE,log.p=TRUE)-log(beta))
  # GS quantile integral; v=z-s centers its normal density.
  # If the lower tail starts beyond 12, rescale around its lower endpoint.
  L <- zb-s
  if (L<=12) {
    gs_factor <- integrate(function(v) dnorm(v)*
      (1+2*Delta-4*Delta*pnorm(s+v,lower.tail=FALSE)/beta)/beta,
      max(-12,L),12,rel.tol=1e-9,abs.tol=1e-12,
      subdivisions=500L)$value
    GS <- x0+exp(logEX+log(gs_factor))
  } else {
    # Generic fallback for unusually small tail probabilities.
    base <- dnorm(L,log=TRUE)
    fac <- integrate(function(t) exp(dnorm(L+t,log=TRUE)-base)*
      (1+2*Delta-4*Delta*pnorm(zb+t,lower.tail=FALSE)/beta)/beta,
      0,12,rel.tol=1e-9,abs.tol=1e-12,subdivisions=500L)$value
    GS <- x0+exp(logEX+base+log(fac))
  }
  # PHT: locate the mode, scale out its height, then integrate nearby.
  # h'(z)=s+(1-r)*Mills(z)-z, so the mode is bracketed below.
  right <- s/min(r,1)+12/sqrt(min(r,1))
  mode <- uniroot(function(z) s+(1-r)*mills(z)-z,
                  c(-12,right),tol=1e-10)$root
  h <- function(z) s*z+(r-1)*pnorm(z,lower.tail=FALSE,log.p=TRUE)+
    dnorm(z,log=TRUE)
  peak <- h(mode)
  width <- 12/sqrt(min(r,1))
  fac <- integrate(function(v) exp(h(mode+v)-peak),-width,width,
                   rel.tol=1e-9,abs.tol=1e-12,subdivisions=500L)$value
  PHT <- x0+exp(m+log(r)+peak+log(fac))
  WT <- x0+exp(logEX+lambda*s)
  setNames(c(VaR,CTE,GS,PHT,WT),measures)
}
truepar <- c(mu,log(sigma))
truth <- risk(truepar)

# Fisher information per conditional, capped observation.
zd <- (ld-mu)/sigma
zu <- (lu-mu)/sigma
trunc_score <- c(mills(zd)/sigma,zd*mills(zd))
score_cap <- c(mills(zu)/sigma,zu*mills(zu))-trunc_score
Info <- matrix(0,2,2)
for (i in 1:2) for (j in 1:2) {
  Info[i,j] <- integrate(function(z) {
    si <- if(i==1) z/sigma-trunc_score[1] else z*z-1-trunc_score[2]
    sj <- if(j==1) z/sigma-trunc_score[1] else z*z-1-trunc_score[2]
    si*sj*exp(dnorm(z,log=TRUE)-logSd)
  },zd,zu,rel.tol=1e-10)$value+atom*score_cap[i]*score_cap[j]
}
if (min(eigen(Info,symmetric=TRUE,only.values=TRUE)$values)<=0)
  stop('Fisher information is not positive definite.')
Vpar <- solve(Info)
grad_risk <- function(h) {
  G <- matrix(0,5,2)
  for (j in 1:2) {
    pp <- pm <- truepar
    pp[j] <- pp[j]+h
    pm[j] <- pm[j]-h
    G[,j] <- (risk(pp)-risk(pm))/(2*h)
  }
  G
}
G <- grad_risk(1e-4)
Gcheck <- grad_risk(5e-5)
if (max(abs(G-Gcheck)/pmax(1,abs(Gcheck)))>1e-4)
  stop('Risk-gradient convergence check failed.')
asym_ml <- data.frame(measure=measures,mean=truth,
  scaled_SD=sqrt(diag(G%*%Vpar%*%t(G))))

# EMP deterministic limits: integrate the influence function, including
# the probability mass at u. Split the grid at the CTE/GS kink.
qcond <- x0+exp(mu+sigma*qnorm(log(beta)+logSd,
                              lower.tail=FALSE,log.p=TRUE))
emp_limits <- function(N) {
  x <- sort(unique(c(seq(d,u,length.out=N+1L),
                     if(qcond<u) qcond else numeric(0))))
  dx <- diff(x)
  trap <- function(v) sum(dx*(head(v,-1L)+tail(v,-1L))/2)
  calc <- function(g,gp) {
    A <- numeric(length(x))
    # Midpoints avoid assigning a finite grid interval to a kink endpoint.
    mid <- (head(x,-1L)+tail(x,-1L))/2
    sm <- exp(logS(mid)-logSd)
    A[-1L] <- cumsum(dx*gp(sm))
    f <- dlnorm(x-x0,mu,sigma)/exp(logSd)
    S <- exp(logS(x)-logSd)
    EA <- trap(A*f)+atom*tail(A,1L)
    variance <- trap(A*A*f)+atom*tail(A,1L)^2-EA^2
    c(mean=d+trap(g(S)),scaled_SD=sqrt(max(0,variance)))
  }
  cte <- calc(function(s) pmin(s/beta,1),
              function(s) ifelse(s<beta,1/beta,0))
  gs <- calc(function(s) {t<-pmin(s/beta,1);t*(1+2*Delta*(1-t))},
             function(s) ifelse(s<beta,(1+2*Delta-4*Delta*s/beta)/beta,0))
  pht <- calc(function(s) s^r,function(s) r*s^(r-1))
  wt <- calc(function(s) pnorm(qnorm(s)+lambda),function(s)
    if(lambda==0) rep(1,length(s)) else exp(-lambda*qnorm(s)-lambda^2/2))
  if (abs(atom-beta)<1e-8) stop('VaR at censoring boundary: nonregular quantile limit.')
  vsd <- if(qcond<u) sqrt(beta*(1-beta))/
    (dlnorm(qcond-x0,mu,sigma)/exp(logSd)) else 0
  data.frame(measure=measures,mean=c(min(qcond,u),cte[1],gs[1],pht[1],wt[1]),
             scaled_SD=c(vsd,cte[2],gs[2],pht[2],wt[2]))
}
asym_emp <- emp_limits(100000L)
echeck <- emp_limits(50000L)
if (max(abs(asym_emp$scaled_SD-echeck$scaled_SD)/
        pmax(1,asym_emp$scaled_SD))>1e-4) stop('EMP grid convergence check failed.')

# Natural coordinates use b=1/sigma^2, NOT the old b=1/(2*sigma^2).
# Conditional U=log((X-x0)/(d-x0)) has density proportional to
# exp(-a*U-b*U^2/2). Its proper closure is b=0,a>0.
H <- cap_log; yd <- ld; delta <- Delta; zb <- qnorm(beta,lower.tail=FALSE)
confidence <- 0.90
zcrit <- qnorm((1+confidence)/2)
score_tolerance <- 1e-6
likelihood_tolerance <- 1e-6
integration_tolerance <- 1e-8
max_t_bound <- 80
D <- new.env(parent=emptyenv()); D$counts <- numeric(0)
record <- function(scope,reason,count=1) {
  key <- paste(scope,reason,sep='|')
  if(!key %in% names(D$counts)) D$counts[key] <- 0
  D$counts[key] <- D$counts[key]+count
}
normalizer <- function(a,b) {
  if(b==0) {
    if(a<=0) stop('Improper limiting density')
    return(list(logZ=-log(a),M1=1/a,M2=2/a^2))
  }
  if(!is.finite(a) || !is.finite(b) || b<0) stop('Invalid natural parameters')
  z <- a/sqrt(b)
  if(z>20) {
    w <- b/a^2
    e <- 1+w*(-2+w*(10+w*(-74+w*(706-8162*w))))
    return(list(logZ=-log(a)-log1p(w*e),M1=e/a,
      M2=(2+w*(-10+w*(74+w*(-706+8162*w))))/a^2))
  }
  lm <- dnorm(z,log=TRUE)-pnorm(z,lower.tail=FALSE,log.p=TRUE)
  excess <- exp(lm)-z
  list(logZ=-log(b)/2-lm,M1=excess/sqrt(b),M2=(1-z*excess)/b)
}
# T1 and T2 include censored exposure H and H^2, respectively.
likelihood <- function(a,b,s) {
  z0 <- normalizer(a,b); zh <- normalizer(a+b*H,b)
  ll <- -a*s$T1-b*s$T2/2-s$m*z0$logZ+s$k*(zh$logZ-z0$logZ)
  ga <- -s$T1+s$n*z0$M1-s$k*zh$M1
  gb <- -s$T2/2+s$n*z0$M2/2-s$k*(zh$M2+2*H*zh$M1)/2
  list(ll=ll,score=c(ga,gb))
}
a_to_t <- function(a,b) {
  if(b==0) return(log(a))
  root <- sqrt(a*a+4*b)
  log(if(a>=0) (a+root)/2 else 2*b/(root-a))
}
fit_model <- function(s,warm=NULL,scope='outer') {
  if(s$m==0 || !is.finite(s$T1) || s$T1<=0) {
    record(scope,'unresolved_all_censored_or_invalid'); return(NULL)
  }
  a0 <- s$m/s$T1
  boundary <- likelihood(a0,0,s)
  starts <- list(c(log(a0),0),c(a_to_t(a0,.5),.5),c(a_to_t(a0,2),2))
  if(!is.null(warm) && warm$type=='finite')
    starts <- c(list(c(a_to_t(warm$a,warm$b),warm$b)),starts)
  candidates <- list()
  bound <- 20
  repeat {
    last <- NULL; cached <- NULL
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
    if(bound>=max_t_bound) break
    bound <- min(max_t_bound,bound*2)
    if(!is.null(best)) starts <- c(list(best$par),starts)
  }
  record(scope,'unresolved_optimizer_or_score'); NULL
}

natural_mills <- function(z) {
  if(z>20) {
    t <- 1/z; e <- t*(1+t^2*(-2+t^2*(10+t^2*(-74+t^2*(706-8162*t^2)))))
    return(list(log=log(z)+log1p(e/z),excess=e,value=z+e))
  }
  lm <- dnorm(z,log=TRUE)-pnorm(z,lower.tail=FALSE,log.p=TRUE)
  list(log=lm,excess=exp(lm)-z,value=exp(lm))
}
# log A(sigma)=coefficient*sigma^2+remainder.
# Separate the leading term to avoid cancellation in ground-up extrapolation.
factor_parts <- function(sigma,m,derivative=FALSE) {
  if(m=='VaR') return(list(coefficient=0,remainder=sigma*zb,prime=zb))
  if(m=='CTE') return(list(coefficient=.5,
    remainder=pnorm(sigma-zb,log.p=TRUE)-log(beta),prime=natural_mills(zb-sigma)$value))
  if(m=='WT') return(list(coefficient=.5,remainder=lambda*sigma,prime=lambda))
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
  if(r==1) return(list(coefficient=.5,remainder=0,prime=0))
  center <- sigma/r; scale <- 1/sqrt(r); shift <- (1-r)*natural_mills(center)$log
  f <- function(t) vapply(t,function(tt)
    exp(-tt^2/2+(1-r)*natural_mills(center+scale*tt)$log-shift),numeric(1))
  I <- integrate(f,-Inf,Inf,rel.tol=integration_tolerance)$value
  prime <- 0
  if(derivative) {
    J <- integrate(function(t) f(t)*vapply(center+scale*t,
      function(z) natural_mills(z)$excess,numeric(1)),-Inf,Inf,rel.tol=integration_tolerance)$value
    prime <- (1-r)*J/(r*I)
  }
  list(coefficient=1/(2*r),remainder=log(r)/2-r*log(2*pi)/2+shift+log(I),prime=prime)
}
critical_slope <- function(m) if(m=='VaR') 0 else if(m=='PHT') 1/(2*r) else .5
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
  q <- yd+(f$coefficient-fit$a)/fit$b+f$remainder
  if(!is.finite(q) || q>log(.Machine$double.xmax))
    return(list(value=NA_real_,gradient=c(NA_real_,NA_real_),logexcess=q))
  excess <- exp(q)
  g <- if(gradient && excess>0) excess*c(-1/fit$b,
    (fit$a-f$coefficient)/fit$b^2-f$prime/(2*fit$b^(3/2))) else c(NA_real_,NA_real_)
  list(value=x0+excess,gradient=g,logexcess=q)
}
all_risks <- function(fit) vapply(measures,function(m)
  tryCatch(risk_value(fit,m)$value,error=function(e) NA_real_),numeric(1))

# Expected information per conditional/censored record in (a,b) coordinates.
# Integrate exact-record score products, then add the censored point mass.
expected_information <- function(fit) {
  a <- fit$a; b <- fit$b
  z0 <- normalizer(a,b); zh <- normalizer(a+b*H,b)
  logpc <- -a*H-b*H^2/2+zh$logZ-z0$logZ
  cens <- c(z0$M1-H-zh$M1,(z0$M2-H^2-2*H*zh$M1-zh$M2)/2)
  sc <- function(t) cbind(z0$M1-t,(z0$M2-t^2)/2)
  I <- matrix(0,2,2)
  for(i in 1:2) for(j in i:2) {
    value <- integrate(function(t) {
      s <- sc(t); exp(-a*t-b*t^2/2-z0$logZ)*s[,i]*s[,j]
    },0,H,rel.tol=integration_tolerance)$value+exp(logpc)*cens[i]*cens[j]
    I[i,j] <- I[j,i] <- value
  }
  I
}
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
fit_ml <- function(y) {
  w <- log((y-x0)/(d-x0)); exact <- y<u
  s <- list(n=length(y),m=sum(exact),k=sum(!exact),T1=sum(w),T2=sum(w*w))
  fit <- fit_model(s)
  if(is.null(fit)) return(list(ok=FALSE,reason='unresolved_fit',type='unresolved',
                             par=c(NA_real_,NA_real_)))
  fit$ok <- TRUE
  fit$reason <- paste0('accepted_',fit$type)
  fit$par <- if(fit$type=='finite') c(ld-fit$a/fit$b,-log(fit$b)/2) else
    c(NA_real_,NA_real_)
  fit
}
# b=0 is an extended conditional MLE, not a finite ground-up LN model.
# Report only path-independent limiting risk values. At the critical slope,
# risk_value returns NA because the extrapolated ground-up limit is ambiguous.
# Standard interior-MLE delta intervals are unavailable for limiting fits.

# Cross-check the new natural-coordinate expected-information calculation
# against the original (mu, log(sigma)) asymptotic delta calculation at truth.
check_fit <- list(type='finite',a=(ld-mu)/sigma^2,b=1/sigma^2)
check_values <- all_risks(check_fit)
stopifnot(max(abs(check_values-truth)/pmax(1,abs(truth)))<1e-6)
check_I <- expected_information(check_fit)
check_R <- chol(check_I)
check_sd <- vapply(measures,function(m) {
  g <- risk_value(check_fit,m,TRUE)$gradient
  sqrt(sum(forwardsolve(t(check_R),g)^2))
},numeric(1))
stopifnot(max(abs(check_sd-asym_ml$scaled_SD)/
              pmax(1,abs(asym_ml$scaled_SD)))<1e-4)
# At b=0 an ordinary two-parameter Wald CI must be unavailable.
stopifnot(delta_ci(list(type='limit',a=1,b=0),50,'VaR')$reason==
            'delta_limit_unavailable')
cat('Expected-information coordinate and limiting-delta checks passed.\n')

# Empirical estimators: table formulas and their relation to eq. (empDM).
# Here b contains the ordered observed losses x*_(i), sampled from F_*:
# Y = min(X,u) conditional on X>d (not the ground-up distribution F).
# For an absolutely continuous distortion, eq. (empDM) assigns weight
# g(1-(i-1)/n) - g(1-i/n) to b[i]. PHT and WT below use these
# increments exactly, for any n.
# CTE and GS retain the original table formulas for integer n*beta.
# For noninteger n*beta, they use exact distortion increments at beta:
# CTE: g(s) = min(s/beta,1).
# GS:  g(s) = q*(1 + 2*Delta*(1-q)), q = min(s/beta,1).
# This includes the fractional tail-boundary weight when n*beta is not
# an integer, without replacing beta by ceiling(n*beta)/n. When n*beta
# is an integer, these reduce to the table's top-k formulas.
# More explicitly, write n*beta = m + a, where m = floor(n*beta)
# and 0 <= a < 1. For CTE, each of the largest m observations has
# weight 1/(n*beta); if a > 0, b[n-m] has weight a/(n*beta).
# All smaller observations have zero weight. GS uses its own distortion
# increment on the same boundary interval, rather than the CTE weight.
# Example: n=10, beta=0.25 gives CTE weights 0.2, 0.4, 0.4 on
# b[8], b[9], b[10], respectively; it does not average the top three.
# If n*beta < 1, both estimators equal the sample maximum.
# Numerical checks (in Python) verified agreement with the old table
# formulas for integer n*beta, unit weight sums and nonnegative weights
# for integer/noninteger cases, and the n*beta < 1 case. R was not
# available for those checks; the full R simulation was not run.
# VaR uses the table's empirical generalized-inverse convention:
# Fhat^{-1}(1-beta) = b[ceiling(n*(1-beta))]. Its step distortion is
# represented by a Stieltjes integral, rather than an ordinary density
# psi(t) dt; the order-statistic formula specifies the quantile convention.
emp_risk <- function(y) {
  b <- sort(y);n <- length(b)
  j <- seq_len(n)
  nbeta <- n*beta
  k <- round(nbeta)
  # Match the previous sample-size tolerance for intended integer n*beta.
  # Use round rather than ceiling here to avoid adding an observation
  # when floating-point arithmetic places an integer slightly above itself.
  # Keep k >= 1 so a very small positive n*beta uses the fractional branch.
  if (k >= 1 && abs(nbeta-k) < 1e-10) {
    # Original table calculations, preserving their arithmetic order.
    i <- (n-k+1L):n
    cte <- mean(tail(b,k))
    gs <- sum(b[i]*(k+2*Delta*(k-2*(n-i)-1))/k^2)
  } else {
    # Noninteger n*beta: use the new exact empirical plug-in weights.
    # For b[j], the quantile interval is ((j-1)/n, j/n]. Its survival
    # endpoints are (n-j+1)/n and (n-j)/n. Divide by beta and cap at 1
    # so intervals below the selected upper tail receive zero weight.
    q_left <- pmin(((n-j+1)/n)/beta,1)
    q_right <- pmin(((n-j)/n)/beta,1)
    w_cte <- q_left-q_right
    # With h(q) = q*(1+2*Delta*(1-q)), the GS increment factors as
    # h(q_left)-h(q_right) = (q_left-q_right)*
    #                       (1+2*Delta*(1-q_left-q_right)).
    # This avoids subtracting two nearly equal h values. For the allowed
    # 0 <= Delta <= 0.5, CTE and GS weights are nonnegative and sum to 1
    # (up to floating-point error); Delta=0 makes GS equal CTE.
    w_gs <- w_cte*(1+2*Delta*(1-q_left-q_right))
    cte <- sum(b*w_cte)
    gs <- sum(b*w_gs)
  }
  c(VaR=b[ceiling(n*(1-beta))],CTE=cte,
    GS=gs,
    PHT=sum(b*((1-(j-1)/n)^r-(1-j/n)^r)),
    WT=sum(b*(pnorm(qnorm(j/n)-lambda)-pnorm(qnorm((j-1)/n)-lambda))))
}
wide <- data.frame(Risk_Measure=rep(measures,each=2L),
                   Estimator=rep(c('ML','EMP'),5L))
for (n in sample_sizes) wide[[paste0('n',n,'_Mean')]] <- NA_real_
wide$nInf_Mean <- as.vector(rbind(asym_ml$mean,asym_emp$mean))
for (n in sample_sizes) wide[[paste0('n',n,'_Scaled_SD')]] <- NA_real_
wide$nInf_Scaled_SD <- as.vector(rbind(asym_ml$scaled_SD,asym_emp$scaled_SD))
diags <- list()
fit_diags <- list()
delta_tables <- list()
optimizer_tables <- list()
save_tables <- function() {
  writeLines('Table 1. Shifted LN; thousands; sqrt(n)*SD; ML ground-up; EMP capped conditional; ML includes finite limiting risks; finite-only retained-risk summaries',output_file)
  write.table(wide,output_file,sep=',',row.names=FALSE,append=TRUE,
              col.names=TRUE,na='',qmethod='double')
  if (length(diags)) {
    cat('\n\nTable 2. Retained estimates and extremes\n',file=output_file,append=TRUE)
    write.table(do.call(rbind,diags),output_file,sep=',',row.names=FALSE,
                append=TRUE,col.names=TRUE,na='',qmethod='double')
    cat('\n\nTable 3. ML fit status counts\n',file=output_file,append=TRUE)
    write.table(do.call(rbind,fit_diags),output_file,sep=',',row.names=FALSE,
                append=TRUE,col.names=TRUE,qmethod='double')
  }
  if(length(delta_tables)) write.csv(do.call(rbind,delta_tables),
    file.path(output_dir,'Shifted_LN_Delta_CI_Summary_v4.csv'),row.names=FALSE,na='')
  if(length(optimizer_tables)) write.csv(do.call(rbind,optimizer_tables),
    file.path(output_dir,'Shifted_LN_Optimizer_Diagnostics_v4.csv'),row.names=FALSE)
}
cat('mu_original:',mu_original,'; internal mu:',mu,'\n')
cat('Cap probability:',atom,'; conditional VaR:',qcond,'\n')
cat('Fisher condition number:',kappa(Info,exact=TRUE),'\n')
print(asym_ml,row.names=FALSE)
print(asym_emp,row.names=FALSE)
save_tables()
set.seed(random_seed)
for (n in sample_sizes) {
  cat('Running n=',n,', B=',B,'\n',sep='');flush.console()
  estimates <- matrix(NA_real_,B,10L)
  pars <- matrix(NA_real_,B,2L)
  reasons <- character(B)
  fit_types <- character(B)
  delta_valid <- delta_covered <- delta_unavailable <- delta_failed <- numeric(5)
  delta_lengths <- numeric(5)
  limiting_infinite <- limiting_ambiguous <- numerical_risk_missing <- numeric(5)
  D$counts <- numeric(0)
  for (b in seq_len(B)) {
    lp <- logSd+log(pmax(runif(n),.Machine$double.xmin))
    z <- qnorm(lp,lower.tail=FALSE,log.p=TRUE)
    y <- pmin(x0+exp(mu+sigma*z),u)
    fit <- fit_ml(y)
    reasons[b] <- fit$reason
    fit_types[b] <- fit$type
    ml <- rep(NA_real_,5L)
    if (fit$ok) {
      if(fit$type=='finite') pars[b,] <- c(fit$par[1L],exp(fit$par[2L]))
      ml <- all_risks(fit)
      if(fit$type=='limit') {
        limiting_infinite <- limiting_infinite+as.integer(is.infinite(ml))
        limiting_ambiguous <- limiting_ambiguous+as.integer(is.na(ml))
      } else numerical_risk_missing <- numerical_risk_missing+as.integer(is.na(ml))
      I <- if(fit$type=='finite') tryCatch(expected_information(fit),
                       error=function(e) NULL) else NULL
      for(j in seq_along(measures)) {
        if(fit$type=='limit') {
          delta_unavailable[j] <- delta_unavailable[j]+1
          next
        }
        dc <- if(is.null(I)) NULL else tryCatch(delta_ci(fit,n,measures[j],I),
                                               error=function(e) NULL)
        if(is.null(dc) || any(!is.finite(dc$ci)) || dc$ci[1]>dc$ci[2]) {
          delta_failed[j] <- delta_failed[j]+1
          next
        }
        delta_valid[j] <- delta_valid[j]+1
        delta_covered[j] <- delta_covered[j]+as.integer(
          dc$ci[1]<=truth[j] && truth[j]<=dc$ci[2])
        delta_lengths[j] <- delta_lengths[j]+diff(dc$ci)
      }
    }
    # EMP is retained regardless of ML success or failure.
    estimates[b,] <- as.vector(rbind(ml,emp_risk(y)))
    if (b%%1000L==0L) {cat('  ',b,'/',B,'\n',sep='');flush.console()}
  }
  retained <- is.finite(estimates)
  avg <- spread <- vmax <- q99 <- rep(NA_real_,10L)
  for (j in 1:10) {
    v <- estimates[retained[,j],j]
    if (length(v)) {
      avg[j] <- mean(v);vmax[j] <- max(v)
      q99[j] <- unname(quantile(v,0.99))
    }
    if (length(v)>=2L) spread[j] <- sqrt(n)*sd(v)
  }
  wide[[paste0('n',n,'_Mean')]] <- avg
  wide[[paste0('n',n,'_Scaled_SD')]] <- spread
  goodpar <- is.finite(pars[,2L])
  diags[[length(diags)+1L]] <- data.frame(n=n,B=B,
    Risk_Measure=wide$Risk_Measure,Estimator=wide$Estimator,
    finite_count=colSums(retained),excluded_count=B-colSums(retained),
    limiting_infinite=as.vector(rbind(limiting_infinite,rep(0,5))),
    limiting_ambiguous=as.vector(rbind(limiting_ambiguous,rep(0,5))),
    numerical_risk_missing=as.vector(rbind(numerical_risk_missing,rep(0,5))),
    finite_q99=q99,finite_max=vmax,
    min_fitted_mu=if(any(goodpar)) min(pars[goodpar,1L]) else NA_real_,
    max_fitted_sigma=if(any(goodpar)) max(pars[goodpar,2L]) else NA_real_)
  counts <- as.data.frame(table(reasons),stringsAsFactors=FALSE)
  names(counts) <- c('fit_status','count')
  fit_diags[[length(fit_diags)+1L]] <- data.frame(n=n,B=B,counts)
  delta_tables[[length(delta_tables)+1L]] <- data.frame(
    n=n,Monte_Carlo_Replications=B,Confidence=confidence,Risk_Measure=measures,
    Original_Finite_Fits=sum(fit_types=='finite'),
    Original_Limit_Fits=sum(fit_types=='limit'),
    Original_Unresolved_Fits=sum(fit_types=='unresolved'),
    Valid=delta_valid,Covered=delta_covered,
    Delta_Limit_Unavailable=delta_unavailable,Delta_Numerical_Failure=delta_failed,
    Coverage_Among_Valid=ifelse(delta_valid>0,delta_covered/delta_valid,NA_real_),
    Covered_Over_Completed=delta_covered/B,
    Mean_Length_Valid_Thousand=ifelse(delta_valid>0,delta_lengths/delta_valid,NA_real_))
  stopifnot(all(delta_valid+delta_unavailable+delta_failed+
                  sum(fit_types=='unresolved')==B))
  if(length(D$counts)) {
    keys <- do.call(rbind,strsplit(names(D$counts),'|',fixed=TRUE))
    optimizer_tables[[length(optimizer_tables)+1L]] <- data.frame(
      n=n,Scope=keys[,1],Reason=keys[,2],Count=unname(D$counts))
  }
  save_tables()
  cat('Saved through n=',n,'\n',sep='')
}
print(wide,row.names=FALSE)
cat('Finished. Output:',output_file,'\n')
# Finite-only MC means and SDs include finite limiting risk values, but exclude
# infinite/ambiguous risks and unresolved computations. They are conditional
# retained-risk summaries, not unconditional moments of the extended estimator.
# Delta intervals are symmetric raw-scale Wald intervals; no endpoint clipping.
# Expected information is evaluated at each finite fit, per conditional censored
# observation, including the point mass at u. The variance divides by n once.
# All numerical values, including interval lengths, are already in thousands.
# Limiting fits do not justify ordinary interior-MLE normal theory. Calibration
# near the boundary must be assessed empirically; 90% coverage is not guaranteed.
# This script performs no parametric bootstrap or profile likelihood calculation.
