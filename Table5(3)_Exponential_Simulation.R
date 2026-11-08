# OUTPUT: ONE CSV beside this R file. Rows: ML first, then EMP, for each risk measure.
# Columns: n=50, n=500, n=1000, n=Inf, each with mean and sqrt(n)*SD.
# n=Inf values are analytical; finite n values are Monte Carlo.
# EXPONENTIAL: simulation + analytic n -> infinity (ML, EMP)
# Base R only. No comparison to published tables. Units: thousands.
# CSV files are saved beside this .R file. Set B <- 1000L for a quick test.
# ===================== EDIT PARAMETERS HERE =====================
# Monetary values use THOUSANDS: x0=1 means 1,000; theta=1 means 1,000.
x0 <- 1
theta <- 1
d <- 4
u <- 14
beta <- 0.10        # VaR/CTE tail probability
# Table 4.1 keeps the fractional VaR boundary when n*beta is noninteger.
# For GS, use 2*delta outside the braces in the EMP row to match its target.
Delta <- 0.25       # GS distortion parameter
r <- 0.75           # PHT distortion parameter
lambda <- 0.50      # Wang transform parameter
sample_sizes <- c(50L, 100L, 500L)
B <- 100000L        # Monte Carlo repetitions
random_seed <- 200L

# CSV SAVE FOLDER: save CSVs BESIDE THIS .R SCRIPT.
# In RStudio, open the saved .R file in the editor and click Source.
# If run from command line, use: Rscript path/to/this_script.R
# If run with source(), use: source("path/to/this_script.R")
script_file <- NULL
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  context <- tryCatch(rstudioapi::getSourceEditorContext(), error=function(e) NULL)
  if (!is.null(context) && nzchar(context$path)) script_file <- context$path
}
if (is.null(script_file)) {
  frames <- sys.frames()
  paths <- lapply(frames, function(frame) frame$ofile)
  paths <- Filter(function(p) is.character(p) && length(p)==1L && nzchar(p), paths)
  if (length(paths)) script_file <- tail(paths, 1L)[[1L]]
}
if (is.null(script_file)) {
  cmd <- commandArgs(trailingOnly=FALSE)
  file_arg <- grep("^--file=", cmd, value=TRUE)
  if (length(file_arg)) script_file <- sub("^--file=", "", file_arg[1L])
}
if (is.null(script_file) || !file.exists(script_file)) {
  stop("Cannot determine this R script's folder. Save and Source the file in RStudio, or run source('full/path/to/file.R').")
}
output_dir <- dirname(normalizePath(script_file, winslash="/", mustWork=TRUE))
# ================================================================
stopifnot(theta > 0, d >= x0, u > d, beta > 0, beta < 1,
          Delta >= 0, r > 0, B >= 2,
          all(sample_sizes >= 2))
delta <- Delta
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir)
output_dir <- normalizePath(output_dir, winslash="/", mustWork=TRUE)
set.seed(random_seed)
outdir <- output_dir
cat('Output directory:', normalizePath(outdir), '\n')


# Assemble one paper-style wide table, without any published values.
make_wide_table <- function(sim, asym_ml, asym_emp, ns) {
  measures <- c("VaR", "CTE", "GS", "PHT", "WT")
  methods <- c("ML", "EMP")
  ans <- data.frame(Risk_Measure=rep(measures, each=2L),
                    Estimator=rep(methods, times=length(measures)),
                    stringsAsFactors=FALSE)
  for (n in ns) {
    means <- sds <- rep(NA_real_, nrow(ans))
    for (i in seq_len(nrow(ans))) {
      sel <- sim$measure == ans$Risk_Measure[i] & sim$estimator == ans$Estimator[i] & sim$n == n
      if (any(sel)) {
        means[i] <- sim$mean[which(sel)[1L]]
        sds[i] <- sim$scaled_SD[which(sel)[1L]]
      }
    }
    ans[[paste0("n", n, "_Mean")]] <- means
    ans[[paste0("n", n, "_Scaled_SD")]] <- sds
  }
  ans[["nInf_Mean"]] <- ans[["nInf_Scaled_SD"]] <- NA_real_
  for (i in seq_len(nrow(ans))) {
    a <- if (ans$Estimator[i] == "ML") asym_ml else asym_emp
    j <- match(ans$Risk_Measure[i], a$measure)
    if (!is.na(j)) {
      ans$nInf_Mean[i] <- a$mean[j]
      ans$nInf_Scaled_SD[i] <- a$scaled_SD[j]
    }
  }
  # Paper-style column orientation: all means first, then all scaled SDs.
  mean_columns <- c(paste0("n", ns, "_Mean"), "nInf_Mean")
  sd_columns <- c(paste0("n", ns, "_Scaled_SD"), "nInf_Scaled_SD")
  ans[, c("Risk_Measure", "Estimator", mean_columns, sd_columns), drop=FALSE]
}

# Wang distortion on survival probabilities: g(s) = Phi(Phi^-1(s)+lambda)
wang <- function(s) {
  z <- numeric(length(s))
  middle <- s > 0 & s < 1
  z[s >= 1] <- 1
  z[middle] <- pnorm(qnorm(s[middle]) + lambda)
  z
}

# Ground-up risk measures: R = x0 + theta * coefficient
# WT coefficient = integral_0^1 [-log(s)] g'(s) ds
wt_coefficient <- integrate(function(z) {
  -pnorm(z, log.p=TRUE) * dnorm(z + lambda)
}, -12, 12, rel.tol=1e-10, subdivisions=500L)$value
coeff <- c(VaR=-log(beta), CTE=1-log(beta),
           GS=1-log(beta)+delta, PHT=1/r, WT=wt_coefficient)
true_risk <- x0 + theta*coeff
p_uncens <- 1-exp(-(u-d)/theta)
# theta_ML asymptotic variance: theta^2 / p_uncens
parametric <- data.frame(estimator='ML', measure=names(coeff),
                         mean=unname(true_risk),
                         scaled_SD=unname(theta*abs(coeff)/sqrt(p_uncens)))
# Analytical ML results retained in memory for combined CSV
cat('Saved Exponential_ML_Asymptotic.csv\n')

# Empirical limit for Y = min(X,u) | X>d, with X exponential shifted by x0.
# On [d,u), S_Y(t)=exp(-(t-d)/theta), density f_Y(t)=S_Y(t)/theta,
# plus a point mass S_Y(u) at u.
# Distortion estimator R_g = d + integral_d^u g(S_Y(t)) dt.
# Influence function IF(y)=A(y)-E[A(Y)], A(y)=integral_d^y g'(S_Y(t))dt.
N <- 50000L
xx <- sort(unique(c(seq(d,u,length.out=N+1L), d-theta*log(beta))))
dx <- diff(xx)
ss <- exp(-(xx-d)/theta)
ff <- ss/theta
atom <- exp(-(u-d)/theta)
trap <- function(v) sum(dx*(head(v,-1L)+tail(v,-1L))/2)
cumtrap <- function(v) c(0,cumsum(dx*(head(v,-1L)+tail(v,-1L))/2))
emp_calc <- function(g,gprime) {
  A <- cumtrap(gprime(ss))
  EA <- trap(A*ff) + atom*tail(A,1)
  EA2 <- trap(A*A*ff) + atom*tail(A,1)^2
  c(mean=d+trap(g(ss)), scaled_SD=sqrt(max(0,EA2-EA^2)))
}
cte <- emp_calc(function(s) pmin(s/beta,1),
                function(s) ifelse(s<beta,1/beta,0))
gs <- emp_calc(function(s) { t<-pmin(s/beta,1); t*(1+2*delta*(1-t)) },
               function(s) ifelse(s<beta,(1+2*delta*(1-2*s/beta))/beta,0))
pht <- emp_calc(function(s) s^r, function(s) r*s^(r-1))
wt <- emp_calc(wang, function(s) {
  z <- exp(-lambda*qnorm(pmin(s,1-1e-15))-lambda^2/2)
  z[s>=1] <- 0
  z
})
# The q_(1-beta) observed quantile is below u for the chosen parameters.
emp_var <- c(mean=d-theta*log(beta), scaled_SD=theta*sqrt((1-beta)/beta))
empirical <- data.frame(estimator='EMP', measure=c('VaR','CTE','GS','PHT','WT'),
                        mean=c(emp_var['mean'],cte['mean'],gs['mean'],pht['mean'],wt['mean']),
                        scaled_SD=c(emp_var['scaled_SD'],cte['scaled_SD'],gs['scaled_SD'],pht['scaled_SD'],wt['scaled_SD']))
# Analytical EMP results retained in memory for combined CSV
cat('Saved Exponential_EMP_Asymptotic.csv\n')

# Simulation: generate conditional exponential above d directly; censor at u.
# Estimator order: ML and EMP for each measure.
# Batch simulation avoids storing all 100,000 x n observations in memory.
measure_names <- names(coeff)
all_rows <- list()
for (n in sample_sizes) {
  cat('Starting n =',n,'with',B,'replications\n')
  accum <- matrix(NA_real_, nrow=B, ncol=10L)
  labels <- as.vector(sapply(measure_names, function(m) paste(m,c('ML','EMP'),sep='_')))
  colnames(accum) <- labels
  batch_size <- 500L
  for (start in seq.int(1L,B,by=batch_size)) {
    stop_idx <- min(B,start+batch_size-1L)
    nb <- stop_idx-start+1L
    raw <- matrix(d+rexp(n*nb,rate=1/theta), nrow=n, ncol=nb)
    # Equivalent to d + Exp(theta); left truncation handled exactly.
    y <- pmin(raw,u)
    y <- apply(y,2,sort)
    if (is.null(dim(y))) y <- matrix(y,nrow=n)
    uncens <- colSums(y<u)
    exposure <- colSums(pmin(y,u)-d)
    th_ml <- ifelse(uncens>0,exposure/uncens,NA_real_)
    risk_ml <- outer(th_ml,coeff,'*')+x0
    # EMP implements eq. (empDM): R[F_EMP] = sum_i y_(i) * w_i,
    # w_i = g(1-(i-1)/n) - g(1-i/n), with y_(i) ascending.
    # For an ordinary spectral density, psi(t) = g'(1-t), so these
    # weights equal the integrals of psi over [(i-1)/n, i/n].
    # PHT and WT below agree directly with their table formulas.
    # For n*beta = m+a, the CTE numerator is the sum of the largest
    # m claims plus a times the VaR boundary claim.
    # Each empirical functional is an L-statistic of sorted observed claims.
    p <- (0:n)/n
    # VaR has a point-mass spectral measure at t=1-beta, rather than
    # an ordinary density psi. The generalized quantile integral gives
    # y_(ceiling(n*(1-beta))), matching the table's empirical VaR.
    w_var <- rep(0,n); w_var[ceiling(n*(1-beta))] <- 1
    # For ascending order, weights use survival probabilities.
    surv_edges <- (n:0)/n
    # Avoid rounding mathematically integral n*beta upward due to
    # floating-point multiplication (e.g., a value just above an integer).
    tail_count <- n*beta
    nearest_count <- round(tail_count)
    if (abs(tail_count-nearest_count) <=
        8*.Machine$double.eps*max(1,abs(tail_count)) && nearest_count >= 1) {
      tail_count <- nearest_count
    }
    # Each order-statistic interval intersects the upper tail in [lo,hi]
    # on the normalized survival scale s/beta. This handles n*beta<1,
    # integer n*beta, and a partially weighted boundary observation.
    tail_lo <- pmin((n-seq_len(n))/tail_count, 1)
    tail_hi <- pmin((n-seq_len(n)+1)/tail_count, 1)
    w_cte <- tail_hi-tail_lo
    # The partially included boundary claim has CTE
    # weight a/(n*beta), where a=n*beta-floor(n*beta); each fully
    # included upper-tail claim has weight 1/(n*beta).
    # Integrate the GS distortion over the same clipped intervals.
    # Algebraically equivalent to H(hi)-H(lo), H(t)=t*(1+2*delta*(1-t)).
    w_gs <- (tail_hi-tail_lo)*(1+2*delta-2*delta*(tail_hi+tail_lo))
    w_pht <- -diff(surv_edges^r)
    w_wt <- -diff(wang(surv_edges))
    emp <- cbind(as.vector(crossprod(w_var,y)),
                 as.vector(crossprod(w_cte,y)),
                 as.vector(crossprod(w_gs,y)),
                 as.vector(crossprod(w_pht,y)),
                 as.vector(crossprod(w_wt,y)))
    ids <- start:stop_idx
    for (j in seq_along(measure_names)) {
      k <- 2L*(j-1L)
      accum[ids,k+1L] <- risk_ml[,j]
      accum[ids,k+2L] <- emp[,j]
    }
  }
  tab <- data.frame(n=n, estimator=rep(c('ML','EMP'),times=5),
                    measure=rep(measure_names,each=2),
                    mean=colMeans(accum,na.rm=TRUE),
                    scaled_SD=apply(accum,2,sd,na.rm=TRUE)*sqrt(n),
                    valid_reps=colSums(is.finite(accum)), total_reps=B,
                    row.names=NULL)
  all_rows[[as.character(n)]] <- tab
  wide <- make_wide_table(do.call(rbind, all_rows), parametric, empirical, sample_sizes)
  write.csv(wide, file.path(outdir, "Exponential_Combined_Results.csv"), row.names=FALSE)
  cat('Saved simulation results through n =',n,'\n')
}
cat('FINISHED. CSV files are in:',normalizePath(outdir),'\n')

cat("\nCSV output folder: ", output_dir, "\n", sep="")
cat("CSV files currently saved:\n")
cat("Combined CSV: ", file.path(outdir, "Exponential_Combined_Results.csv"), "\n", sep="")
