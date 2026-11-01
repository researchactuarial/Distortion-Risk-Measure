# OUTPUT: ONE CSV with Table 1 (simulation/asymptotic results) and Table 2 (boundary diagnostics).
# Columns: all means (sample_sizes and Inf), then all scaled SDs.
# n=Inf values are analytical; finite n values are Monte Carlo.
# PARETO TABLE 5.1 - REVIEWED VERSION (BASE R)
# IMPORTANT: ML estimates ground-up X; EMP estimates observed Y=min(X,u)|X>d.
# These are DIFFERENT estimands, not competing estimators of one quantity.
# Monte Carlo moments for nonfinite functionals are reported FINITE-ONLY.
# In particular, even the conditional PHT moments may be mathematically infinite
# because alpha_hat can approach 1/r continuously from above.
# Analytical values are calculated BEFORE Monte Carlo simulation begins.
# The single combined CSV is updated after EACH sample-size batch.
# CSV files will be saved in the SAME folder as this .R file.
# ===================== EDIT PARAMETERS HERE =====================
# Monetary values use THOUSANDS: x0=1 means 1,000; d=4 means 4,000.
x0 <- 1
alpha <- 2
d <- 4
u <- 14
beta <- 0.10        # VaR/CTE tail probability
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
stopifnot(x0 > 0, alpha > 1, d >= x0, u > d, beta > 0, beta < 1,
          d * beta^(-1/alpha) < u,
          Delta >= 0, r > 0, r*alpha > 1, lambda >= 0, B >= 2,
          all(sample_sizes >= 2), all(sample_sizes*beta >= 1))
delta <- Delta
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir)
output_dir <- normalizePath(output_dir, winslash="/", mustWork=TRUE)

cat('Output directory:', output_dir, '\n')
cat('Writing analytical CSV files first...\n')


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

cat("\n=== ANALYTICAL ML n -> infinity ===\n")
local({
  # Analytic n -> infinity Pareto Type I, manuscript Table 5.1
  # No simulation; no packages. Units: thousands.
  p_uncensored <- 1 - (d/u)^alpha
  sd_mle <- function(derivative) abs(derivative) * alpha / sqrt(p_uncensored)
  VaR <- x0 * beta^(-1/alpha)
  CTE <- VaR * alpha/(alpha-1)
  GS <- CTE * (1 + 2*delta/(2*alpha-1))
  PHT <- x0 * (1 + 1/(r*alpha-1))
  # WT integral over Z~N(0,1); log survival avoids numerical underflow.
  wt <- function(a) {
    f <- function(z) x0 * exp(-pnorm(z,lower.tail=FALSE,log.p=TRUE)/a +
                               lambda*z-lambda^2/2 + dnorm(z,log=TRUE))
    integrate(f, lower=-12, upper=12, rel.tol=1e-9,
              subdivisions=500L)$value
  }
  WT <- wt(alpha)
  # Derivatives of ground-up risk functionals with respect to alpha.
  d_VaR <- VaR * log(beta)/alpha^2
  d_CTE <- CTE * (log(beta)/alpha^2 - 1/(alpha*(alpha-1)))
  d_GS <- GS * (log(beta)/alpha^2 - 1/(alpha*(alpha-1)) -
                   4*delta/((2*alpha-1)*(2*alpha-1+2*delta)))
  d_PHT <- -x0*r/(r*alpha-1)^2
  h <- 1e-4
  d_WT <- (wt(alpha+h)-wt(alpha-h))/(2*h)
  measure <- c('VaR','CTE','GS','PHT','WT')
  mean_ML <- c(VaR,CTE,GS,PHT,WT)
  scaled_SD_ML <- sd_mle(c(d_VaR,d_CTE,d_GS,d_PHT,d_WT))
  # Empirical VaR: conditional Pareto quantile and quantile CLT.
  mean_EMP_VaR <- d*beta^(-1/alpha)
  scaled_SD_EMP_VaR <- (d/alpha)*sqrt((1-beta)/beta^(1+2/alpha))
  # Empirical PHT: exact asymptotic variance via influence function for
  # right-censored observations is more involved; do not substitute a formula
  # for the uncensored distribution. Other EMP results need independent
  # derivation / validation before being marked as calculated.
  result <- data.frame(measure, ML_mean=mean_ML,
                       ML_scaled_SD=scaled_SD_ML)
  cat('Analytical ML n=infinity values (thousands):\n')
  print(result, digits=8, row.names=FALSE)
  cat('\nEmpirical VaR: asymptotic mean =',mean_EMP_VaR,
      '; scaled SD =',scaled_SD_EMP_VaR,'\n')
  pareto_asym_ml <<- data.frame(measure=result$measure, mean=result$ML_mean, scaled_SD=result$ML_scaled_SD)
})

cat("\n=== ANALYTICAL EMP n -> infinity ===\n")
local({
  # Deterministic verification of Pareto EMP n -> infinity entries in Table 5.1.
  # Uses only base R. All values are in thousands. No Monte Carlo simulation.
  # X | X>d follows Pareto I with d=4, alpha=2; observations Y=min(X,u).
  # The empirical distortion functional of Y has limit
  # R_g = d + integral_d^u g(S(t)) dt, S(t)=(d/t)^alpha.
  # Its influence function is IF(y) = integral_d^u g'(S(t))
  #                              * (I(y>t)-S(t)) dt.
  # Thus Var(IF(Y)) = E[A(Y)^2]-E[A(Y)]^2, with
  # A(y)=integral_d^y g'(S(t)) dt, and Y has density
  # alpha*d^alpha*y^(-alpha-1) on (d,u), plus mass S(u) at u.
  
  N <- 200000L  # deterministic trapezoidal grid; increase to check convergence
  x <- seq(d,u,length.out=N+1L)
  h <- (u-d)/N
  S <- (d/x)^alpha
  f <- alpha*d^alpha/x^(alpha+1)
  atom <- (d/u)^alpha
  trap <- function(z) h*(sum(z)-0.5*(z[1]+z[length(z)]))
  cumtrap <- function(z) c(0,cumsum(h*(z[-1L]+z[-length(z)])/2))
  # For CTE and GS, include beta-quantile as an exact integration node
  # to avoid quadrature artifacts at the kink in g.
  q <- d*beta^(-1/alpha)
  x <- sort(unique(c(x,q)))
  step <- diff(x)
  S <- (d/x)^alpha
  f <- alpha*d^alpha/x^(alpha+1)
  trap <- function(z) sum(step*(z[-1L]+z[-length(z)])/2)
  cumtrap <- function(z) c(0,cumsum(step*(z[-1L]+z[-length(z)])/2))
  
  calc <- function(g, gp) {
    v <- gp(S)
    A <- cumtrap(v)
    EA <- trap(A*f)+atom*tail(A,1)
    EA2 <- trap(A^2*f)+atom*tail(A,1)^2
    c(mean=d+trap(g(S)), scaled_sd=sqrt(max(0,EA2-EA^2)))
  }
  
  # Conditional top-beta CTE distortion.
  cte <- calc(function(s) pmin(s/beta,1),
              function(s) ifelse(s<beta,1/beta,0))
  # Gini shortfall: g(s)=t[1+2 Delta(1-t)] below beta, 1 above beta.
  gs <- calc(function(s) {t <- pmin(s/beta,1); t*(1+2*Delta*(1-t))},
             function(s) ifelse(s<beta,(1+2*Delta*(1-2*s/beta))/beta,0))
  pht <- calc(function(s) s^r,function(s) r*s^(r-1))
  # Wang transform g(s)=Phi(Phi^-1(s)+lambda); derivative phi(z+lambda)/phi(z).
  # The first grid point S(d)=1 gives qnorm(1)=Inf and derivative zero.
  # R evaluates this endpoint as an improper endpoint limit; replace safely.
  # Recalculate WT with a finite limiting endpoint derivative g'(1)=0.
  wt_gp <- function(s) {
    out <- exp(-lambda*qnorm(pmin(s,1-1e-15))-lambda^2/2)
    out[s>=1] <- 0
    out
  }
  wt <- calc(function(s) pnorm(qnorm(s)+lambda),wt_gp)
  var_mean <- d*beta^(-1/alpha)
  var_sd <- (d/alpha)*sqrt((1-beta)/beta^(1+2/alpha))
  results <- data.frame(measure=c('VaR','CTE','GS','PHT','WT'),
                        mean=c(var_mean,cte['mean'],gs['mean'],pht['mean'],wt['mean']),
                        scaled_sd=c(var_sd,cte['scaled_sd'],gs['scaled_sd'],pht['scaled_sd'],wt['scaled_sd']))
  print(results, row.names=FALSE, digits=9)
  pareto_asym_emp <<- data.frame(measure=results$measure, mean=results$mean, scaled_SD=results$scaled_sd)
  cat('\nNo simulation was used. n=infinity SD is sqrt(Var(IF(Y))).\n')
})

cat("\n=== FINITE-SAMPLE SIMULATION ===\n")
local({
  # Monte Carlo estimates for the Pareto model.
  # Run: Rscript Pareto_Table5_1_reproduction_audit.R
  # No external R packages required. Results saved alongside this script.
  # Replications and sample sizes are set above; units in thousands.
  set.seed(random_seed)
  ns <- sample_sizes
  # Work in thousands throughout.
  # Conditional Pareto simulation X | X>d, followed by censoring at u.
  sample_obs <- function(n) pmin(d * runif(n)^(-1/alpha), u)
  # Exact closed form via 1D numerical integration in probability space.
  # Weighted distortion WT = E[ Q(U) * exp(lambda*Phi^{-1}(U)-lambda^2/2) ].
  # For lambda >= 0, the ground-up WT is infinite when a <= 1.
  wt_pareto <- function(a) {
    if (!is.finite(a) || a <= 0) return(NA_real_)
    if (a <= 1) return(Inf)
    f <- function(z) {
      p <- pnorm(z)
      exp(log(x0) - log(pmax(1-p, .Machine$double.xmin))/a +
            lambda*z - lambda^2/2) * dnorm(z)
    }
    tryCatch(integrate(f, -12, 10, rel.tol=1e-6, subdivisions=200L)$value,
             error=function(e) NA_real_)
  }
  # Fast midpoint quadrature for finite WT; check against adaptive integration
  # at the true parameter and several admissible nearby parameters.
  zgrid <- seq(-9, 9, length.out=2401L)
  zw <- (zgrid[-1L]+zgrid[-length(zgrid)])/2
  wgt <- diff(zgrid) * dnorm(zw) * exp(lambda*zw-lambda^2/2)
  logsurv <- pnorm(zw, lower.tail=FALSE, log.p=TRUE)
  wt_fast <- function(a) {
    if (!is.finite(a) || a <= 0) return(NA_real_)
    if (a <= 1) return(Inf)
    val <- x0*sum(wgt*exp(-logsurv/a))
    if (!is.finite(val)) return(Inf)
    val
  }
  # Quadrature audit: only finite, well-separated shape parameters.
  wt_test_shapes <- unique(c(alpha, max(1.2, alpha*0.8), alpha*1.2))
  wt_check <- data.frame(alpha=wt_test_shapes,
                         adaptive=vapply(wt_test_shapes, wt_pareto, numeric(1)),
                         midpoint=vapply(wt_test_shapes, wt_fast, numeric(1)))
  wt_check$absolute_error <- abs(wt_check$adaptive - wt_check$midpoint)
  print(wt_check, row.names=FALSE)
  if (any(!is.finite(wt_check$absolute_error)) ||
      any(wt_check$absolute_error > 1e-3 * pmax(1, abs(wt_check$adaptive)))) {
    warning("WT fixed-grid quadrature failed its spot check; increase grid accuracy")
  }
  # MLE under left truncation and right censoring; all observed claims >= d.
  one_rep <- function(n) {
    b <- sort(sample_obs(n))
    unc <- b < u
    a <- sum(unc) / sum(log(b/d))
    # If all claims censored, MLE is at alpha=0 boundary.
    if (!is.finite(a) || a <= 0) return(list(a=a, values=rep(NA_real_,10)))
    # 5 ML and 5 empirical risk measures (same order as published table).
    var <- x0*beta^(-1/a)
    cte <- if (a>1) var*a/(a-1) else Inf
    gs <- if (a>1) cte*(1+2*delta/(2*a-1)) else Inf
    pht <- if (r*a>1) x0*(1+1/(r*a-1)) else Inf
    wt <- wt_fast(a)
    # Exact EMP distortion plug-in estimators (eq. empDM), with b sorted:
    #   sum_i b[i] * {g(1-(i-1)/n) - g(1-i/n)}.
    # Table 4.1 uses k=ceiling(n*beta) for CTE and GS, effectively using
    # beta_n=k/n. The weights below retain the original beta even when
    # n*beta is noninteger, giving the boundary claim a fractional weight.
    # They agree with the table when n*beta is an integer (as for the
    # current beta=0.10 and n=50,100,500). No rounding of beta is needed.
    # For CTE this is equivalent to
    #   (sum of largest m claims + theta*b[n-m])/(n*beta),
    # m=floor(n*beta), theta=n*beta-m; the sum is empty when m=0.
    # CTE distortion: g(s)=min(s/beta,1).
    # GS distortion: g(s)=t*(1+2*delta*(1-t)), t=min(s/beta,1).
    # VaR is the empirical quantile b[ceiling(n*(1-beta))]. Its distortion
    # has a jump, so its quantile representation uses a Stieltjes integral
    # with the matching quantile convention, not a density psi(t).
    emp_var <- b[ceiling(n*(1-beta))]
    j <- seq_len(n)
    nbeta <- n*beta
    # Treat values within floating-point roundoff of an integer as integer.
    # The bounds prevent a tiny positive nbeta from being rounded to zero.
    k_nearest <- round(nbeta)
    nbeta_is_integer <- k_nearest >= 1 && k_nearest <= n &&
      abs(nbeta-k_nearest) <= 8*.Machine$double.eps*max(1,abs(nbeta))
    if (nbeta_is_integer) {
      # Original attached code's CTE/GS formulas and evaluation order.
      # Using round rather than ceiling avoids an extra claim when nbeta
      # exceeds an integer only because of floating-point roundoff.
      k <- as.integer(k_nearest)
      top <- tail(b,k)
      emp_cte <- mean(top)
      i <- (n-k+1):n
      emp_gs <- sum(b[i]*(k+2*delta*(k-2*(n-i)-1))/k^2)
    } else {
      # Noninteger nbeta: exact distortion increments at the specified beta.
      # If nbeta=m+theta, CTE gives the largest m claims weight 1/nbeta
      # and the boundary claim b[n-m] weight theta/nbeta.
      # For nbeta<1, both CTE and GS put all weight on the largest claim.
      t_left <- pmin((n-j+1)/nbeta, 1)
      t_right <- pmin((n-j)/nbeta, 1)
      emp_cte <- sum(b*(t_left-t_right))
      # Factored GS difference: g(l)-g(r) =
      # (l-r)*(1+2*delta-2*delta*(l+r)).
      emp_gs <- sum(b*(t_left-t_right)*
                    (1+2*delta-2*delta*(t_left+t_right)))
    }
    # Formula checks in Python covered integer and fractional nbeta,
    # nbeta<1, beta=1, and the integer-branch roundoff tolerance.
    # R was unavailable; the full R simulation was not rerun.
    emp_pht <- sum(b*((1-(j-1)/n)^r-(1-j/n)^r))
    emp_wt <- sum(b*(pnorm(qnorm(j/n)-lambda)-pnorm(qnorm((j-1)/n)-lambda)))
    list(a=a,values=c(var,emp_var,cte,emp_cte,gs,emp_gs,pht,emp_pht,wt,emp_wt))
  }
  labels <- c('VaR_ML','VaR_EMP','CTE_ML','CTE_EMP','GS_ML','GS_EMP',
              'PHT_ML','PHT_EMP','WT_ML','WT_EMP')
  results <- list()
  boundary_results <- list()
  output_file <- file.path(output_dir, 'Pareto_Combined_Results.csv')
  diagnostics_file <- file.path(output_dir, 'Pareto_PHT_Near_Boundary_Diagnostics.csv')
  pht_diagnostics <- list()
  for (h in seq_along(ns)) {
    n <- ns[h]
    cat('Running n =', n, 'with', B, 'replicates; output:', output_file, '\n')
    flush.console()
    estimates <- matrix(NA_real_, nrow=B, ncol=length(labels))
    shape_est <- numeric(B)
    for (rep in seq_len(B)) {
      out <- one_rep(n)
      shape_est[rep] <- out$a
      estimates[rep,] <- out$values
    }
    # Count nonfinite ML risk measures due to the fitted Pareto shape.
    # CTE and GS require alpha_hat > 1; PHT requires r*alpha_hat > 1.
    # NA/NaN/nonpositive alpha_hat is counted as invalid for all three.
    bad_shape <- !is.finite(shape_est) | shape_est <= 0
    cte_gs_bad <- bad_shape | shape_est <= 1
    pht_bad <- bad_shape | r*shape_est <= 1
    cte_gs_count <- sum(cte_gs_bad)
    pht_count <- sum(pht_bad)
    # Finite estimates close to the pole can dominate reported PHT moments.
    finite_pht <- estimates[is.finite(estimates[,7L]), 7L]
    valid_shapes <- shape_est[is.finite(shape_est) & shape_est > 1/r]
    pht_diagnostics[[h]] <- data.frame(
      n=n, B=B, boundary=1/r, violation_count=pht_count,
      violation_percent=100*pht_count/B,
      finite_count=length(finite_pht),
      count_within_0_01_above_boundary=sum(valid_shapes < 1/r+0.01),
      count_within_0_05_above_boundary=sum(valid_shapes < 1/r+0.05),
      min_positive_distance_to_boundary=if(length(valid_shapes)) min(valid_shapes-1/r) else NA_real_,
      max_finite_PHT=if(length(finite_pht)) max(finite_pht) else NA_real_,
      PHT_finite_q99=if(length(finite_pht)) unname(quantile(finite_pht,.99)) else NA_real_,
      PHT_finite_q999=if(length(finite_pht)) unname(quantile(finite_pht,.999)) else NA_real_)
    # CTE and GS share the same boundary event; calculate once and display
    # separately for readability. Rates are proportions and percentages.
    boundary_results[[h]] <- data.frame(
      Risk_Measure=c("CTE", "GS", "PHT"),
      Violation_Condition=c(
        "alpha_hat <= 1 (or invalid alpha_hat)",
        "alpha_hat <= 1 (or invalid alpha_hat)",
        paste0("r * alpha_hat <= 1; r = ", format(r, trim=TRUE),
               "; alpha_hat <= ", format(1/r, digits=8, trim=TRUE),
               " (or invalid alpha_hat)")
      ),
      Sample_Size=rep(n, 3L),
      Replications=rep(B, 3L),
      Violation_Count=c(cte_gs_count, cte_gs_count, pht_count),
      Violation_Rate=c(cte_gs_count, cte_gs_count, pht_count)/B,
      Violation_Percent=100*c(cte_gs_count, cte_gs_count, pht_count)/B,
      stringsAsFactors=FALSE
    )
    # Descriptive FINITE-ONLY Monte Carlo summaries. These are not estimates
    # of unconditional finite moments. PHT moments can also be infinite even
    # after conditioning on alpha_hat > 1/r (pole at the boundary).
    finite_n <- colSums(is.finite(estimates))
    conditional_means <- sapply(seq_along(labels), function(j) {
      v <- estimates[,j]; if (!any(is.finite(v))) NA_real_ else mean(v[is.finite(v)])
    })
    conditional_sds <- sapply(seq_along(labels), function(j) {
      v <- estimates[,j]; if (sum(is.finite(v))<2) NA_real_ else sd(v[is.finite(v)])*sqrt(n)
    })
    results[[h]] <- data.frame(n=n, measure=labels,
       simulated_mean_finite_only=conditional_means,
       simulated_scaled_sd_finite_only=conditional_sds,
       finite_replicates=finite_n, total_replicates=B,
       pht_boundary_count=pht_count, pht_boundary_rate=pht_count/B)
    partial <- do.call(rbind, results)
    partial$estimator <- sub(".*_(ML|EMP)$", "\\1", partial$measure)
    partial$measure <- sub("_(ML|EMP)$", "", partial$measure)
    names(partial)[names(partial)=="simulated_mean_finite_only"] <- "mean"
    names(partial)[names(partial)=="simulated_scaled_sd_finite_only"] <- "scaled_SD"
    wide <- make_wide_table(partial, pareto_asym_ml, pareto_asym_emp, sample_sizes)
    # Save both tables into ONE CSV, with a heading before each table.
    # Table 2 is placed below Table 1 in the same file.
    boundary_wide <- do.call(rbind, boundary_results)
    boundary_wide <- boundary_wide[order(match(boundary_wide$Risk_Measure,
                                              c("CTE", "GS", "PHT")),
                                         boundary_wide$Sample_Size), , drop=FALSE]
    writeLines(paste0(
      "Table 1. Pareto simulation and asymptotic results; ",
      "finite-n moments are finite-only descriptive summaries; ",
      "scaled SD = sqrt(n) * SD; ML targets ground-up X, ",
      "EMP targets observed capped conditional Y"), output_file)
    write.table(wide, file=output_file, sep=",", row.names=FALSE,
                col.names=TRUE, append=TRUE, na="", qmethod="double")
    cat("\n\nTable 2. Pareto ML finiteness-boundary violation rates\n",
        file=output_file, append=TRUE)
    write.table(boundary_wide, file=output_file, sep=",", row.names=FALSE,
                col.names=TRUE, append=TRUE, qmethod="double")
    write.csv(do.call(rbind, pht_diagnostics), diagnostics_file, row.names=FALSE)
    cat('Saved both tables and boundary diagnostics through n =', n, '; CTE/GS:', cte_gs_count, '; PHT:', pht_count, 'out of', B, '\n')
    flush.console()
  }
  res <- do.call(rbind,results)
  # Combined CSV already written after each sample-size batch.
  print(res, row.names=FALSE, digits=5)
  cat("\nBoundary rates (rates are proportions; percent columns are percentages):\n")
  print(do.call(rbind, boundary_results), row.names=FALSE)
  cat('NOTE: finite-only PHT moments are descriptive and may not converge\n',
      'to finite population moments because of the pole at alpha_hat=1/r.\n')
})


cat('\nFINISHED. CSV files saved in: ', output_dir, '\n')
print(file.path(output_dir, "Pareto_Combined_Results.csv"))

cat("\nCSV output folder: ", output_dir, "\n", sep="")
cat("CSV files currently saved:\n")
cat("Combined CSV: ", file.path(output_dir, "Pareto_Combined_Results.csv"), "\n", sep="")
cat("Table 2 is below Table 1 in the combined CSV.\n")

cat("PHT near-boundary diagnostics: ", file.path(output_dir, "Pareto_PHT_Near_Boundary_Diagnostics.csv"), "\n", sep="")
