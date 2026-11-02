# SHIFTED EXPONENTIAL ONLY: truncated/censored Monte Carlo study.
# Ground-up: X=x0+E, E~Exponential(scale=theta).
# Observe n claims conditional on X>=d; right-censor at censor_at.
# Ground-up risks: VaR, CTE, GS, PHT, WT.
# 90% CIs: delta, parametric percentile bootstrap, likelihood ratio.
# Closed-form MLE; vectorized bootstrap; LR ratios cached by exact-claim count.
# Delta variance matches the supplied table: theta^2*A^2/(n*p(theta)),
# p(theta)=1-exp(-(censor_at-d)/theta); x0 is fixed.
# R execution was unavailable during authoring. No simulation results supplied.

# 1. STUDY SETTINGS
# Production uses the full Monte Carlo budget; validation is a small pipeline
# check, not a study from which to draw statistical conclusions. N_MC counts
# outer datasets per sample size, whereas B counts bootstrap datasets per fit.
MODE <- 'production'  # 'validation' for the first pipeline check
N_MC <- if (MODE=='production') 10000L else 10L
B <- if (MODE=='production') 500L else 30L
# Each n is the number of observed claims above the deductible, including
# claims whose values are only known to exceed the censoring threshold.
sample_sizes <- c(50L,100L,500L)
# x0 is the known ground-up location. theta_true is the true exponential
# scale (and mean excess); d is left truncation, not a subtracted payment.
# Claims at censor_at are recorded as censored, not as exact observations.
x0 <- 1000
theta_true <- 1000  # SCALE; exponential rate is 1/theta_true
d <- 4000
censor_at <- 14000
# A central two-sided confidence level is used for all three methods.
# beta is the upper-tail probability for VaR/CTE/GS (VaR level=1-beta).
# delta adjusts GS; r is the proportional-hazards distortion parameter;
# lambda is the Wang-transform parameter. All remain fixed during simulation.
confidence <- .90
beta <- .10; delta <- .25; r <- .75; lambda <- .50
# Require at least ceiling(B * bootstrap_required_fraction) usable bootstrap
# estimates before returning endpoints. At 1, any failed draw fails that CI.
# Blocking limits memory without changing the number of attempted draws.
bootstrap_required_fraction <- 1  # .979 means >=979 valid when B=1000
bootstrap_block <- 250L          # bound temporary matrix memory
# Each outer replication later receives its own deterministic seed.
# Checkpoints are saved periodically and at the final replication.
# resume reads an existing checkpoint only if its entire configuration agrees.
seed <- 200L
checkpoint_every <- 100L
resume <- TRUE
RUN_CHECKS <- TRUE
# Version identifies the output/checkpoint schema. All risk factors and
# interval vectors must follow the same ordering as measures.
version <- 'Exp_publication_metrics_v5'
delta_information <- 'expected_Fisher_per_observation'
measures <- c('VaR','CTE','GS','PHT','WT')
methods <- c('Delta','Bootstrap_percentile','Profile_likelihood')

# 2. OUTPUT LOCATION AND INPUT CHECKS
# Locate this file under Rscript, source(), or RStudio, in that order.
# If no usable script path is found, write outputs beneath getwd().
# The optional RStudio package is used only for path discovery.
script <- NULL
args <- grep('^--file=',commandArgs(FALSE),value=TRUE)
if(length(args)) script <- sub('^--file=','',args[1])
if(is.null(script)) for(fr in rev(sys.frames())) {
  if(!is.null(fr$ofile) && file.exists(fr$ofile)) {script <- fr$ofile; break}
}
if(is.null(script) && requireNamespace('rstudioapi',quietly=TRUE)) {
  script <- tryCatch(if(rstudioapi::isAvailable())
    rstudioapi::getSourceEditorContext()$path else NULL,error=function(e) NULL)
}
folder <- if(length(script)==1 && !is.na(script) && nzchar(script) && file.exists(script))
  dirname(normalizePath(script)) else getwd()
# Separate output folders by mode, distortion settings and schema version.
# Other settings are protected on resume by the saved configuration check.
outdir <- file.path(folder,sprintf('MC_Exp_%s_r%.2f_lambda%.2f_%s',MODE,r,lambda,version))
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
# Reject unsupported parameter values before starting the expensive loop.
stopifnot(MODE %in% c('production','validation'),N_MC>0,B>=20,
          all(sample_sizes>=2),theta_true>0,x0>0,d>x0,censor_at>d,
          beta>0,beta<1,delta>=0,delta<.5,r>0,r<=1,lambda>0,
          confidence>0,confidence<1,bootstrap_block>=1,
          bootstrap_required_fraction>0,bootstrap_required_fraction<=1)
# H is the observable excess range above d. Memorylessness removes x0 and
# d from the conditional excess distribution; censoring still acts at H.
# zcrit gives central normal endpoints; cutoff gives the one-parameter LR cutoff.
H <- censor_at-d
zcrit <- qnorm((1+confidence)/2)
cutoff <- qchisq(confidence,1)

# 3. GROUND-UP RISK FUNCTIONS AND DELTA VARIANCES
# All five targets have positive slope A with respect to theta, so intervals
# for theta can be mapped to risks without reversing endpoints.
# These targets describe the original X distribution, not X conditional on X>=d.
# Every ground-up risk is R=x0+A*theta.
# WT factor E[-log(1-Phi(Z))], Z~Normal(lambda,1).
# Compute the Wang slope once. A log survival probability avoids cancellation
# in 1-Phi(z) in the upper tail; dnorm weights Z with mean lambda and SD 1.
# Integration tolerances concern numerical quadrature, not Monte Carlo error.
wt_factor <- integrate(function(z) {
  -pnorm(z,lower.tail=FALSE,log.p=TRUE)*dnorm(z,mean=lambda)
},-Inf,Inf,rel.tol=1e-9,subdivisions=1000)$value
# Named slopes encode the closed-form exponential risk expressions.
# PHT denotes the proportional-hazards transform and WT the Wang transform.
A <- c(VaR=-log(beta),CTE=1-log(beta),GS=1-log(beta)+delta,
       PHT=1/r,WT=wt_factor)
risk <- function(theta) x0+A*theta
# Expected information per observation: I_theta=p(theta)/theta^2.
# The table gives theta^2*A^2/p(theta), the variance coefficient for
# sqrt(n)*(Rhat-R); divide this coefficient by n for the actual CI.
delta_table_components <- function(theta) {
  # p(theta) is the probability of an exact observation conditional on entry.
  # -expm1(-H/theta) evaluates 1-exp(-H/theta) accurately for small H/theta.
  # The delta method multiplies inverse information by (dR/dtheta)^2=A^2.
  probability <- -expm1(-H/theta)
  coefficient <- theta^2*A^2/probability
  list(probability=probability,coefficient=coefficient)
}
# Evaluate population targets once for subsequent coverage indicators.
# The following CSVs document the true targets, censoring probability and
# variance coefficients independently of any estimated interval.
true_risks <- risk(theta_true)
stopifnot(all(is.finite(A)),all(A>0))
write.csv(data.frame(Risk_Measure=measures,True_Risk=unname(true_risks)),
          file.path(outdir,'True_Ground_Up_Risks.csv'),row.names=FALSE)
write.csv(data.frame(Conditional_Censoring_Probability=exp(-H/theta_true)),
          file.path(outdir,'True_Censoring_Probability.csv'),row.names=FALSE)
true_table <- delta_table_components(theta_true)
table_audit <- data.frame(Risk_Measure=measures,True_Risk=unname(true_risks),
  Risk_Factor=unname(A),Uncensored_Probability=true_table$probability,
  Variance_Coefficient=unname(true_table$coefficient))
# The audit SEs use the true theta and are theoretical reference values.
# Actual delta intervals below plug in the fitted theta for each dataset.
for(nn in sample_sizes) table_audit[[paste0('Asymptotic_SE_n',nn)]] <-
  sqrt(unname(true_table$coefficient)/nn)
write.csv(table_audit,file.path(outdir,'True_Risks_And_Delta_Table.csv'),row.names=FALSE)


# 4. CONDITIONAL DATA GENERATION AND SUFFICIENT STATISTICS
# The exponential memoryless property gives X-d ~ Exp(scale=theta) after entry.
# A censored claim contributes exposure H and no exact-event log(theta) term.
# Retain only m, T and k: individual claims are unnecessary for these fits.
# Conditional on X>=d: X=d+E, E exponential with scale theta.
# m=number exact; T=sum(min(E,H)), including censored exposure H.
# ell(theta)=-m*log(theta)-T/theta, up to constants.
# theta_hat=T/m. Information for eta=log(theta) at MLE equals m.
simulate_statistics <- function(theta,n) {
  # Generate excesses directly from the conditional distribution; this avoids
  # rejection sampling of ground-up claims below the truncation threshold.
  e <- rexp(n,rate=1/theta)
  list(m=sum(e<H),T=sum(pmin(e,H)),k=sum(e>=H))
}

# 5. LIKELIHOOD-RATIO INTERVALS
# The fitted scale is positive when m>0. In log-ratio coordinates the LR
# boundary no longer depends on T, so all datasets with the same m reuse roots.
# The returned pair contains multiplicative lower and upper scale ratios.
# LR depends only on m and t=log(theta/theta_hat).
# D(t)=2*m*(t+exp(-t)-1). Cache ratios for each observed m.
lr_cache <- new.env(parent=emptyenv())
lr_ratios <- function(m) {
  if(m<=0) return(c(NA_real_,NA_real_))
  key <- as.character(m)
  if(exists(key,envir=lr_cache,inherits=FALSE)) return(get(key,envir=lr_cache))
  # Subtract the chi-square cutoff to turn the two CI boundaries into roots.
  # expm1 improves accuracy near t=0, where exp(-t)-1 involves cancellation.
  f <- function(t) 2*m*(t+expm1(-t))-cutoff
  # Expand brackets away from zero until both LR boundaries are enclosed.
  # Solve each side separately, then exponentiate back from log-ratio units.
  left <- -1; right <- 1
  while(f(left)<0) left <- left*2
  while(f(right)<0) right <- right*2
  ans <- exp(c(uniroot(f,c(left,0),tol=1e-10)$root,
               uniroot(f,c(0,right),tol=1e-10)$root))
  assign(key,ans,envir=lr_cache)
  ans
}

# 6. PARAMETRIC PERCENTILE BOOTSTRAP
# Resample from the fitted conditional/censored model using the original n
# and H. This reproduces both exact observations and censoring under the fit.
# Bootstrap fit failures remain visible; they are never replaced by new draws.
# B conditional/censored bootstrap samples in memory-limited matrix blocks.
# Exactly B samples are attempted; failed samples are not replaced.
bootstrap_theta <- function(theta_hat,n) {
  # Initialize missing estimates so unsuccessful fits remain explicitly absent.
  estimates <- rep(NA_real_,B)
  for(first in seq.int(1L,B,by=bootstrap_block)) {
    # Each matrix column is one bootstrap dataset; the last block may be smaller.
    last <- min(B,first+bootstrap_block-1L)
    nb <- last-first+1L
    e <- matrix(rexp(n*nb,rate=1/theta_hat),nrow=n,ncol=nb)
    # Column sums compute the same sufficient statistics as the outer simulation.
    m <- colSums(e<H)
    T <- colSums(pmin(e,H))
    values <- rep(NA_real_,nb)
    # No exact claims means no finite scale MLE. Also reject invalid exposure
    # before division; a later check rejects nonfinite or nonpositive estimates.
    good <- m>0 & is.finite(T) & T>0
    values[good] <- T[good]/m[good]
    estimates[first:last] <- values
  }
  # Percentiles use only usable fits, subject to the prespecified threshold.
  # R quantile type=1 is the inverse empirical CDF (an observed order statistic).
  # If the threshold is missed, both endpoints stay NA and the outer CI fails.
  good <- is.finite(estimates) & estimates>0
  ci <- c(NA_real_,NA_real_)
  if(sum(good)>=ceiling(B*bootstrap_required_fraction)) {
    ci <- unname(quantile(estimates[good],
      c((1-confidence)/2,(1+confidence)/2),type=1))
  }
  list(ci=ci,valid=sum(good),failed=sum(!good))
}

# 7. PRE-SIMULATION CONSISTENCY CHECKS
# These deterministic checks verify formulas and implementations, not empirical
# coverage. stopifnot halts the run if an identity fails its numerical tolerance.
# Executable consistency checks run before a long study.
run_checks <- function() {
  # Use a small fixed censored dataset to check that both LR endpoints attain
  # the intended deviance cutoff when evaluated in the original likelihood.
  theta <- 1200; e <- c(100,500,3000,H,H)
  m <- sum(e<H); T <- sum(e); estimate <- T/m
  ll <- function(th) -m*log(th)-T/th
  ratio <- lr_ratios(m)
  stopifnot(abs(2*(ll(estimate)-ll(estimate*ratio[1]))-cutoff)<1e-6,
            abs(2*(ll(estimate)-ll(estimate*ratio[2]))-cutoff)<1e-6)
  # A symmetric difference in log(theta) checks dR/dlog(theta)=A*theta.
  for(j in seq_along(A)) {
    step <- 1e-5
    numerical <- (risk(theta*exp(step))[j]-risk(theta*exp(-step))[j])/(2*step)
    stopifnot(abs(numerical-A[j]*theta)/(1+abs(numerical))<1e-6)
  }
  # Check scalar/matrix sufficient statistics on identical observations.
  mat <- matrix(c(100,H,1000,500,200,H),nrow=3)
  stopifnot(all(colSums(mat<H)==c(2,2)),
            all(colSums(pmin(mat,H))==c(H+1100,H+700)))
  # Check explicit table entries, including signs that disappear on squaring.
  for(th in c(500,1000,3000)) {
    tab <- delta_table_components(th)
    pp <- -expm1(-H/th)
    expected <- th^2/pp*c(VaR=log(beta)^2,CTE=(log(beta)-1)^2,
      GS=(log(beta)-1-delta)^2,PHT=r^(-2),WT=wt_factor^2)
    stopifnot(all(abs(tab$coefficient-expected)/(1+abs(expected))<1e-12))
  }
  # Independently verify C_WT using its survival-integral definition.
  independent_wt <- integrate(function(v)
    exp(pnorm(qnorm(-v,log.p=TRUE)+lambda,log.p=TRUE)),0,Inf,
    rel.tol=1e-8,subdivisions=1000)$value
  stopifnot(abs(wt_factor-independent_wt)/(1+abs(independent_wt))<1e-7)
  cat('LR, gradients, sufficient statistics and expected-information table checks passed.\n')
}
if(RUN_CHECKS) run_checks()

# 8. PER-METHOD, PER-RISK ACCUMULATORS
# There is one row per risk/method combination. Valid includes finite and
# infinite computed lengths. Covered counts coverage within that valid set;
# Finite_Covered counts coverage only within the finite-length subset.
# Undefined, reversed and indeterminate lengths are mutually exclusive failures.
# Outer fit failures are counted separately because no method is attempted.
new_accumulator <- function() {
  a <- expand.grid(Risk_Measure=measures,Method=methods,stringsAsFactors=FALSE)
  a$Valid <- a$Covered <- a$Finite_Length_Count <- a$Infinite_Length_Count <- 0L
  a$Finite_Covered <- a$Undefined_Count <- a$Reversed_Count <- a$Indeterminate_Length_Count <- 0L
  # List columns retain exact finite lengths for medians and tail quantiles.
  # WT diagnostics retain only the ten largest finite intervals per row.
  # These lists stay in checkpoints but are removed before scalar CSV export.
  a$Finite_Lengths <- I(rep(list(numeric(0)),nrow(a)))
  a$WT_Diagnostics <- I(rep(list(data.frame()),nrow(a)))
  a$Sum_Finite_Length <- 0
  a
}
# Classify one method's interval vector and update its risk-specific rows.
# Optional fit context identifies unusually long WT intervals for later review.
# Checks are ordered: missing endpoints, reversed bounds, undefined difference,
# then finite versus infinite length. Infinity alone does not imply failure.
update_accumulator <- function(a,method,lower,upper,replication=NA_integer_,
                               theta_hat=NA_real_,m=NA_integer_,T=NA_real_) {
  for(j in seq_along(measures)) {
    i <- which(a$Method==method & a$Risk_Measure==measures[j])
    ci <- c(lower[j],upper[j])
    # anyNA catches both NA and NaN. Equal infinite endpoints pass the ordering
    # check but yield NaN length and are classified as indeterminate below.
    if(anyNA(ci)) {a$Undefined_Count[i] <- a$Undefined_Count[i]+1L; next}
    if(ci[1]>ci[2]) {a$Reversed_Count[i] <- a$Reversed_Count[i]+1L; next}
    length <- ci[2]-ci[1]
    if(is.na(length)) {a$Indeterminate_Length_Count[i] <- a$Indeterminate_Length_Count[i]+1L; next}
    # Endpoint inclusion is closed: a truth exactly on either bound is covered.
    # This valid-set counter includes an infinite interval if it contains the truth.
    a$Valid[i] <- a$Valid[i]+1L
    a$Covered[i] <- a$Covered[i]+as.integer(ci[1]<=true_risks[j] && true_risks[j]<=ci[2])
    # Only finite lengths enter precision summaries and finite-valid coverage.
    # An infinite difference can also arise from finite-endpoint numeric overflow.
    if(is.finite(length)) {
      a$Finite_Length_Count[i] <- a$Finite_Length_Count[i]+1L
      a$Sum_Finite_Length[i] <- a$Sum_Finite_Length[i]+length
      a$Finite_Covered[i] <- a$Finite_Covered[i]+as.integer(ci[1]<=true_risks[j] && true_risks[j]<=ci[2])
      a$Finite_Lengths[[i]] <- c(a$Finite_Lengths[[i]],length)
      # Keep the largest ten finite WT lengths in descending order after each
      # update. Replication plus n/method identifies the originating outer dataset;
      # Theta_Hat, Exact_Claims and Exposure help diagnose unstable fitted scales.
      if(measures[j]=='WT') {
        rows <- rbind(a$WT_Diagnostics[[i]],data.frame(Replication=replication,
          Lower=ci[1],Upper=ci[2],Length=length,Theta_Hat=theta_hat,
          Exact_Claims=m,Exposure=T))
        a$WT_Diagnostics[[i]] <- head(rows[order(rows$Length,decreasing=TRUE),],10L)
      }
    } else a$Infinite_Length_Count[i] <- a$Infinite_Length_Count[i]+1L
  }
  a
}
# 9. REPRODUCIBILITY RECORD AND CHECKPOINT WRITES
# Record parameters and, when discoverable, the exact script text. Because
# source text is compared on resume, even comment edits invalidate an existing
# checkpoint configuration: use a fresh output folder for an edited script.
source_text <- if(length(script)==1 && file.exists(script)) readLines(script,warn=FALSE) else NULL
configuration <- list(version=version,MODE=MODE,N_MC=N_MC,B=B,sample_sizes=sample_sizes,
  x0=x0,theta_true=theta_true,d=d,censor_at=censor_at,confidence=confidence,
  beta=beta,delta=delta,r=r,lambda=lambda,bootstrap_block=bootstrap_block,
  bootstrap_required_fraction=bootstrap_required_fraction,seed=seed,source=source_text,
  delta_information=delta_information)
saveRDS(configuration,file.path(outdir,'Simulation_Configuration.rds'))

# Write a temporary RDS before rotating the old checkpoint to a backup.
# If replacement fails, attempt to restore the previous checkpoint and stop.
# On restart, the backup is used only when the main checkpoint is absent.
save_checkpoint <- function(state,path) {
  tmp <- paste0(path,'.tmp'); backup <- paste0(path,'.bak')
  saveRDS(state,tmp)
  if(file.exists(backup)) unlink(backup)
  if(file.exists(path) && !file.rename(path,backup)) stop('Cannot back up checkpoint')
  if(!file.rename(tmp,path)) {
    if(file.exists(backup)) file.rename(backup,path)
    stop('Cannot replace checkpoint')
  }
  if(file.exists(backup)) unlink(backup)
}
# 10. MONTE CARLO LOOP FOR ONE SAMPLE SIZE
# cell_id determines the seed offset; sample-size ordering therefore matters.
# State contains completed counts, fit/bootstrap diagnostics and accumulators.
# Returning scalar rows lets the caller combine all sample-size results.
run_cell <- function(n,cell_id) {
  path <- file.path(outdir,sprintf('checkpoint_Exp_n%d.rds',n))
  backup <- paste0(path,'.bak')
  state <- list(done=0L,fit_failures=0L,censored_total=0,
                bootstrap_attempted=0,bootstrap_valid=0,bootstrap_failed=0,
                acc=new_accumulator(),configuration=configuration)
  # Require exact configuration agreement to prevent combining incompatible
  # runs. No saved RNG state is needed because each replication resets its seed.
  if(resume && (file.exists(path) || file.exists(backup))) {
    state <- readRDS(if(file.exists(path)) path else backup)
    if(!identical(state$configuration,configuration))
      stop('Checkpoint configuration differs; use a new output folder')
  }
  # Continue from the first unfinished replication, or skip an already complete
  # cell. All methods use the same outer dataset for paired comparison.
  if(state$done<N_MC) for(replication in seq.int(state$done+1L,N_MC)) {
    # Separate seed streams by sample-size cell and replication. Resuming the
    # same configuration regenerates each remaining dataset/bootstrap sequence.
    set.seed(seed+cell_id*1000000L+replication)
    s <- simulate_statistics(theta_true,n)
    state$censored_total <- state$censored_total+s$k
    # If no finite fitted scale can be constructed, count one outer fit failure.
    # No method's interval counter is updated; final Failed includes this case.
    if(s$m==0L || !is.finite(s$T) || s$T<=0) {
      state$fit_failures <- state$fit_failures+1L
    } else {
      # Closed-form MLE from censored exposure; map it to all five point estimates.
      theta_hat <- s$T/s$m
      point <- risk(theta_hat)
      # Expected Fisher information, matching the supplied table.
      # This replaces the previous observed-information denominator sqrt(m).
      tab <- delta_table_components(theta_hat)
      # Convert the asymptotic variance coefficient into an SE by dividing by n.
      # Raw-scale Wald limits may fall below x0; they are deliberately not clipped.
      se <- sqrt(tab$coefficient/n)
      state$acc <- update_accumulator(state$acc,'Delta',point-zcrit*se,point+zcrit*se,replication,theta_hat,s$m,s$T)
      # Map the LR scale endpoints through increasing risk functions.
      ratios <- lr_ratios(s$m)
      state$acc <- update_accumulator(state$acc,'Profile_likelihood',
        risk(theta_hat*ratios[1]),risk(theta_hat*ratios[2]),replication,theta_hat,s$m,s$T)
      # Bootstrap diagnostics count attempted individual bootstrap fits, not risk
      # rows. A single successful scale CI supplies all five risk intervals.
      boot <- bootstrap_theta(theta_hat,n)
      state$bootstrap_attempted <- state$bootstrap_attempted+B
      state$bootstrap_valid <- state$bootstrap_valid+boot$valid
      state$bootstrap_failed <- state$bootstrap_failed+boot$failed
      # Risks are increasing affine functions of theta, so empirical type=1
      # percentile endpoints map exactly from bootstrap theta endpoints.
      state$acc <- update_accumulator(state$acc,'Bootstrap_percentile',
                                     risk(boot$ci[1]),risk(boot$ci[2]),replication,theta_hat,s$m,s$T)
    }
    # Mark the outer dataset complete even if it failed, so all-replication
    # percentages retain failures in their denominator. Save progress periodically.
    state$done <- replication
    if(replication %% checkpoint_every==0 || replication==N_MC) {
      save_checkpoint(state,path)
      cat('Shifted exponential n=',n,': ',replication,'/',N_MC,'\n',sep='')
      flush.console()
    }
  }
  # 11. CELL SUMMARIES AND DIAGNOSTIC EXPORTS
  # All length summaries condition on finite valid intervals unless explicitly
  # named otherwise. No usable finite lengths gives NA, never zero.
  a <- state$acc
  # Apply a scalar statistic to each retained length vector. type=7 quantiles
  # use R's interpolated default, distinct from bootstrap endpoint type=1.
  finite_stat <- function(f) vapply(a$Finite_Lengths,function(x)
    if(length(x)) unname(f(x)) else NA_real_,numeric(1))
  a$Median_Finite_Length <- finite_stat(median)
  a$P95_Finite_Length <- finite_stat(function(x) quantile(x,.95,type=7))
  a$P99_Finite_Length <- finite_stat(function(x) quantile(x,.99,type=7))
  a$Maximum_Finite_Length <- finite_stat(max)
  # Export at most ten finite WT extremes for each method and sample size.
  # No file is written for a method with no finite WT intervals.
  for(i in which(a$Risk_Measure=='WT')) {
    rows <- a$WT_Diagnostics[[i]]
    if(nrow(rows)) {
      rows$n <- n; rows$Method <- a$Method[i]
      write.csv(rows,file.path(outdir,sprintf('Largest_WT_Intervals_n%d_%s.csv',n,a$Method[i])),row.names=FALSE)
    }
  }
  # Remove nested raw-data columns from the summary copy only; state$acc still
  # retains them in the checkpoint for reproducibility and exact quantiles.
  a$Finite_Lengths <- NULL; a$WT_Diagnostics <- NULL
  a$Model <- 'Shifted_Exponential'; a$n <- n; a$Completed <- state$done
  a$Failed <- state$done-a$Valid; a$Fit_Failures <- state$fit_failures
  # Three coverage denominators answer different questions:
  # (1) all completed replications, with construction failures counted as misses;
  # (2) all valid intervals, including infinite lengths;
  # (3) finite valid intervals only, used as the primary conditional measure.
  # Conditional coverage can be selective; always read it with failure percentages.
  a$Coverage_Success_Over_All <- a$Covered/state$done
  a$Coverage_Among_All_Valid <- ifelse(a$Valid>0,a$Covered/a$Valid,NA_real_)
  a$Coverage_Among_Finite_Valid <- ifelse(a$Finite_Length_Count>0,
    a$Finite_Covered/a$Finite_Length_Count,NA_real_)
  # Binomial plug-in Monte Carlo SE conditional on the observed finite count.
  # This measures simulation uncertainty in coverage, not a risk-estimate SE.
  # It is NA when that count is zero; plug-in SE can be zero at coverage 0 or 1.
  a$MCSE_Coverage_Among_Finite_Valid <- sqrt(a$Coverage_Among_Finite_Valid*
    (1-a$Coverage_Among_Finite_Valid)/a$Finite_Length_Count)
  # Report percentages on a 0-100 scale using all completed outer datasets.
  # Coverage proportions above remain on the 0-1 scale.
  a$Percent_Infinite_Intervals <- 100*a$Infinite_Length_Count/state$done
  a$Percent_Undefined_Or_Numerical_Failures <- 100*a$Failed/state$done
  # Accounting identities ensure every completed dataset is classified once
  # for each risk/method and finite/infinite valid categories reconcile.
  stopifnot(all(a$Valid==a$Finite_Length_Count+a$Infinite_Length_Count),
    all(a$Failed==state$fit_failures+a$Undefined_Count+a$Reversed_Count+a$Indeterminate_Length_Count))
  # The extended-real valid-set mean is Inf if any valid length is infinite.
  # The finite-only mean omits those intervals and can differ greatly from the
  # median under skewness. The all-dataset mean is NA if any CI is undefined.
  a$Mean_Length_Valid <- ifelse(a$Valid==0,NA_real_,
    ifelse(a$Infinite_Length_Count>0,Inf,a$Sum_Finite_Length/a$Valid))
  a$Mean_Finite_Length <- ifelse(a$Finite_Length_Count>0,
    a$Sum_Finite_Length/a$Finite_Length_Count,NA_real_)
  a$Mean_Length_All <- ifelse(a$Failed==0,a$Mean_Length_Valid,NA_real_)
  a$Mean_Length_Valid_Thousand <- a$Mean_Length_Valid/1000
  # Average censoring over all generated claims, including outer fit failures.
  # Per-cell summaries include all scalar counters for an exclusion audit.
  a$Observed_Censoring_Fraction <- state$censored_total/(state$done*n)
  write.csv(a,file.path(outdir,sprintf('Summary_Exp_n%d.csv',n)),row.names=FALSE,na='NA')
  # Bootstrap totals cover only outer datasets with successful original fits;
  # Bootstrap_Failed counts resample fits, not failed outer bootstrap intervals.
  diagnostics <- data.frame(n=n,Completed=state$done,Outer_Fit_Failures=state$fit_failures,
    Bootstrap_Attempted=state$bootstrap_attempted,Bootstrap_Valid=state$bootstrap_valid,
    Bootstrap_Failed=state$bootstrap_failed)
  write.csv(diagnostics,file.path(outdir,sprintf('Diagnostics_Exp_n%d.csv',n)),row.names=FALSE)
  a
}

# 12. COMBINE SAMPLE SIZES AND WRITE PUBLICATION TABLES
# After each cell, refresh the combined detailed summary so completed work
# remains available even if a later sample-size run is interrupted.
summaries <- vector('list',length(sample_sizes))
for(i in seq_along(sample_sizes)) {
  summaries[[i]] <- run_cell(sample_sizes[i],i)
  write.csv(do.call(rbind,summaries[seq_len(i)]),
    file.path(outdir,'Monte_Carlo_Summary.csv'),row.names=FALSE,na='NA')
}
summary <- do.call(rbind,summaries)
# Primary table: conditional finite-valid coverage with MCSE, typical finite
# precision (median), and frequencies of infinite or undefined intervals.
# The finite count makes the coverage denominator directly inspectable.
main_columns <- c('Model','n','Risk_Measure','Method','Completed','Finite_Length_Count',
  'Coverage_Among_Finite_Valid','MCSE_Coverage_Among_Finite_Valid',
  'Median_Finite_Length','Percent_Infinite_Intervals','Percent_Undefined_Or_Numerical_Failures')
write.csv(summary[,main_columns],file.path(outdir,'Main_Results.csv'),row.names=FALSE,na='NA')
# Supplementary table: average and extreme finite precision plus detailed
# exclusion counts and alternative coverage denominators for interpretation.
write.csv(summary[,c('Model','n','Risk_Measure','Method','Mean_Finite_Length',
  'P95_Finite_Length','P99_Finite_Length','Maximum_Finite_Length','Coverage_Among_All_Valid',
  'Coverage_Success_Over_All','Fit_Failures','Undefined_Count','Reversed_Count',
  'Indeterminate_Length_Count','Finite_Length_Count','Infinite_Length_Count','Failed')],
  file.path(outdir,'Supplementary_Results.csv'),row.names=FALSE,na='NA')
# 13. RESULTS GUIDE AND CONSOLE DISPLAY
# Write interpretation notes alongside outputs so denominator conventions,
# quantile definitions and limitations travel with the exported tables.
writeLines(c(
  'Shifted exponential only. theta is SCALE, x0 fixed.',
  sprintf('N_MC=%d per n; B=%d; confidence=%.2f; r=%.2f; lambda=%.2f.',N_MC,B,confidence,r,lambda),
  'n counts conditional observations above d, including right-censored records.',
  'Risks are ground-up; truncation and censoring enter the likelihood.',
  'MLE=T/m; T includes exposure H for every censored record.',
  'Delta intervals use expected Fisher information I_theta=p(theta)/theta^2 per observation.',
  'p(theta)=1-exp(-(censor_at-d)/theta), evaluated at theta_hat.',
  'Delta variance=theta_hat^2*A^2/(n*p(theta_hat)); it is not divided by m.',
  'Delta intervals remain raw-scale Wald intervals; lower endpoints are not clipped.',
  'True_Risks_And_Delta_Table.csv exports the variance coefficients and standard errors.',
  'Bootstrap intervals use empirical type=1 percentile quantiles; no failed draws are replaced.',
  'Profile intervals use asymptotic chi-square(1) cutoffs.',
  'Mean_Length_Valid_Thousand=Mean_Length_Valid/1000.',
  'Main_Results.csv reports coverage among finite valid intervals, its binomial Monte Carlo SE, median finite length, and failure/infinite percentages.',
  'Median finite length describes typical precision; mean finite length describes average precision and can be dominated by rare enormous intervals.',
  'Coverage_Among_Finite_Valid excludes infinite intervals and undefined/numerically failed intervals; denominator=Finite_Length_Count.',
  'MCSE uses sqrt(p*(1-p)/Finite_Length_Count); conditional on the observed finite valid count.',
  'Coverage_Among_All_Valid includes ordered infinite-length intervals; denominator=Valid.',
  'Infinite and failure percentages use all completed Monte Carlo replications as denominator.',
  'Failed includes outer fit failures, missing/NaN endpoints, reversed endpoints, and indeterminate lengths.',
  'Infinite_Length_Count includes ordered infinite endpoints and finite-endpoint subtraction overflow.',
  'Supplementary_Results.csv reports mean, P95/P99 (type=7), maximum finite length, and detailed exclusion counts.',
  'Largest_WT_Intervals files retain up to 10 largest finite WT intervals per n/method with replication, endpoints and fit statistics.',
  'Exact finite lengths are retained in checkpoints for median and quantile calculation; v5 uses a new checkpoint folder.',
  'No finite valid intervals implies NA for finite-length summaries and conditional coverage.',
  'For this shifted exponential model all risks are affine in theta; WT need not exhibit exceptional skewness.', 
  'Coverage_Success_Over_All counts failed intervals as noncoverage.',
  'Mean_Length_All is undefined (NA) if any outer interval failed.',
  'Diagnostic bootstrap totals count fits, not risk-measure rows.',
  'No numerical optimizer is needed. No uncertainty in fixed x0 is included.'
),file.path(outdir,'README_Results.txt'))
# Display the primary results without changing full-precision CSV values.
# All unscaled interval lengths use the same monetary units as x0 and theta;
# only columns explicitly named Thousand are divided by 1000.
print(summary[,main_columns],row.names=FALSE)
cat('Results saved in:',outdir,'\n')
