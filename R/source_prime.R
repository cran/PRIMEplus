
get_iargs <- function(op) {

  len <- 11
  ret <- c(op$DEBUG, op$N, op$J, op$num_rep_EM, op$num_rand, op$num_rep,
           op$n_failure, op$FUNCTION, op$n_trt, op$print, op$betaMat_nr)
  if (length(ret) != len) stop("ERROR") 
  ret

}

get_dargs <- function(op) {

  len <- 2
  ret <- c(op$tstar, op$eps_EM)
  if (length(ret) != len) stop("ERROR") 
  ret

}

getHazard <- function(time, treatment, event_status, t.fail.o=NULL) {

  control         <- coxph.control()
  control$timefix <- FALSE
  if (!length(t.fail.o)) t.fail.o <- get_tfailo(time, event_status)$t.fail.o
  n.failure       <- length(t.fail.o)
  t.all           <- c(0, t.fail.o)
  t.diff          <- t.all[-1]-t.all[-length(t.all)]
  fit             <- coxph(Surv(time, event_status) ~ treatment, control=control)
  ss              <- survfit(fit)
  cum.hazard      <- unique(-log(ss$surv))

  if (length(cum.hazard) != n.failure) {
    tmp <- cum.hazard != 0
    cum.hazard <- cum.hazard[tmp]
    if (length(cum.hazard) != n.failure) stop("ERROR computing hazard")
  }

  hazard.dis   = c(cum.hazard[1], cum.hazard[-1] - cum.hazard[-length(cum.hazard)])
  hazard       = hazard.dis / t.diff

  list(hazard=hazard, t.all=t.all)

} # END: getHazard

get_tfailo <- function(X, event_status) {

  xvec          <- X[event_status == 1]
  t.failure     <- unique(xvec)
  tmp           <- order(t.failure)
  t.fail.o      <- t.failure[tmp]
  delta.l       <- as.numeric(table(xvec))

  list(t.fail.o=t.fail.o, delta.l=delta.l)

} # END: get_tfailo

#####################################################
# !!! Data must be ordered before calling EM.main !!!
#####################################################
call_EM_C <- function(X, trt, event, t_all, 
               lambda0, delta_l, beta0, P, op) {

  iargs         <- get_iargs(op)
  dargs         <- get_dargs(op)
  N             <- op$N
  J             <- op$J
  nfail         <- op$n_failure
  ret_conv      <- as.integer(-9999)
  ret_beta      <- as.numeric(rep(-9999, J))
  ret_lambda    <- as.numeric(rep(-9999, nfail))
  ret_W         <- as.numeric(rep(-9999, N*J))
  ret_loglik    <- as.numeric(rep(-9999, 3))

  tmp <- .C("C_EM", as.integer(iargs), as.numeric(dargs), as.numeric(X), as.integer(event),
            as.integer(trt), as.numeric(t_all), as.numeric(lambda0), as.integer(delta_l), 
            as.numeric(beta0), as.numeric(P),
            ret_conv=as.integer(ret_conv), ret_lambda=as.numeric(ret_lambda), 
            ret_beta=as.numeric(ret_beta), ret_W=as.numeric(ret_W), ret_loglik=as.numeric(ret_loglik))

  ret_conv  <- tmp$ret_conv == 1
  if (!ret_conv) {
    ret_loglik[] <- NA
  } else {
    ret_loglik   <- tmp$ret_loglik
  }

  W              <- matrix(tmp$ret_W, nrow=N, ncol=J, byrow=TRUE)
  result         <- list(converged=ret_conv, beta=tmp$ret_beta[1:(J-1)], probResponder=W[, -J, drop=FALSE],
                         lambda=tmp$ret_lambda, loglike=ret_loglik[1], 
                         loglike.marg=ret_loglik[2])

  result

} 

getOp <- function(trt, lambda0, effect_p, t1, stopTol, maxiter, print, betaMat_nr, num_rand=1000) {

  nfail <- length(lambda0)
  ntrt  <- sum(trt %in% 1)
  ret   <- list(DEBUG=0, N=length(trt), J=length(effect_p), num_rep_EM=maxiter,
                num_rand=num_rand, num_rep=1000,
                n_failure=nfail, FUNCTION=0, n_trt=ntrt, print=print,
                eps_EM=stopTol, tstar=t1, betaMat_nr=betaMat_nr)
  ret

}

getProbResponderDF <- function(W, data, time.var, trt.var, status.var, id.var, ord) {

  nr <- nrow(data)
  if (nr != nrow(W)) stop("ERROR")
  if (length(id.var)) {
    ids <- data[, id.var, drop=TRUE]
  } else {
    id.var <- "_ID_" 
    ids    <- 1:nr
  }
  vars  <- c(time.var, trt.var, status.var)
  Wvars <- paste0("ProbGroup_", 1:ncol(W)) 
 
  # data is the original order, W is ordered with treated subjects first
  ret   <- data.frame(ids, data[, vars, drop=FALSE], W[ord, , drop=FALSE], stringsAsFactors=FALSE)
  colnames(ret) <- c(id.var, vars, Wvars)
  if (length(ord) != nr) stop("ERROR 2")
 
  ret
}

getNinitBeta <- function(beta0) {

  ret <- nrow(beta0)
  if (is.null(ret)) ret <- 0
  ret

}

PRIMEplus.EM <- function(data, effect_p, beta0, time.var="X", trt.var="trt", status.var="event_status", id.var="id", 
                         t1=1, lambda0=NULL, stopTol=1e-4, maxiter=100000, print=0) {

  n        <- length(effect_p)
  effect_p <- check_effect_p(effect_p)
  beta0    <- check_beta(beta0, n)
  check_t1(t1)
  check_stopTol(stopTol)
  check_maxiter(maxiter)
  check_print(print)

  tmp      <- setUpAndCheckData(data, time.var, trt.var, status.var, lambda0, effect_p, beta0)
  X        <- tmp$time
  trt      <- tmp$treatment
  event    <- tmp$event_status
  lambda0  <- tmp$lambda0
  t_all    <- c(0, tmp$t.fail)
  delta_l  <- tmp$delta.l
  ord      <- tmp$order
  check_lambda0(lambda0, event, X)
  betaMat_nr <- getNinitBeta(beta0)

  op  <- getOp(trt, lambda0, effect_p, t1, stopTol, maxiter, print, betaMat_nr)
  tmp <- call_EM_C(X, trt, event, t_all, lambda0, delta_l, beta0, effect_p, op)
  tmp$probResponder <- getProbResponderDF(tmp$probResponder, data, time.var, trt.var, status.var, id.var, ord)
  
  tmp

}

call_LRT_C <- function(X, trt, event, t_all, 
               lambda0, delta_l, beta0, P, op) {

  iargs         <- get_iargs(op)
  dargs         <- get_dargs(op)
  N             <- op$N
  J             <- op$J
  nfail         <- op$n_failure
  
  ret_conv   <- -9999
  ret_beta   <- rep(-9999, J)
  ret_lambda <- rep(-9999, nfail)
  ret_W      <- rep(-9999, N*J)
  ret_loglik <- rep(-9999, 3)
  ret_LRT    <- -9999

  tmp <- .C("C_LRT", as.integer(iargs), as.numeric(dargs), as.numeric(X), as.integer(event),
            as.integer(trt), as.numeric(t_all), as.numeric(lambda0), as.integer(delta_l), 
            as.numeric(beta0), as.numeric(P),
            ret_conv=as.integer(ret_conv), ret_lambda=as.numeric(ret_lambda), 
            ret_beta=as.numeric(ret_beta), ret_W=as.numeric(ret_W), ret_loglik=as.numeric(ret_loglik),
            ret_LRT=as.numeric(ret_LRT))

  ret_conv  <- tmp$ret_conv == 1
  if (!ret_conv) {
    ret_loglik[] <- NA
    ret_LRT      <- NA
    pvalue       <- NA
  } else {
    ret_loglik   <- tmp$ret_loglik
    ret_LRT      <- tmp$ret_LRT
    pvalue       <- pchisq(ret_LRT, df=J-1, lower.tail=FALSE)
  }

  W              <- matrix(tmp$ret_W, nrow=N, ncol=J, byrow=TRUE)
  result         <- list(converged=ret_conv, beta=tmp$ret_beta[1:(J-1)], probResponder=W[, -J, drop=FALSE],
                         lambda=tmp$ret_lambda, loglike=ret_loglik[1], 
                         loglike.marg=ret_loglik[2], loglike.marg.0=ret_loglik[3],
                         LRT=ret_LRT, pvalue=pvalue)

  result

} 

PRIMEplus.LRT <- function(data, effect_p, beta0, time.var="X", trt.var="trt", status.var="event_status", id.var="id", 
                          t1=1, lambda0=NULL, stopTol=1e-4, maxiter=100000, print=0) {

  n        <- length(effect_p)
  effect_p <- check_effect_p(effect_p)
  beta0    <- check_beta(beta0, n) 
  check_t1(t1)
  check_stopTol(stopTol)
  check_maxiter(maxiter)
  check_print(print)

  tmp      <- setUpAndCheckData(data, time.var, trt.var, status.var, lambda0, effect_p, beta0)
  X        <- tmp$time
  trt      <- tmp$treatment
  event    <- tmp$event_status
  lambda0  <- tmp$lambda0
  t_all    <- c(0, tmp$t.fail)
  delta_l  <- tmp$delta.l
  ord      <- tmp$order
  check_lambda0(lambda0, event, X)
  betaMat_nr <- getNinitBeta(beta0)

  op  <- getOp(trt, lambda0, effect_p, t1, stopTol, maxiter, print, betaMat_nr)
  tmp <- call_LRT_C(X, trt, event, t_all, lambda0, delta_l, beta0, effect_p, op)
  tmp$probResponder <- getProbResponderDF(tmp$probResponder, data, time.var, trt.var, status.var, id.var, ord)
  
  tmp

}

call_ReRandLRT_C <- function(X, trt, event, t_all, 
               lambda0, delta_l, beta0, P, op) {

  iargs          <- get_iargs(op)
  dargs          <- get_dargs(op)
  N              <- op$num_rand
  ret_sumNrand   <- -9999
  ret_p_LRT      <- -9999
  ret_p_loglikm  <- -9999
  ret_LRT_vec    <- rep(-9999, N)
  tmp <- .C("C_ReRandLRT", as.integer(iargs), as.numeric(dargs), as.numeric(X), as.integer(event),
            as.integer(trt), as.numeric(t_all), as.numeric(lambda0), as.integer(delta_l), 
            as.numeric(beta0), as.numeric(P),
            ret_sumNrand=as.integer(ret_sumNrand), ret_p_LRT=as.numeric(ret_p_LRT), 
            ret_p_loglikm=as.numeric(ret_p_loglikm), ret_LRT_vec=as.numeric(ret_LRT_vec))

  M <- tmp$ret_sumNrand
  if (M < N) {
    msg <- paste0("Randomization results are based on ", M, " randomizations")
    warning(msg)
  }
  if (M < 1) {
    ret_sumNrand  <- 0
    ret_p_LRT     <- NA
    ret_p_loglikm <- NA
  } else {
    ret_sumNrand  <- M
    ret_p_LRT     <- tmp$ret_p_LRT
    ret_p_loglikm <- tmp$ret_p_loglikm
  }

  result <- list(pvalue.LRT=ret_p_LRT, pvalue.loglike.marg=ret_p_loglikm, n.randomizations=M)

  result

} 

PRIMEplus.ReRand.LRT <- function(data, effect_p, beta0, time.var="X", trt.var="trt", status.var="event_status", id.var="id", 
                                 t1=1, lambda0=NULL, stopTol=1e-4, maxiter=100000, print=0, num_rand=1000) {

  n        <- length(effect_p)
  effect_p <- check_effect_p(effect_p)
  beta0    <- check_beta(beta0, n) 
  check_t1(t1)
  check_stopTol(stopTol)
  check_maxiter(maxiter)
  check_print(print)
  check_num_rand(num_rand)

  tmp      <- setUpAndCheckData(data, time.var, trt.var, status.var, lambda0, effect_p, beta0)
  X        <- tmp$time
  trt      <- tmp$treatment
  event    <- tmp$event_status
  lambda0  <- tmp$lambda0
  t_all    <- c(0, tmp$t.fail)
  delta_l  <- tmp$delta.l
  check_lambda0(lambda0, event, X)
  betaMat_nr <- getNinitBeta(beta0)

  op  <- getOp(trt, lambda0, effect_p, t1, stopTol, maxiter, print, betaMat_nr, num_rand=num_rand)
  tmp <- call_ReRandLRT_C(X, trt, event, t_all, lambda0, delta_l, beta0, effect_p, op)
  
  tmp

}

PRIMEplus.Power_main <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=1000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=1000, print=0, powerL=-1, powerU=-1,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  # effect_p and HR must be appended when this function is called, note they are assumed to not be
  #   appended in generate_data
  neff      <- length(effect_p)
  effect_p0 <- effect_p[1:neff]
  HR0       <- HR
  beta0     <- log(HR)
  p.val.all <- rep(NA, nsim) 
  nok       <- 0
  pflag     <- (powerL > 0) && (powerU > 0)
  oper      <- "="
  prt2      <- 0
  for (i in 1:nsim)
  {  
    if (print) cat(paste("Simulation ", i, "\n", sep=""))
    # 1-sum(effect_p) and 1 are added to effect_p and HR respectively in generate_data
    dat      <- generate_data(N=nmax, rand_ratio=rand_ratio, effect_p=effect_p0, 
                              enroll_rate=enroll_rate, lambda1=lambda1, HR=HR0, tau=tau, t1=t1)
    tmp      <- setUpAndCheckData(dat, "X", "trt", "event_status", NULL, NULL, NULL)
    X        <- tmp$time
    trt      <- tmp$treatment
    event    <- tmp$event_status
    lambda0  <- tmp$lambda0
    t_all    <- c(0, tmp$t.fail)
    delta_l  <- tmp$delta.l
    ord      <- tmp$order
    check_lambda0(lambda0, event, X)
    betaMat_nr <- 0 #getNinitBeta(beta0)
    tmp        <- checkTrtStatus(trt, event, min.sample.size=min.sample.size, 
                                 min.n.event=min.n.event, min.per.trt=min.per.trt)

    # Appended versions of effect_p, beta0, HR
    op  <- getOp(trt, lambda0, effect_p, t1, stopTol, maxiter, prt2, betaMat_nr, num_rand=num_rand)
    ret <- try(call_ReRandLRT_C(X, trt, event, t_all, lambda0, delta_l, beta0, effect_p, op), silent=FALSE)
    if (!("try-error" %in% class(ret))) {
      p.val.all[i] <- ret$pvalue.LRT
      nok          <- nok + 1
    }

    # compute min and max power
    if (pflag) {
      M        <- nsim - i
      tmp      <- sum(p.val.all <= alpha, na.rm=TRUE) 
      denom    <- M + nok
      minPower <- tmp/denom
      maxPower <- (tmp + M)/denom
      if (minPower > powerL) {
        if (M) oper <- ">"
        break
      } else if (maxPower < powerU) {
        if (M) oper <- "<"
        break
      } 
    }  
  }
  m     <- sum(is.na(p.val.all))
  n     <- nsim - m
  if (m && !pflag) warning(paste("power based only on ", n, " simulated datasets", sep=""))
  if (!n) stop("ERROR: power could not be estimated")
  if (!pflag) minPower <- mean(as.numeric(p.val.all <= alpha), na.rm=TRUE) 
  
  list(power=minPower, n=n, oper=oper)

} # END: PRIMEplus.Power_main


PRIMEplus.Power <- function(nmax=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=10000, print=0,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  check_nmax(nmax)
  check_rand_ratio(rand_ratio)
  check_enroll_rate(enroll_rate)
  check_lambda1(lambda1)
  n        <- length(effect_p)
  effect_p <- check_effect_p(effect_p)
  HR       <- check_HR(HR, n)
  check_tau(tau)
  check_t1(t1)
  check_stopTol(stopTol)
  check_maxiter(maxiter)
  check_print(print)
  check_alpha(alpha)
  check_num_rand(num_rand)
  check_nsim(nsim)
  check_min.sample.size(min.sample.size)
  check_min.n.event(min.n.event)
  check_min.per.trt(min.per.trt)

  tmp <- PRIMEplus.Power_main(nmax=nmax, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=print,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
  
  ret <- list(power=tmp$power, n.simulations=tmp$n)
  ret

} # END: PRIMEplus.Power

PRIMEplus.SampleSize <- function(power=0.8, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
                       lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
                       alpha=0.05, num_rand=1000, nsim=10000, min.N=100, max.N=700, 
                       tol.power=0.01, tol.N=1, print=1,
                       min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  check_power(power)
  check_rand_ratio(rand_ratio)
  check_enroll_rate(enroll_rate)
  check_lambda1(lambda1)
  n        <- length(effect_p)
  effect_p <- check_effect_p(effect_p)
  HR       <- check_HR(HR, n)
  check_tau(tau)
  check_t1(t1)
  check_stopTol(stopTol)
  check_maxiter(maxiter)
  check_print(print)
  check_alpha(alpha)
  check_num_rand(num_rand)
  check_nsim(nsim)
  check_min.sample.size(min.sample.size)
  check_min.n.event(min.n.event)
  check_min.per.trt(min.per.trt)
  check_min.N(min.N)
  check_max.N(max.N, min.N)
  check_tol.power(tol.power)
  check_tol.N(tol.N)

  PWRL <- power + tol.power + 1e-6
  PWRU <- power - tol.power - 1e-6

  # Test left endpoint
  tmp <- PRIMEplus.Power_main(nmax=min.N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=0,
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
  pwr1 <- tmp$power
  if (print) cat(paste("N = ", min.N, " power ", tmp$oper, " ", pwr1, "\n", sep=""))
  if (abs(power - pwr1) <= tol.power)  return(list(sampleSize=min.N, power=pwr1))
  if (pwr1 > power) {
    warning("The desired power is less than min.N]")
    return(list(sampleSize=min.N, power=pwr1)) 
  }

  # Test midpoint
  Nmid <- floor((min.N + max.N)/2)
  tmp  <- PRIMEplus.Power_main(nmax=Nmid, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=0,
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
  pwrm <- tmp$power
  if (print) cat(paste("N = ", Nmid, " power ", tmp$oper, " ",  pwrm, "\n", sep=""))
  if (abs(power - pwrm) <= tol.power)  return(list(sampleSize=Nmid, power=pwrm))

  step <- 100
  if (pwrm > power) {
    N0    <- min.N
    N1    <- Nmid
    diff1 <- abs(pwr1 - power)
    diff2 <- abs(pwrm - power)
    pwr2  <- pwrm
  } else {
    # Test right endpoint
    tmp  <- PRIMEplus.Power_main(nmax=max.N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=0,
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    pwr2 <- tmp$power
    if (print) cat(paste("N = ", max.N, " power ", tmp$oper, " ", pwr2, "\n", sep=""))
    if (abs(power - pwr2) <= tol.power)  return(list(sampleSize=max.N, power=pwr2))
    if (pwr2 < power) {
      warning("The desired power is greater than max.N")
      return(list(sampleSize=max.N, power=pwr2))
    }
    N0    <- Nmid
    N1    <- max.N
    diff1 <- abs(pwrm - power)
    diff2 <- abs(pwr2 - power)
    pwr1  <- pwrm
  }
 
  if ((pwr1 == 0) || (pwr2 == 1)) {
    N      <- floor((N0 + N1)/2)
  } else if (diff1 < diff2) {
    N      <- min(N0 + step, N1)
  } else {
    N      <- floor((N0 + N1)/2)
  }
  iter    <- 0
  minDiff <- 1e100
  minN    <- 0
  minPwr  <- -1
  while (1) {
    tmp  <- PRIMEplus.Power_main(nmax=N, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=0,
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    pwr  <- tmp$power
    oper <- tmp$oper
    if (print) cat(paste("N = ", N, " power ", oper, " ",  pwr, "\n", sep=""))
    diff <- abs(power - pwr)
    if (diff < minDiff) {
      minDiff <- diff
      minN    <- N
      minPwr  <- pwr
    }
    if ((diff <= tol.power) || (abs(N1 - N0) <= tol.N)) break
    if (pwr < power) {
      N0 <- N
    } else {
      N1 <- N
    }
    N  <- floor((N0 + N1)/2) 
  }
  if (oper != "=") {
    tmp  <- PRIMEplus.Power_main(nmax=minN, rand_ratio=rand_ratio, effect_p=effect_p, 
             enroll_rate=enroll_rate, lambda1=lambda1, HR=HR, tau=tau, t1=t1, maxiter=maxiter,
             stopTol=stopTol, alpha=alpha, num_rand=num_rand, nsim=nsim, print=0,
             powerL=PWRL, powerU=PWRU,
             min.sample.size=min.sample.size, min.n.event=min.n.event, min.per.trt=min.per.trt)
    minPwr <- tmp$power
    if (print) cat(paste("N = ", N, " power ", tmp$oper, " ",  minPwr, "\n", sep=""))
  }
  
  list(sampleSize=minN, power=minPwr)

} # END: PRIMEplus.SampleSize



