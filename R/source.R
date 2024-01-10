generate_data <- function(N=500,
                          rand_ratio=0.5,
                          effect_p=0.6,
                          enroll_rate=380*0.25/6,
                          lambda1=0.117,
                          HR=0.5,
                          tau=12*5,
                          t1=1) 
{
  
  # output data
  if (length(effect_p) != length(HR)) stop("ERROR: length(effect_p) != length(HR)")
  psum  <- sum(effect_p)
  if (psum > 1) stop("ERROR: sum(effect_p) > 1")
  p.vec  <- c(effect_p, 1-psum)
  HR.vec <- c(HR, 1)
  M      <- length(p.vec)
  dat    <- as.data.frame(matrix(0, nrow=N, ncol=8+M))
  Dnms   <- paste0("D", 1:M)
  names(dat) <- c("id", "trt", Dnms, "tau", "enroll_time", "time_to_event", "event_status", "X", "t1")
  
  dat$id = 1:N
  
  # treatment allocation: 1 = experimental arm
  # sort by treatment assignment status to make the responders only arise among treated subjects 
  # set.seed(seed)
  dat$trt <- sort(rbinom(n = N, size = 1, prob=rand_ratio), decreasing = TRUE)
  n_trt = sum(dat$trt == 1)
  n_cnt = sum(dat$trt == 0)
  
  # simulate the response statue among the experimental arm with proportion effect_P: constraint: proportion of responders = effect_p
  D.trt  = t(rmultinom(n = n_trt, size = 1, prob = p.vec))
  D.cnt  = t(replicate(n_cnt, c(rep(0, M-1),1)))
  for (i in 1:M) {
    dat[[Dnms[i]]] <- rbind(D.trt, D.cnt)[,i]
  }

  # total study duration 
  dat$tau = tau
  
  # enrollment follows expected trajectory from Poisson process with memoryless property
  ## waiting time between two consective subjects follows an Exponential distribution with enrollment rate
  ## arrival time for each subject is converted from the waitig by taking the cumulative sum of the waiting times
  waiting_times <- rexp(n = N, rate=enroll_rate)
  dat$enroll_time <- cumsum(waiting_times)
  
  # t*
  dat$t1 <- t1
  
  # time to event
  n_nresp <- sum(dat$Z == 0) # number of non-responders among the experimental and control groups
  tvec    <- c(0, dat$t1[1])
  for (i in 1:M) {
    vec                    <- dat[, Dnms[i], drop=TRUE]
    tmp                    <- vec == 1
    n_resp.D               <- sum(tmp)
    dat$time_to_event[tmp] <- rpexp(n_resp.D, rate=c(lambda1, HR.vec[i]*lambda1), t=tvec)
  }  

  # build in the administrative censoring 
  dat$event_status <- ifelse(dat$time_to_event <= tau - dat$enroll_time, 1, 0)
  
  # observational time
  dat$X <- ifelse(dat$time_to_event <= tau - dat$enroll_time, dat$time_to_event, tau - dat$enroll_time)
  tmp   <- dat$X <= 0
  if (any(tmp)) {
    warning("Some observational times are <= 0. Try increasing tau or changing other parameters.")  
    dat <- dat[!tmp, , drop=FALSE]
  }
  row.names(dat) <- NULL
  
  dat
}

checkTrtStatus <- function(trt, status, min.sample.size=50, min.n.event=5, min.per.trt=0.25) {

  ret <- 0

  N <- length(trt)
  if (N < min.sample.size) {
    msg <- paste("WARNING: Sample size < ", min.sample.size, ". This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }
  M     <- ceiling(N*min.per.trt)
  ntrt0 <- sum(trt == 0)
  ntrt1 <- N - ntrt0
  if (!ntrt0) stop("ERROR: all subjects have treatment=1")
  if (!ntrt1) stop("ERROR: all subjects have treatment=0")
  if (ntrt0 < M) {
    msg <- paste("WARNING: Too few subjects with treatment=0. This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }  
  if (ntrt1 < M) {
    msg <- paste("WARNING: Too few subjects with treatment=1. This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }  
  Nev <- sum(status == 1)
  if (!Nev) stop("ERROR: all subjects have event_status=0")
  if (Nev < min.n.event) {
    msg <- paste("WARNING: Too few subjects with event_status=1. This can lead to poor results.", sep="")
    warning(msg)
    ret <- 1
  }

  ret

}

getProbResponder <- function(P, N, n_trt) {

  J   <- length(P)
  ret <- matrix(data=0, nrow=N, ncol=J)
  ret[1:n_trt, ] <- matrix(data=rep(P, times=n_trt), nrow=n_trt, ncol=J, byrow=TRUE)
  ret

}

setUpAndCheckData <- function(data, time.var, trt.var, status.var, lambda0, effect_p, beta0) {

  if ( (!is.data.frame(data)) && (!is.matrix(data)) ) stop("ERROR: data must be a data frame or matrix")
  if (!nrow(data)) stop("data contains no rows")
  if (!ncol(data)) stop("data contains no columns")

  cx <- colnames(data)
  check_variable(time.var, "time.var", cx) 
  check_variable(trt.var, "trt.var", cx) 
  check_variable(status.var, "status.var", cx) 

  # data must be ordered by treatment
  ord          <- order(data[, trt.var], decreasing = TRUE)
  data         <- data[ord, , drop=FALSE]
  treatment    <- data[, trt.var]
  event_status <- data[, status.var]
  time         <- data[, time.var]
  n            <- length(time)
  if (length(treatment) != n) stop("ERROR: length(treatment) != length(time)")
  if (length(event_status) != n) stop("ERROR: length(event_status) != length(time)")
  tmp <- (!is.finite(time)) | (time < 0)
  tmp[is.na(tmp)] <- TRUE
  if (any(tmp)) stop("ERROR: with time vector")
  if (!all(treatment %in% 0:1)) stop("ERROR: with treatment vector")
  if (!all(event_status %in% 0:1)) stop("ERROR: with event_status vector")
  
  #if (!length(effect_p)) { 
  #  check_group.var(data, group.var)
  #  grpvec   <- data[, group.var, drop=TRUE]
  #  effect_p <- estimate_P(grpvec) 
  #  beta0    <- estimate_beta(data, time.var, trt.var, status.var, group.var)
  #  n        <- length(effect_p)
  #  effect_p <- check_effect_p(effect_p)
  #  beta0    <- check_beta(beta0, n)
  #}

  tmp     <- get_tfailo(time, event_status)
  t.fail  <- tmp$t.fail.o
  delta.l <- tmp$delta.l
  if (!length(lambda0)) lambda0 <- getHazard(time, treatment, event_status, t.fail)$hazard

  n_trt <- sum(treatment)
  #if (is.null(probResponder)) probResponder <- getProbResponder(effect_p, nrow(data), n_trt)
  #if (nrow(probResponder) != nrow(data)) stop("ERROR: nrow(probResponder) != nrow(data)")

  # Return the reverse order
  list(time=time, treatment=treatment, event_status=event_status, order=order(ord),
       lambda0=lambda0, t.fail=t.fail, delta.l=delta.l,
       effect_p=effect_p, beta0=beta0)

} # END: setUpAndCheckData

check_numeric <- function(obj, nm, min=0, max=1, inc.min=1, inc.max=1, len=Inf) {

  n <- length(obj)
  if (!n) stop(paste0("ERROR: ", nm,  " must be specified"))
  if (is.finite(len) && (n != len)) stop(paste0("ERROR: ", nm,  " must have length ", len))
  if (!is.numeric(obj)) stop(paste0("ERROR: ", nm,  " must be numeric"))
  if (any(!is.finite(obj))) {
    if (len == 1) {
      stop(paste0("ERROR: ", nm,  " must be a finite numeric value"))
    } else {
      stop(paste0("ERROR: ", nm,  " must have finite numeric values"))
    }
  }

  tmp <- TRUE
  if (inc.min) {
    tmp <- tmp & (obj >= min) 
    op1 <- " >= "
  } else {
    tmp <- tmp & (obj > min) 
    op1 <- " > "
  }
  if (inc.max) {
    tmp <- tmp & (obj <= max) 
    op2 <- " <= "
  } else {
    tmp <- tmp & (obj < max) 
    op2 <- " < "
  }
  if (!all(tmp)) stop(paste0("ERROR: ", nm,  " must be", op1, min, " and", op2, max))
  NULL

}

check_numeric_matrix <- function(obj, nm, min.nr=0, max.nr=Inf, min.nc=1, max.nc=Inf) {

  n <- length(obj)
  if (!n) stop(paste0("ERROR: ", nm,  " must be specified"))
  if (!is.numeric(obj)) stop(paste0("ERROR: ", nm,  " must be numeric"))
  if (!is.vector(obj) && !is.matrix(obj)) stop(paste0("ERROR: ", nm,  " must be a numeric vector or matrix"))
  nr <- nrow(obj)
  if (is.null(nr)) nr <- 0
  nc <- ncol(obj)
  if (is.null(nc)) nc <- 0
  if (nr && min.nr && (nr < min.nr)) stop(paste0("ERROR: ", nm,  " must have at least ", min.nr, " rows"))
  if (nr && is.finite(max.nr) && (nr > max.nr)) stop(paste0("ERROR: ", nm,  " cannot have more than ", max.nr, " rows"))
  if (nc && min.nc && (nc < min.nc)) stop(paste0("ERROR: ", nm,  " must have at least ", min.nc, " columns"))
  if (nc && is.finite(max.nc) && (nc > max.nc)) stop(paste0("ERROR: ", nm,  " cannot have more than ", max.nc, " columns"))
  if (!is.matrix(obj)) {
    if (n != min.nc) stop(paste0("ERROR: ", nm,  " must be a vector of length ", min.nc, " or a matrix with ", min.nc, " columns."))
  }
  if (any(!is.finite(obj))) stop(paste0("ERROR: ", nm,  " must have finite numeric values"))
 
  NULL

}

check_integer <- function(obj, nm, min=0, max=1, inc.min=1, inc.max=1, len=Inf) {

  check_numeric(obj, nm, min=min, max=max, inc.min=inc.min, inc.max=inc.max)
  if (obj != as.integer(obj)) stop(paste0("ERROR: ", nm, " must be an integer"))

  NULL

}

check_effect_p <- function(obj) {

  check_numeric(obj, "effect_p", min=0, max=1, inc.min=0, inc.max=0, len=Inf)
  s <- sum(obj)
  if (s >= 1) stop(paste0("ERROR: sum(effect_p) must be < 1"))
  ret <- c(obj, 1-s)
  ret

}

check_t1 <- function(obj) {
  check_numeric(obj, "t1", min=0, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_lambda0 <- function(obj, status, time) {
  if (is.null(obj)) return(NULL)
  tmp   <- status %in% 1
  nfail <- length(unique(time[tmp])) 
  check_numeric(obj, "lambda0", min=-Inf, max=Inf, inc.min=0, inc.max=0, len=nfail)
  NULL
}

check_stopTol <- function(obj) {
  check_numeric(obj, "stopTol", min=0, max=Inf, inc.min=0, inc.max=0, len=1)
  NULL
}

check_maxiter <- function(obj) {
  check_integer(obj, "maxiter", min=1, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_print <- function(obj) {
  check_integer(obj, "print", min=0, max=2, inc.min=1, inc.max=1, len=1)
  NULL
}

check_num_rand <- function(obj) {
  check_integer(obj, "num_rand", min=1, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_isString <- function(obj) {
  if ((length(obj) == 1) && (is.character(obj))) {
    ret <- TRUE
  } else {
    ret <- FALSE
  }
  ret
}

check_variable <- function(var, name, cx) {

  if (!check_isString(var)) stop(paste0("ERROR: ", name, " must be a string"))
  if (!(var %in% cx)) stop(paste0("ERROR: ", name, " must be a variable in data"))
  NULL

}

get_random_beta <- function(len, min=-1, max=1) {

  while (1) {
    ret <- runif(len, min=min, max=max)
    if (any(ret == 0)) next
    if (any(duplicated(ret))) next
    break
  }
  ret
}

check_beta <- function(obj, len) {

  if (length(obj)) {
    check_numeric_matrix(obj, "beta0", min.nr=0, max.nr=Inf, min.nc=len, max.nc=len) 
    #check_numeric(obj, "beta0", min=-Inf, max=Inf, inc.min=0, inc.max=0, len=len) 
    tmp <- obj %in% 0
    if (any(tmp)) stop("ERROR: elements of beta0 cannot be 0.")
    if (!is.matrix(obj)) {
      tmp <- duplicated(obj)
      if (any(tmp)) stop("ERROR: beta0 cannot have duplicated values.")
      ret <- c(obj, 0)
    } else {
      nc  <- ncol(obj)
      tmp <- rep(FALSE, nrow(obj))
      if (nc > 1) {
        for (i in 1:(nc-1)) {
          for (j in (i+1):nc) tmp <- tmp | (obj[, i] == obj[, j])
        }
        if (any(tmp)) stop("ERROR: beta0 cannot have duplicated values within each row.")
      }
      ret <- cbind(obj, 0)
    }
  } else {
    ret <- c(get_random_beta(len, min=-1, max=1), 0)
  }
  ret

}

check_nmax <- function(obj) {
  check_integer(obj, "nmax", min=10, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_rand_ratio <- function(obj) {
  check_numeric(obj, "rand_ratio", min=0, max=1, inc.min=0, inc.max=0, len=1)
  NULL
}

check_enroll_rate <- function(obj) {
  check_numeric(obj, "enroll_rate", min=0, max=Inf, inc.min=0, inc.max=0, len=1)
  NULL
}

check_lambda1 <- function(obj) {
  check_numeric(obj, "lambda1", min=-Inf, max=Inf, inc.min=0, inc.max=0, len=1)
  NULL
}

check_HR <- function(obj, len) {
  if (len != length(obj)) stop("ERROR: length(HR) != length(effect_p)")
  check_numeric(obj, "HR", min=0, max=Inf, inc.min=0, inc.max=0, len=len)
  c(obj, 1)
}

check_tau <- function(obj) {
  check_numeric(obj, "tau", min=0, max=Inf, inc.min=0, inc.max=0, len=1)
  NULL
}

check_alpha <- function(obj) {
  check_numeric(obj, "alpha", min=0, max=1, inc.min=0, inc.max=0, len=1)
  NULL
}

check_num_rand <- function(obj) {
  check_integer(obj, "num_rand", min=1, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_nsim <- function(obj) {
  check_integer(obj, "nsim", min=1, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_min.sample.size <- function(obj) {
  check_integer(obj, "min.sample.size", min=10, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_min.n.event <- function(obj) {
  check_integer(obj, "min.n.event", min=1, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_min.per.trt <- function(obj) {
  check_numeric(obj, "min.per.trt", min=0, max=1, inc.min=0, inc.max=0, len=1)
  NULL
}

check_power <- function(obj) {
  check_numeric(obj, "power", min=0, max=1, inc.min=0, inc.max=0, len=1)
  NULL
}

check_tol.power <- function(obj) {
  check_numeric(obj, "tol.power", min=0, max=1, inc.min=0, inc.max=0, len=1)
  NULL
}

check_min.N <- function(obj) {
  check_integer(obj, "min.N", min=10, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_max.N <- function(obj, min.N) {
  check_integer(obj, "max.N", min=min.N, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}

check_tol.N <- function(obj, min.N) {
  check_integer(obj, "tol.N", min=1, max=Inf, inc.min=1, inc.max=0, len=1)
  NULL
}


check_group.var <- function(data, v) {

  check_variable(v, "group.var", colnames(data)) 
  vec  <- as.numeric(data[, v, drop=TRUE])
  uvec <- unique(vec)
  m    <- length(uvec)
  if (!all(vec %in% 0:(m-1))) stop("ERROR: group.var must be coded as 0, 1, ..., with 0=controls")
  NULL
}