
\name{PRIMEplus.ReRand.LRT}
\alias{PRIMEplus.ReRand.LRT}
\title{Randomization Test}
\description{
   PRIMEplus randomization likelihood ratio test                           
}
\usage{ PRIMEplus.ReRand.LRT(data, effect_p, beta0, time.var="X", trt.var="trt", 
                             status.var="event_status", id.var="id", t1=1, lambda0=NULL, 
                             stopTol=1e-4, maxiter=100000, print=0, num_rand=1000)
}
\arguments{
  \item{data}{Data frame or matrix containing a time-to-event variable (\code{time.var}),
              a treatment variable (\code{trt.var}), 
              and a censoring variable (\code{status.var}).}
  \item{effect_p}{Vector of proportions for groups of responders in the treatment arm at baseline (see details).}
  \item{beta0}{Vector or matrix of initial estimates for the log-hazard ratios (see details).}
  \item{time.var}{Time-to-event variable name in \code{data}. The default is "X".}
  \item{trt.var}{Binary treatment variable name in \code{data} coded as 
                 0 for controls and 1 for subjects that received treatment.}
  \item{status.var}{Name of the binary censoring variable in \code{data} coded
                   as 0 for censored subjects and 1 for subjects that
                   experienced an event.}
  \item{id.var}{NULL or subject id variable in \code{data}. The default is "id".}
  \item{t1}{Delayed duration. The default is 1.}
  \item{lambda0}{NULL or vector of initial estimates for the baseline hazards
                 corrsponding to the ordered event times.}
  \item{stopTol}{Stopping tolerance. The default is 1e-4.}
  \item{maxiter}{Maximum number of iterations. The default is 100000.}
  \item{print}{0-2 to print information. Larger values will print more information.
               The default is 0.}
  \item{num_rand}{The number of randomizations. The default is 1000.}
}
\details{
  The EM algorithm is sensitive to the initial values of the log-hazard ratios (\code{beta0}), so different
  initial estimates should be tried to ensure the maximum log-likelihood is obtained. 
  Thus, \code{beta0} can be a vector or matrix, where in the case of a matrix, each row corresponds to a different set
  of intial estimates. Each set of initial estimates must contain distinct non-zero values.
  The length and order of \code{effect_p} must be the same as the columns of \code{beta0}. Both of these 
  should contain information only for the groups of responders. For example, if there are
  full responders and partial responders, then \code{effect_p} would be a vector of length two,
  and \code{beta0} would be a vector of length two or a matrix with two columns.
  The initial values are only used for the observed data; each randomization uses the MLE estimates
  as initial estimates.
}
\value{ A list containing the objects:
  \tabular{ll}{
    Name                \tab Description \cr
    pvalue.LRT          \tab p-value of the randomization test based on the likelihood ratio test \cr
    pvalue.loglike.marg \tab p-value of the randomization test based on the marginal likelihood \cr
    n.randomizations    \tab the number of randomizations that the p-values are based on
  }
}
%\references{}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="PRIMEplus")
  effp  <- c(0.3, 0.3)
  beta0 <- c(log(0.4), log(0.6))
  ret   <- PRIMEplus.ReRand.LRT(data, effp, beta0, num_rand=100)
  ret
}
\keyword{ permutation } 
\seealso{ \code{\link{PRIMEplus.LRT}} }
