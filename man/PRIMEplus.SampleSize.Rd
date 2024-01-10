
\name{PRIMEplus.SampleSize}
\alias{PRIMEplus.SampleSize}
\title{Sample Size}
\description{
   Compute the sample size for a given power                   
}
\usage{ 
PRIMEplus.SampleSize(power=0.8, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
              lambda1=0.117, HR=0.5, tau=12*5, t1=1, maxiter=100000, stopTol=1e-4,
              alpha=0.05, num_rand=1000, nsim=10000, min.N=100, max.N=700,
              tol.power=0.01, tol.N=1, print=1,
              min.sample.size=50, min.n.event=5, min.per.trt=0.25)
}
\arguments{
  \item{power}{The desired power. The default is 0.8.}
  \item{rand_ratio}{Allocation ratio}
  \item{effect_p}{Vector for proportion of responders in the treatment arm at baseline (see details)}
  \item{enroll_rate}{Enrollment rate in subjects per month}
  \item{lambda1}{Baseline hazard in terms of months} 
  \item{HR}{Vector of hazard ratios for responders against controls (see details)}
  \item{tau}{Total study duration}
  \item{t1}{Delayed duration}
  \item{maxiter}{Maximum number of iterations in the EM algorithm. The default is 100000.}
  \item{stopTol}{Stopping tolerance in the EM algorithm. The default is 1e-4.}
  \item{alpha}{Significance level. The default is 0.05.}
  \item{num_rand}{Number of replications in the re-randomization test. The default is 1000.}
  \item{nsim}{Number of simulations in computing power (\bold{see Details}). 
            The default is 10000.}
  \item{min.N}{Lower bound for the sample size. The default is 100.}
  \item{max.N}{Upper bound for the sample size. The default is 700.}
  \item{tol.power}{Stopping tolerance for the power. The default is 0.01.}
  \item{tol.N}{Stopping tolerance for the sample size. The default is 1.}
  \item{print}{0 or 1 to print information. The default is 1.}
  \item{min.sample.size}{Minimum sample size. The default is 50.}
  \item{min.n.event}{Minimum number of events. The default is 5.}
  \item{min.per.trt}{Minimum proportion of controls and treated subjects. The default is 0.25.}
}
\details{
  The length and order of \code{effect_p} must be the same as \code{HR}. Both of these vectors
  should contain information only for the groups of responders. For example, if there are
  full responders and partial responders, then \code{effect_p} and \code{HR} would be vectors 
  of length two.

  This uses a bisection method to estimate the sample size. At each iteration,
  the estimated power \code{power_est} is computed using \code{\link{PRIMEplus.Power}} 
  for a given sample size holding all other parameters fixed.
  The algorithm terminates when \code{abs(power - power_est) <= tol.power} or
  when the length of the estimated interval containing the sample size is
  less than or equal to \code{tol.N}.
  
  \bold{NOTE:} \cr
  \bold{It is important to note that the power for a given sample size is estimated by 
  running a simulation. Thus, by setting a different seed, a different result may 
  be returned. Therefore, to ensure a more precise estimated sample size, set the
  option \code{nsim} to a large value and/or run this function several times by
  setting different seeds and examine the distribution of returned sample sizes.}

}
\value{ A list containing the sample size and power.
}
%\references{}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
%\examples{}
\keyword{ sample size } 
\seealso{ \code{\link{PRIMEplus.Power}} }
