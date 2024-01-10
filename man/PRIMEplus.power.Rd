
\name{PRIMEplus.Power}
\alias{PRIMEplus.Power}
\title{Power}
\description{
   Compute the power using LRT Re-randomization test                   
}
\usage{ 
PRIMEplus.Power(nmax=500, rand_ratio=0.5, effect_p=0.6, 
    enroll_rate=380*0.25/6, lambda1=0.117, HR=0.5, tau=12*5, t1=1, 
    maxiter=100000, stopTol=1e-4, alpha=0.05, num_rand=1000, nsim=10000, 
    print=0, min.sample.size=50, min.n.event=5, min.per.trt=0.25)
}
\arguments{
  \item{nmax}{Sample size}
  \item{rand_ratio}{Probability of assignment to treatment arm}
  \item{effect_p}{Vector for proportion of responders in the treatment arm at baseline (see details)}
  \item{enroll_rate}{Enrollment rate in subjects per month}
  \item{lambda1}{Baseline hazard in terms of months} 
  \item{HR}{Vector of hazard ratios for responders against controls (see details)}
  \item{tau}{Total study duration}
  \item{t1}{Delayed duration in months}
  \item{maxiter}{Maximum number of iterations in the EM algorithm. The default is 100000.}
  \item{stopTol}{Stopping tolerance in the EM algorithm. The default is 1e-4.}
  \item{alpha}{Significance level. The default is 0.05.}
  \item{num_rand}{The number of replications in the re-randomization test. The default is 1000.}
  \item{nsim}{The number of simulations. The default is 1000.}
  \item{print}{0 or 1 to print information. The default is 0.}
  \item{min.sample.size}{Minimum sample size. The default is 50.}
  \item{min.n.event}{Minimum number of events. The default is 5.}
  \item{min.per.trt}{Minimum proportion of controls and treated subjects. The default is 0.25.}
}
\details{
  The length and order of \code{effect_p} must be the same as \code{HR}. Both of these vectors
  should contain information only for the groups of responders. For example, if there are
  full responders and partial responders, then \code{effect_p} and \code{HR} would be vectors 
  of length two.

  For each simulation, a simulated data set is created from the 
  \code{\link{generate_data}} function and then an estimated p-value is computed 
  by calling \code{\link{PRIMEplus.ReRand.LRT}}.
  The power is calculated as the proportion of iterations whose estimated p-value
  was less than or equal to \code{alpha}.
}
\value{ A list containing the power and the number of simulated datasets used
        in the calculation.
}
%\references{}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
%\examples{}
\keyword{ power} 
\seealso{ \code{\link{PRIMEplus.ReRand.LRT}} }
