
\name{generate_data}
\alias{generate_data}
\title{Simulated data}
\description{
   Generate simulated data                          
}
\usage{ 
generate_data(N=500, rand_ratio=0.5, effect_p=0.6, enroll_rate=380*0.25/6, 
              lambda1=0.117, HR=0.5, tau=12*5, t1=1)
}
\arguments{
  \item{N}{Maximum sample size}
  \item{rand_ratio}{Allocation ratio}
  \item{effect_p}{Vector for proportion of responders in the treatment arm at baseline (see details)}
  \item{enroll_rate}{Enrollment rate in subjects per month}
  \item{lambda1}{Baseline hazard in terms of months} 
  \item{HR}{Vector of hazard ratios for responders against controls (see details)}
  \item{tau}{Total study duration}
  \item{t1}{Delayed duration}
}
\details{
  The length and order of \code{effect_p} must be the same as \code{HR}. Both of these vectors
  should contain information only for the groups of responders. For example, if there are
  full responders and partial responders, then \code{effect_p} and \code{HR} would be vectors 
  of length two.
}
\value{ A data frame with columns:
  \tabular{ll}{
    Name          \tab Description                                  \cr
    id            \tab id variable                                  \cr
    trt           \tab treatment allocation: 1 = treatment arm      \cr
    D1            \tab patient's response status for group 1        \cr
    D2            \tab patient's response status for group 2        \cr
    Dm            \tab patient's response status for non-responders \cr
    tau           \tab total study duration                         \cr
    enroll_time   \tab patients' enrollment times                   \cr
    time_to_event \tab patients' event times                        \cr
    event_status  \tab censoring indicator                          \cr
    X             \tab observational time                           \cr
    t1            \tab delayed duration                                              
  }
}
%\references{}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov>, Yongsoek Park <yongpark@pitt.edu> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data <- generate_data()
  data[1:5, ]
}
\keyword{ simulation } 
%\seealso{ \code{\link{getHazard}}, \code{\link{EM.nonparametric}} }
