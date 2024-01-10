
\name{getHazard}
\alias{getHazard}
\title{Compute initial estimates for the baseline hazard}
\description{
   Calls the coxph function to compute initial estimates for the baseline hazard                        
}
\usage{ getHazard(time, treatment, event_status, t.fail.o=NULL)
}
\arguments{
  \item{time}{Vector of times.}
  \item{treatment}{Binary vector of treatments (1=subject received treatment).}
  \item{event_status}{Binary vector of event status (1=subject experienced an event).}
  \item{t.fail.o}{NULL or vector of event times.}
}
%\details{}
\value{ List containing the baseline hazards ordered by the event times.
}
%\references{}
\author{
Zhenzhen Xu <Zhenzhen.Xu@fda.hhs.gov> and Bin Zhu <bin.zhu@nih.gov>
}
\examples{
  data(data, package="PRIMEplus")
  lambda0 <- getHazard(data[, "X"], data[, "trt"], data[, "event_status"])$hazard
  lambda0[1:10]
}
\keyword{ hazard } 
\seealso{ \code{\link{PRIMEplus.EM}} }
