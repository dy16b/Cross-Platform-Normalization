\name{MatchMixeR}
\alias{MM}

\title{MatchMixeR: A cross-platform normalization procedure for matched samples.}

\description{
  Matched sample based cross-platform normalization function using mixed effect regression
}
\usage{MM(Xmat, Ymat)}

\arguments{
  \item{X}{gene expression level matrix on platform X (to be transformed)}

  \item{Y}{gene expression level matrix on platform Y ( not to be transformed)}
}

\details{
  This function fits mixed effect regression model.
}

\value{8 lists including betahat, betamat, Yhat}

\author{Xing Qiu}

\examples{
fit_mm <- MM(gpl96,gpl570)
}

\keyword{MatchMixeR}
\keyword{MM}