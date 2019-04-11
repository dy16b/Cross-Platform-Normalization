\name{MatchMixeR}
\alias{flmer}

\title{MatchMixeR: A cross-platform normalization procedure for matched samples.}

\description{
  Fast linear mixed effects regression (flmer); a moment-based method for the MM Model
}

\usage{flmer(Xmat, Ymat)}

\arguments{
  \item{X}{gene expression level matrix on platform X (to be transformed)}

  \item{Y}{gene expression level matrix on platform Y ( not to be transformed)}
}

\details{
  This function fits the mixed effect regression model (MM) using fast moment based approaches.
}

\value{8 lists including betahat, betamat, Yhat}

\author{Xing Qiu}

\examples{
fit_flmer <- flmer(gpl96,gpl570)
}

\keyword{MatchMixeR}
\keyword{flmer}
