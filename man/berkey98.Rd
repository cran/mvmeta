\name{berkey98}
\alias{berkey98}
\docType{data}

\title{ Five Published Trials on Periodontal Disease }

\description{
The dataset contains the results of 5 published trials comparing surgical and non-surgical treatments for medium-severity periodontal disease, one year after treatment. The 2 estimated outcomes are average improvements (surgical minus non-surgical, in mm) in probing depth (PD) and attachment level (AL).
}

\usage{data(berkey98)}

\format{
  A data frame with 5 observations on the following 8 variables.
  \describe{
    \item{\code{trial} }{ number of trial.}
    \item{\code{pubyear} }{publication year of the trial.}
    \item{\code{npat} }{ number of patients included in the trial.}
    \item{\code{PD}}{ estimated improvement of surgical vs. non-surgical treatments in probing depth (mm).}
    \item{\code{AL}}{ estimated improvement of surgical vs. non-surgical treatments in attachment level (mm).}
    \item{\code{var_PD}}{ variance of the estimated outcome for \code{PD}.}
    \item{\code{cov_PD_AL}}{ variance of the estimated outcome for \code{AL}.}
    \item{\code{var_AL}}{ covariance of the estimated outcomes for \code{PD} and \code{AL}.}
  }
}

\source{
Berkey C. S., Hoaglin D. C., Antczak-Bouckoms A., Mosteller F, Colditz, G. A. Meta-analysis of multiple outcomes by regression with random effects. \emph{Statistics in Medicine}. 1998; \bold{17}:2537-2550.
}

\keyword{datasets}