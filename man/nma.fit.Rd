% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nma.fit.R
\name{nma.fit}
\alias{nma.fit}
\title{Assess Model Fit}
\usage{
nma.fit(nma, plot.pD = TRUE, plot.DIC = TRUE, plot.Dres = TRUE, ...)
}
\arguments{
\item{nma}{A \code{BUGSnetRun} object produced by running \code{nma.run()}.}

\item{plot.pD}{Whether to include pD on the plot. Default is TRUE.}

\item{plot.DIC}{Whether to include DIC on the plot. Default is TRUE.}

\item{plot.Dres}{Whether to include Dres on the plot. Default is TRUE.}

\item{...}{Graphical arguments such as main=, ylab=, and xlab= may be passed as in \code{plot()}. These arguments will only effect the
leverage plot.}
}
\value{
\code{DIC} - A number indicating the Deviance Information Criteria. The DIC is calculated as the sum of \code{Dres} and \code{pD}.
A larger DIC is indicative of a worse model fit.

\code{leverage} - A vector with one value per study arm indicating the leverage of each data point (study arm). Leverage is defined as \code{pmdev} minus the
deviance at the posterior mean of the fitted values.

\code{w} - A vector with one value per study arm. The magnitude of \code{w} represents the data point's contribution to the posterior mean deviance of the
model and is simply the square root of \code{pmdev}. The sign indicates whether the data point is being over (negative sign) or under (positive sign) estimated by
the model and is calculated as the sign of the difference of the observed outcome minus the predicted outcome.

\code{pmdev} - A vector with one value per study arm representing the posterior mean residual deviance for each data point (study arm).

\code{Dres} - The posterior mean of the residual deviance.

\code{pD} - The effective number of parameters, calculated as the sum of the leverages.
}
\description{
Computes the Deviance Information Criteria (DIC) and produces a leverage plot (as in the NICE Technical Support Document 2)
for a given model. These can be used to assess and compare the fit of different models (i.e fixed vs random effects, consistency vs
inconsistency). \code{nma.fit} also produces a plot comparing the leverage of each data point against their contribution to
the total posterior deviance. Points lying outside the purple dotted line are generally identified as contributing to the model's poor fit.
Points with high leverage are influencial i.e. they have a stong influence on the estimates.
}
\examples{
data(diabetes.sim)

diabetes.slr <- data.prep(
  arm.data = diabetes.sim, 
  varname.t = "Treatment", 
  varname.s = "Study"
)

#Random effects, consistency model.
#Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
diabetes.re.c <- nma.model(
  data = diabetes.slr,
  outcome = "diabetes", 
  N = "n",
  reference = "Placebo",
  family = "binomial",
  link = "cloglog",
  effects = "random",
  type = "consistency",
  time = "followup"
)
 
diabetes.re.c.res <- nma.run(
  model = diabetes.re.c,
  n.adapt = 100,
  n.burnin = 0,
  n.iter = 100
)

#Fixed effects, consistency model.
#Binomial family, cloglog link. This implies that the scale will be the Hazard Ratio.
diabetes.fe.c <- nma.model(
  data = diabetes.slr,
  outcome = "diabetes", 
  N = "n",
  reference = "Placebo",
  family = "binomial",
  link = "cloglog",
  effects = "fixed",
  type="consistency",
  time="followup"
)
 
diabetes.fe.c.res <- nma.run(
  model = diabetes.fe.c,
  n.adapt = 100,
  n.burnin = 0,
  n.iter = 100
)  

#Compare fixed vs random effects via leverage plots and DIC 
par(mfrow=c(1,2))
nma.fit(diabetes.fe.c.res, main = "Fixed Effects Model")
nma.fit(diabetes.re.c.res, main= "Random Effects Model")
}
\seealso{
\code{\link{nma.run}}
}
