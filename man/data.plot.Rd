% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.plot.R
\name{data.plot}
\alias{data.plot}
\title{Patient Characteristic Plot}
\usage{
data.plot(
  data,
  covariate,
  half.length = NULL,
  by = "study",
  avg.hline = TRUE,
  fill.str = NULL,
  text.size = 20,
  coord.flip = FALSE
)
}
\arguments{
\item{data}{A \code{BUGSnetData} object produced by \code{data.prep()}}

\item{covariate}{A string indicating the name of the patient characteristic to be plotted}

\item{half.length}{A string indicating how to calculate the half-length of error bars (optional)}

\item{by}{If by="study" then data from arms will be grouped by study/trial. If by="treatment" then
bar graph is grouped by treatment.}

\item{avg.hline}{If TRUE, adds overall average line to plot. Default is TRUE.}

\item{fill.str}{An optional string indicating the variable to categorize measurements. For instance,
some studies report the mean treatment effect and others may report the median treatment effect. If there
is a variable in \code{data} called "type.measure" indicating whether the mean or median is
reported, setting fill.str="type.measure" would colour all studies reporting the mean as red, and all
the studies reporting the median as turquoise.}

\item{text.size}{Font size of the text. Default is 20.}

\item{coord.flip}{If TRUE, coordinates will be flipped to display the plot in portrait mode. Default is FALSE.}
}
\description{
Plots a particular patient characteristic by study or by treatment. Useful for assessing
differences in potential effect modifiers.
}
\examples{
data(diabetes.sim)

diabetes.slr <- data.prep(arm.data = diabetes.sim, 
varname.t = "Treatment", 
varname.s = "Study")

# Example containing a fill.str, an overall average, and no error

data.plot(data = diabetes.slr,
             covariate = "age", 
             fill.str="age_type",
             by = "study")
             
# Example containing no fill.str, no overall average, but contains errorbars

data.plot(data = diabetes.slr,
             covariate = "age", 
             half.length = "age_SD",
             avg.hline=FALSE,
             by = "study")
}
\seealso{
\code{\link{data.prep}}
}
