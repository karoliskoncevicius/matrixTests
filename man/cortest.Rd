% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor.R
\name{cortest}
\alias{cortest}
\alias{cor_pearson}
\title{correlation}
\usage{
cor_pearson(x, y, alternative = "two.sided", conf.level = 0.95)
}
\arguments{
\item{x}{numeric matrix.}

\item{y}{numeric matrix for the second group of observations.}

\item{alternative}{alternative hypothesis to use for each row of x.
A single string or a vector of length nrow(x).
Must be one of "two.sided" (default), "greater" or "less".}

\item{conf.level}{confidence levels used for the confidence intervals.
A single number or a numeric vector of length nrow(x).
All values must be in the range of [0;1].}
}
\value{
a data.frame where each row contains the results of a correlation
test performed on the corresponding row of x. The columns will vary
depending on the type of test performed.
}
\description{
Performs a correlation test on each row of a the input matrix.
}
\details{
Functions to perform various correlation tests for rows of matrices.
Main arguments and results were intentionally matched to the \code{cor.test()}
function from default stats package.

\code{cor_pearson()} - test for Pearson correlation coefficient.
Same as \code{cor.test(x,y,method="pearson")}
}
\examples{
X <- t(iris[iris$Species=="setosa",1:4])
Y <- t(iris[iris$Species=="virginica",1:4])
cor_pearson(X, Y)

}
\seealso{
\code{cor.test()}
}
\author{
Karolis Koncevičius
}