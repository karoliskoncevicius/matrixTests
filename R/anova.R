#' One Sample t-test
#'
#' Performs a t-test on each row of a matrix.
#'
#' Functions to perform one sample and two sample t-tests for rows of matrices.
#' Main arguments and results were intentionally matched to the \code{t.test()}
#' function from default stats package. Other arguments were split into separate
#' functions:
#'
#' \code{mat_ttest_single()} - t-test for mean of a single group. Same as \code{t.test(x)}
#'
#' \code{mat_ttest_equalvar()} - groups have equal variance. Same as \code{t.test(x, y, var.equal=TRUE)}
#'
#' \code{mat_ttest_welch()} - Welch approximation. Same as \code{t.test(x, y)}
#'
#' \code{mat_ttest_paired()} - paired t-test. Same as \code{t.test(x, y, paired=TRUE)}
#'
#' @name anova
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param mu true values of the means for the null hypothesis.
#' A single number or numeric vector of length nrow(x).
#' @param alternative alternative hypothesis to use for each row of x.
#' A single string or a vector of length nrow(x).
#' Must be one of "two.sided" (default), "greater" or "less".
#' @param conf.level confidence levels used for the confidence intervals.
#' A single number or a numeric vector of length nrow(x).
#' All values must be in the range of [0;1].
#'
#' @return a data.frame where each row contains the results of a t.test
#' performed on the corresponding row of x. The columns will vary depending on
#' the type of test performed.
#'
#' @seealso \code{t.test()}
#'
#' @examples
#' X <- t(iris[iris$Species=="setosa",1:4])
#' Y <- t(iris[iris$Species=="virginica",1:4])
#' mat_ttest_welch(X, Y)
#'
#' # same row using different confidence levels
#' mat_ttest_equalvar(X[c(1,1,1),], Y[c(1,1,1),], conf.level=c(0.9, 0.95, 0.99))
#'
#' @author Karolis Koncevičius
#' @export
mat_anova_oneway <- function(x, groups) {

  if(!is.null(x) && is.vector(x))
    x <- matrix(x, nrow=1)

  if(!is.null(x) && is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  M <- rowMeans(x, na.rm=TRUE)
  withinScatter <- 0
  betweenScatter <- 0
  for(g in unique(groups)) {
    gMeans <- rowMeans(x[,groups==g, drop=FALSE], na.rm=TRUE)
    groupScatter   <- rowSums((x[,groups==g, drop=FALSE]-gMeans)^2, na.rm=TRUE)
    withinScatter  <- withinScatter + groupScatter
    betweenScatter <- betweenScatter + rowSums(!is.na(x[,groups==g, drop=FALSE]))*(gMeans-M)^2
  }

  N <- length(unique(groups))
  n <- rowSums(!is.na(x))
  F <- (betweenScatter/(N-1)) / (withinScatter/(n-N))
  p <- pf(F, N-1, n-N, lower.tail=FALSE)

  data.frame(p.value=p, df1=N-1, df2=n-N, F=F)
}

