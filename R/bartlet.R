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
bartlett <- function(x, groups) {

  if(!is.null(x) && is.vector(x))
    x <- matrix(x, nrow=1)

  if(!is.null(x) && is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  if(!is.null(groups))
    groups <- as.character(groups)

  assert_character_vec_length(groups, ncol(x))

  bad <- is.na(groups)
  if(any(bad)) {
    warning(sum(bad), " columns skipped due to missing group information")
    x      <- x[,!is.na(groups)]
    groups <- groups[!is.na(groups)]
  }

  nGroups  <- numeric(nrow(x))
  nSamples <- numeric(nrow(x))
  tooSmall <- logical(nrow(x))
  for(g in unique(groups)) {
    nGroupObs <- rowSums(!is.na(x[,groups==g, drop=FALSE]))
    tooSmall[nGroupObs < 2] <- TRUE
    nSamples  <- nSamples + nGroupObs
    nGroups   <- nGroups + ifelse(nGroupObs==0, 0, 1)
    vtot      <- rVars(x[,groups==g, drop=FALSE]) * nGroupObs
  }

  if(any(nGroups < 2))
    warning(sum(nGroups < 2), " rows only had observations for one group")

  if(any(tooSmall))
    warning(sum(tooSmall), " rows had groups with less than 2 observations")

  vtot <- vtot/(nSamples - nGroups)
  df   <- nGroups-1
  ksq  <- (n.total * log(v.total) - sum(n * log(v))) /
           (1 + (sum(1/n) - 1/n.total)/(3 * (k - 1)))

  p <- pchisqksq, df, lower.tail=FALSE)

  data.frame(sum.sq.treatment=betweenScatter, sum.sq.residuals=withinScatter,
             mean.sq.treatment=betweenScatter/(nGroups-1),
             mean.sq.residuals=withinScatter/(nSamples-nGroups),
             df.treatment=nGroups-1, df.residuals=nSamples-nGroups,
             F.statistic=F, p.value=p
             )
}

bartlett_real <- function(x, g, ...) {
  g <- factor(g[OK])
  k <- nlevels(g)
  if (k < 2)
    stop("all observations are in the same group")
  x <- split(x, g)
  n <- sapply(x, "length") - 1
  if (any(n <= 0))
    stop("there must be at least 2 observations in each group")
  v <- sapply(x, "var")
  n.total <- sum(n)
  v.total <- sum(n * v)/n.total
  STATISTIC <- ((n.total * log(v.total) - sum(n * log(v)))/(1 + (sum(1/n) - 1/n.total)/(3 * (k - 1))))
  PARAMETER <- k - 1
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
  names(STATISTIC) <- "Bartlett's K-squared"
  names(PARAMETER) <- "df"
  RVAL <- list(statistic = STATISTIC, parameter = PARAMETER, p.value = PVAL, data.name = DNAME, method = "Bartlett test of homogeneity of variances")
}

