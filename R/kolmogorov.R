#' Kolmogorov-Smirnov test
#'
#' Performs a Kolmogorov-Smirnov test on each row/column of the input matrix.
#'
#' Function to perform two sample Kolmogorov-Smirnov test on rows/columns of
#' matrices. Main arguments and results were intentionally matched to the
#' \code{ks.test()} function from default stats package.
#'
#' Results should be the same as running \code{ks.test(x, y)} on every row (or
#' column) of \code{x} and \code{y}.
#'
#' By default if ‘exact’ argument is set to 'NA', exact p-values are computed
#' if the product of 'x' and 'y' sample sizes is less than 10000. Otherwise,
#' asymptotic distributions are used.
#'
#' Alternative hypothesis setting specifies null and alternative hypotheses.
#' The possible values of 'two sided', 'less', and 'greater'.
#' 'two sided' sets the null hypothesis for the distributions of 'x' being equal to the distribution 'y'.
#' 'less' sets the null hypothesis for the distribution of x not being less than the distribution of y.
#' 'greater' sets the null hypothesis for the distribution of x not being greater than the distribution of y.
#' See \code{help(ks.test)} for more details.
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with values for each observation.
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param exact logical or NA (default) indicator whether an exact p-value should be computed (see Details).
#' A single value or a logical vector with values for each observation.
#'
#' @return a data.frame where each row contains the results of a Kolmogorov-Smirnov test
#' performed on the corresponding row/column of x and y.
#' Each row contains the following information (in order):\cr
#' 1. obs.x - number of x observations\cr
#' 2. obs.y - number of y observations\cr
#' 3. obs.tot - total number of observations\cr
#' 5. statistic - Wilcoxon test statistic\cr
#' 6. pvalue - p-value\cr
#' 8. alternative - chosen alternative hypothesis\cr
#' 9. exact - indicates if exact p-value was computed\cr
#'
#' @seealso \code{ks.test()}
#'
#' @examples
#' X <- iris[iris$Species=="setosa", 1:4]
#' Y <- iris[iris$Species=="virginica", 1:4]
#' col_kolmogorovsmirnov_twosample(X, Y)
#'
#' # same column using different alternative hypotheses
#' col_kolmogorovsmirnov_twosample(X[,c(1,1,1)], Y[,c(1,1,1)], alternative=c("t", "g", "l"))
#'
#' @author Karolis Koncevičius
#' @name kolmogorov
#' @export
row_kolmogorovsmirnov_twosample <- function(x, y, alternative="two.sided", exact=NA) {
  is.null(x)
  is.null(y)

  if (is.vector(x))
    x <- matrix(x, nrow=1)
  if(is.vector(y))
    y <- matrix(y, nrow=1)

  if (is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)
  if(is.data.frame(y) && all(sapply(y, is.numeric)))
    y <- data.matrix(y)

  assert_numeric_mat_or_vec(x)
  assert_numeric_mat_or_vec(y)

  if(nrow(y)==1 & nrow(x)>1) {
    y <- matrix(y, nrow=nrow(x), ncol=ncol(y), byrow=TRUE)
  }

  assert_equal_nrow(x, y)

  if(length(alternative)==1)
    alternative <- rep.int(alternative, nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(exact)==1)
    exact <- rep.int(exact, nrow(x))
  assert_logical_vec_length(exact, 1, nrow(x))
  assert_all_in_set(exact, c(TRUE, FALSE, NA))


  nxs <- as.numeric(ncol(x) - matrixStats::rowCounts(x, value=NA))
  nys <- as.numeric(ncol(y) - matrixStats::rowCounts(y, value=NA))

  naexact <- is.na(exact)
  exact[naexact] <- (nxs[naexact] * nys[naexact]) < 10000

  w <- cbind(x, y)
  z <- col(w)
  z[is.na(w)] <- NA

  # order w and z by row
  ord <- order(row(w), w)
  w   <- matrix(w[ord], nrow=nrow(w), byrow=TRUE)
  z   <- matrix(z[ord], nrow=nrow(w), byrow=TRUE)

  inds <- z <= ncol(x)
  z[which(inds)]  <- rep(1/nxs, ncol(inds))[which(inds)]    # replace values by row
  z[which(!inds)] <- rep(-1/nys, ncol(inds))[which(!inds)]  # replace values by row
  # z <- rowCumsums(z)  # NOTE: precision problems which cause bugs
  z <- t(apply(z, 1, cumsum))

  w[is.nan(w)] <- NA  # so that NaNs will not show ties error
  diffs   <- matrixStats::rowDiffs(w)
  ties    <- cbind(diffs == 0 | is.nan(diffs), rep(FALSE, nrow(diffs)))  # NOTE: is.nan needed to show ties with infinities
  hasties <- matrixStats::rowAnys(ties)
  z[which(ties)] <- NA

  statistic <- rep.int(NA_real_, nrow(x))
  p         <- rep.int(NA_real_, nrow(x))

  inds <- alternative == "two.sided" & nxs > 0 & nys > 0  # NOTE: because psmirnove gives an error otherwise
  if(any(inds)) {
    statistic[inds] <- matrixStats::rowMaxs(abs(z[inds,,drop=FALSE]), na.rm=TRUE)
    p[inds]         <- mapply(stats::psmirnov, statistic[inds], asplit(cbind(nxs[inds], nys[inds]), 1),
                              asplit(w[inds,,drop=FALSE], 1), exact[inds], two.sided = TRUE, lower.tail = FALSE)
  }

  inds <- alternative == "greater" & nxs > 0 & nys > 0
  if(any(inds)) {
    statistic[inds] <- matrixStats::rowMaxs(z[inds,,drop=FALSE], na.rm=TRUE)
    p[inds]         <- mapply(stats::psmirnov, statistic[inds], asplit(cbind(nxs[inds], nys[inds]), 1),
                              asplit(w[inds,,drop=FALSE], 1), exact[inds], two.sided = FALSE, lower.tail = FALSE)
  }

  inds <- alternative == "less" & nxs > 0 & nys > 0
  if(any(inds)) {
    statistic[inds] <- -matrixStats::rowMins(z[inds,,drop=FALSE], na.rm=TRUE)
    p[inds]         <- mapply(stats::psmirnov, statistic[inds], asplit(cbind(nxs[inds], nys[inds]), 1),
                              asplit(w[inds,,drop=FALSE], 1), exact[inds], two.sided = FALSE, lower.tail = FALSE)
  }

  # NOTE: not sure if necessary
  # p <- pmin(1, pmax(0, p))


  w1 <- nxs < 1
  showWarning(w1, 'kolmogorovsmirnov_twosample', 'had less than 1 "x" observation')

  w2 <- nys < 1
  showWarning(w2, 'kolmogorovsmirnov_twosample', 'had less than 1 "y" observation')

  w3 <- !exact & hasties
  showWarning(w3, 'kolmogorovsmirnov_twosample', 'had ties: asymptotic p-values will be approximate')

  statistic[w1 | w2] <- NA
  p[w1 | w2]         <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxs+nys, statistic=statistic,
             pvalue=p, alternative=alternative, exact=exact,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname kolmogorov
#' @export
col_kolmogorovsmirnov_twosample <- function(x, y, alternative="two.sided", exact=NA) {
  row_kolmogorovsmirnov_twosample(t(x), t(y), alternative=alternative, exact=exact)
}
