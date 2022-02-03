#' Anderson-Darling test
#'
#' Performs Anderson-Darling goodness of fit test for normality.
#'
#' \code{row_andersondarling(x)} - Anderson-Darling test on rows.
#' \code{col_andersondarling(x)} - Anderson-Darling test on columns.
#'
#' Results should be the same as running \code{nortest::ad.test(x)}
#' on every row (or column) of \code{x}
#'
#' @param x numeric matrix.
#'
#' @return a data.frame where each row contains the results of Anderson-Darling
#' test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs - number of observations\cr
#' 5. statistic - test statistic\cr
#' 6. pvalue - p-value
#'
#' @seealso \code{shapiro.test()}
#'
#' @examples
#' col_andersondarling(iris[,1:4])
#' row_andersondarling(t(iris[,1:4]))
#'
#' @author Karolis Koncevičius
#' @name andersondarling
#' @export
row_andersondarling <- function(x) {
  is.null(x)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)


  xi <- matrix(x[order(row(x), x)], ncol=ncol(x), byrow=TRUE)  # order by row (increasing)
  xd <- matrix(x[order(row(x), -x)], ncol=ncol(x), byrow=TRUE) # order by row (decreasing)
  n <- ncol(xi) - matrixStats::rowCounts(x, value=NA)

  m  <- rowMeans(x, na.rm=TRUE) # means and sds are the same for x, xi, and xd
  s  <- sqrt(rowVars(x, na.rm=TRUE))
  lg <- pnorm((xi-m)/s, log.p=TRUE) + pnorm(-(xd-m)/s, log.p=TRUE)

  h <- (2 * col(x) - 1) * lg
  A <- -n - rowMeans(h, na.rm=TRUE)
  AA <- (1 + 0.75/n + 2.25/n^2) * A

  p <- rep(3.7e-24, nrow(x))
  inds <- AA < 0.2
  p[inds] <- 1 - exp(-13.436 + 101.14 * AA[inds] - 223.73 * AA[inds]^2)
  inds <- AA >=0.2 & AA < 0.34
  p[inds] <- 1 - exp(-8.318 + 42.796 * AA[inds] - 59.938 * AA[inds]^2)
  inds <- AA >=0.34 & AA < 0.6
  p[inds] <- exp(0.9177 - 4.279 * AA[inds] - 1.38 * AA[inds]^2)
  inds <- AA >=0.6 & AA < 10
  p[inds] <- exp(1.2937 - 5.709 * AA[inds] + 0.0186 * AA[inds]^2)


  w1 <- n < 8
  showWarning(w1, 'had less than 8 total observations')

  w2 <- !w1 & s <= 0
  showWarning(w2, 'had essentially constant values')

  A[w1 | w2] <- NA
  p[w1 | w2] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs=n, statistic=A, pvalue=p, row.names=rnames)
}

#' @rdname andersondarling
#' @export
col_andersondarling <- function(x) {
  row_andersondarling(t(x))
}
