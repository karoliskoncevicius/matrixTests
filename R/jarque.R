#' Jarque-Bera test
#'
#' Performs a Jarque-Bera goodness of fit test for normality.
#'
#' \code{row_jarquebera(x)} - Jarque-Bera test on rows.
#' \code{col_jarquebera(x)} - Jarque-Bera test on columns.
#'
#' Results should be the same as running \code{moments::jarque.test(x)}
#' on every row (or column) of \code{x}
#'
#' @param x numeric matrix.
#'
#' @return a data.frame where each row contains the results of a Jarque-Bera
#' test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs - number of observations\cr
#' 2. skewness - skewness\cr
#' 3. kurtosis - kurtosis\cr
#' 4. df - degrees of freedom\cr
#' 5. statistic - chi-squared statistic\cr
#' 6. pvalue - p-value
#'
#' @seealso \code{shapiro.test()}
#'
#' @examples
#' col_jarquebera(iris[,1:4])
#' row_jarquebera(t(iris[,1:4]))
#'
#' @author Karolis Koncevičius
#' @name jarquebera
#' @export
row_jarquebera <- function(x) {
  is.null(x)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)


  n    <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  m0   <- x - rowMeans(x, na.rm=TRUE)
  m2   <- rowMeans(m0*m0, na.rm=TRUE)
  m3   <- rowMeans(m0^3, na.rm=TRUE)
  m4   <- rowMeans(m0^4, na.rm=TRUE)
  skew <- (m3 / m2^(1.5))
  kurt <- (m4 / (m2*m2))

  df   <- rep.int(2, length(n))
  ksq  <- n * skew^2/6 + n * (kurt-3)^2/24
  p    <- 1 - stats::pchisq(ksq, df=df)


  w1 <- n < 2
  showWarning(w1, 'had less than 2 total observations')

  w2 <- !w1 & matrixStats::rowAlls(m0, value=0, na.rm=TRUE)
  showWarning(w2, 'had essentially constant values')

  df[w1 | w2]  <- NA
  ksq[w1 | w2] <- NA
  p[w1 | w2]   <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs=n, skewness=skew, kurtosis=kurt, df=df,
             statistic=ksq, pvalue=p, row.names=rnames
             )
}

#' @rdname jarquebera
#' @export
col_jarquebera <- function(x) {
  row_jarquebera(t(x))
}
