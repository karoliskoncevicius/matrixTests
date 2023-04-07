#' Correlation
#'
#' Performs a correlation test on each row/column of a the input matrix.
#'
#' Functions to perform various correlation tests for rows/columns of matrices.
#' Main arguments and results were intentionally matched to the \code{cor.test()}
#' function from default stats package.
#'
#' \code{row_cor_pearson(x, y)} - test for Pearson correlation on rows.
#' \code{col_cor_pearson(x, y)} - test for Pearson correlation on columns.
#'
#' Results should be the same as running \code{cor.test(x, y, method="pearson")}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with value for each observation.
#' Must be one of "two.sided" (default), "greater" or "less".
#' @param conf.level confidence levels used for the confidence intervals.
#' A single number or a numeric vector with value for each observation.
#' All values must be in the range of [0;1] or NA.
#'
#' @return a data.frame where each row contains the results of a correlation
#' test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.paired - number of paired observations (present in x and y)\cr
#' 2. cor - estimated correlation coefficient\cr
#' 3. df - degrees of freedom\cr
#' 4. statistic - t statistic\cr
#' 5. pvalue - p-value\cr
#' 6. conf.low - lower confidence interval\cr
#' 7. conf.high - higher confidence interval\cr
#' 8. alternative - chosen alternative hypothesis\cr
#' 9. cor.null - correlation of the null hypothesis (=0)\cr
#' 10. conf.level - chosen confidence level
#'
#' @note
#' For a marked increase in computation speed turn off the calculation of
#' confidence interval by setting \code{conf.level} to NA.
#'
#' @seealso \code{cor.test()}
#'
#' @examples
#' X <- iris[iris$Species=="setosa",1:4]
#' Y <- iris[iris$Species=="virginica",1:4]
#' col_cor_pearson(X, Y)
#' row_cor_pearson(t(X), t(Y))
#'
#' @author Karolis Koncevičius
#' @name cortest
#' @export
row_cor_pearson <- function(x, y, alternative="two.sided", conf.level=0.95) {
  is.null(x)
  is.null(y)

  if(is.vector(x))
    x <- matrix(x, nrow=1)
  if(is.vector(y))
    y <- matrix(y, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)
  if(is.data.frame(y) && all(sapply(y, is.numeric)))
    y <- data.matrix(y)

  assert_numeric_mat_or_vec(x)
  assert_numeric_mat_or_vec(y)

  if(nrow(y)==1 & nrow(x)>1) {
    y <- matrix(y, nrow=nrow(x), ncol=ncol(y), byrow=TRUE)
  }

  assert_equal_nrow(x, y)
  assert_equal_ncol(x, y)

  if(length(alternative)==1)
    alternative <- rep.int(alternative, nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(all(is.na(conf.level)))
    conf.level[] <- NA_real_
  if(length(conf.level)==1)
    conf.level <- rep.int(conf.level, nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_closed_interval(conf.level, 0, 1, na.allow=TRUE)


  mu <- rep.int(0, nrow(x)) # can't be changed because different test should be used in that case.

  isna <- is.na(x+y)
  x[isna] <- NA
  y[isna] <- NA

  ns <- ncol(x) - matrixStats::rowCounts(isna)

  x  <- x - rowMeans(x, na.rm=TRUE)
  y  <- y - rowMeans(y, na.rm=TRUE)
  sx <- sqrt(rowSums(x * x, na.rm=TRUE) / (ns-1))
  sy <- sqrt(rowSums(y * y, na.rm=TRUE) / (ns-1))

  rs <- rowSums(x * y, na.rm=TRUE) / (sx*sy*(ns-1))
  rs[abs(rs - 1) < .Machine$double.eps^0.5] <- 1  # if not different from 1 use 1
  rs[abs(rs + 1) < .Machine$double.eps^0.5] <- -1 # if not different from -1 use -1
  df <- ns-2

  pres <- do_pearson(rs, df, alternative, conf.level)

  w1 <- ns < 3
  showWarning(w1, 'cor_pearson', 'had less than 3 complete observations')

  w2 <- !w1 & sx==0
  showWarning(w2, 'cor_pearson', 'had zero standard deviation in x')

  w3 <- !w1 & sy==0
  showWarning(w3, 'cor_pearson', 'had zero standard deviation in y')

  w4 <- !w2 & !w3 & ns == 3
  showWarning(w4, 'cor_pearson', 'had exactly 3 complete observations: no confidence intervals produced')

  w5 <- !w1 & abs(rs)==1
  showWarning(w5, 'cor_pearson', 'had essentially perfect fit: results might be unreliable for small sample sizes')

  pres[w4, 3:4]       <- NA
  df[w1 | w2 | w3]    <- NA
  pres[w1 | w2 | w3,] <- NA
  rs[w1 | w2 | w3]    <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.paired=ns, cor=rs, df=df, statistic=pres[,1],
             pvalue=pres[,2], conf.low=pres[,3], conf.high=pres[,4],
             alternative=alternative, cor.null=mu, conf.level=conf.level,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname cortest
#' @export
col_cor_pearson <- function(x, y, alternative="two.sided", conf.level=0.95) {
  row_cor_pearson(t(x), t(y), alternative=alternative, conf.level=conf.level)
}

