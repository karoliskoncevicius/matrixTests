#' F Variance test
#'
#' Performs the F test of equality of variances for two normal populations on
#' each row/column of the two input matrices.
#'
#' NA values are always ommited.
#'
#' \code{row_f_var(x, y)} - F-test for variance on rows.
#' \code{col_f_var(x, y)} - F-test for variance on columns.
#'
#' Results should be the same as as running \code{var.test(x, y)}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param null - hypothesized 'x' and 'y' variance ratio.
#' A single number or numeric vector with values for each observation.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with values for each observation.
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param conf.level confidence levels used for the confidence intervals.
#' A single number or a numeric vector with values for each observation.
#' All values must be in the range of [0:1].
#'
#' @return a data.frame where each row contains the results of the F variance
#' test performed on the corresponding row/column of x and y.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.x - number of x observations\cr
#' 2. obs.y - number of y observations\cr
#' 3. obs.tot - total number of observations\cr
#' 4. var.x - variance of x\cr
#' 5. var.y - variance of y\cr
#' 6. var.ratio - x/y variance ratio\cr
#' 7. df.num - numerator degrees of freedom\cr
#' 8. df.denom - denominator degrees of freedom\cr
#' 9. statistic - F statistic\cr
#' 10 pvalue - p-value\cr
#' 11. conf.low - lower bound of the confidence interval\cr
#' 12. conf.high - higher bound of the confidence interval\cr
#' 13. ratio.null - variance ratio of the null hypothesis\cr
#' 14. alternative - chosen alternative hypothesis\cr
#' 15. conf.level - chosen confidence level
#'
#' @seealso \code{var.test()}
#'
#' @examples
#' X <- iris[iris$Species=="setosa",1:4]
#' Y <- iris[iris$Species=="virginica",1:4]
#' col_f_var(X, Y)
#'
#' @author Karolis Koncevičius
#' @name fvar
#' @export
row_f_var <- function(x, y, null=1, alternative="two.sided", conf.level=0.95) {
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

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep(null, length.out=nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_open_interval(null, 0, Inf)

  if(length(conf.level)==1)
    conf.level <- rep(conf.level, length.out=nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_closed_interval(conf.level, 0, 1)


  hasinfx <- is.infinite(x)
  x[hasinfx] <- NA
  hasinfx <- rowSums(hasinfx) > 0

  hasinfy <- is.infinite(y)
  y[hasinfy] <- NA
  hasinfy <- rowSums(hasinfy) > 0

  nxs  <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  nys  <- ncol(y) - matrixStats::rowCounts(y, value=NA)
  nxys <- nxs + nys

  vxs <- rowVars(x, n=nxs, na.rm=TRUE)
  vys <- rowVars(y, n=nys, na.rm=TRUE)

  estimate  <- vxs/vys

  dfx <- nxs - 1
  dfy <- nys - 1

  fres <- do_ftest(estimate, null, alternative, dfx, dfy, conf.level)


  w1 <- hasinfx
  showWarning(w1, 'had infinite "x" observations that were removed')

  w2 <- hasinfy
  showWarning(w2, 'had infinite "y" observations that were removed')

  w3 <- nxs <= 1
  showWarning(w3, 'had less than 1 "x" observation')

  w4 <- nys <= 1
  showWarning(w4, 'had less than 1 "y" observation')

  w5 <- vxs == 0
  showWarning(w5, 'had zero variance in "x"')

  w6 <- vys == 0
  showWarning(w6, 'had zero variance in "y"')


  dfx[w3 | w4 | (w5 & w6)]   <- NA
  dfy[w3 | w4 | (w5 & w6)]   <- NA
  fres[w3 | w4 | (w5 & w6),] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxys, var.x=vxs, var.y=vys,
             var.ratio=estimate, df.num=dfx, df.denom=dfy, statistic=fres[,1],
             pvalue=fres[,2], conf.low=fres[,3], conf.high=fres[,4],
             ratio.null=null, alternative=alternative, conf.level=conf.level,
             row.names=rnames, stringsAsFactors=FALSE
             )
}

#' @rdname fvar
#' @export
col_f_var <- function(x, y, null=1, alternative="two.sided", conf.level=0.95) {
  row_f_var(t(x), t(y), null, alternative, conf.level)
}

