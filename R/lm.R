#' F-test for nested linear regression models
#'
#' Performs an F-test comparing two nested linear regression models on each row/column of the input matrix.
#'
#' \code{row_lm_f} - F-test for linear regression on rows.
#' \code{col_lm_f} - F-test for linear regression on columns.
#'
#' @param x numeric matrix.
#' @param m - a model matrix for a linear regression model to be tested.
#' @param null - a null model which the original model will be compared against (default = intercept-only).
#'
#' @return a data.frame where each row contains the results of an F-test
#' performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs - total number of observations\cr
#' 2. betaN - estimated coefficients (as many as there are column in m)\cr
#' 3. rsquared.model - R-squared of the full model\cr
#' 4. rsquared.null - R-squared of the null model\cr
#' 5. df.model - model terms degrees of freedom\cr
#' 6. df.residual - residual degrees of freedom\cr
#' 7. statistic - F statistic\cr
#' 8. pvalue - p-value\cr
#'
#' @seealso \code{lm()}, \code{anova()}
#'
#' @examples
#' X   <- t(iris[,-5])
#' mod <- model.matrix(~ iris$Species)
#' row_lm_f(X, mod)
#'
#' # with specified null model
#' X <- mtcars[,c("mpg", "disp", "qsec")]
#' mod  <- model.matrix(~ mtcars$drat + mtcars$cyl + mtcars$wt + mtcars$vs)
#' mod0 <- model.matrix(~               mtcars$cyl + mtcars$wt + mtcars$vs)
#' col_lm_f(X, mod, mod0)
#'
#' @author Karolis Koncevičius
#' @name linearmodel
#' @export
row_lm_f <- function(x, m, null=stats::model.matrix(~ 1, data=data.frame(seq_len(nrow(m))))) {
  is.null(x)
  is.null(m)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)
  assert_numeric_mat(m)
  assert_numeric_mat(null)
  assert_ncol_equal_nrow(x, m)
  assert_equal_nrow(m, null)
  assert_unique_colnames(m)
  assert_unique_colnames(null)
  assert_nested_model(null, m)
  assert_all_in_open_interval(m, -Inf, Inf)

  hasinfx <- is.infinite(x)
  x[hasinfx] <- NA
  hasinfx <- rowSums(hasinfx) > 0


  nobs <- rowSums(!is.na(x))

  res1 <- do_regression(x, m)
  res0 <- do_regression(x, null)

  df    <- res1$stats$dfmod - res0$stats$dfmod
  dfres <- res1$stats$dfres

  statistic <- (res0$stats$ssres - res1$stats$ssres) / df
  statistic <- statistic / (res1$stats$ssres / dfres)
  p <- stats::pf(statistic, df, dfres, lower.tail=FALSE)


  # TODO: gather all the warnings
  w1 <- hasinfx
  showWarning(w1, 'had infinite observations that were removed')


  # TODO: decide what to do about beta values
  rownames(res1$betas) <- paste0("beta.", 1:nrow(res1$betas)-1)

  # TODO: add partial r-squared (eta squared?) to the output
  # TODO: maybe also return RSS values

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs=nobs, t(res1$betas),
             rsquared.model=res1$stats$rsq, rsquared.null=res0$stats$rsq,
             df.model=df, df.residual=dfres, statistic=statistic, pvalue=p,
             row.names=rnames
             )
}


#' @rdname linearmodel
#' @export
col_lm_f <- function(x, m, null=stats::model.matrix(~ 1, data=data.frame(seq_len(nrow(m))))) {
  row_lm_f(t(x), m, null)
}

