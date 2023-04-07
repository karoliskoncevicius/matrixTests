#' Cosinor
#'
#' Performs a Cosinor test for periodicity on each row/column of the input matrix.
#'
#' \code{row_cosinor} - cosinor test on rows.
#' \code{col_cosinor} - cosinor test on columns.
#'
#' @param x numeric matrix.
#' @param t a vector specifying time variable for each observation of x.
#' @param period oscillation period in the units of \code{t} (default = 24, suitable when inspecting diurnal rhythms with hourly data).
#'
#' @return a data.frame where each row contains the results of a cosinor test
#' performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs - total number of observations\cr
#' 2. mesor - "Midline Estimating Statistic Of Rhythm" - the average value around which the variable oscillates\cr
#' 3. amplitude - difference between mesor and the peak of the rhythm\cr
#' 4. acrophase - time when rhythm reaches its peak\cr
#' 5. rsquared - R-squared\cr
#' 6. df.model - model terms degrees of freedom\cr
#' 7. df.residual - residual degrees of freedom\cr
#' 8. statistic - F statistic for the omnibus test against intercept-only model\cr
#' 9. pvalue - p-value\cr
#' 10. period - the period used within the model\cr
#'
#' @seealso \code{\link[cosinor]{cosinor.lm}}
#'
#' @examples
#  # sinus wave with Gaussian noise
#' wave <- sin(2*pi*1:24/24) + rnorm(24)
#' row_cosinor(wave, 1:24, 24)
#'
#' @author Karolis Koncevičius
#' @name cosinor
#' @export
row_cosinor <- function(x, t, period=24) {
  is.null(x)
  is.null(t)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)
  assert_numeric_vec_length(t, ncol(x))
  assert_numeric_vec_length(period, 1)
  assert_all_in_open_interval(period, 0, Inf)

  if(anyNA(t)) {
    bad <- is.na(t)
    x   <- x[,!bad, drop=FALSE]
    t   <- t[!bad]
    warning(sum(bad), ' columns dropped due to missing time information')
  }

  period <- rep.int(period, min(1, nrow(x)))

  hasinfx <- is.infinite(x)
  x[hasinfx] <- NA
  hasinfx <- rowSums(hasinfx) > 0


  nobs  <- rowSums(!is.na(x))

  B <- matrix(nrow = ncol(x), ncol = 3)
  B[,1] <- rep.int(1, ncol(x))
  B[,2] <- if(length(period) > 0) sinpi(2*t/period) else NA
  B[,3] <- if(length(period) > 0) cospi(2*t/period) else NA

  res <- do_regression(x, B)

  mesor     <- res$betas[1,]
  phase     <- atan2(res$betas[2,], res$betas[3,])
  acrophase <- (phase / (2*pi) * period + period) %% period
  amplitude <- sqrt(res$betas[2,]^2 + res$betas[3,]^2)


  w1 <- hasinfx
  showWarning(w1, 'cosinor', 'had infinite observations that were removed')

  w2 <- nobs < 3
  showWarning(w2, 'cosinor', 'had less than 3 complete observations: no p-values produced, amplitude and acrophase will be unreliable')

  w3 <- nobs == 3
  showWarning(w3, 'cosinor', 'had exactly 3 complete observations: no p-values produced')

  w4 <- !w2 & !w3 & res$stats$sstot == 0
  showWarning(w4, 'cosinor', 'had essentially constant values')

  w5 <- !w2 & res$stats$dfmod == 0
  showWarning(w5, 'cosinor', 'had only 1 unique timepoint within the specified period: no p-values produced, amplitude and acrophase will be unreliable')

  w6 <- !w2 & res$stats$dfmod == 1
  showWarning(w6, 'cosinor', 'had only 2 unique timepoints within the specified period: amplitude and acrophase will be unreliable')

  w7 <- !w2 & !w3 & !w4 & res$stats$rsq == 1
  showWarning(w7, 'cosinor', 'had essentially perfect fit')

  res$stats[w2 | w3 | w4 | w5, c("dfmod","dfres","f","p")] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs=nobs, mesor=mesor, amplitude=amplitude, acrophase=acrophase,
             rsquared=res$stats$rsq, df.model=res$stats$dfmod,
             df.residual=res$stats$dfres, statistic=res$stats$f,
             pvalue=res$stats$p, period=period, row.names=rnames
             )
}


#' @rdname cosinor
#' @export
col_cosinor <- function(x, t, period=24) {
  row_cosinor(t(x), t, period)
}

