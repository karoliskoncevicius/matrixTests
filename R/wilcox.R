#' Wilcoxon test
#'
#' Performs a Wilcoxon test on each row/column of the input matrix.
#'
#' Functions to perform one sample and two sample Wilcoxon tests on rows/columns of matrices.
#' Main arguments and results were intentionally matched to the \code{wilcox.test()}
#' function from default stats package. Other arguments were split into separate
#' functions:
#'
#' \code{row_wilcoxon_onesample(x)} - one sample Wilcoxon test on rows.
#' \code{col_wilcoxon_onesample(x)} - one sample Wilcoxon test on columns.
#'
#' Results should be the same as running \code{wilcox.test(x)}
#' on every row (or column) of \code{x}.
#'
#' \code{row_wilcoxon_twosample(x, y)} - two sample Wilcoxon test on rows.
#' \code{col_wilcoxon_twosample(x, y)} - two sample Wilcoxon test on columns.
#'
#' Results should be the same as running \code{wilcox.test(x, y)}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' \code{row_wilcoxon_paired(x, y)} - two sample paired Wilcoxon test on rows.
#' \code{col_wilcoxon_paired(x, y)} - two sample paired Wilcoxon test on columns.
#'
#' Results should be the same as running \code{wilcox.test(x, y, paired=TRUE)}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' By default if ‘exact’ argument is set to 'NA', exact p-values are computed
#' only if both 'x' and 'y' contain less than 50 values and there are no
#' ties. Single sample and paired tests have additional requirement of not
#' having zeroe values (values equal to null hypothesis location argument 'mu').
#' Otherwise, a normal approximation is used. Be wary of using 'exact=TRUE' on
#' large sample sizes as computations can take a very long time.
#'
#' 'correct' argument controls the continuity correction of p-values but only
#' when exact p-values cannot be computed and normal approximation is used.
#' For cases where exact p-values are returned 'correct' is switched to FALSE.
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param null true values of the location shift for the null hypothesis.
#' A single number or numeric vector with values for each observation.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with values for each observation.
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param exact logical or NA (default) indicator whether an exact p-value
#' should be computed (see Details).
#' A single value or a logical vector with values for each observation.
#' @param correct logical indicator whether continuity correction should be
#' applied in the cases where p-values are obtained using normal approximation.
#' A single value or logical vector with values for each observation.
#'
#' @return a data.frame where each row contains the results of a wilcoxon test
#' performed on the corresponding row/column of x.
#' The columns will vary depending on the type of test performed.\cr\cr
#' They will contain a subset of the following information:\cr
#' 1. obs.x - number of x observations\cr
#' 2. obs.y - number of y observations\cr
#' 3. obs.tot - total number of observations\cr
#' 4. obs.paired - number of paired observations (present in x and y)\cr
#' 5. statistic - Wilcoxon test statistic\cr
#' 6. pvalue - p-value\cr
#' 7. location.null - location shift of the null hypothesis\cr
#' 8. alternative - chosen alternative hypothesis\cr
#' 9. exact - indicates if exact p-value was computed\cr
#' 10. correct - indicates if continuity correction was performed
#'
#' @note Confidence interval and pseudo-median calculations are not implemented.
#'
#' @seealso \code{wilcox.test()}
#'
#' @examples
#' X <- iris[iris$Species=="setosa",1:4]
#' Y <- iris[iris$Species=="virginica",1:4]
#' col_wilcoxon_twosample(X, Y)
#'
#' # same row using different alternative hypotheses
#' col_wilcoxon_twosample(X[,c(1,1,1)], Y[,c(1,1,1)], alternative=c("t", "g", "l"))
#'
#' @author Karolis Koncevičius
#' @name wilcoxon
#' @export
row_wilcoxon_twosample <- function(x, y, null=0, alternative="two.sided",
                                   exact=NA, correct=TRUE
                                   ) {
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
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep(null, length.out=nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_open_interval(null, -Inf, Inf)

  if(length(exact)==1)
    exact <- rep(exact, length.out=nrow(x))
  assert_logical_vec_length(exact, 1, nrow(x))
  assert_all_in_set(exact, c(TRUE, FALSE, NA))

  if(length(correct)==1)
    correct <- rep(correct, length.out=nrow(x))
  assert_logical_vec_length(correct, 1, nrow(x))
  assert_all_in_set(correct, c(TRUE, FALSE))


  nxs  <- rep.int(ncol(x), nrow(x)) - matrixStats::rowCounts(is.na(x))
  nys  <- rep.int(ncol(y), nrow(y)) - matrixStats::rowCounts(is.na(y))

  naexact <- is.na(exact)
  exact[naexact] <- (nxs[naexact] < 50) & (nys[naexact] < 50)

  r <- matrixStats::rowRanks(cbind(x - null, y), ties.method="average")

  statistic <- rowSums(r[,seq_len(ncol(x)),drop=FALSE], na.rm=TRUE) - nxs * (nxs + 1)*0.5

  nties   <- rowTies(r)
  hasties <- rowSums(nties>0) > 0

  wres <- rep(NA_integer_, nrow(x))
  inds <- exact & !hasties
  wres[inds]  <- do_wilcox_2_exact(statistic[inds], nxs[inds], nys[inds], alternative[inds])
  wres[!inds] <- do_wilcox_2_approx(statistic[!inds], nxs[!inds], nys[!inds], alternative[!inds],
                                  nties[!inds,,drop=FALSE], correct[!inds])


  w1 <- nxs < 1
  showWarning(w1, 'had less than 1 "x" observation')

  w2 <- nys < 1
  showWarning(w2, 'had less than 1 "y" observation')

  w3 <- exact & hasties
  showWarning(w3, 'had ties: cannot compute exact p-values with ties')

  statistic[w1 | w2] <- NA
  wres[w1 | w2]      <- NA

  exact <- exact & !hasties
  correct <- correct & !exact


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxs+nys, statistic=statistic,
             pvalue=wres, location.null=null, alternative=alternative,
             exact=exact, corrected=correct,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname wilcoxon
#' @export
col_wilcoxon_twosample <- function(x, y, null=0, alternative="two.sided",
                                   exact=NA, correct=TRUE) {
  row_wilcoxon_twosample(t(x), t(y), null=null, alternative=alternative,
                         exact=exact, correct=correct
                         )
}


#' @rdname wilcoxon
#' @export
row_wilcoxon_onesample <- function(x, null=0, alternative="two.sided",
                                   exact=NA, correct=TRUE
                                   ) {
  is.null(x)

  if (is.vector(x))
    x <- matrix(x, nrow=1)

  if (is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep(null, length.out=nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_open_interval(null, -Inf, Inf)

  if(length(exact)==1)
    exact <- rep(exact, length.out=nrow(x))
  assert_logical_vec_length(exact, 1, nrow(x))
  assert_all_in_set(exact, c(TRUE, FALSE, NA))

  if(length(correct)==1)
    correct <- rep(correct, length.out=nrow(x))
  assert_logical_vec_length(correct, 1, nrow(x))
  assert_all_in_set(correct, c(TRUE, FALSE))


  x <- x - null

  haszeroes <- x==0
  x[haszeroes] <- NA
  haszeroes <- rowSums(haszeroes, na.rm=TRUE) > 0

  nxs <- rep.int(ncol(x), nrow(x)) - matrixStats::rowCounts(is.na(x))

  naexact <- is.na(exact)
  exact[naexact] <- nxs[naexact] < 50

  r <- matrixStats::rowRanks(abs(x), ties.method="average")
  rtmp <- r
  rtmp[x<=0] <- NA
  statistic <- rowSums(rtmp, na.rm=TRUE)

  nties   <- rowTies(r)
  hasties <- rowSums(nties>0) > 0

  wres <- rep(NA_integer_, nrow(x))
  inds <- exact & !hasties & !haszeroes
  wres[inds]  <- do_wilcox_1_exact(statistic[inds], nxs[inds], alternative[inds])
  wres[!inds] <- do_wilcox_1_approx(statistic[!inds], nxs[!inds], alternative[!inds],
                                    nties[!inds,,drop=FALSE], correct[!inds])


  w1 <- haszeroes
  showWarning(w1, 'had observations with "x" equal "null" that were removed')

  w2 <- nxs < 1
  showWarning(w2, 'had less than 1 "x" observation')

  w3 <- exact & haszeroes
  showWarning(w3, 'had zeroes: cannot compute exact p-values with zeroes')

  w4 <- !w3 & exact & hasties
  showWarning(w4, 'had ties: cannot compute exact p-values with ties')

  statistic[w2] <- NA
  wres[w2]      <- NA

  exact <- exact & !hasties & !haszeroes
  correct <- correct & !exact


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs=nxs, statistic=statistic, pvalue=wres, location.null=null,
             alternative=alternative, exact=exact, corrected=correct,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname wilcoxon
#' @export
col_wilcoxon_onesample <- function(x, null=0, alternative="two.sided",
                                   exact=NA, correct=TRUE) {
  row_wilcoxon_onesample(t(x), null=null, alternative=alternative,
                         exact=exact, correct=correct
                         )
}


#' @rdname wilcoxon
#' @export
row_wilcoxon_paired <- function(x, y, null=0, alternative="two.sided",
                                exact=NA, correct=TRUE
                                ) {
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
  assert_equal_ncol(x, y)

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep(null, length.out=nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_open_interval(null, -Inf, Inf)

  if(length(exact)==1)
    exact <- rep(exact, length.out=nrow(x))
  assert_logical_vec_length(exact, 1, nrow(x))
  assert_all_in_set(exact, c(TRUE, FALSE, NA))

  if(length(correct)==1)
    correct <- rep(correct, length.out=nrow(x))
  assert_logical_vec_length(correct, 1, nrow(x))
  assert_all_in_set(correct, c(TRUE, FALSE))


  xy <- x - y
  xy <- xy - null

  haszeroes <- xy==0
  xy[haszeroes] <- NA
  haszeroes <- rowSums(haszeroes, na.rm=TRUE) > 0

  nxs  <- rep.int(ncol(x), nrow(x)) - matrixStats::rowCounts(is.na(x))
  nys  <- rep.int(ncol(y), nrow(y)) - matrixStats::rowCounts(is.na(y))
  nxys <- rep.int(ncol(xy), nrow(xy)) - matrixStats::rowCounts(is.na(xy))

  naexact <- is.na(exact)
  exact[naexact] <- nxys[naexact] < 50

  r <- matrixStats::rowRanks(abs(xy), ties.method="average")
  rtmp <- r
  rtmp[xy<=0] <- NA
  statistic <- rowSums(rtmp, na.rm=TRUE)

  nties   <- rowTies(r)
  hasties <- rowSums(nties>0) > 0

  wres <- rep(NA_integer_, nrow(x))
  inds <- exact & !hasties & !haszeroes
  wres[inds]  <- do_wilcox_1_exact(statistic[inds], nxys[inds], alternative[inds])
  wres[!inds] <- do_wilcox_1_approx(statistic[!inds], nxys[!inds], alternative[!inds],
                                    nties[!inds,,drop=FALSE], correct[!inds])


  w1 <- haszeroes
  showWarning(w1, 'had observations with "x-y" equal "null" that were removed')

  w2 <- nxys < 1
  showWarning(w2, 'had less than 1 paired observation')

  w3 <- exact & haszeroes
  showWarning(w3, 'had zeroes: cannot compute exact p-values with zeroes')

  w4 <- !w3 & exact & hasties
  showWarning(w4, 'had ties: cannot compute exact p-values with ties')

  statistic[w2] <- NA
  wres[w2]      <- NA

  exact <- exact & !hasties & !haszeroes
  correct <- correct & !exact


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.paired=nxys, statistic=statistic,
             pvalue=wres, location.null=null, alternative=alternative,
             exact=exact, corrected=correct,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname wilcoxon
#' @export
col_wilcoxon_paired <- function(x, y, null=0, alternative="two.sided",
                                   exact=NA, correct=TRUE) {
  row_wilcoxon_paired(t(x), t(y), null=null, alternative=alternative,
                      exact=exact, correct=correct
                      )
}


