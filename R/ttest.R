#' t-test
#'
#' Performs a t-test on each row/column of the input matrix.
#'
#' Functions to perform one sample and two sample t-tests for rows/columns of matrices.
#' Main arguments and results were intentionally matched to the \code{t.test()}
#' function from default stats package. Other arguments were split into separate
#' functions:
#'
#' \code{row_t_onesample(x)} - one sample t-test on rows.
#' \code{col_t_onesample(x)} - one sample t-test on columns.
#'
#' Results should be the same as running \code{t.test(x)}
#' on every row (or column) of \code{x}.
#'
#' \code{row_t_equalvar(x, y)} - two sample equal variance t-test on rows.
#' \code{col_t_equalvar(x, y)} - two sample equal variance t-test on columns.
#'
#' Results should be the same as running \code{t.test(x, y, var.equal=TRUE)}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' \code{row_t_welch(x, y)} - two sample t-test with Welch correction on rows.
#' \code{col_t_welch(x, y)} - two sample t-test with Welch correction on columns.
#'
#' Results should be the same as running \code{t.test(x, y)}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' \code{row_t_paired(x, y)} - two sample paired t-test on rows.
#' \code{col_t_paired(x, y)} - two sample paired t-test on columns.
#'
#' Results should be the same as running \code{t.test(x, y, paired=TRUE)}
#' on every row (or column) of \code{x} and \code{y}.
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param null true values of the means for the null hypothesis.
#' A single number or numeric vector with values for each observation.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with values for each observation.
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param conf.level confidence levels used for the confidence intervals.
#' A single number or a numeric vector with values for each observation.
#' All values must be in the range of [0:1] or NA.
#'
#' @return a data.frame where each row contains the results of a t.test
#' performed on the corresponding row/column of x.
#' The columns will vary depending on the type of test performed.\cr\cr
#' They will contain a subset of the following information:\cr
#' 1. obs.x - number of x observations\cr
#' 2. obs.y - number of y observations\cr
#' 3. obs.tot - total number of observations\cr
#' 4. obs.paired - number of paired observations (present in x and y)\cr
#' 5. mean.x - mean estiamte of x\cr
#' 6. mean.y - mean estiamte of y\cr
#' 7. mean.diff - mean estiamte of x-y difference\cr
#' 8. var.x - variance estiamte of x\cr
#' 9. var.y - variance estiamte of y\cr
#' 10. var.diff - variance estiamte of x-y difference\cr
#' 11. var.pooled - pooled variance estimate of x and y\cr
#' 12. stderr - standard error\cr
#' 13. df - degrees of freedom\cr
#' 14. statistic - t statistic\cr
#' 15. pvalue - p-value\cr
#' 16. conf.low - lower bound of the confidence interval\cr
#' 17. conf.high - higher bound of the confidence interval\cr
#' 18. mean.null - mean of the null hypothesis\cr
#' 19. alternative - chosen alternative hypothesis\cr
#' 20. conf.level - chosen confidence level
#'
#' @note
#' For a marked increase in computation speed turn off the calculation of
#' confidence interval by setting \code{conf.level} to NA.
#'
#' @seealso \code{t.test()}
#'
#' @examples
#' X <- iris[iris$Species=="setosa",1:4]
#' Y <- iris[iris$Species=="virginica",1:4]
#' col_t_welch(X, Y)
#'
#' # same row using different confidence levels
#' col_t_equalvar(X[,c(1,1,1)], Y[,c(1,1,1)], conf.level=c(0.9, 0.95, 0.99))
#'
#' @author Karolis Koncevičius
#' @name ttest
#' @export
row_t_equalvar <- function(x, y, null=0, alternative="two.sided", conf.level=0.95) {
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
    alternative <- rep.int(alternative, nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep.int(null, nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_closed_interval(null, -Inf, Inf)

  if(all(is.na(conf.level)))
    conf.level[] <- NA_real_
  if(length(conf.level)==1)
    conf.level <- rep.int(conf.level, nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_closed_interval(conf.level, 0, 1, na.allow=TRUE)


  mxs  <- rowMeans(x, na.rm=TRUE)
  mys  <- rowMeans(y, na.rm=TRUE)
  mxys <- mxs - mys

  nxs  <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  nys  <- ncol(y) - matrixStats::rowCounts(y, value=NA)
  nxys <- nxs + nys

  vxs <- rowVars(x, n=nxs, m=mxs, na.rm=TRUE)
  vys <- rowVars(y, n=nys, m=mys, na.rm=TRUE)

  dfs  <- nxs + nys - 2

  vs <- rep.int(0, nrow(x))
  inds <- nxs > 1
  vs[inds] <- vs[inds] + (nxs[inds]-1) * vxs[inds]
  inds <- nys > 1
  vs[inds] <- vs[inds] + (nys[inds]-1) * vys[inds]
  vs <- vs/dfs

  stders <- sqrt(vs * (1/nxs + 1/nys))

  tres <- do_ttest(mxys, null, stders, alternative, dfs, conf.level)


  w1 <- nxys < 3
  showWarning(w1, 't_equalvar', 'had less than 3 total observations')

  w2 <- !w1 & nxs < 1
  showWarning(w2, 't_equalvar', 'had zero "x" observations')

  w3 <- !w1 & nys < 1
  showWarning(w3, 't_equalvar', 'had zero "y" observations')

  w4 <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), abs(mys))
  showWarning(w4, 't_equalvar', 'had essentially constant values')


  stders[w1 | w2 | w3 | w4] <- NA
  dfs[w1 | w2 | w3 | w4]    <- NA
  tres[w1 | w2 | w3 | w4,]  <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxys, mean.x=mxs, mean.y=mys,
             mean.diff=mxys, var.x=vxs, var.y=vys, var.pooled=vs,
             stderr=stders, df=dfs, statistic=tres[,1], pvalue=tres[,2],
             conf.low=tres[,3], conf.high=tres[,4], mean.null=null,
             alternative=alternative, conf.level=conf.level,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname ttest
#' @export
col_t_equalvar <- function(x, y, null=0, alternative="two.sided", conf.level=0.95) {
  row_t_equalvar(t(x), t(y), null=null, alternative=alternative, conf.level=conf.level)
}


#' @rdname ttest
#' @export
row_t_welch <- function(x, y, null=0, alternative="two.sided", conf.level=0.95) {
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
    alternative <- rep.int(alternative, nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep.int(null, nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_closed_interval(null, -Inf, Inf)

  if(all(is.na(conf.level)))
    conf.level[] <- NA_real_
  if(length(conf.level)==1)
    conf.level <- rep.int(conf.level, nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_closed_interval(conf.level, 0, 1, na.allow=TRUE)


  mxs  <- rowMeans(x, na.rm=TRUE)
  mys  <- rowMeans(y, na.rm=TRUE)
  mxys <- mxs - mys

  nxs  <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  nys  <- ncol(y) - matrixStats::rowCounts(y, value=NA)
  nxys <- nxs + nys

  vxs <- rowVars(x, n=nxs, m=mxs, na.rm=TRUE)
  vys <- rowVars(y, n=nys, m=mys, na.rm=TRUE)

  stderxs <- vxs/nxs
  stderys <- vys/nys
  stders  <- stderxs + stderys
  dfs     <- stders*stders / (stderxs*stderxs/(nxs - 1) + stderys*stderys/(nys - 1))
  stders  <- sqrt(stders)

  tres <- do_ttest(mxys, null, stders, alternative, dfs, conf.level)


  w1 <- nxs < 2
  showWarning(w1, 't_welch', 'had less than 2 "x" observations')

  w2 <- !w1 & nys < 2
  showWarning(w2, 't_welch', 'had less than 2 "y" observations')

  w3 <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), abs(mys))
  showWarning(w3, 't_welch', 'had essentially constant values')

  stders[w1 | w2 | w3] <- NA
  dfs[w1 | w2 | w3]     <- NA
  tres[w1 | w2 | w3,]   <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxys, mean.x=mxs, mean.y=mys,
             mean.diff=mxys, var.x=vxs, var.y=vys, stderr=stders, df=dfs,
             statistic=tres[,1], pvalue=tres[,2], conf.low=tres[,3],
             conf.high=tres[,4], mean.null=null, alternative=alternative,
             conf.level=conf.level, stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname ttest
#' @export
col_t_welch <- function(x, y, null=0, alternative="two.sided", conf.level=0.95) {
  row_t_welch(t(x), t(y), null=null, alternative=alternative, conf.level=conf.level)
}


#' @rdname ttest
#' @export
row_t_onesample <- function(x, null=0, alternative="two.sided", conf.level=0.95) {
  is.null(x)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)


  if(length(alternative)==1)
    alternative <- rep.int(alternative, nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(null)==1)
    null <- rep.int(null, nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_closed_interval(null, -Inf, Inf)

  if(all(is.na(conf.level)))
    conf.level[] <- NA_real_
  if(length(conf.level)==1)
    conf.level <- rep.int(conf.level, nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_closed_interval(conf.level, 0, 1, na.allow=TRUE)


  nxs <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  mxs <- rowMeans(x, na.rm=TRUE)
  vxs <- rowVars(x, n=nxs, m=mxs, na.rm=TRUE)

  stders <- sqrt(vxs/nxs)
  dfs <- nxs-1

  tres <- do_ttest(mxs, null, stders, alternative, dfs, conf.level)


  w1 <- nxs < 2
  showWarning(w1, 't_onesample', 'had less than 2 "x" observations')

  w2 <- !w1 & stders <= 10 * .Machine$double.eps * abs(mxs)
  showWarning(w2, 't_onesample', 'had essentially constant values')


  stders[w1 | w2]  <- NA
  dfs[w1 | w2]     <- NA
  tres[w1 | w2,]   <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs=nxs, mean=mxs, var=vxs, stderr=stders, df=dfs,
             statistic=tres[,1], pvalue=tres[,2], conf.low=tres[,3],
             conf.high=tres[,4], mean.null=null, alternative=alternative,
             conf.level=conf.level, stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname ttest
#' @export
col_t_onesample <- function(x, null=0, alternative="two.sided", conf.level=0.95) {
  row_t_onesample(t(x), null=null, alternative=alternative, conf.level=conf.level)
}


#' @rdname ttest
#' @export
row_t_paired <- function(x, y, null=0, alternative="two.sided", conf.level=0.95) {
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

  if(length(null)==1)
    null <- rep.int(null, nrow(x))
  assert_numeric_vec_length(null, 1, nrow(x))
  assert_all_in_closed_interval(null, -Inf, Inf)

  if(all(is.na(conf.level)))
    conf.level[] <- NA_real_
  if(length(conf.level)==1)
    conf.level <- rep.int(conf.level, nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_closed_interval(conf.level, 0, 1, na.allow=TRUE)


  xy <- x-y

  mxs  <- rowMeans(x, na.rm=TRUE)
  mys  <- rowMeans(y, na.rm=TRUE)
  mxys <- rowMeans(xy, na.rm=TRUE)

  nxs  <- ncol(x)  - matrixStats::rowCounts(x, value=NA)
  nys  <- ncol(y)  - matrixStats::rowCounts(y, value=NA)
  nxys <- ncol(xy) - matrixStats::rowCounts(xy, value=NA)

  vxs  <- rowVars(x, n=nxs, m=mxs, na.rm=TRUE)
  vys  <- rowVars(y, n=nys, m=mys, na.rm=TRUE)
  vxys <- rowVars(xy, n=nxys, m=mxys, na.rm=TRUE)

  stders <- sqrt(vxys/nxys)
  dfs <- nxys-1

  tres <- do_ttest(mxys, null, stders, alternative, dfs, conf.level)


  w1 <- nxys < 2
  showWarning(w1, 't_paired', 'had less than 2 paired observations')

  w2 <- stders <= 10 * .Machine$double.eps * abs(mxys)
  showWarning(w2, 't_paired', 'had essentially constant values')

  stders[w1 | w2] <- NA
  dfs[w1 | w2]    <- NA
  tres[w1 | w2,]  <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.paired=nxys, mean.x=mxs, mean.y=mys,
             mean.diff=mxys, var.x=vxs, var.y=vys, var.diff=vxys,
             stderr=stders, df=dfs, statistic=tres[,1], pvalue=tres[,2],
             conf.low=tres[,3], conf.high=tres[,4], mean.null=null,
             alternative=alternative, conf.level=conf.level,
             stringsAsFactors=FALSE, row.names=rnames
             )
}

#' @rdname ttest
#' @export
col_t_paired <- function(x, y, null=0, alternative="two.sided", conf.level=0.95) {
  row_t_paired(t(x), t(y), null=null, alternative=alternative, conf.level=conf.level)
}
