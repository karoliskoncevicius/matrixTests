#' t-test
#'
#' Performs a t-test on each row of a the input matrix.
#'
#' Functions to perform one sample and two sample t-tests for rows of matrices.
#' Main arguments and results were intentionally matched to the \code{t.test()}
#' function from default stats package. Other arguments were split into separate
#' functions:
#'
#' \code{ttest_onegroup()} - t-test for mean of a single group. Same as \code{t.test(x)}
#'
#' \code{ttest_equalvar()} - groups have equal variance. Same as \code{t.test(x, y, var.equal=TRUE)}
#'
#' \code{ttest_welch()} - Welch approximation. Same as \code{t.test(x, y)}
#'
#' \code{ttest_paired()} - paired t-test. Same as \code{t.test(x, y, paired=TRUE)}
#'
#' @name ttest
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
#' ttest_welch(X, Y)
#'
#' # same row using different confidence levels
#' ttest_equalvar(X[c(1,1,1),], Y[c(1,1,1),], conf.level=c(0.9, 0.95, 0.99))
#'
#' @author Karolis Koncevičius
#' @export
ttest_onegroup <- function(x, alternative="two.sided", mu=0, conf.level=0.95) {
  force(x)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(mu)==1)
    mu <- rep(mu, length.out=nrow(x))
  assert_numeric_vec_length(mu, 1, nrow(x))
  assert_all_in_range(mu, -Inf, Inf)

  if(length(conf.level)==1)
    conf.level <- rep(conf.level, length.out=nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_range(conf.level, 0, 1)

  mxs <- rowMeans(x, na.rm=TRUE)

  nxs <- rowSums(!is.na(x))

  bad <- nxs < 2
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows had less than 2 "x" observations'))

  vxs    <- ifelse(bad, NA, rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1))
  dfs    <- ifelse(bad, NA, nxs-1)
  stders <- ifelse(bad, NA, sqrt(vxs/nxs))

  bad <- stders <= 10 * .Machine$double.eps * abs(mxs)
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows were essentially constant'))

  dfs[bad]      <- NA
  stders[bad]   <- NA


  tres <- do_ttest(mxs, mu, stders, alternative, dfs, conf.level)

  data.frame(mean.x=mxs, var.x=vxs, obs.x=nxs, t.statistic=tres$t.statistic,
             p.value=tres$p.value, ci.low=tres$ci.low, ci.high=tres$ci.high,
             stderr=stders, df=dfs, mean.null=mu, conf.level=conf.level,
             alternative=alternative, stringsAsFactors=FALSE,
             row.names=rownames(x)
             )
}



#' @export
#' @rdname ttest
ttest_equalvar <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
  force(x)
  force(y)

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
  assert_equal_nrow(x, y)

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(mu)==1)
    mu <- rep(mu, length.out=nrow(x))
  assert_numeric_vec_length(mu, 1, nrow(x))
  assert_all_in_range(mu, -Inf, Inf)

  if(length(conf.level)==1)
    conf.level <- rep(conf.level, length.out=nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_range(conf.level, 0, 1)


  mxs  <- rowMeans(x, na.rm=TRUE)
  mys  <- rowMeans(y, na.rm=TRUE)
  mxys <- mxs - mys

  nxs  <- rowSums(!is.na(x))
  nys  <- rowSums(!is.na(y))
  nxys <- nxs + nys

  bad <- (nxs+nys) < 3
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows had less than 3 total observations'))

  dfs <- ifelse(bad, NA, nxs + nys - 2)

  bad <- (nxs+nys) > 2 & nys < 1
  if(any(bad, na.rm=TRUE)) {
    warning(paste0(sum(bad), ' of the rows had zero "y" observations'))
  }

  vys <- ifelse(bad, NA, rowSums((y-mys)^2, na.rm=TRUE) / (nys-1))
  dfs[bad] <- NA

  bad <- (nxs+nys) > 2 & nxs < 1
  if(any(bad, na.rm=TRUE)) {
    warning(paste0(sum(bad), ' of the rows had zero "x" observations'))
  }

  vxs <- ifelse(bad, NA, rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1))
  dfs[bad] <- NA


  vs   <- rep(0, nrow(x))
  vs   <- ifelse(nxs > 1, vs + (nxs-1) * vxs, vs)
  vs   <- ifelse(nys > 1, vs + (nys-1) * vys, vs)
  vs   <- vs/dfs
  stders <- sqrt(vs * (1/nxs + 1/nys))


  bad <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), abs(mys))
  if(any(bad, na.rm=TRUE)) {
    warning(paste0(sum(bad), ' of the rows were essentially constant'))
  }

  dfs[bad]      <- NA
  stders[bad]   <- NA


  tres <- do_ttest(mxys, mu, stders, alternative, dfs, conf.level)

  data.frame(mean.x=mxs, mean.y=mys, mean.diff=mxys, var.x=vxs, var.y=vys,
             var.pooled=vs, obs.x=nxs, obs.y=nys, obs.tot=nxys,
             t.statistic=tres$t.statistic, p.value=tres$p.value,
             ci.low=tres$ci.low, ci.high=tres$ci.high, stderr=stders, df=dfs,
             mean.null=mu, conf.level=conf.level, alternative=alternative,
             stringsAsFactors=FALSE,
             row.names=rownames(x)
             )
}


#' @export
#' @rdname ttest
ttest_welch <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
  force(x)
  force(y)

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
  assert_equal_nrow(x, y)

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(mu)==1)
    mu <- rep(mu, length.out=nrow(x))
  assert_numeric_vec_length(mu, 1, nrow(x))
  assert_all_in_range(mu, -Inf, Inf)

  if(length(conf.level)==1)
    conf.level <- rep(conf.level, length.out=nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_range(conf.level, 0, 1)


  mxs <- rowMeans(x, na.rm=TRUE)
  mys <- rowMeans(y, na.rm=TRUE)
  mxys    <- mxs - mys

  nxs <- rowSums(!is.na(x))
  nys <- rowSums(!is.na(y))
  nxys    <- nxs + nys

  bad <- nys < 2
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows had less than 2 "y" observations'))

  vys <- ifelse(bad, NA, rowSums((y-mys)^2, na.rm=TRUE) / (nys-1))

  bad <- nxs < 2
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows had less than 2 "x" observations'))

  vxs <- ifelse(bad, NA, rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1))


  stderxs <- sqrt(vxs/nxs)
  stderys <- sqrt(vys/nys)
  stders  <- sqrt(stderxs^2 + stderys^2)
  dfs     <- stders^4/(stderxs^4/(nxs - 1) + stderys^4/(nys - 1))


  bad <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), abs(mys))
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows were essentially constant'))

  dfs[bad]    <- NA
  stders[bad] <- NA


  tres <- do_ttest(mxys, mu, stders, alternative, dfs, conf.level)

  data.frame(mean.x=mxs, mean.y=mys, mean.diff=mxys, var.x=vxs, var.y=vys,
             obs.x=nxs, obs.y=nys, obs.tot=nxys, t.statistic=tres$t.statistic,
             p.value=tres$p.value, ci.low=tres$ci.low, ci.high=tres$ci.high,
             stderr=stders, df=dfs, mean.null=mu, conf.level=conf.level,
             alternative=alternative, stringsAsFactors=FALSE,
             row.names=rownames(x)
             )
}

#' @export
#' @rdname ttest
ttest_paired <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
  force(x)
  force(y)

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
  assert_equal_nrow(x, y)
  assert_equal_ncol(x, y)

  if(length(alternative)==1)
    alternative <- rep(alternative, length.out=nrow(x))
  assert_character_vec_length(alternative, 1, nrow(x))

  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
  assert_all_in_set(alternative, choices)

  if(length(mu)==1)
    mu <- rep(mu, length.out=nrow(x))
  assert_numeric_vec_length(mu, 1, nrow(x))
  assert_all_in_range(mu, -Inf, Inf)

  if(length(conf.level)==1)
    conf.level <- rep(conf.level, length.out=nrow(x))
  assert_numeric_vec_length(conf.level, 1, nrow(x))
  assert_all_in_range(conf.level, 0, 1)


  xy <- x-y

  mxs  <- rowMeans(x, na.rm=TRUE)
  mys  <- rowMeans(y, na.rm=TRUE)
  mxys <- rowMeans(xy, na.rm=TRUE)

  nxs  <- rowSums(!is.na(x))
  nys  <- rowSums(!is.na(y))
  nxys <- rowSums(!is.na(xy))

  dfs <- nxys-1

  bad <- nxys < 2
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows had less than 2 paired observations'))

  dfs  <- ifelse(bad, NA, nxys-1)
  vxys <- ifelse(bad, NA, rowSums((xy-mxys)^2, na.rm=TRUE) / (nxys-1))

  vxs <- ifelse(nxs < 2, NA, rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1))
  vys <- ifelse(nys < 2, NA, rowSums((y-mys)^2, na.rm=TRUE) / (nys-1))

  stders <- sqrt(vxys/nxys)

  bad <- stders <= 10 * .Machine$double.eps * abs(mxys)
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows were essentially constant'))

  dfs[bad]    <- NA
  stders[bad] <- NA


  tres <- do_ttest(mxys, mu, stders, alternative, dfs, conf.level)

  data.frame(mean.x=mxs, mean.y=mys, mean.diff=mxys, var.x=vxs, var.y=vys,
             var.diff=vxys, obs.x=nxs, obs.y=nys, obs.pair=nxys,
             t.statistic=tres$t.statistic, p.value=tres$p.value,
             ci.low=tres$ci.low, ci.high=tres$ci.high, stderr=stders, df=dfs,
             mean.null=mu, conf.level=conf.level, alternative=alternative,
             stringsAsFactors=FALSE, row.names=rownames(x)
             )
}

