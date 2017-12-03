#' t-test
#'
#' Performs a t-test on each row/column of a the input matrix.
#'
#' Functions to perform one sample and two sample t-tests for rows/columns of matrices.
#' Main arguments and results were intentionally matched to the \code{t.test()}
#' function from default stats package. Other arguments were split into separate
#' functions:
#'
#' \code{row.t.onesample()}, \code{col.t.onesample()}
#' - t-test for mean of a single group. Same as \code{t.test(x)}
#'
#' \code{row.t.equalvar()}, \code{col.t.equalvar()}
#' - t-test for groups with equal variance. Same as \code{t.test(x, y, var.equal=TRUE)}
#'
#' \code{row.t.welch()}, \code{col.t.welch()}
#' - t.test with Welch approximation. Same as \code{t.test(x, y)}
#'
#' \code{row.t.paired()}, \code{col.t.paired()}
#' - paired t-test. Same as \code{t.test(x, y, paired=TRUE)}
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param mu true values of the means for the null hypothesis.
#' A single number or numeric vector of length nrow(x)/ncol(x).
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector of length nrow(x)/ncol(x).
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param conf.level confidence levels used for the confidence intervals.
#' A single number or a numeric vector of length nrow(x)/ncol(x).
#' All values must be in the range of [0;1].
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
#' 9. var.x - variance estiamte of x\cr
#' 9. var.y - variance estiamte of y\cr
#' 10. var.diff - variance estiamte of x-y difference\cr
#' 11. var.pooled - pooled variance estimate of x and y\cr
#' 12. stderr - standard error\cr
#' 13. df - degrees of freedom\cr
#' 14. statistic.t - t statistic\cr
#' 15. p.value - p-value\cr
#' 16. ci.low - lower bound of the confidence interval\cr
#' 17. ci.high - higher bound of the confidence interval\cr
#' 18. alternative - chosen alternative hypothesis\cr
#' 19. correlation.null - mean of the null hypothesis\cr
#' 20. conf.level - chosen confidence level
#'
#' @seealso \code{t.test()}
#'
#' @examples
#' X <- iris[iris$Species=="setosa",1:4]
#' Y <- iris[iris$Species=="virginica",1:4]
#' col.t.welch(X, Y)
#'
#' # same row using different confidence levels
#' col.t.equalvar(X[,c(1,1,1)], Y[,c(1,1,1)], conf.level=c(0.9, 0.95, 0.99))
#'
#' @author Karolis Koncevičius
#' @name ttest
#' @export
row.t.onesample <- function(x, alternative="two.sided", mu=0, conf.level=0.95) {
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
  nxs <- matrixStats::rowCounts(!is.na(x))
  vxs <- rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1)
  dfs <- nxs-1
  stders <- sqrt(vxs/nxs)


  tres <- do_ttest(mxs, mu, stders, alternative, dfs, conf.level)


  w1 <- nxs < 2
  showWarning(w1, 'had less than 2 "x" observations')

  w2 <- !w1 & stders <= 10 * .Machine$double.eps * abs(mxs)
  showWarning(w2, 'were essentially constant')


  tres[w1 | w2,] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, mean.x=mxs, var.x=vxs, stderr=stders, df=dfs,
             statistic.t=tres[,1], p.value=tres[,2], ci.low=tres[,3],
             ci.high=tres[,4], alternative=alternative, mean.null=mu,
             conf.level=conf.level, stringsAsFactors=FALSE,
             row.names=rnames
             )
}

#' @rdname ttest
#' @export
col.t.onesample <- function(x, alternative="two.sided", mu=0, conf.level=0.95) {
  row.t.onesample(t(x), alternative=alternative, mu=mu, conf.level=conf.level)
}

#' @rdname ttest
#' @export
row.t.equalvar <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
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

  nxs  <- matrixStats::rowCounts(!is.na(x))
  nys  <- matrixStats::rowCounts(!is.na(y))
  nxys <- nxs + nys

  vxs <- rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1)
  vys <- rowSums((y-mys)^2, na.rm=TRUE) / (nys-1)

  dfs <- nxs + nys - 2

  vs <- rep(0, nrow(x))
  vs <- ifelse(nxs > 1, vs + (nxs-1) * vxs, vs)
  vs <- ifelse(nys > 1, vs + (nys-1) * vys, vs)
  vs <- vs/dfs
  stders <- sqrt(vs * (1/nxs + 1/nys))

  tres <- do_ttest(mxys, mu, stders, alternative, dfs, conf.level)


  w1 <- nxys < 3
  showWarning(w1, 'had less than 3 total observations')

  w2 <- !w1 & nxs < 1
  showWarning(w2, 'had zero "x" observations')

  w3 <- !w1 & nys < 1
  showWarning(w3, 'had zero "y" observations')

  w4 <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), abs(mys))
  showWarning(w4, 'were essentially constant')


  tres[w1 | w2 | w3 | w4,] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxys, mean.x=mxs, mean.y=mys,
             mean.diff=mxys, var.x=vxs, var.y=vys, var.pooled=vs,
             stderr=stders, df=dfs, statistic.t=tres[,1], p.value=tres[,2],
             ci.low=tres[,3], ci.high=tres[,4], alternative=alternative,
             mean.null=mu, conf.level=conf.level, stringsAsFactors=FALSE,
             row.names=rnames
             )
}


#' @rdname ttest
#' @export
col.t.equalvar <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
  row.t.equalvar(t(x), t(y), alternative=alternative, mu=mu, conf.level=conf.level)
}


#' @rdname ttest
#' @export
row.t.welch <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
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

  nxs  <- matrixStats::rowCounts(!is.na(x))
  nys  <- matrixStats::rowCounts(!is.na(y))
  nxys <- nxs + nys

  vxs <- rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1)
  vys <- rowSums((y-mys)^2, na.rm=TRUE) / (nys-1)

  stderxs <- sqrt(vxs/nxs)
  stderys <- sqrt(vys/nys)
  stders  <- sqrt(stderxs^2 + stderys^2)
  dfs     <- stders^4/(stderxs^4/(nxs - 1) + stderys^4/(nys - 1))

  tres <- do_ttest(mxys, mu, stders, alternative, dfs, conf.level)


  w1 <- nxs < 2
  showWarning(w1, 'had less than 2 "x" observations')

  w2 <- !w1 & nys < 2
  showWarning(w2, 'had less than 2 "y" observations')

  w3 <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), abs(mys))
  showWarning(w3, 'were essentially constant')

  tres[w1 | w2 | w3,] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.tot=nxys, mean.x=mxs, mean.y=mys,
             mean.diff=mxys, var.x=vxs, var.y=vys, stderr=stders, df=dfs,
             statistic.t=tres[,1], p.value=tres[,2], ci.low=tres[,3],
             ci.high=tres[,4], alternative=alternative, mean.null=mu,
             conf.level=conf.level, stringsAsFactors=FALSE,
             row.names=rnames
             )
}


#' @rdname ttest
#' @export
col.t.welch <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
  row.t.welch(t(x), t(y), alternative=alternative, mu=mu, conf.level=conf.level)
}


#' @rdname ttest
#' @export
row.t.paired <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
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

  nxs  <- matrixStats::rowCounts(!is.na(x))
  nys  <- matrixStats::rowCounts(!is.na(y))
  nxys <- matrixStats::rowCounts(!is.na(xy))

  vxs <- rowSums((x-mxs)^2, na.rm=TRUE) / (nxs-1)
  vxs[nxs < 2] <- NA
  vys <- rowSums((y-mys)^2, na.rm=TRUE) / (nys-1)
  vys[nys < 2] <- NA

  vxys <- rowSums((xy-mxys)^2, na.rm=TRUE) / (nxys-1)

  stders <- sqrt(vxys/nxys)
  dfs <- nxys-1

  tres <- do_ttest(mxys, mu, stders, alternative, dfs, conf.level)


  w1 <- nxys < 2
  showWarning(w1, 'had less than 2 paired observations')

  w2 <- stders <= 10 * .Machine$double.eps * abs(mxys)
  showWarning(w2, 'were essentially constant')

  tres[w1 | w2,] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.x=nxs, obs.y=nys, obs.paired=nxys, mean.x=mxs, mean.y=mys,
             mean.diff=mxys, var.x=vxs, var.y=vys, var.diff=vxys,
             stderr=stders, df=dfs, statistic.t=tres[,1], p.value=tres[,2],
             ci.low=tres[,3], ci.high=tres[,4], alternative=alternative,
             mean.null=mu, conf.level=conf.level, stringsAsFactors=FALSE,
             row.names=rnames
             )
}


#' @rdname ttest
#' @export
col.t.paired <- function(x, y, alternative="two.sided", mu=0, conf.level=0.95) {
  row.t.paired(t(x), t(y), alternative=alternative, mu=mu, conf.level=conf.level)
}
