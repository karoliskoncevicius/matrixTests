#' Levene test
#'
#' Levene's test and Brown-Forsythe test for equality of variances
#' between groups on each row/column of the input matrix.
#'
#' NA values are always ommited.
#' If values are missing for a whole group - that group is discarded.
#'
#' \code{row_levene(x, g)} - Levene's test on rows.
#' \code{col_levene(x, g)} - Levene's test on columns.
#'
#' \code{row_brownforsythe(x, g)} - Brown-Forsythe test on rows.
#' \code{col_brownforsythe(x, g)} - Brown-Forsythe test on columns.
#'
#' @note Difference between Levene's test and Brown-Forsythe test is that
#' the Brown-Forsythe test uses the median instead of the mean in computing the
#' spread within each group. Many software implementations use the name
#' "Levene's test" for both variants.
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.
#'
#' @return a data.frame where each row contains the results of the Levene's test
#' performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 3. df.between - between group (treatment) degrees of freedom\cr
#' 4. df.within - within group (residual) degrees of freedom\cr
#' 5. statistic - F statistic\cr
#' 6. pvalue - p.value
#'
#' @seealso \code{\link[car]{leveneTest}}
#'
#' @examples
#' col_levene(iris[,1:4], iris$Species)
#' row_brownforsythe(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name levene
#' @export
row_levene <- function(x, g) {
  is.null(x)
  is.null(g)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_vec_length(g, ncol(x))

  if(anyNA(g)) {
    bad <- is.na(g)
    x   <- x[,!bad, drop=FALSE]
    g   <- g[!bad]
    warning(sum(bad), ' columns dropped due to missing group information')
  }

  g  <- as.character(g)
  gs <- unique(g)

  hasinfx <- is.infinite(x)
  x[hasinfx] <- NA
  hasinfx <- rowSums(hasinfx) > 0

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(gs))
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(gs)) {
    tmpx <- x[,g==gs[i], drop=FALSE]
    tmpx <- abs(tmpx - rowMeans(tmpx, na.rm=TRUE))
    x[,g==gs[i]] <- tmpx
    nPerGroup[,i] <- ncol(tmpx) - matrixStats::rowCounts(tmpx, value=NA)
    mPerGroup[,i] <- rowMeans(tmpx, na.rm=TRUE)
    vPerGroup[,i] <- rowVars(tmpx, n=nPerGroup[,i], m=mPerGroup[,i], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)
  M <- rowMeans(x, na.rm=TRUE)

  betweenScatter <- rowSums(nPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  withinScatter  <- rowSums((nPerGroup-1) * vPerGroup, na.rm=TRUE)

  dft <- nGroups-1
  dfr <- nSamples-nGroups

  F <- (betweenScatter/dft) / (withinScatter/dfr)
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)


  w1 <- hasinfx
  showWarning(w1, 'levene', 'had infinite observations that were removed')

  w2 <- nGroups < 2
  showWarning(w2, 'levene', 'had less than 2 groups with enough observations')

  w3 <- !w2 & all(nPerGroup < 3)
  showWarning(w3, 'levene', 'had no groups with at least 3 observations')

  w4 <- !w2 & !w3 & withinScatter==0
  showWarning(w4, 'levene', 'had zero within group variance of absolute residuals from the mean')

  w5 <- !w2 & !w3 & !w4 & withinScatter <= .Machine$double.eps
  showWarning(w5, 'levene', 'had essentially constant absolute residuals from the mean: results might be unreliable')

  dft[w2 | w3 | w4] <- NA
  dfr[w2 | w3 | w4] <- NA
  F[w2 | w3 | w4] <- NA
  p[w2 | w3 | w4] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups,
             df.between=dft, df.within=dfr, statistic=F, pvalue=p,
             row.names=rnames
             )

}

#' @rdname levene
#' @export
col_levene <- function(x, g) {
  row_levene(t(x), g)
}


#' @rdname levene
#' @export
row_brownforsythe <- function(x, g) {
  is.null(x)
  is.null(g)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_vec_length(g, ncol(x))


  if(anyNA(g)) {
    bad <- is.na(g)
    x   <- x[,!bad, drop=FALSE]
    g   <- g[!bad]
    warning(sum(bad), ' columns dropped due to missing group information')
  }

  g  <- as.character(g)
  gs <- unique(g)

  hasinfx <- is.infinite(x)
  x[hasinfx] <- NA
  hasinfx <- rowSums(hasinfx) > 0

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(gs))
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(gs)) {
    tmpx <- x[,g==gs[i], drop=FALSE]
    tmpx <- abs(tmpx - matrixStats::rowMedians(tmpx, na.rm=TRUE))
    x[,g==gs[i]] <- tmpx
    nPerGroup[,i] <- ncol(tmpx) - matrixStats::rowCounts(tmpx, value=NA)
    mPerGroup[,i] <- rowMeans(tmpx, na.rm=TRUE)
    vPerGroup[,i] <- rowVars(tmpx, n=nPerGroup[,i], m=mPerGroup[,i], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)
  M <- rowMeans(x, na.rm=TRUE)

  betweenScatter <- rowSums(nPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  withinScatter  <- rowSums((nPerGroup-1) * vPerGroup, na.rm=TRUE)

  dft <- nGroups-1
  dfr <- nSamples-nGroups

  F <- (betweenScatter/dft) / (withinScatter/dfr)
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)

  w1 <- hasinfx
  showWarning(w1, 'brownforsythe', 'had infinite observations that were removed')

  w2 <- nGroups < 2
  showWarning(w2, 'brownforsythe', 'had less than 2 groups with enough observations')

  w3 <- !w2 & all(nPerGroup < 3)
  showWarning(w3, 'brownforsythe', 'had no groups with at least 3 observations')

  w4 <- !w2 & !w3 & withinScatter==0
  showWarning(w4, 'brownforsythe', 'had zero within group variance of absolute residuals from the median')

  w5 <- !w2 & !w3 & !w4 & withinScatter <= .Machine$double.eps
  showWarning(w5, 'brownforsythe', 'had essentially constant absolute residuals from the median: results might be unreliable')

  dft[w2 | w3 | w4] <- NA
  dfr[w2 | w3 | w4] <- NA
  F[w2 | w3 | w4]   <- NA
  p[w2 | w3 | w4]   <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups,
             df.between=dft, df.within=dfr, statistic=F, pvalue=p,
             row.names=rnames
             )

}

#' @rdname levene
#' @export
col_brownforsythe <- function(x, g) {
  row_brownforsythe(t(x), g)
}

