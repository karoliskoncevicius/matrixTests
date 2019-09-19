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
#' the Brown–Forsythe test uses the median instead of the mean in computing the
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
  force(x)
  force(g)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_vec_length(g, ncol(x))

  bad <- is.na(g)
  if(any(bad)) {
    warning(sum(bad), ' columns dropped due to missing group information')
    x <- x[,!bad, drop=FALSE]
    g <- g[!bad]
  }

  g <- as.character(g)

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(g)))
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(unique(g))) {
    tmpx <- x[,g==unique(g)[i], drop=FALSE]
    tmpx <- abs(tmpx - rowMeans(tmpx, na.rm=TRUE))
    x[,g==unique(g)[i]] <- tmpx
    mPerGroup[,i] <- rowMeans(tmpx, na.rm=TRUE)
    nPerGroup[,i] <- rep.int(ncol(tmpx), nrow(tmpx)) - matrixStats::rowCounts(is.na(tmpx))
    vPerGroup[,i] <- rowSums((tmpx-mPerGroup[,i])^2, na.rm=TRUE) / (nPerGroup[,i]-1)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)
  M <- rowMeans(x, na.rm=TRUE)

  betweenScatter <- rowSums(nPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  withinScatter  <- rowSums((nPerGroup-1) * vPerGroup, na.rm=TRUE)
  withinScatter[withinScatter <= 10 * .Machine$double.eps] <- 0

  dft <- nGroups-1
  dfr <- nSamples-nGroups

  F <- (betweenScatter/dft) / (withinScatter/dfr)
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & all(nPerGroup < 3)
  showWarning(w2, 'had no groups with more than 2 observations')

  w3 <- !w1 & !w2 & withinScatter==0
  showWarning(w3, 'had zero within group variance of absolute residuals from the mean')

  F[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3] <- NA

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
  force(x)
  force(g)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_vec_length(g, ncol(x))

  bad <- is.na(g)
  if(any(bad)) {
    warning(sum(bad), ' columns dropped due to missing group information')
    x <- x[,!bad, drop=FALSE]
    g <- g[!bad]
  }

  g <- as.character(g)

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(g)))
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(unique(g))) {
    tmpx <- x[,g==unique(g)[i], drop=FALSE]
    tmpx <- abs(tmpx - matrixStats::rowMedians(tmpx, na.rm=TRUE))
    x[,g==unique(g)[i]] <- tmpx
    mPerGroup[,i] <- rowMeans(tmpx, na.rm=TRUE)
    nPerGroup[,i] <- rep.int(ncol(tmpx), nrow(tmpx)) - matrixStats::rowCounts(is.na(tmpx))
    vPerGroup[,i] <- rowSums((tmpx-mPerGroup[,i])^2, na.rm=TRUE) / (nPerGroup[,i]-1)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)
  M <- rowMeans(x, na.rm=TRUE)

  betweenScatter <- rowSums(nPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  withinScatter  <- rowSums((nPerGroup-1) * vPerGroup, na.rm=TRUE)
  withinScatter[withinScatter <= 10 * .Machine$double.eps] <- 0

  dft <- nGroups-1
  dfr <- nSamples-nGroups

  F <- (betweenScatter/dft) / (withinScatter/dfr)
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & all(nPerGroup < 3)
  showWarning(w2, 'had no groups with more than 2 observations')

  w3 <- !w1 & !w2 & withinScatter==0
  showWarning(w3, 'had zero within group variance of absolute residuals from the median')

  F[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3] <- NA

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

