#' ONEWAY ANOVA
#'
#' Performs an analysis of variance tests on each row/column of the input matrix.
#'
#' Functions to perform ONEWAY ANOVA analysis for rows/columns of matrices.
#'
#' \code{row_oneway_equalvar(x, g)} - oneway ANOVA on rows.
#' \code{col_oneway_equalvar(x, g)} - oneway ANOVA on columns.
#'
#' Results should be the same as running \code{aov(x ~ g)}
#' on every row (or column) of \code{x}
#'
#' \code{row_oneway_welch(x, g)} - oneway ANOVA with Welch correction on rows.
#' \code{col_oneway_welch(x, g)} - oneway ANOVA with Welch correction on columns.
#'
#' Results should be the same as running \code{oneway.test(x, g, var.equal=FALSE)}
#' on every row (or column) of \code{x}
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.

#' @return a data.frame where each row contains the results of an oneway anova
#' test performed on the corresponding row/column of x.
#' The columns will vary depending on the type of test performed.\cr\cr
#' They will contain a subset of the following information:\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 3. sumsq.between - between group (treatment) sum of squares\cr
#' 4. sumsq.within - within group (residual) sum of squares\cr
#' 5. meansq.between - between group mean squares\cr
#' 6. meansq.within - within group mean squares\cr
#' 7. df.between - between group (treatment) degrees of freedom\cr
#' 8. df.within - within group (residual) degrees of freedom\cr
#' 9. statistic - F statistic\cr
#' 10. pvalue - p.value
#'
#' @seealso \code{aov()}, \code{oneway.test()}
#'
#' @examples
#' col_oneway_welch(iris[,1:4], iris$Species)
#' row_oneway_equalvar(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name oneway
#' @export
row_oneway_equalvar <- function(x, g) {
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
    nPerGroup[,i] <- matrixStats::rowCounts(!is.na(tmpx))
    mPerGroup[,i] <- rowMeans(tmpx, na.rm=TRUE)
    vPerGroup[,i] <- rowSums((tmpx-mPerGroup[,i])^2, na.rm=TRUE) / (nPerGroup[,i]-1)
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


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & nGroups==nSamples
  showWarning(w2, 'had one observation per group')

  w3 <- !w1 & !w2 & withinScatter==0 & betweenScatter==0
  showWarning(w3, 'had essentially constant values')

  w4 <- !w1 & !w2 & withinScatter==0 & betweenScatter!=0
  showWarning(w4, 'had zero within group variance: result might be unreliable')

  F[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups,
             sumsq.between=betweenScatter, sumsq.within=withinScatter,
             meansq.between=betweenScatter/dft, meansq.within=withinScatter/dfr,
             df.between=dft, df.within=dfr, statistic=F, pvalue=p,
             row.names=rnames
             )
}

#' @rdname oneway
#' @export
col_oneway_equalvar <- function(x, g) {
  row_oneway_equalvar(t(x), g)
}

#' @rdname oneway
#' @export
row_oneway_welch <- function(x, g) {
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
    nPerGroup[,i] <- matrixStats::rowCounts(!is.na(tmpx))
    mPerGroup[,i] <- rowMeans(tmpx, na.rm=TRUE)
    vPerGroup[,i] <- rowSums((tmpx-mPerGroup[,i])^2, na.rm=TRUE) / (nPerGroup[,i]-1)
  }
  mPerGroup[nPerGroup<2] <- NA
  vPerGroup[nPerGroup<2] <- NA
  nPerGroup[nPerGroup<2] <- NA

  nSamples  <- rowSums(nPerGroup, na.rm=TRUE)
  nGroups   <- matrixStats::rowCounts(nPerGroup!=0, na.rm=TRUE)
  wPerGroup <- nPerGroup/vPerGroup
  wTot      <- rowSums(wPerGroup, na.rm=TRUE)
  M         <- rowSums((wPerGroup * mPerGroup)/wTot, na.rm=TRUE)

  betweenScatter <- rowSums(wPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  tmp <- rowSums((1 - wPerGroup/wTot)^2/(nPerGroup - 1), na.rm=TRUE) / (nGroups*nGroups - 1)
  dfr <- 1/(3*tmp)
  dft <- nGroups-1

  F <- betweenScatter / (dft * (1 + 2 * (nGroups-2) * tmp))
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & nGroups < length(unique(g))
  showWarning(w2, 'had groups with less than 2 observations: those groups were removed')

  w3 <- !w1 & rowSums(vPerGroup!=0, na.rm=TRUE)==0
  showWarning(w3, 'had zero variance in all of the groups')

  w4 <- !w1 & !w3 & rowSums(vPerGroup==0, na.rm=TRUE) > 0
  showWarning(w4, 'had groups with zero variance: result might be unreliable')


  F[w1 | w3] <- NA
  p[w1 | w3] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, df.between=dft,
             df.within=dfr, statistic=F, pvalue=p,
             row.names=rnames
             )
}

#' @rdname oneway
#' @export
col_oneway_welch <- function(x, g) {
  row_oneway_welch(t(x), g)
}
