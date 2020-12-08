#' Fligner-Killeen test
#'
#' Performs the Fligner-Killeen test of homogeneity of variances (with median
#' centering of the groups) on each row/column of the input matrix.
#'
#' NA values are always ommited. If values are missing for a whole group - that
#' group is discarded. Groups with only one observation are also discarded.
#'
#' \code{row_flignerkilleen(x, g)} - Fligner-Killeen test on rows.
#' \code{col_flignerkilleen(x, g)} - Fligner-Killeen test on columns.
#'
#' Results should be the same as as running \code{fligner.test(x, g)}
#' on every row (or column) of \code{x}.
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.

#' @return a data.frame where each row contains the results of the
#' Fligner-Killeen test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 3. df - degrees of freedom\cr
#' 4. statistic - squared statistic\cr
#' 5. pvalue - p-value
#'
#' @seealso \code{fligner.test()}
#'
#' @examples
#' col_flignerkilleen(iris[,1:4], iris$Species)
#' row_flignerkilleen(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name fligner
#' @export
row_flignerkilleen <- function(x, g) {
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

  g <- as.character(g)


  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(g)))
  for(i in seq_along(unique(g))) {
    inds <- g==unique(g)[i]
    tmpx <- x[,inds, drop=FALSE]
    nPerGroup[,i] <- rep.int(ncol(tmpx), nrow(tmpx)) - matrixStats::rowCounts(is.na(tmpx))
    x[,inds] <- tmpx - matrixStats::rowMedians(tmpx, na.rm=TRUE)
  }

  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)
  nSamples <- rowSums(nPerGroup, na.rm=TRUE)

  a <- stats::qnorm((1 + matrixStats::rowRanks(abs(x), ties.method="average") / (nSamples+1)) / 2)
  a <- matrix(a, nrow=nrow(x), ncol=ncol(x))

  mPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(g)))
  for(i in seq_along(unique(g))) {
    mPerGroup[,i] <- rowSums(a[,g==unique(g)[i], drop=FALSE], na.rm=TRUE)
  }

  df   <- nGroups-1
  stat <- rowSums(mPerGroup*mPerGroup / nPerGroup, na.rm=TRUE)
  stat <- (stat - nSamples * rowMeans(a, na.rm=TRUE)^2) / rowVars(a, na.rm=TRUE)
  p    <- stats::pchisq(stat, df, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & matrixStats::rowAlls(nPerGroup < 2)
  showWarning(w2, 'had one observation per group')

  w3 <- !w1 & !w2 & matrixStats::rowAlls(x==0 | is.na(x))
  showWarning(w3, 'had zero variance in all of the groups')

  df[w1 | w2 | w3]   <- NA
  stat[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3]    <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, df=df, statistic=stat,
             pvalue=p, row.names=rnames
             )
}

#' @rdname fligner
#' @export
col_flignerkilleen <- function(x, g) {
  row_flignerkilleen(t(x), g)
}

