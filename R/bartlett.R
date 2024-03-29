#' Bartlett test
#'
#' Performs the Bartlett's test of homogeneity of variances on each row/column
#' of the input matrix.
#'
#' NA values are always ommited. If values are missing for a whole group - that
#' group is discarded. Groups with only one observation are also discarded.
#'
#' \code{row_bartlett(x, g)} - Bartlet's test on rows.
#' \code{col_bartlett(x, g)} - Bartlet's test on columns.
#'
#' Results should be the same as as running \code{bartlett.test(x, g)}
#' on every row (or column) of \code{x}.
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.
#'
#' @return a data.frame where each row contains the results of the bartlett test
#' performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 3. var.pooled - pooled variance estimate\cr
#' 4. df - degrees of freedom\cr
#' 5. statistic - chi-squared statistic\cr
#' 6. pvalue - p-value
#'
#' @seealso \code{bartlett.test()}
#'
#' @examples
#' col_bartlett(iris[,1:4], iris$Species)
#' row_bartlett(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name bartlett
#' @export
row_bartlett <- function(x, g) {
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

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(gs))
  vPerGroup <- nPerGroup
  for(i in seq_along(gs)) {
    tmpx <- x[,g==gs[i], drop=FALSE]
    nPerGroup[,i] <- ncol(tmpx) - matrixStats::rowCounts(tmpx, value=NA)
    vPerGroup[,i] <- rowVars(tmpx, n=nPerGroup[,i], na.rm=TRUE)
  }
  nPerGroup[nPerGroup < 2] <- NA # drop group with less than 2 observations
  nGroups <- ncol(nPerGroup) - matrixStats::rowCounts(nPerGroup, value=NA)

  nSamples <- rowSums(nPerGroup, na.rm=TRUE)
  vtot <- rowSums(vPerGroup*(nPerGroup-1), na.rm=TRUE) / (nSamples - nGroups)

  df  <- nGroups-1
  ksq <- ((nSamples-nGroups) * log(vtot) - rowSums((nPerGroup-1) * log(vPerGroup), na.rm=TRUE)) /
           (1 + (rowSums(1/(nPerGroup-1), na.rm=TRUE) - 1/(nSamples-nGroups)) / (3 * df))
  p <- stats::pchisq(ksq, df, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'bartlett', 'had less than 2 groups with enough observations')

  w2 <- !w1 & nGroups < length(gs)
  showWarning(w2, 'bartlett', 'had groups with less than 2 observations: those groups were removed')

  w3 <- !w1 & vtot==0
  showWarning(w3, 'bartlett', 'had zero variance in all of the groups')

  w4 <- !w1 & !w3 & rowSums(vPerGroup==0, na.rm=TRUE) > 0
  showWarning(w4, 'bartlett', 'had groups with zero variance: result might be unreliable')

  df[w1 | w3]  <- NA
  ksq[w1 | w3] <- NA
  p[w1 | w3]   <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, var.pooled=vtot, df=df,
             statistic=ksq, pvalue=p, row.names=rnames
             )
}

#' @rdname bartlett
#' @export
col_bartlett <- function(x, g) {
  row_bartlett(t(x), g)
}

