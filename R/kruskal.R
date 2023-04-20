#' Kruskal-Wallis rank sum test
#'
#' Performs a Kruskal-Wallis rank sum test on each row/column of the input matrix.
#'
#' \code{row_kruskalwallis(x, g)} - Kruskal Wallis test on rows.
#' \code{col_kruskalwallis(x, g)} - Kruskal Wallis test on columns.
#'
#' Results should be the same as running \code{kruskal.test(x, g)}
#' on every row (or column) of \code{x}
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.
#'
#' @return a data.frame where each row contains the results of a Kruskal-Wallis
#' test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 4. df - degrees of freedom\cr
#' 5. statistic - chi-squared statistic\cr
#' 6. pvalue - p.value
#'
#' @seealso \code{kruskal.test()}
#'
#' @examples
#' col_kruskalwallis(iris[,1:4], iris$Species)
#' row_kruskalwallis(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name kruskalwallis
#' @export
row_kruskalwallis <- function(x, g) {
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

  ranks <- matrixStats::rowRanks(x, ties.method="average")

  ties <- rowRankTies(ranks)

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(gs))
  rPerGroup <- nPerGroup
  for(i in seq_along(gs)) {
    inds <- g == gs[i]
    tmpx <- x[,inds, drop=FALSE]
    nPerGroup[,i] <- ncol(tmpx) - matrixStats::rowCounts(tmpx, value=NA)
    rPerGroup[,i] <- rowSums(ranks[,inds, drop=FALSE], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)

  df <- nGroups - 1

  st0 <- rowSums(rPerGroup*rPerGroup/nPerGroup, na.rm=TRUE)
  st1 <- 12*st0 / (nSamples * (nSamples + 1)) - 3 * (nSamples + 1)
  st2 <- 1 - rowSums(ties^3 - ties) / (nSamples^3 - nSamples)
  stat <- st1/st2

  p <- stats::pchisq(stat, df, lower.tail=FALSE)


  w1 <- nSamples < 2
  showWarning(w1, 'kruskalwallis', 'had less than 2 total observations')

  w2 <- !w1 & nGroups < 2
  showWarning(w2, 'kruskalwallis', 'had less than 2 groups with enough observations')

  w3 <- !w1 & !w2 & matrixStats::rowMaxs(ties)==nSamples
  showWarning(w3, 'kruskalwallis', 'had essentially constant values')

  df[w1 | w2 | w3]   <- NA
  stat[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3]    <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, df=df, statistic=stat,
             pvalue=p, row.names=rnames
             )
}

#' @rdname kruskalwallis
#' @export
col_kruskalwallis <- function(x, g) {
  row_kruskalwallis(t(x), g)
}
