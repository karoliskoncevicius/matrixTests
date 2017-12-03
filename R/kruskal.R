#' Kruskal-Wallis Rank Sum Test
#'
#' Performs a Kruskal-Wallis rank sum test on each row/column of the input matrix.
#'
#' \code{row.kruskalwallis} - Kruskal Wallis test on rows.
#' \code{col.kruskalwallis} - Kruskal Wallis test on columns.
#' Same as as \code{kruskal.test()}
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.

#' @return a data.frame where each row contains the results of a Kruskal-Wallis
#' test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.group - number of groups\cr
#' 4. df - degrees of freedom\cr
#' 5. statistic.chsq - chi-squared statistic\cr
#' 6. p.value - p.value
#'
#' @seealso \code{kruskal.test()}
#'
#' @examples
#' col.kruskalwallis(iris[,1:4], iris$Species)
#' row.kruskalwallis(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name kruskalwallis
#' @export
row.kruskalwallis <- function(x, g) {
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

  ranks <- matrixStats::rowRanks(x, ties.method="average")

  ties <- rowTables(ranks)

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(g)))
  rPerGroup <- nPerGroup
  for(i in seq_along(unique(g))) {
    group <- unique(g)[i]
    nPerGroup[,i] <- matrixStats::rowCounts(!is.na(x[,g==group, drop=FALSE]))
    rPerGroup[,i] <- rowSums(ranks[,g==group, drop=FALSE], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)

  st0 <- rowSums(rPerGroup^2/nPerGroup, na.rm=TRUE)
  st1 <- 12*st0 / (nSamples * (nSamples + 1)) - 3 * (nSamples + 1)
  st2 <- 1 - rowSums(ties^3 - ties) / (nSamples^3 - nSamples)
  stat <- st1/st2
  df <- nGroups - 1

  p <- stats::pchisq(stat, df, lower.tail=FALSE)


  w1 <- nSamples < 2
  showWarning(w1, 'had less than 2 total observations')

  w2 <- !w1 & nGroups < 2
  showWarning(w2, 'had less than 2 groups with enough observations')

  w3 <- !w1 & !w2 & matrixStats::rowCounts(ties!=0)==1
  showWarning(w3, 'were essentially constant')

  stat[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3]    <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, df=df, statistic.chsq=stat,
             p.value=p, row.names=rnames
             )
}

#' @rdname kruskalwallis
#' @export
col.kruskalwallis <- function(x, g) {
  row.kruskalwallis(t(x), g)
}
