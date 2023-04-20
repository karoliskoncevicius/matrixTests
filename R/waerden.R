#' Van der Waerden test
#'
#' Performs van der Waerden test on each row/column of the input matrix.
#'
#' \code{row_waerden(x, g)} - van der Waerden test on rows.
#' \code{col_waerden(x, g)} - van det Waerden test on columns.
#'
#' @param x numeric matrix.
#' @param g a vector specifying group membership for each observation of x.
#'
#' @return a data.frame where each row contains the results of van det Waerden
#' test performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 3. df - degrees of freedome\cr
#' 4. statistic - van det Waerden chi-squared statistic\cr
#' 5. pvalue - p.value
#'
#' @seealso \code{\link[PMCMRplus]{vanWaerdenTest}}, \code{row_oneway_equalvar}, \code{row_kruskalwallis}
#'
#' @examples
#' col_waerden(iris[,1:4], iris$Species)
#' row_waerden(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name waerden
#' @export
row_waerden <- function(x, g) {
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

  r <- matrixStats::rowRanks(x, ties.method="average")
  n <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  z <- stats::qnorm(r/(n+1))
  z <- matrix(z, nrow=nrow(x), ncol=ncol(x))

  nPerGroup <- matrix(numeric(), nrow=nrow(z), ncol=length(gs))
  sPerGroup <- nPerGroup
  for(i in seq_along(gs)) {
    tmpx <- z[,g==gs[i], drop=FALSE]
    nPerGroup[,i] <- ncol(tmpx) - matrixStats::rowCounts(tmpx, value=NA)
    sPerGroup[,i] <- rowSums(tmpx, na.rm=TRUE)
  }

  nGroups <- matrixStats::rowCounts(nPerGroup!=0)

  s2   <- rowSums(z^2, na.rm=TRUE) / (n - 1)
  stat <- rowSums(sPerGroup^2 / nPerGroup, na.rm=TRUE) / s2

  df <- nGroups - 1
  p  <- stats::pchisq(stat, df, lower.tail = FALSE)


  w1 <- n < 2
  showWarning(w1, 'waerden', 'had less than 2 total observations')

  w2 <- !w1 & nGroups < 2
  showWarning(w2, 'waerden', 'had less than 2 groups with enough observations')

  w3 <- !w1 & !w2 & s2==0
  showWarning(w3, 'waerden', 'had essentially constant values')

  df[w1 | w2 | w3]   <- NA
  stat[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3]    <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=n, obs.groups=nGroups, df=df, statistic=stat, pvalue=p,
             row.names=rnames
             )
}

#' @rdname waerden
#' @export
col_waerden <- function(x, g) {
  row_waerden(t(x), g)
}

