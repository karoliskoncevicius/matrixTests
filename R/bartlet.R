#' Bartlett test
#'
#' Performs the Bartlett's test of homogeneity of variances on each row/column
#' of the input matrix.
#'
#' NA values are always ommited. If values are missing for a whole group - that
#' group is discarded. Groups with only one observation are also discarded.
#'
#' \code{row.bartlett(x, groups)} - Bartlet's test on rows.
#' \code{col.bartlett(x, groups)} - Bartlet's test on columns.
#' Same as \code{bartlett.test(x,  groups)}
#'
#' @param x numeric matrix.
#' @param groups a vector specifying group membership for each observation of x.

#' @return a data.frame where each row contains the results of the bartlett test
#' performed on the corresponding row/column of x.\cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.tot - total number of observations\cr
#' 2. obs.groups - number of groups\cr
#' 3. var.pooled - pooled variance estimate\cr
#' 4. df - degrees of freedom\cr
#' 5. statistic.chsq - chi-squared statistic\cr
#' 6. p.value - p.value
#'
#' @seealso \code{bartlett.test()}
#'
#' @examples
#' col.bartlett(iris[,1:4], iris$Species)
#' row.bartlett(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name bartlett
#' @export
row.bartlett <- function(x, groups) {
  force(x)
  force(groups)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_vec_length(groups, ncol(x))


  bad <- is.na(groups)
  if(any(bad)) {
    warning(sum(bad), ' columns dropped due to missing group information')
    x      <- x[,!bad, drop=FALSE]
    groups <- groups[!bad]
  }

  groups <- as.character(groups)

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(groups)))
  vPerGroup <- nPerGroup
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    nPerGroup[,i] <- matrixStats::rowCounts(!is.na(x[,groups==g, drop=FALSE]))
    vPerGroup[,i] <- rowVars(x[,groups==g, drop=FALSE], na.rm=TRUE)
  }
  nPerGroup[nPerGroup < 2] <- NA # drop groups with less than 2 observations
  nGroups <- matrixStats::rowCounts(!is.na(nPerGroup))

  nSamples <- rowSums(nPerGroup, na.rm=TRUE)
  vtot <- rowSums(vPerGroup*(nPerGroup-1), na.rm=TRUE) / (nSamples - nGroups)
  df   <- nGroups-1

  ksq  <- ((nSamples-nGroups) * log(vtot) - rowSums((nPerGroup-1) * log(vPerGroup), na.rm=TRUE)) /
           (1 + (rowSums(1/(nPerGroup-1), na.rm=TRUE) - 1/(nSamples-nGroups)) / (3 * df))
  p <- stats::pchisq(ksq, df, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & nGroups < length(unique(groups))
  showWarning(w2, 'had groups with less than 2 observations: those groups were removed')

  w3 <- !w1 & vtot==0 & nGroups!=0
  showWarning(w3, 'had zero variance in all of the groups')

  w4 <- !w1 & !w3 & rowSums(vPerGroup==0, na.rm=TRUE) > 0
  showWarning(w4, 'had groups with zero variance: result might be unreliable')

  ksq[w1 | w3] <- NA
  p[w1 | w3] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, var.pooled=vtot, df=df,
             statistic.chsq=ksq, p.value=p, row.names=rnames
             )
}

#' @rdname bartlett
#' @export
col.bartlett <- function(x, groups) {
  row.bartlett(t(x), groups)
}

