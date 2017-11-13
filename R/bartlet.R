#' Bartlett test
#'
#' Performs the Bartlett's test of homogeneity of variances on each row of the
#' input matrix.
#'
#' Function to perform Bartlett's tests for homogeneity of variances in each of
#' the groups on the rows of a matrix.
#'
#' NA values are always ommited. If values are missing for a whole group - that
#' group is discarded. Groups with only one observation are also discarded.
#'
#' \code{bartlett(x, groups)} - Bartlet's test.
#' Same as \code{bartlett.test(x,  groups)}
#'
#' @name bartlett
#'
#' @param x numeric matrix.
#' @param groups a vector specifying group membership for each column of x.

#' @return a data.frame where each row contains the results of the bartlett test
#' performed on the corresponding row of x.
#'
#' @seealso \code{bartlett.test()}
#'
#' @examples
#' bartlett(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @export
bartlett <- function(x, groups) {

  if(!is.null(x) && is.vector(x))
    x <- matrix(x, nrow=1)

  if(!is.null(x) && is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  if(!is.null(groups))
    groups <- as.character(groups)

  assert_character_vec_length(groups, ncol(x))


  bad <- is.na(groups)
  if(any(bad)) {
    warning(sum(bad), " columns dropped due to missing group information")
    x      <- x[,!bad, drop=FALSE]
    groups <- groups[!bad]
  }

  nPerGroup <- matrix(nrow=nrow(x), ncol=length(unique(groups)),
                      dimnames=list(rownames(x))
                      )
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    nPerGroup[,i] <- rowSums(!is.na(x[,groups==g, drop=FALSE]))
  }
  nPerGroup[nPerGroup < 2] <- NA
  nGroups <- rowSums(!is.na(nPerGroup))

  bad <- nGroups < 2
  if(any(bad))
    warning(sum(bad), " rows had less than 2 groups with enough observations")

  bad <- nGroups < length(unique(groups))
  if(any(bad))
    warning(sum(bad), " rows had groups with less than 2 observations")


  vPerGroup <- matrix(nrow=nrow(x), ncol=length(unique(groups)))
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    vPerGroup[,i] <- rowVars(x[,groups==g, drop=FALSE], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup, na.rm=TRUE)
  vtot <- rowSums(vPerGroup*(nPerGroup-1), na.rm=TRUE) / (nSamples - nGroups)
  df   <- nGroups-1
  ksq  <- ((nSamples-nGroups) * log(vtot) - rowSums((nPerGroup-1) * log(vPerGroup), na.rm=TRUE)) /
           (1 + (rowSums(1/(nPerGroup-1), na.rm=TRUE) - 1/(nSamples-nGroups)) / (3 * df))

  p <- pchisq(ksq, df, lower.tail=FALSE)

  data.frame(var.tot=vtot, obs.tot=nSamples, obs.groups=nGroups,
             ksq.statistic=ksq, p.value=p, df=df
             )
}

