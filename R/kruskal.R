#' Kruskal-Wallis Rank Sum Test
#'
#' Performs a Kruskal-Wallis rank sum test on each row of the input matrix.
#'
#' \code{kruskalwallis} - sam as as \code{kruskal.test()}
#'
#' @param x numeric matrix.
#' @param groups a vector specifying group membership for each column of x.

#' @return a data.frame where each row contains the results of a Kruskal-Wallis
#' test performed on the corresponding row of x.
#'
#' @seealso \code{kruskal.test()}
#'
#' @examples
#' kruskalwallis(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @export
kruskalwallis <- function(x, groups) {
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

  ranks <- matrix(numeric(), nrow=nrow(x), ncol=ncol(x))
  ranks[] <- t(apply(x, 1, rank))
  ranks[is.na(x)] <- NA

  ties <- rowTables(ranks)

  nPerGroup <- matrix(numeric(), nrow=nrow(x), ncol=length(unique(groups)))
  rPerGroup <- nPerGroup
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    nPerGroup[,i] <- rowSums(!is.na(x[,groups==g, drop=FALSE]))
    rPerGroup[,i] <- rowSums(ranks[,groups==g, drop=FALSE], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- rowSums(nPerGroup!=0)

  st0 <- rowSums(rPerGroup^2/nPerGroup, na.rm=TRUE)
  st1 <- 12*st0 / (nSamples * (nSamples + 1)) - 3 * (nSamples + 1)
  st2 <- 1 - rowSums(ties^3 - ties) / (nSamples^3 - nSamples)
  stat <- st1/st2
  df <- nGroups - 1

  bad <- nSamples < 2
  if(any(bad, na.rm=TRUE))
    warning(paste0(sum(bad), ' of the rows had less than 2 total observations'))

  stat[bad] <- NA
  df[bad]   <- NA


  bad <- nGroups < 2 & nSamples > 1
  if(any(bad)) {
    warning(sum(bad), ' of the rows had less than 2 groups with enough observations')
  }

  stat[bad] <- NA
  df[bad]   <- NA

  p <- pchisq(stat, df, lower.tail=FALSE)

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, df=df, chsq.statistic=stat,
             p.value=p, row.names=rnames
             )
}

