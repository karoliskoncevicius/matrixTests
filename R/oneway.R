#' ONEWAY ANOVA
#'
#' Performs an analysis of variance tests on each row of the input matrix.
#'
#' Functions to perform ONEWAY ANOVA analysis for rows of matrices.
#'
#' \code{oneway_equalvar} - one-way anova. Same as \code{aov(x ~ groups)}
#' \code{oneway_welch} _ one-way anova with Welch correction for variances.
#' Same as \code{oneway.test(var.equal=FALSE)}
#'
#' @name oneway
#'
#' @param x numeric matrix.
#' @param groups a vector specifying group membership for each column of x.

#' @return a data.frame where each row contains the results of an oneway anova
#' test performed on the corresponding row of x.
#'
#' @seealso \code{aov()}, \code{oneway.test()}
#'
#' @examples
#' oneway_equalvar(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @export
oneway_equalvar <- function(x, groups) {

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

  M <- rowMeans(x, na.rm=TRUE)
  withinScatter  <- numeric(nrow(x))
  betweenScatter <- numeric(nrow(x))
  nGroups        <- numeric(nrow(x))
  nSamples       <- numeric(nrow(x))
  for(g in unique(groups)) {
    nGroupObs      <- rowSums(!is.na(x[,groups==g, drop=FALSE]))
    nSamples       <- nSamples + nGroupObs
    nGroups        <- nGroups + ifelse(nGroupObs==0, 0, 1)
    gMeans <- rowMeans(x[,groups==g, drop=FALSE], na.rm=TRUE)
    gMeans <- ifelse(nGroups==0, 0, gMeans)
    groupScatter   <- rowSums((x[,groups==g, drop=FALSE]-gMeans)^2, na.rm=TRUE)
    withinScatter  <- withinScatter + groupScatter
    betweenScatter <- betweenScatter + nGroupObs*(gMeans-M)^2
  }

  F <- (betweenScatter/(nGroups-1)) / (withinScatter/(nSamples-nGroups))
  p <- pf(F, nGroups-1, nSamples-nGroups, lower.tail=FALSE)

  bad <- nGroups < 2
  if(any(bad)) {
    warning(sum(bad), " rows had less than 2 groups with enough observations")
  }

  bad <- nGroups==nSamples
  if(any(bad)) {
    warning(sum(bad), " rows had one observation per group")
  }

  bad <- withinScatter==0 & nGroups!=nSamples & nGroups > 1
  if(any(bad)) {
    warning(sum(bad), " rows had essentially perfect fit")
  }

  data.frame(sum.sq.treatment=betweenScatter, sum.sq.residuals=withinScatter,
             mean.sq.treatment=betweenScatter/(nGroups-1),
             mean.sq.residuals=withinScatter/(nSamples-nGroups),
             obs.tot=nSamples, obs.groups=nGroups,
             df.treatment=nGroups-1, df.residuals=nSamples-nGroups,
             F.statistic=F, p.value=p
             )
}

