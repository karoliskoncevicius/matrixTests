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

  nPerGroup <- matrix(nrow=nrow(x), ncol=length(unique(groups)))
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    vPerGroup[,i] <- rowVars(x[,groups==g, drop=FALSE], na.rm=TRUE)
    mPerGroup[,i] <- rowMeans(x[,groups==g, drop=FALSE], na.rm=TRUE)
    nPerGroup[,i] <- rowSums(!is.na(x[,groups==g, drop=FALSE]))
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- rowSums(nPerGroup!=0)
  M <- rowMeans(x, na.rm=TRUE)

  betweenScatter <- rowSums(nPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  withinScatter  <- rowSums((nPerGroup-1) * vPerGroup, na.rm=TRUE)

  F <- (betweenScatter/(nGroups-1)) / (withinScatter/(nSamples-nGroups))
  p <- pf(F, nGroups-1, nSamples-nGroups, lower.tail=FALSE)

  bad <- nGroups < 2
  if(any(bad)) {
    warning(sum(bad), ' of the rows had less than 2 groups with enough observations')
  }

  bad <- nGroups==nSamples & nGroups > 1
  if(any(bad)) {
    warning(sum(bad), ' of the rows had one observation per group')
  }

  bad <- withinScatter==0 & nGroups!=nSamples & nGroups > 1
  if(any(bad)) {
    warning(sum(bad), ' of the rows had essentially perfect fit')
  }

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(sum.sq.treatment=betweenScatter, sum.sq.residuals=withinScatter,
             mean.sq.treatment=betweenScatter/(nGroups-1),
             mean.sq.residuals=withinScatter/(nSamples-nGroups),
             obs.tot=nSamples, obs.groups=nGroups,
             df.treatment=nGroups-1, df.residuals=nSamples-nGroups,
             F.statistic=F, p.value=p,
             row.names=rnames
             )
}

# TODO: finalize
oneway_welch <- function(x, groups) {

  if(!is.null(x) && is.vector(x))
    x <- matrix(x, nrow=1)

  if(!is.null(x) && is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  if(!is.null(groups))
    groups <- as.character(groups)

  assert_vec_length(groups, ncol(x))

  bad <- is.na(groups)
  if(any(bad)) {
    warning(sum(bad), ' columns dropped due to missing group information')
    x      <- x[,!bad, drop=FALSE]
    groups <- groups[!bad]
  }

  nPerGroup <- matrix(nrow=nrow(x), ncol=length(unique(groups)))
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    vPerGroup[,i] <- rowVars(x[,groups==g, drop=FALSE], na.rm=TRUE)
    mPerGroup[,i] <- rowMeans(x[,groups==g, drop=FALSE], na.rm=TRUE)
    nPerGroup[,i] <- rowSums(!is.na(x[,groups==g, drop=FALSE]))
  }

  nSamples  <- rowSums(nPerGroup)
  nGroups   <- rowSums(nPerGroup!=0)
  wPerGroup <- nPerGroup/vPerGroup
  wTot      <- rowSums(wPerGroup)
  M         <- rowSums((wPerGroup * mPerGroup)/wTot)

  betweenScatter <- rowSums(wPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  tmp <- rowSums((1 - wPerGroup/wTot)^2/(nPerGroup - 1)) / (nGroups^2 - 1)

  F <- betweenScatter / ((nGroups-1) * (1 + 2 * (nGroups-2) * tmp))
  p <- pf(F, nGroups-1, 1/(3*tmp), lower.tail=FALSE)

  bad <- nGroups < 2
  if(any(bad)) {
    warning(sum(bad), ' of the rows had less than 2 groups with enough observations')
  }

  bad <- nGroups==nSamples
  if(any(bad)) {
    warning(sum(bad), ' of the rows had one observation per group')
  }

  bad <- withinScatter==0 & nGroups!=nSamples & nGroups > 1
  if(any(bad)) {
    warning(sum(bad), ' of the rows had essentially perfect fit')
  }

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(sum.sq.treatment=betweenScatter, sum.sq.residuals=withinScatter,
             mean.sq.treatment=betweenScatter/(nGroups-1),
             mean.sq.residuals=withinScatter/(nSamples-nGroups),
             obs.tot=nSamples, obs.groups=nGroups,
             df.treatment=nGroups-1, df.residuals=1/(3*tmp),
             F.statistic=F, p.value=p,
             row.names=rnames
             )
}
