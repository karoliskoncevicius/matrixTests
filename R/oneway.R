#' ONEWAY ANOVA
#'
#' Performs an analysis of variance tests on each row of the input matrix.
#'
#' Functions to perform ONEWAY ANOVA analysis for rows of matrices.
#'
#' \code{row.oneway.equalvar} - one-way anova. Same as \code{aov(x ~ groups)}
#' \code{row.oneway.welch} _ one-way anova with Welch correction for variances.
#' Same as \code{oneway.test(var.equal=FALSE)}
#'
#' @param x numeric matrix.
#' @param groups a vector specifying group membership for each column of x.

#' @return a data.frame where each row contains the results of an oneway anova
#' test performed on the corresponding row of x.
#'
#' @seealso \code{aov()}, \code{oneway.test()}
#'
#' @examples
#' row.oneway.equalvar(t(iris[,1:4]), iris$Species)
#'
#' @author Karolis Koncevičius
#' @name oneway
#' @export
row.oneway.equalvar <- function(x, groups) {
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
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    nPerGroup[,i] <- matrixStats::rowCounts(!is.na(x[,groups==g, drop=FALSE]))
    mPerGroup[,i] <- rowMeans(x[,groups==g, drop=FALSE], na.rm=TRUE)
    vPerGroup[,i] <- rowVars(x[,groups==g, drop=FALSE], na.rm=TRUE)
  }

  nSamples <- rowSums(nPerGroup)
  nGroups  <- matrixStats::rowCounts(nPerGroup!=0)
  M <- rowMeans(x, na.rm=TRUE)

  betweenScatter <- rowSums(nPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  withinScatter  <- rowSums((nPerGroup-1) * vPerGroup, na.rm=TRUE)

  dft <- nGroups-1
  dfr <- nSamples-nGroups

  F <- (betweenScatter/dft) / (withinScatter/dfr)
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & nGroups==nSamples
  showWarning(w2, 'had one observation per group')

  w3 <- !w1 & !w2 & withinScatter==0 & betweenScatter==0
  showWarning(w3, 'had essentially constant values')

  w4 <- !w1 & !w2 & withinScatter==0 & betweenScatter!=0
  showWarning(w4, 'had zero within group variance: result might be unreliable')

  F[w1 | w2 | w3] <- NA
  p[w1 | w2 | w3] <- NA


  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(sum.sq.treatment=betweenScatter, sum.sq.residuals=withinScatter,
             mean.sq.treatment=betweenScatter/dft,
             mean.sq.residuals=withinScatter/dfr,
             obs.tot=nSamples, obs.groups=nGroups,
             df.treatment=dft, df.residuals=dfr,
             statistic.F=F, p.value=p,
             row.names=rnames
             )
}

#' @export
#' @rdname oneway
row.oneway.welch <- function(x, groups) {
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
  mPerGroup <- vPerGroup <- nPerGroup
  for(i in seq_along(unique(groups))) {
    g <- unique(groups)[i]
    nPerGroup[,i] <- matrixStats::rowCounts(!is.na(x[,groups==g, drop=FALSE]))
    mPerGroup[,i] <- rowMeans(x[,groups==g, drop=FALSE], na.rm=TRUE)
    vPerGroup[,i] <- rowVars(x[,groups==g, drop=FALSE], na.rm=TRUE)
  }
  mPerGroup[nPerGroup<2] <- NA
  vPerGroup[nPerGroup<2] <- NA
  nPerGroup[nPerGroup<2] <- NA

  nSamples  <- rowSums(nPerGroup, na.rm=TRUE)
  nGroups   <- matrixStats::rowCounts(nPerGroup!=0, na.rm=TRUE)
  wPerGroup <- nPerGroup/vPerGroup
  wTot      <- rowSums(wPerGroup, na.rm=TRUE)
  M         <- rowSums((wPerGroup * mPerGroup)/wTot, na.rm=TRUE)

  betweenScatter <- rowSums(wPerGroup * (mPerGroup-M)^2, na.rm=TRUE)
  tmp <- rowSums((1 - wPerGroup/wTot)^2/(nPerGroup - 1), na.rm=TRUE) / (nGroups^2 - 1)
  dfr <- 1/(3*tmp)
  dft <- nGroups-1

  F <- betweenScatter / (dft * (1 + 2 * (nGroups-2) * tmp))
  p <- stats::pf(F, dft, dfr, lower.tail=FALSE)


  w1 <- nGroups < 2
  showWarning(w1, 'had less than 2 groups with enough observations')

  w2 <- !w1 & nGroups < length(unique(groups))
  showWarning(w2, 'had groups with less than 2 observations: those groups were removed')

  w3 <- !w1 & rowSums(vPerGroup!=0, na.rm=TRUE)==0
  showWarning(w3, 'had zero variance in all of the groups')

  w4 <- !w1 & !w3 & rowSums(vPerGroup==0, na.rm=TRUE) > 0
  showWarning(w4, 'had groups with zero variance: result might be unreliable')


  F[w1 | w3] <- NA
  p[w1 | w3] <- NA

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.tot=nSamples, obs.groups=nGroups, df.treatment=dft,
             df.residuals=dfr, statistic.F=F, p.value=p,
             row.names=rnames
             )
}

