#' iEVORA
#'
#' Epigenetic Variable Outliers for cancer Risk prediction Analysis
#'
#' Measures differential variability between two groups. The algorithm has
#' 2 steps: detecting difference in variance (Bartlett's test) and detecting
#' difference in means(t-test). The second step is done to regularize the
#' variability test which is overly sensitive to single outliers.
#'
#' By default the result is considered significant if variability test produces
#' a significant p-value (below selected threshold) after FDR correction and
#' t-test returns a significant p-value without using the FDR correction.
#'
#' The algorithm is aimed at large DNA methylation data sets, where one may wish
#' to find features (mostly CpGs) which differ between two normal cellular
#' phenotypes, but with one of these phenotypes representing cells which are at
#' risk of neoplastic transformation.
#'
#' @param x numeric matrix
#' @param groups a vector specifying group membership for each observation of x.
#' Must contain two unique groups one labeled "1" and another "0".
#' If the vector is neither numeric nor logical the group appearing first is
#' labeled "0" and the remaining one as "1".
#' @param cutT cutoff threshold for the raw p-value of the t-test step.
#' (default 0.05)
#' @param cutBfdr cutoff threshold for the FDR-corrected p-value of the
#' Bartlett's test step. (default 0.001)
#'
#' @return a data.frame where each row contains result of the iEVORA algorithm
#' for the corresponding row/column of x. \cr\cr
#' Each row contains the following information (in order):\cr
#' 1. obs.0 - number of observations in 0 group\cr
#' 2. obs.1 - number of observations in 1 group\cr
#' 3. obs.tot - number of total observations\cr
#' 4. mean.0 - mean of the 0 group \cr
#' 5. mean.1 - mean of the 1 group \cr
#' 6. mean.diff - mean difference (group1 - group0)\cr
#' 7. var.0 - variance of the 0 group \cr
#' 8. var.1 - variance of the 1 group \cr
#' 9. var.log2.ratio - log ratio of variances log(var1/var0) \cr
#' 10. statistic.t - t.statistic of the t-test step \cr
#' 11. p.value.t - raw p-value of the t-test step \cr
#' 12. statistic.bt - chsq.statistic of the bartlett test step \cr
#' 13. p.value.bt - raw p-value of the Bartlett's test step \cr
#' 14. q.value.bt - fdr-adjusted p-value of the Bartlett's test step \cr
#' 15. significant - indicator showing if the result was significant \cr
#' 16. rank - rank of the significant results (ordered by t.test p-value)
#'
#' @seealso \code{row.bartlett}, \code{row.t.welch}
#'
#' @examples
#' # perform iEVORA on iris dataset for setosa against all other groups
#' col.ievora(iris[,1:4], iris$Species=="setosa")
#'
#' @references Andrew E Teschendorff et.al. DNA methylation outliers in normal
#' breast tissue identify field defects that are enriched in cancer.
#' Nature Communications 7, 10478 (2016) doi:10.1038/ncomms10478
#'
#' @author Karolis Koncevičius
#' @name ievora
#' @export
row.ievora <- function(x, groups, cutT=0.05, cutBfdr=0.001) {
  force(x)
  force(groups)

  if(is.vector(x))
    x <- matrix(x, nrow=1)

  if(is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_vec_length(groups, ncol(x))
  assert_max_number_of_levels(groups, 2)

  bad <- is.na(groups)
  if(any(bad)) {
    warning(sum(bad), ' columns dropped due to missing group information')
    x      <- x[,!bad, drop=FALSE]
    groups <- groups[!bad]
  }

  if(is.logical(groups)) {
    groups <- as.numeric(groups)
  }

  if(!is.numeric(groups) | !(all(groups %in% c(0,1)))) {
    groups <- match(groups, unique(groups))-1
  }

  assert_numeric_vec_length(cutT,  1)
  assert_numeric_vec_length(cutBfdr,  1)
  assert_all_in_range(cutT, 0, 1)
  assert_all_in_range(cutBfdr, 0, 1)

  tres <- row.t.welch(x[,groups==1, drop=FALSE], x[,groups==0, drop=FALSE])
  bres <- row.bartlett(x, groups)

  brq <- stats::p.adjust(bres$p.value, "fdr")
  isSig <- brq < cutBfdr & tres$p.value < cutT
  isSig[is.na(isSig)] <- FALSE
  rank  <- rep(NA, length(isSig))
  rank[isSig] <- rank(tres$p.value[isSig], ties.method="first")

  var0 <- rowVars(x[,groups==0, drop=FALSE], na.rm=TRUE)
  var1 <- rowVars(x[,groups==1, drop=FALSE], na.rm=TRUE)
  logR <- log2(var1/var0)

  rnames <- rownames(x)
  if(!is.null(rnames)) rnames <- make.unique(rnames)
  data.frame(obs.0=tres$obs.y, obs.1=tres$obs.x, obs.tot=tres$obs.tot,
             mean.0=tres$mean.y, mean.1=tres$mean.x, mean.diff=tres$mean.diff,
             var.0=var0, var.1=var1, var.log2.ratio=logR,
             statistic.t=tres$statistic.t, p.value.t=tres$p.value,
             statistic.bt=bres$statistic.chsq, p.value.bt=bres$p.value,
             q.value.bt=brq, significant=isSig,
             rank=rank, row.names=rnames
             )
}

#' @rdname ievora
#' @export
col.ievora <- function(x, groups, cutT=0.05, cutBfdr=0.001) {
  row.ievora(t(x), groups, cutT=cutT, cutBfdr=cutBfdr)
}

