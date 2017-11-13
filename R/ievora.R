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
#' @param groups a vector specifying group membership for each column of x.
#' Must contain two unique groups one labeled "1" and another "0". Otherwise
#' the group appearing first is labeled "0" and the remaining one as "1".
#' @param cutT cutoff threshold for the raw p-value of the t-test step.
#' (default 0.05)
#' @param cutBfdr cutoff threshold for the FDR-corrected p-value of the
#' Bartlett's test step. (default 0.001)
#'
#' @return a data.frame where each row contains result of the iEVORA algorithm
#' for the corresponding row of x. \cr
#' Each row contains the following information (in order): \cr
#' 1. mean of the 0 group. \cr
#' 2. mean of the 1 group. \cr
#' 3. variance of the 0 group. \cr
#' 4. variance of the 1 group. \cr
#' 5. log ratio of variances log(var1/var0). \cr
#' 6. t.statistic. \cr
#' 7. raw p-value of the t-test step. \cr
#' 8. raw p-value of the Bartlett's test step. \cr
#' 9. fdr-adjusted p-value of the Bartlett's test step. \cr
#' 10. indicator showing whether or not the result was significant. \cr
#' 11. rank of the significant results (ordered by t.test p-value)
#'
#' @seealso \code{bartlett}, \code{ttest_welch}
#'
#' @examples
#' # perform ievora on iris dataset for setosa against all other groups
#' ievora(t(iris[,1:4]), iris$Species=="setosa")
#'
#' @references Andrew E Teschendorff et.al. DNA methylation outliers in normal
#' breast tissue identify field defects that are enriched in cancer.
#' Nature Communications 7, 10478 (2016) doi:10.1038/ncomms10478
#'
#' @author Karolis Koncevičius
#' @export
ievora <- function(x, groups, cutT=0.05, cutBfdr=0.001) {

  if(!is.null(x) && is.vector(x))
    x <- matrix(x, nrow=1)

  if(!is.null(x) && is.data.frame(x) && all(sapply(x, is.numeric)))
    x <- data.matrix(x)

  assert_numeric_mat_or_vec(x)

  assert_number_of_levels(groups, 2)

  if(!is.null(groups) && is.logical(groups)) {
    groups <- as.numeric(groups)
  }

  if(!is.numeric(groups) | !(all(groups %in% c(0,1)))) {
    groups <- ifelse(groups==unique(na.omit(groups))[1], 0, 1)
  }

  assert_numeric_vec_length(groups, ncol(x))

  bad <- is.na(groups)
  if(any(bad)) {
    warning(sum(bad), " columns dropped due to missing group information")
    x      <- x[,!bad, drop=FALSE]
    groups <- groups[!bad]
  }


  tres <- ttest_welch(x[,groups==0, drop=FALSE], x[,groups==1, drop=FALSE])
  bres <- bartlett(x, groups)

  brq <- p.adjust(bres$p.value, "fdr")
  isSig <- brq < cutBfdr & tres$p.value < cutT
  rank  <- rep(NA, length(isSig))
  rank[isSig] <- rank(tres$p.value[isSig], ties.method="first")

  var0 <- rowVars(x[,groups==0, drop=FALSE])
  var1 <- rowVars(x[,groups==1, drop=FALSE])
  logR <- log2(var1/var0)

  data.frame(mean.0=tres$mean.x, mean.1=tres$mean.y, var.0=var0, var.1=var1,
             logR=logR, t.statistic=-tres$t.statistic, tt.p.value=tres$p.value,
             bt.p.value=bres$p.value, bt.q.value=brq, significant=isSig,
             rank=rank, row.names=rownames(tres)
             )
}

