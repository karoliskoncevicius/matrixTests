rowVars <- function(x, n=NULL, m=NULL, na.rm=FALSE) {
  if(is.null(n))
    n <- ncol(x) - matrixStats::rowCounts(x, value=NA)
  if(is.null(m))
    m <- rowMeans(x, na.rm=na.rm)
  res <- rowSums((x-m)^2, na.rm=na.rm) / (n-1)
  res[n <= 1] <- NA
  res[!is.finite(m)] <- NaN
  res
}

rowRankTies <- function(r) {
  if(storage.mode(r) != "integer") {
    storage.mode(r) <- "integer"
  }
  r <- t(r)
  res <- vector("list", ncol(r))
  for(i in seq_len(ncol(r))) {
    res[[i]] <- tabulate(r[,i], nrow(r))
  }
  if(length(res) == 0) {  # needed to work with zero-row inputs
    res <- list(integer())
  }
  do.call(rbind, res)
}

showWarning <- function(isWarning, err) {
  if(any(isWarning, na.rm=TRUE)) {
    parentFun <- deparse(as.list(sys.call(-1))[[1]])
    grandFun  <- as.list(sys.call(-2))
    if(length(grandFun) > 0) {
      grandFun <- deparse(grandFun[[1]])
      if(grandFun %in% getNamespaceExports("matrixTests")) {
        parentFun <- grandFun
      }
    }
    pref <- "row"
    if(grepl("^col_", parentFun)) pref <- "column"
    n <- sum(isWarning, na.rm=TRUE)
    i <- match(TRUE, isWarning)
    err <- paste0(parentFun, ": ", n, ' of the ', pref, 's ', err, ".",
                  '\nFirst occurrence at ', pref, ' ', i
                  )
    warning(err, call.=FALSE)
  }
}
