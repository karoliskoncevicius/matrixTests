rowVars <- function(x, na.rm=FALSE) {
  M <- rowMeans(x, na.rm=na.rm)
  n <- rowSums(!is.na(x))-1
  ifelse(n > 0, rowSums((x-M)^2, na.rm=na.rm) / n, NA)
}

rowTables <- function(x) {
  lvl <- sort(unique(as.numeric(x)))
  res <- matrix(numeric(), nrow=nrow(x), ncol=length(lvl))
  colnames(res) <- lvl
  for(i in seq_along(lvl)) {
    res[,i] <- rowSums(x==lvl[i], na.rm=TRUE)
  }
  res
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
