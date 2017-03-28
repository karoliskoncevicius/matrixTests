rowVars <- function(x, na.rm=FALSE) {
  M <- rowMeans(x, na.rm=na.rm)
  n <- rowSums(!is.na(x))-1
  ifelse(n > 0, rowSums((x-M)^2, na.rm=na.rm) / n, NA)
}

