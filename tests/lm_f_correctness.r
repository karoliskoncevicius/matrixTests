library(matrixTests)

#--- functions -----------------------------------------------------------------

base_lm <- function(mat, mod1, mod0) {
  stopifnot(ncol(mat) == nrow(mod1))
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  n <- rsq1 <- rsq0 <- dfm <- dfr <- f <- p <- numeric(nrow(mat))
  betas <- matrix(nrow=nrow(mat), ncol=ncol(mod1))
  colnames(betas) <- paste0("beta.", 1:ncol(betas)-1)

  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    lm0 <- lm(vec ~ 0 + mod0[!bad,])
    lm1 <- lm(vec ~ 0 + mod1[!bad,])
    res <- anova(lm0, lm1)

    n[i]      <- length(vec)
    betas[i,] <- coefficients(lm1)
    rsq0[i]   <- 1 - (res$RSS[1] / sum((X[i,] - mean(X[i,]))^2))
    rsq1[i]   <- 1 - (res$RSS[2] / sum((X[i,] - mean(X[i,]))^2))
    dfm[i]    <- res$Df[2]
    dfr[i]    <- res$Res.Df[2]
    f[i]      <- res$F[2]
    p[i]      <- res[2,6]
  }

  data.frame(obs=n, betas, rsquared.model=rsq1, rsquared.null=rsq0,
             df.model=dfm, df.residual=dfr, statistic=f, pvalue=p
             )
}


#--- montecarlo ----------------------------------------------------------------

# random 100 rows with various covariates
X  <- matrix(rnorm(1000), ncol=10)
m  <- model.matrix(~ rnorm(10) + sample(1:10) + runif(10, -1, 1) + sample(c(0,1), 10, replace=TRUE))
m0 <- m[,c(1,sample(2:ncol(m), 3))]
res1 <- base_lm(X, m, m0)
res2 <- row_lm_f(X, m, m0 )
stopifnot(all.equal(res1, res2))

