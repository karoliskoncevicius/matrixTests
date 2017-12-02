context("correctness of ievora")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_ievora <- function(mat, groups, tCut=0.05, bCut=0.001) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat <- mat[,!bad, drop=FALSE]
  groups <- groups[!bad]

  doDV <- function(tmp.v,pheno.v) {
    co.idx <- which(pheno.v==0)
    ca.idx <- which(pheno.v==1)
    bt.o <- bartlett.test(x=tmp.v,g=pheno.v)
    pv <- bt.o$p.value
    logR <- log2(var(tmp.v[ca.idx])/var(tmp.v[co.idx]))
    avCA <- mean(tmp.v[ca.idx])
    avCO <- mean(tmp.v[co.idx])
    out.v <- c(logR,pv,avCA,avCO)
    names(out.v) <- c("log(V1/V0)","P(BT)","Av1","Av0")
    return(out.v)
  }

  doTT <- function(tmp.v,pheno.v){
    tt.o <- t.test(tmp.v ~ pheno.v)
    out.v <- c(-tt.o$stat,tt.o$p.val)
    names(out.v) <- c("t","P")
    return(out.v)
  }

  m0 <- m1 <- v0 <- v1 <- n0 <- n1 <- vlr <- ts <- tp <- bp <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    naInds <- is.na(mat[i,])
    vec <- mat[i,!naInds]
    grp <- groups[!naInds]
    statDVC <- doDV(vec, grp)
    statDMC <- doTT(vec, grp)
    m0[i] <- statDVC[4]
    m1[i] <- statDVC[3]
    v0[i] <- var(vec[grp==0])
    v1[i] <- var(vec[grp==1])
    n0[i] <- sum(grp==0)
    n1[i] <- sum(grp==1)
    vlr[i] <- statDVC[1]
    ts[i]  <- statDMC[1]
    tp[i]  <- statDMC[2]
    bp[i]  <- statDVC[2]
  }

  bq  <- p.adjust(bp, "fdr")
  sig <- bp < bCut & tp < tCut
  idx <- rank(tp[sig], ties.method="first")
  rnk <- rep(NA, length(sig))
  rnk[sig] <- idx

  data.frame(mean.0=m0, mean.1=m1, var.0=v0, var.1=v1, obs.0=n0, obs.1=n1,
             var.log2.ratio=vlr, statistic.t=ts, tt.p.value=tp, bt.p.value=bp,
             bt.q.value=bq, significant=sig, rank=rnk,
             row.names=rownames(mat)
             )
}

################################################################################
#################### TEST CONSISTENCY WITH ORIGINAL VERSION ####################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X1 <- matrix(rnorm(10000,0,0.5), ncol=50)
  X2 <- matrix(rnorm(10000,0,1), ncol=50)
  X <- cbind(X1, X2)
  X[sample(length(X),500)] <- NA
  rownames(X) <- 1:nrow(X)
  groups <- rep(c(0,1), each=50)

  t1 <- base_ievora(X, groups, 0.05, 0.001)
  t2 <- row.ievora(X, groups, 0.05, 0.001)

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  x <- rnorm(12, sd=0.000001); g <- rep(c("a","b"), each=6)
  expect_equal(base_ievora(x, g=="a"), row.ievora(x, g=="a"))
})

test_that("minumum allowed sample sizes give equal results", {
  x <- rnorm(4); g <- c(0,1,0,1)
  expect_equal(base_ievora(x, g), row.ievora(x, g))
  expect_equal(base_ievora(c(x, NA), c(g, 0)), row.ievora(c(x, NA), c(g, 0)))
})

test_that("first group is treated as 0 group", {
  x <- rnorm(12); g <- rep(c("a","b"), each=6)
  expect_equal(row.ievora(x, g), row.ievora(x, as.numeric(g=="b")))
  expect_equal(row.ievora(x, g), row.ievora(x, g=="b"))
})

test_that("constant values give equal results", {
  x <- c(1,1,2,3); g <- c(0,0,1,1)
  expect_equal(base_ievora(x, g), suppressWarnings(row.ievora(x, g)))
})

