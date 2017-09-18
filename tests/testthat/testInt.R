
test_that("Acceptable Intercept Fit Absolute Error",
{
  Fit.Int  <- acml.linear(y=odsInt$Y,
                          x=as.matrix(cbind(1, odsInt[,c("time","snp","snptime","confounder")])),
                          z=as.matrix(cbind(1, odsInt$time)),
                          id=odsInt$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-2.569621, 9.992718),
                          SampProb=c(1, 0.1228, 1),
                          w.function="intercept")

  Fit.Int2  <- acml.lmem(y=odsInt$Y,
                          x=as.matrix(cbind(1, odsInt[,c("time","snp","snptime","confounder")])),
                          z=as.matrix(cbind(1, odsInt$time)),
                          id=odsInt$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-2.569621, 9.992718),
                          SampProb=c(1, 0.1228, 1),
                          w.function="intercept")

    expect_true( all(abs(Fit.Int$est - Fit.Int2$est) < 1e-10) )
})
  