test_that("Acceptable Slope Fit Absolute Error",
{
  Fit.Slp  <- acml.linear(y=odsSlp$Y,
                          x=as.matrix(cbind(1, odsSlp[,c("time","snp","snptime","confounder")])),
                          z=as.matrix(cbind(1, odsSlp$time)),
                          id=odsSlp$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-0.7488912,  3.4557775),
                          SampProb=c(1, 0.1228, 1),
                          w.function="slope")

  Fit.Slp2  <- acml.lmem(y=odsSlp$Y,
                          x=as.matrix(cbind(1, odsSlp[,c("time","snp","snptime","confounder")])),
                          z=as.matrix(cbind(1, odsSlp$time)),
                          id=odsSlp$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-0.7488912,  3.4557775),
                          SampProb=c(1, 0.1228, 1),
                          w.function="slope")

  expect_true( all(abs(Fit.Slp$est - Fit.Slp$est) < 1e-10) )
})
  