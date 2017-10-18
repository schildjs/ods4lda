
test_that("Acceptable Intercept Fit Absolute Error",
{
  Fit.Int  <- acml.linear(y=odsInt$Y,
                          x=as.matrix(cbind(1, odsInt[,c("time","snp","confounder", "snptime")])),
                          z=as.matrix(cbind(1, odsInt$time)),
                          id=odsInt$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-2.569621, 9.992718),
                          SampProb=c(1, 0.1228, 1),
                          w.function="intercept")
  
  
  odsInt$SampProbi <- 1
  Fit.Int2  <- acml.lmem(Y~time*snp+confounder,
                          ~time,
                          data=odsInt,
                          id=id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-2.569621, 9.992718),
                          SampProb=c(1, 0.1228, 1),
                          SampProbiWL =SampProbi,
                          w.function="intercept")

  expect_true(length(Fit.Int$Ests) == 9)
  
  expect_true( all(abs(Fit.Int$Ests - Fit.Int2$Ests) < 1e-10) )
    
  expect_true( all(diag(abs(Fit.Int$covar - Fit.Int2$covar)) < 1e-8) )
  
  expect_true( all(diag(Fit.Int2$covar) > 0) )
})
  