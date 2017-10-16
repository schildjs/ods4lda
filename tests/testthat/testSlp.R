test_that("Acceptable Slope Fit Absolute Error",
{
  Fit.Slp  <- acml.linear(y=odsSlp$Y,
                          x=as.matrix(cbind(1, odsSlp[,c("time","snp","confounder","snptime")])),
                          z=as.matrix(cbind(1, odsSlp$time)),
                          id=odsSlp$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-0.7488912,  3.4557775),
                          SampProb=c(1, 0.1228, 1),
                          w.function="slope")
  odsSlp$SampProbi <- 1
  Fit.Slp2  <- acml.lmem(Y~time*snp+confounder,
                          ~time,
                          data=odsSlp,
                          id=id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-0.7488912,  3.4557775),
                          SampProb=c(1, 0.1228, 1),
                          SampProbiWL =SampProbi,
                          w.function="slope")

  expect_true(length(Fit.Slp$Ests) == 9)
  
  expect_true( all(abs(Fit.Slp$Ests - Fit.Slp2$Ests) < 1e-10) )
    
  expect_true( all(diag(abs(Fit.Slp$cov - Fit.Slp2$cov)) < 1e-8) )
})
  