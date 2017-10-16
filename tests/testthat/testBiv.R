
test_that("Acceptable Bivariate Fit Absolute Error",
{
  Fit.Biv  <- acml.linear(y=odsBiv$Y,
                          x=as.matrix(cbind(1, odsBiv[,c("time","snp","confounder","snptime")])),
                          z=as.matrix(cbind(1, odsBiv$time)),
                          id=odsBiv$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768),
                          SampProb=c(0.122807, 1),
                          w.function="bivar")
  
  odsBiv$SampProbi <- 1
  Fit.Biv2 <- acml.lmem(Y~time*snp+confounder,
                          ~time,
                          data=odsBiv,
                        id=id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768),
                        SampProb=c(0.122807, 1),
                        SampProbiWL=SampProbi,
                        w.function="bivar")
  
  expect_true(length(Fit.Biv$Ests) == 9)

  expect_true( all(abs(Fit.Biv$Ests - Fit.Biv2$Ests) < 1e-10) )
  
  expect_true( all(diag(abs(Fit.Biv$cov - Fit.Biv2$cov)) < 1e-8) )
})
  