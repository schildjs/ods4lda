
test_that("Acceptable Bivariate Fit Absolute Error",
{
  Fit.Biv  <- acml.linear(y=odsBiv$Y,
                          x=as.matrix(cbind(1, odsBiv[,c("time","snp","snptime","confounder")])),
                          z=as.matrix(cbind(1, odsBiv$time)),
                          id=odsBiv$id,
                          InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                          ProfileCol=NA,
                          cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768),
                          SampProb=c(0.122807, 1),
                          w.function="bivar")

  Fit.Biv2 <- acml.lmem(y=odsBiv$Y,
                        x=as.matrix(cbind(1, odsBiv[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsBiv$time)),
                        id=odsBiv$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768),
                        SampProb=c(0.122807, 1),
                        w.function="bivar")

  expect_true( all(abs(Fit.Biv$est - Fit.Biv2$est) < 1e-10) )
})
  