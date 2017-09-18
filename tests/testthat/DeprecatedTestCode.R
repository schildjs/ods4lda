## Run the code for one realization from simulation scenarion e) in Schildcrout et al (2015) with the exception that we
## estimate the parameter associated with C_i here.  In the paper it was not estimated.
##
## Expect Fit.Biv to take substantially longer to fit than the others
setwd("~/rsch/OutDepSamp/Software/TestODS4LDA")
rm(list=ls())
library(ODS4LDA)
odsInt <- read.csv("odsInt.csv")
odsSlp <- read.csv("odsSlp.csv")
odsBiv <- read.csv("odsBiv.csv")

date()
Fit.Int  <- acml.linear(y=odsInt$Y,
                        x=as.matrix(cbind(1, odsInt[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsInt$time)),
                        id=odsInt$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-2.569621, 9.992718),
                        SampProb=c(1, 0.1228, 1),
                        w.function="intercept")
date()
Fit.Int2  <- acml.lmem(y=odsInt$Y,
                        x=as.matrix(cbind(1, odsInt[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsInt$time)),
                        id=odsInt$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-2.569621, 9.992718),
                        SampProb=c(1, 0.1228, 1),
                        w.function="intercept")
date()
Fit.Slp  <- acml.linear(y=odsSlp$Y,
                        x=as.matrix(cbind(1, odsSlp[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsSlp$time)),
                        id=odsSlp$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-0.7488912,  3.4557775),
                        SampProb=c(1, 0.1228, 1),
                        w.function="slope")
date()
Fit.Slp2  <- acml.lmem(y=odsSlp$Y,
                        x=as.matrix(cbind(1, odsSlp[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsSlp$time)),
                        id=odsSlp$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-0.7488912,  3.4557775),
                        SampProb=c(1, 0.1228, 1),
                        w.function="slope")
date()
Fit.Biv  <- acml.linear(y=odsBiv$Y,
                        x=as.matrix(cbind(1, odsBiv[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsBiv$time)),
                        id=odsBiv$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768),
                        SampProb=c(0.122807, 1),
                        w.function="bivar")
date()
Fit.Biv2  <- acml.lmem(y=odsBiv$Y,
                        x=as.matrix(cbind(1, odsBiv[,c("time","snp","snptime","confounder")])),
                        z=as.matrix(cbind(1, odsBiv$time)),
                        id=odsBiv$id,
                        InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                        ProfileCol=NA,
                        cutpoints=c(-4.413225, 11.935188, -1.390172, 4.084768),
                        SampProb=c(0.122807, 1),
                        w.function="bivar")
date()

Fit.Int
Fit.Slp
Fit.Biv


library(MASS)
GenPopnData <- function(N = 1000, n = 11, beta = c(1, 0.25, 0, 0.25),
                        sig.b0 = 0.25, sig.b1 = 0.25, rho = 0, sig.e = 0.5,
                        pr.conf=0.25, pr.grp.0 = 0.1, pr.grp.1 = 0.2, n.low=11, n.high=11){
    id       <- rep(1:N, each=n)
    conf.i   <- rbinom(N, 1, pr.conf)
    pr.grp.i <- ifelse(conf.i==1, pr.grp.1, pr.grp.0)
    grp.i    <- rbinom(N,1,pr.grp.i)
    conf     <- rep(conf.i, each=n)
    grp      <- rep(grp.i, each=n)
    #grp      <- rep(rbinom(N, 1, prob.grp), each=n)
    time     <- rep(seq(-1,1, length.out=n), N)
    cov.mat  <- matrix(c(sig.b0^2, rho*sig.b0*sig.b1, rho*sig.b0*sig.b1, sig.b1^2),2,2)
    bi       <- mvrnorm(N, mu=c(0,0), Sigma=cov.mat)
    b        <- cbind(rep(bi[,1], each=n),rep(bi[,2], each=n))
    error    <- rnorm(N*n, 0, sig.e)


    X <- as.matrix(cbind(1, time, grp, time*grp, conf))
    Z <- X[,c(1:2)]
    Y <- X %*% beta + Z[,1]*b[,1] + Z[,2]*b[,2] + error

    ## Induce MCAR dropout
    n.obs <- rep(sample(rep(c(n.low:n.high),2), N, replace=TRUE), each=n)
    obs.num <- rep(c(1:n), N)
    X <- X[obs.num<= n.obs,]
    Z <- Z[obs.num<= n.obs,]
    Y <- Y[obs.num<= n.obs]
    id <- id[obs.num<= n.obs]
    list(id=id, X=X, Y=Y, Z=Z, ni=c(unlist(tapply(Y,id,length))),
         N=N, n=n, beta=beta,sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e, pr.grp.1=pr.grp.1, pr.grp.0=pr.grp.0)
}

LinRegFn <- function(data){ X <- cbind(1, data[,3])
Y <- data$Y
solve(t(X)%*%X) %*% t(X) %*% Y}

est.quants <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, pr.grp.0, pr.grp.1, pr.conf, quant, n.low=11, n.high=11){
    d        <- GenPopnData(N=N, n=n, beta=beta, sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e,
                            pr.conf=pr.conf, pr.grp.0 = pr.grp.0, pr.grp.1 = pr.grp.1)
    data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
    out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
    out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
    r        <- rbind( quant,
                       quantile(out[,1], probs=quant),
                       quantile(out[,2], probs=quant),
                       quantile(out[,3], probs=quant))
    r}

## Do a search to find the quantiles that correspond to the central rectangle that contains 60 and 80 percent of
## the subject specific intercepts and slopes.  This is not necessary if slope and intercepts are independent
## but with unequal followup they were positiviely correlated.  Searches for the smallest 'rectangle' defined
## by quantiles that contains 60 and 80 percent of the data
est.bivar.lims <- function(N, n, beta, sig.b0, sig.b1, rho, sig.e, pr.grp.0, pr.grp.1, pr.conf, quants, n.low=11, n.high=11){
    d        <- GenPopnData(N=N, n=n, beta=beta, sig.b0=sig.b0, sig.b1=sig.b1, rho=rho, sig.e=sig.e,
                            pr.conf=pr.conf, pr.grp.0 = pr.grp.0, pr.grp.1 = pr.grp.1)
    #d        <- GenPopnData(N, n, beta, sig.b0, sig.b1, rho, sig.e, prob.grp, n.low=n.low, n.high=n.high)
    data.tmp <- data.frame(id=d$id, Y=d$Y, X.time=d$X[,2])
    out      <- matrix(unlist(lapply(split(data.tmp, data.tmp$id), LinRegFn)), byrow=TRUE, ncol=2)
    out      <- cbind(c(tapply(data.tmp$Y, data.tmp$id, mean)), out)
    print(cor(out[,2], out[,3]))

    q1 <- .99
    Del <- 1
    while (Del>0.001){ q1 <- q1-.00025
    Del <- abs(mean(out[,2] > quantile(out[,2], probs=1-q1) &
                        out[,2] < quantile(out[,2], probs=q1) &
                        out[,3] > quantile(out[,3], probs=1-q1) &
                        out[,3] < quantile(out[,3], probs=q1)) - quants[1])
    #print(c(q1,Del))
    q1}
    q2 <- .99
    Del <- 1
    while (Del>0.001){ q2 <- q2-.00025
    Del <- abs(mean(out[,2] > quantile(out[,2], probs=1-q2) &
                        out[,2] < quantile(out[,2], probs=q2) &
                        out[,3] > quantile(out[,3], probs=1-q2) &
                        out[,3] < quantile(out[,3], probs=q2)) - quants[2])
    q2}

    rbind( c( quantile(out[,2], probs=1-q1), quantile(out[,2], probs=q1), quantile(out[,3], probs=1-q1), quantile(out[,3], probs=q1)),
           c( quantile(out[,2], probs=1-q2), quantile(out[,2], probs=q2), quantile(out[,3], probs=1-q2), quantile(out[,3], probs=q2)))

}

ODS.Sampling.Bivar <- function(dat,                      ## a list generated from the GenPopnData() function
                               PopnQuantsBivar,          ## a matrix from est.bivar.lims function
                               PopnPropInRectangle,      ## Proportion of subjects in the central rectangle
                               TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
                               SamplingStrategy){        ## Options are "IndepODS" and "DepODS"

    NCohort                    <- length(unique(dat$id))
    NStratumThry               <- round(NCohort*c(PopnPropInRectangle, (1-PopnPropInRectangle)))
    SampProbThry               <- TargetNSampledPerStratum / NStratumThry

    Lims <- (PopnPropInRectangle==.6)*PopnQuantsBivar[1,] + (PopnPropInRectangle==.8)*PopnQuantsBivar[2,]

    uid <- unique(dat$id)
    ni  <- c(unlist(tapply(dat$id, dat$id, length)))
    SampVar <- NULL
    for(i in uid){ yi      <- dat$Y[dat$id==i]
    xi      <- dat$X[dat$id==i,1:2]
    SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
    }
    print(sum(SampVar[,2]>Lims[1] & SampVar[,2]<Lims[2]))
    print(sum(SampVar[,3]>Lims[3] & SampVar[,3]<Lims[4]))

    SampStratum  <- ifelse(SampVar[,2]>Lims[1] & SampVar[,2]<Lims[2] &
                               SampVar[,3]>Lims[3] & SampVar[,3]<Lims[4], 1,2)
    print(table(SampStratum))

    NperStratum  <- unlist(tapply(uid, SampStratum, length))

    SampProbiThry <- ifelse(SampStratum==1, SampProbThry[1],
                            ifelse(SampStratum==2, SampProbThry[2], NA))

    Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
    SampProbObs <- c(tapply(Sampled, SampStratum, mean))

    SampProbiObs  <- ifelse(SampStratum==1, SampProbObs[1],
                            ifelse(SampStratum==2, SampProbObs[2], NA))

    ## Independent Sampling
    if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1]

    TheSample <- dat$id %in% InODS
    X.ods     <- dat$X[TheSample,]
    Y.ods     <- dat$Y[TheSample]
    Z.ods     <- dat$Z[TheSample,]
    id.ods    <- dat$id[TheSample]
    if (SamplingStrategy=="IndepODS"){
        SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
        SampProb.ods  <- SampProbThry
    }
    #if (SamplingStrategy=="DepODS"){
    #	   SampProbi.ods <- rep(SampProbiObs, ni)[TheSample]
    #	   SampProb.ods  <- SampProbObs
    #}
    SampStrat.ods <- SampStratum[InODS]
    Qi <- SampVar[InODS,2:3]
    dup.id <- duplicated(dat$id)
    dat.univariate <- dat$X[!dup.id,]
    dat.univ.ods   <- dat.univariate[InODS,]

    list(#X=X.ods, Y=Y.ods, Z=Z.ods, id=id.ods,
        X=dat$X, Y=dat$Y, Z=dat$Z, id=dat$id,
        SampProb=SampProb.ods, SampProbi=SampProbi.ods,
        N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0,
        sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
        prob.grp=dat$prob.grp,w.function="bivar",
        #cutpoint=c(C1,C2),
        SampStratum=SampStrat.ods, Qi=Qi, dat.univ.ods=dat.univ.ods,
        cutpoint=Lims, InSample=TheSample)
}

ODS.Sampling <- function(dat,                      ## a list generated from the GenPopnData() function
                         PopnQuants,               ## a matrix from est.quants function
                         w.function,               ## Response summary to sample on ("mean","intercept", or "slope")
                         quants,                   ## Population quantiles to define the theoretical sampling strata
                         TargetNSampledPerStratum, ## Theoretical (and possibly observed) number sampled per stratum
                         SamplingStrategy,         ## Options are "IndepODS" and "DepODS"
                         Univariate=FALSE){

    NCohort                    <- length(unique(dat$id))
    NStratumThry               <- round(NCohort*c(quants[1], quants[2]-quants[1], 1-quants[2]))
    print(TargetNSampledPerStratum)
    print(NStratumThry)
    SampProbThry               <- TargetNSampledPerStratum / NStratumThry
    C1 <- ifelse(w.function=="mean",      PopnQuants[2,match(quants[1], PopnQuants[1,])],
                 ifelse(w.function=="intercept", PopnQuants[3,match(quants[1], PopnQuants[1,])],
                        ifelse(w.function=="slope",     PopnQuants[4,match(quants[1], PopnQuants[1,])])))
    C2 <- ifelse(w.function=="mean",      PopnQuants[2,match(quants[2], PopnQuants[1,])],
                 ifelse(w.function=="intercept", PopnQuants[3,match(quants[2], PopnQuants[1,])],
                        ifelse(w.function=="slope",     PopnQuants[4,match(quants[2], PopnQuants[1,])])))

    #     C1 <- ifelse(quants[1]==.1   & w.function=="mean",       PopnQuants[1,1],
    #           ifelse(quants[1]==.1   & w.function=="intercept",  PopnQuants[2,1],
    #           ifelse(quants[1]==.1   & w.function=="slope",      PopnQuants[3,1],
    #           ifelse(quants[1]==.125   & w.function=="mean",       PopnQuants[1,2],
    #           ifelse(quants[1]==.125   & w.function=="intercept",  PopnQuants[2,2],
    #           ifelse(quants[1]==.125   & w.function=="slope",      PopnQuants[3,2],
    #           ifelse(quants[1]==.2   & w.function=="mean",       PopnQuants[1,3],
    #           ifelse(quants[1]==.2   & w.function=="intercept",  PopnQuants[2,3],
    #           ifelse(quants[1]==.2   & w.function=="slope",      PopnQuants[3,3])))))))))
    #     C2 <- ifelse(quants[2]==.8   & w.function=="mean",       PopnQuants[1,4],
    #           ifelse(quants[2]==.8   & w.function=="intercept",  PopnQuants[2,4],
    #           ifelse(quants[2]==.8   & w.function=="slope",      PopnQuants[3,4],
    #           ifelse(quants[2]==.875   & w.function=="mean",       PopnQuants[1,5],
    #           ifelse(quants[2]==.875   & w.function=="intercept",  PopnQuants[2,5],
    #           ifelse(quants[2]==.875   & w.function=="slope",      PopnQuants[3,5],
    #           ifelse(quants[2]==.9   & w.function=="mean",       PopnQuants[1,6],
    #           ifelse(quants[2]==.9   & w.function=="intercept",  PopnQuants[2,6],
    #           ifelse(quants[2]==.9   & w.function=="slope",      PopnQuants[3,6])))))))))


    uid <- unique(dat$id)
    ni  <- c(unlist(tapply(dat$id, dat$id, length)))
    SampVar <- NULL
    for(i in uid){ yi      <- dat$Y[dat$id==i]
    xi      <- dat$X[dat$id==i,1:2]
    SampVar <- rbind(SampVar, c(mean(yi), solve(t(xi)%*%xi) %*% t(xi) %*% yi))
    }
    SampVar <- (w.function=="mean")*SampVar[,1] +
        (w.function=="intercept")*SampVar[,2] +
        (w.function=="slope")*SampVar[,3]

    SampStratum  <- ifelse(SampVar<C1, 1,
                           ifelse(SampVar<C2, 2, 3))
    NperStratum  <- unlist(tapply(uid, SampStratum, length))

    SampProbiThry <- ifelse(SampVar<C1, SampProbThry[1],
                            ifelse(SampVar<C2, SampProbThry[2], SampProbThry[3]))

    Sampled     <- rbinom(length(SampProbiThry), 1, SampProbiThry)
    SampProbObs <- c(tapply(Sampled, SampStratum, mean))

    SampProbiObs  <- ifelse(SampVar<C1, SampProbObs[1],
                            ifelse(SampVar<C2, SampProbObs[2], SampProbObs[3]))
    #print(rbind(SampProbThry, SampProbObs))
    #print(cbind(SampProbiThry,SampProbiObs))


    ## Independent Sampling
    if (SamplingStrategy=="IndepODS") InODS <- uid[ Sampled==1]

    TheSample <- dat$id %in% InODS
    X.ods     <- dat$X[TheSample,]
    Y.ods     <- dat$Y[TheSample]
    Z.ods     <- dat$Z[TheSample,]
    id.ods    <- dat$id[TheSample]
    if (SamplingStrategy=="IndepODS"){
        SampProbi.ods <- rep(SampProbiThry, ni)[TheSample]
        SampProb.ods  <- SampProbThry
    }
    SampStrat.ods <- SampStratum[InODS]
    Qi <- SampVar[InODS]
    dup.id <- duplicated(dat$id)
    dat.univariate <- dat$X[!dup.id,]
    dat.univ.ods   <- dat.univariate[InODS,]

    SampProb  <- SampProbThry
    SampProbi <- rep(SampProbiThry, each=ni)
    Qi        <- SampVar

    list(X=dat$X, Y=dat$Y, Z=dat$Z, id=dat$id,
         SampProb=SampProb, SampProbi=SampProbi,
         N=dat$N, n=dat$n, beta=dat$beta, sig.b0=dat$sig.b0,
         sig.b1=dat$sig.b1, rho=dat$rho, sig.e=dat$sig.e,
         prob.grp=dat$prob.grp, w.function=w.function,
         cutpoint=c(C1,C2), SampStratum=SampStratum, Qi=Qi, InSample=TheSample)

}
#########################################################
#########################################################
#########################################################
## note that we really do not need quants, PopnQuants, and w.function but to make the fitting function work
Random.Sampling <- function(d, quants=c(.1, .9), PopnQuants, w.function="mean", n=225){
    s <- sample(unique(d$id), n)
    TheSample <- d$id %in% s
    X.rand     <- d$X
    Y.rand     <- d$Y
    Z.rand     <- d$Z
    id.rand    <- d$id

    C1 <- PopnQuants[2,match(quants[1], PopnQuants[1,])]
    C2 <- PopnQuants[2,match(quants[2], PopnQuants[1,])]

    list(X=X.rand, Y=Y.rand, Z=Z.rand, id=id.rand, InSample=TheSample,
         N=d$N, n=d$n, beta=d$beta, sig.b0=d$sig.b0,
         sig.b1=d$sig.b1, rho=d$rho, sig.e=d$sig.e,
         prob.grp=d$prob.grp,
         ## output below is only needed for the fitting function to run.  They are not used.
         SampProb=c(1,1,1), cutpoint=c(C1,C2), SampProbi=rep(1, length(Y.rand)), w.function=w.function)
}





Fit  <- acml.lmem(y=odsSlp$Y,
                       x=as.matrix(cbind(1, odsSlp[,c("time","snp","snptime","confounder")])),
                       z=as.matrix(cbind(1, odsSlp$time)),
                       id=odsSlp$id,
                       InitVals=c(5,  1, -2.5,  0.75,  0,  1.6094379,  0.2231436, -0.5108256,  1.6094379),
                       ProfileCol=NA,
                       cutpoints=c(-0.7488912,  3.4557775),
                       SampProb=c(1, 0.1228, 1),
                       w.function="slope")
date()

params <- Fit$Ests
vcovs <- Fit$covar

Xe             <- x1.insample[,3]
Xo             <- x1.insample[,5]
Xe.Xo.Seq1.mod <- glm(Xe ~ Xo, family=binomial)
print(paste(date(), "Start MI", sep="-----"))
Ests.mi <- Var.mi <- NULL

Fit  <- ACML.LME(y=y.insample, x=x.insample, z=z.insample, id=id.insample, InitVals=inits, ProfileCol=ProfileCol.s,
                 cutpoints=CutPoint.s, SampProb=SampProb.s, w.function=w.function.s)


## Get estimates and variance-covariance matrix
if (is.na(ProfileCol.s)){ n.par <- length(Fit$Ests)
params <- Fit$Ests
vcovs <- Fit$covar
for (m1 in 2:n.par){ for (m2 in 1:(m1-1)){ #print(c(m1,m2))
    vcovs[m1,m2] <- vcovs[m2,m1] }}
}else{ n.par <- length(Fit$Ests)
params <- c(Fit$Ests[1:(ProfileCol.s-1)], 0, Fit$Ests[ProfileCol.s:n.par])
vcovs  <- cbind(Fit$covar[,1:(ProfileCol.s-1) ], 0, Fit$covar[,ProfileCol.s:n.par])
vcovs  <- rbind(vcovs[1:(ProfileCol.s-1), ], 0, vcovs[ProfileCol.s:n.par,])
for (m1 in 2:(n.par+length(ProfileCol.s))){ for (m2 in 1:(m1-1)){ #print(c(m1,m2))
    vcovs[m1,m2] <- vcovs[m2,m1] }}
}

Xe             <- x1.insample[,3]
Xo             <- x1.insample[,5]
Xe.Xo.Seq1.mod <- glm(Xe ~ Xo, family=binomial)
print(paste(date(), "Start MI", sep="-----"))
Ests.mi <- Var.mi <- NULL
for (j in 1:n.imp){ #print(j)

    ## Draw a sample from the parameter estimate distribution for the model Xe | Xo, S=1
    Xe.Xo.Seq1.MI.params <- rmvnorm(1, Xe.Xo.Seq1.mod$coef, summary(Xe.Xo.Seq1.mod)$cov.unscaled)
    prXeEq1.Xo <- expit(cbind(1,x1.outsample[,5]) %*% t(Xe.Xo.Seq1.MI.params))

    ## Draw a sample from the parameter estimate distribution for the model Y | Xe, Xo, S=1
    Y.Xe.Xo.Seq1.MI.params <- rmvnorm(1, params, vcovs)

    ## Condtl log likelihood if Xe=0: log(pr(y | Xo, Xe=0, S=1))
    x.outsample.0 <- x.outsample
    x.outsample.0[,3] <- 0
    x.outsample.0[,4] <- x.outsample.0[,2] * x.outsample.0[,3]
    prY.Xo.S1.Xe0 <- LogLikeAndScore(params= Y.Xe.Xo.Seq1.MI.params, y=y.outsample, x=x.outsample.0, z=z.outsample, id=id.outsample,
                                     w.function=w.function.s, cutpoints=CutPoint.s, SampProb = SampProb.s, SampProbi=rep(1, length(y.outsample)),
                                     ##this w.function, cutpoints, SampProb, and SampProbi dont really factor into the calculation, but need tospecify them for the function to work
                                     ProfileCol=ProfileCol.s, Keep.liC=TRUE)
    ## Condtl log likelihood if Xe=1: log(pr(y | Xo, Xe=1, S=1))
    x.outsample.1 <- x.outsample
    x.outsample.1[,3] <- 1
    x.outsample.1[,4] <- x.outsample.1[,2] * x.outsample.1[,3]
    prY.Xo.S1.Xe1  <- LogLikeAndScore(params= Y.Xe.Xo.Seq1.MI.params, y=y.outsample, x=x.outsample.1, z=z.outsample, id=id.outsample,
                                      w.function=w.function.s, cutpoints=CutPoint.s, SampProb = SampProb.s, SampProbi=rep(1, length(y.outsample)),
                                      ##this w.function, cutpoints, SampProb, and SampProbi dont really factor into the calculation, but need tospecify them for the function to work
                                      ProfileCol=ProfileCol.s, Keep.liC=TRUE)

    odds <- exp(prY.Xo.S1.Xe1)*prXeEq1.Xo/(exp(prY.Xo.S1.Xe0)*(1-prXeEq1.Xo))
    probs <- Odds2Prob(odds)

    Impute.Xe.1 <- rbinom(length(probs), 1, probs)
    y.mi  <- Sampled.InAndOut$Y
    x.mi  <- Sampled.InAndOut$X
    z.mi  <- Sampled.InAndOut$Z
    id.mi <- Sampled.InAndOut$id

    for (j in id1.outsample){ x.mi[id.mi == j,3] <- Impute.Xe.1[id1.outsample==j] }
    x.mi[,4] <- x.mi[,2]*x.mi[,3]

    lme.mi  <- lmer(y.mi ~ x.mi[,2]+x.mi[,3]+x.mi[,4]+x.mi[,5]+(z.mi[,2] | id.mi), REML=FALSE)
    Ests.mi <- rbind(Ests.mi, lme.mi@fixef )
    Var.mi  <- rbind(Var.mi,diag(vcov(lme.mi)) )
    #    ACML.mi  <- ACML.LME(y=y.mi, x=x.mi, z=z.mi, id=id.mi, InitVals=ACML$Ests)
    #    Ests.mi <- rbind(Ests.mi,  ACML.mi$Ests)
    #    Var.mi <- rbind(Var.mi, diag(ACML.mi$covar))
}


