summary.lmem <- summary.lmem2 <-
    function(object, ...) {

    oo      <- object
    L.coef <- length(oo$coefficients)
    L.beta <- L.coef-4
    m.i     <- seq(L.beta)
    a.i     <- seq(L.coef)
    vc.i    <- a.i[-m.i]
    mean.oo        <- with(oo, data.frame(estimate = coefficients[m.i], mod.se = sqrt(diag(covariance)[m.i])))
    mean.oo$X2     <- with(mean.oo, (estimate/mod.se)^2)
    mean.oo$pX2    <- with(mean.oo, pchisq(X2, df=1,lower.tail=FALSE))
    names(mean.oo) = c('Estimate','Model SE', 'Chi Square','Pr(>Chi)')

    assoc.oo <- with(oo, data.frame(estimate = coefficients[vc.i], mod.se = sqrt(diag(covariance)[vc.i])))
    assoc.oo$X2  <- with(assoc.oo, (estimate/mod.se)^2)
    assoc.oo$pX2 <- with(assoc.oo, pchisq(X2, df=1,lower.tail=FALSE))
    names(assoc.oo) = c('Estimate','Model SE','Chi Square','Pr(>Chi)')

    if (length(unique(attr(oo,'args')[["Weights"]])) !=1){

        mean.oo <- with(oo, data.frame(estimate = beta, rob.se = sqrt(diag(robcov)[m.i])))
        mean.oo$X2  <- with(mean.oo, (estimate/rob.se)^2)
        mean.oo$pX2 <- with(mean.oo, pchisq(X2, df=1,lower.tail=FALSE))
        names(mean.oo) = c('Estimate','Robust SE', 'Chi Square','Pr(>Chi)')

        assoc.oo <- with(oo, data.frame(estimate = alpha, rob.se = sqrt(diag(robcov)[a.i])))
        assoc.oo$X2  <- with(assoc.oo, (estimate/rob.se)^2)
        assoc.oo$pX2 <- with(assoc.oo, pchisq(X2, df=1,lower.tail=FALSE))
        names(assoc.oo) = c('Estimate','Robust SE','Chi Square','Pr(>Chi)')
        warning('When performing a weighted likelihood analysis (by specifying non-constant weights), robust standard errors are reported. Model based standard errors will not be correct and should not be used.')
        }

        out = list(class = class(oo), call = oo$call, #control=object$control,
                   #info=object$info_stats,
                   mean.table = mean.oo, assoc.table=assoc.oo,
                   Log.likelihood=oo$logLik,
                   code = oo$control[["code"]],
                   n.iter = oo$control[["n.iter"]],
                   n.clusters = oo$control[["n.clusters"]],
                   n.obs = oo$control[["n.obs"]],
                   max.cluster.size = oo$control[["max.cluster.size"]])
        class(out) = 'summary.acml.lmem'
        out
    }
