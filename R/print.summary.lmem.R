print.summary.lmem <-
    function(x,...) {
        cat('\nClass:\n',x$class,'\n',sep='')
        cat("\nCall:\n", paste(x$call, sep="\n", collapse="\n"), "\n", sep="")
        #cat('\nInformation Criterion:\n')
        #print.default(format(x$info, ...), print.gap = 2L, quote = FALSE)
        cat("\nMarginal Mean Parameters:\n")
        printCoefmat(x$mean.table,signif.stars = FALSE)
        cat('\n')
        cat("Association Parameters:\n")
        printCoefmat(x$assoc.table,signif.stars = FALSE)
        cat('\n')
        cat('Number of clusters:            ',x$n.clusters,'\n')
        cat('Number of observations:        ',x$n.obs,'\n')
        cat('Maximum cluster size:          ',x$max.cluster.size,'\n')
        cat('Convergence status (nlm code): ',x$code,'\n')
        cat('Number of iterations:          ',x$n.iter)
        cat('\n')
    }
