
disco <- function(freqs) {
    
    if(!is.loaded("disco_chisq")) {
        
        lib.file <- file.path(paste("QuantCombine", .Platform$dynlib.ext,sep=""))
        dyn.load(lib.file)
        cat("Loaded ",lib.file,"\n")
    }

    labels <- sort(unique(as.numeric(names(freqs))))
    n.labels <- length(labels)
    #res = rep(0.0,9*nrow(freqs))
    
    output = .C("disco_chisq",
    freqs=as.integer(as.vector(t(freqs))),
    n.genes=as.integer(nrow(freqs)),
    labels=as.integer(labels),
    n.labels=as.integer(n.labels),
    res=as.double(rep(0.0,9*nrow(freqs)))
    )

    lapply(0:(nrow(freqs)-1),function(g){list("chisq"=output$res[1+(g*9)],"p.value"=output$res[2+(g*9)],"beta1"=c(output$res[(3+g*9):(6+g*9)]),"beta2"=c(output$res[(7+g*9):(9+g*9)]))}) 
}
