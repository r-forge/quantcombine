
discoChiSq <- function(freqs) {

    ##check that labels present
    names = names(freqs)
    if(is.null(names)) {
        stop("Input freqs table must have column labels, which are the scores calculated by the getQuantileScores function, e.g. -2,-1,0,1,2 if you used n.quantiles=4 as the fourth argument to getQuantileScores. If you've read a text or CSV file into R and there are X's or V's in front of the numbers, or the negative signs are replaced by periods, this is okay, since discoChiSq will fix these.")
    }
    
    ##check that freqs is a data frame; convert to if otherwise
    if(!is.data.frame(freqs)) {

        if(is.matrix(freqs)) {
            freqs.df <- data.frame(freqs,row.names=dimnames(freqs)[[1]])
            names(freqs.df) <- dimnames(freqs)[[2]]
            freqs <- freqs.df
        }
        else {
            stop("Please specify a data frame or matrix for the freqs argument.")
        }
    }


    
    ##check which columns have data and which have IDs
    f1 = function(x) {is.factor(freqs[1,x])}
    freqs.IDs = sapply(1:ncol(freqs),f1)
    
    ##store IDs
    freqs.IDcols = data.frame(freqs[,freqs.IDs],row.names=row.names(freqs))
    names(freqs.IDcols) = names(freqs)[freqs.IDs]
    
    ##extract only freqs
    freqs.only = freqs[,!freqs.IDs]
    names(freqs.only) = names(freqs)[!freqs.IDs]
    
    ##if loaded freqs from file, labels will likely have X's or V's at front and . for -
    names = names(freqs.only)


    f2 = function(i){
        parts = strsplit(names[i],"\\.")[[1]]

        if(length(parts) == 2) {

            if(nchar(parts[1]) == 0) {
                return(-1 * as.integer(parts[2]))
            }
            else {
                return(as.integer(parts[1]))
            }
        }
        else {
            if(length(parts) == 3) {
                return(-1 * as.integer(parts[2]))
            }
            else {
                return(as.integer(names[i]))
            }
        }
    }
    
    
    if(length(grep("X",names))>0){
        names = sapply(1:length(names),function(i){strsplit(names[i],"X")[[1]][2]})

        if(length(grep("\\.",names))>0) {
            names = sapply(1:length(names),f2)
            
            #names = sub("\\.","-",names)
        }
    }
    if(length(grep("V",names))>0){
        names = sapply(1:length(names),function(i){strsplit(names[i],"V")[[1]][2]})

        if(length(grep("\\.",names))>0) {
            names = sapply(1:length(names),f2)
            #names = sub("\\.","-",names)
        }
    }
        
    ##shouldn't be any NAs if freqs generated from the package
    if(any(is.na(freqs.only))) {
        stop("Missing values are present in the freqs data. Missing values can be present in the raw expression data fed to get.quantile.scores, but not here. If your data was generated by get.quantile.freqs, and you get this message, this is a bug. Please contact the QuantCombine package maintainer.")
    }
    

        
    ###END OF INPUT CHECKS###

    if(!is.loaded("disco_chisq")) {
        
        lib.file <- file.path(paste("QuantCombine", .Platform$dynlib.ext,sep=""))
        dyn.load(lib.file)
        cat("Loaded ",lib.file,"\n")
    }

    labels <- sort(unique(as.numeric(names)))
    n.labels <- length(labels)

    
    output = .C("disco_chisq",
    freqs=as.integer(as.vector(t(freqs.only))),
    n.genes=as.integer(nrow(freqs.only)),
    labels=as.integer(labels),
    n.labels=as.integer(n.labels),
    res=as.double(rep(0.0,9*nrow(freqs.only))),
    PACKAGE="QuantCombine"
    )

    results = data.frame(
               t(sapply(0:(nrow(freqs.only)-1),
                      function(g){
                          c(
                               "p.value"=output$res[(1+g*9)],
                               "chisq"=output$res[(2+g*9)],
                               "beta1-1"=output$res[(3+g*9)],
                               "beta1-2"=output$res[(4+g*9)],
                               "beta1-3"=output$res[(5+g*9)],
                               "beta1-4"=output$res[(6+g*9)],
                               "beta2-1"=output$res[(7+g*9)],
                               "beta2-2"=output$res[(8+g*9)],
                               "beta2-3"=output$res[(9+g*9)]
                               )
                      }
                      )
                 ),
               row.names=row.names(freqs.only)
               )

    ##add IDs to end of results
    results.ncol = ncol(results)
    freqs.IDs.sum = sum(freqs.IDs==TRUE)
    
    if(freqs.IDs.sum > 0) {
        for(i in seq(1,freqs.IDs.sum)) {

            results = cbind(results,freqs.IDcols[,i])
            n = names(freqs.IDcols)[i]
            attributes(results)$names[ncol(results)] = n
        }
    }

    return(results)
}
