
get.quantile.scores.C <- function(exprs,grp1,grp2,n.quantiles=4) {

    ##check that n.quantiles is an integer
    is.wholenumber <- function(x,tol=.Machine$double.eps^0.5) abs(x - round(x)) < tol
    if(!is.wholenumber(n.quantiles)) {
        stop("Please specify a whole number value for n.quantiles.")
    }

    ##check that n.quantiles is positive
    if(n.quantiles < 0) {
        stop("Please specify a value > 0 for n.quantiles.")
    }
    
    ##check that exprs is a data frame; convert to if otherwise
    if(!is.data.frame(exprs)) {

        if(is.matrix(exprs)) {
            exprs.df <- data.frame(exprs,row.names=dimnames(exprs)[[1]])
            names(exprs.df) <- dimnames(exprs)[[2]]
            exprs <- exprs.df
        }
        else {
            ##can it be a vector and have a name?
            stop("Please specify a data frame or matrix for the exprs argument.")
        }
    }
    
    ##check that grp1 and grp2 are numeric
    ##if(!is.numeric(grp1)) {

        ##maybe they were stored as factors by mistake
    ##    grp1.n = as.numeric(grp1)
    ##    if (any(is.na(grp1.n))) {
    ##        stop("Please specify grp1 as a numeric index or index range, e.g. 1:4.")
    ##    }
    ##    else {
    ##        grp1 = grp1.n
    ##    }
    ##}
    ##if(!is.numeric(grp2)) {
        
        ##maybe they were stored as factors by mistake
    ##    grp2.n = as.numeric(grp2)
    ##    if (any(is.na(grp2.n))) {
    ##        stop("Please specify grp2 as a numeric index or index range, e.g. 5:10.")
    ##    }
    ##    else {
    ##        grp2 = grp2.n
    ##    }
    ##}
    
    ##check that grp1 and grp2 are separate, and grp2 follows grp1
    ##if(max(grp1) > min(grp2)) {
    ##    exprs <- cbind(exprs[,grp1],exprs[,grp2])
        
        ##update grp1 and grp2 indices - only lengths needed below
        ##grp1 = seq(from=1,to=length(grp1))
        ##grp2 = seq(from=(length(grp1)+1),to=(length(grp1)+length(grp2)))

    ##    warn("NOTE: Rearranging samples in output so grp1 and grp2 are separate, and grp2 follows grp1.")
    ##}

    exprs <- try(cbind(exprs[,grp1],exprs[,grp2]))
    if(class(exprs) == "try-error") {
        stop("Error in 'cbind(exprs[,grp1],exprs[,grp2])'. Please be sure that grp1 and grp2 are valid column labels/indices.")
    }
    
    ##check for NAs in data
    if(any(is.na(exprs))) {
        stop("NOTE: Missing values were found in the exprs data. This is okay; these will be given a value of -1000 by get.quantile.scores, and will be ignored by get.quantile.freqs.")
    }

    ##check that values are all numbers?... not likely that they aren't

    
    ###END OF INPUT CHECKS###
    
    if(!is.loaded("get_quantile_scores")) {
        
        lib.file <- file.path(paste("QuantCombine", .Platform$dynlib.ext,sep=""))
        dyn.load(lib.file)
        cat("Loaded ",lib.file,"\n")
    }


    n.labels <- n.quantiles + 1
    seq <- seq(length.out = n.labels) #1,2,3,...
    labels <- seq - median(seq) #shift to -2,-1,0,...
    
    ## if num.labels is odd, will get half values from median subtraction, so add 0.5
    if(round(labels[1]) != labels[1]) {
        labels <- labels + 0.5	
    }
    ##print(labels)
    quantile.unit <- 1/n.labels
    quantile.thresholds <- sapply(1:n.quantiles,function(x){quantile.unit * x})

    grp1.length <- length(grp1)
    grp2.length <- length(grp2)
    n.scores <- grp1.length + grp2.length
    
    #should i use ncol(data)

    output <- .C("get_quantile_scores",
                 exprs = as.double(t(exprs)),
                 n.genes = as.integer(nrow(exprs)),
                 n.grp1 = as.integer(grp1.length),
                 n.grp2 = as.integer(grp2.length),
                 labels = as.integer(labels),
                 n.labels = as.integer(n.labels),
                 quantile.thresholds = as.double(quantile.thresholds),
                 scores = as.integer(rep(0,n.scores*nrow(exprs))),
                 PACKAGE="QuantCombine"
                 )

    df <- data.frame(matrix(output$scores,ncol=n.scores,byrow=TRUE),row.names = row.names(exprs))
    names(df) <- names(exprs)
    return(df)
}
