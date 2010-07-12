
get.quantile.freqs.C <- function(scores,grp1,grp2) {


    ##check that scores is a data frame; convert to if otherwise
    if(!is.data.frame(scores)) {

        if(is.matrix(scores)) {
            scores.df <- data.frame(scores,row.names=dimnames(scores)[[1]])
            names(scores.df) <- dimnames(scores)[[2]]
            scores <- scores.df
        }
        else {
            ##can it be a vector and have a name?
            stop("Please specify a data frame or matrix for the scores argument.")
        }
    }
   
    ##extract only specified columns
    scores <- try(cbind(scores[,grp1],scores[,grp2]))
    if(class(scores) == "try-error") {
        stop("Error in 'cbind(scores[,grp1],scores[,grp2])'. Please be sure that grp1 and grp2 are valid column labels/indices.")
    }
    
    ##check for NAs in data
    if(any(is.na(scores))) {
        stop("NAs not allowed...")
    }


    
    ###END OF INPUT CHECKS###
    
    if(!is.loaded("get_quantile_freqs")) {
        
        lib.file <- file.path(paste("QuantCombine", .Platform$dynlib.ext,sep=""))
        dyn.load(lib.file)
        cat("Loaded ",lib.file,"\n")
    }

    labels.range <- range(scores[1,]) #each gene will always have a max and min value
    labels <- seq(from = labels.range[1], to = labels.range[2])
    n.labels = length(labels)
    n.freqs = 2 * n.labels
    
    output <- .C("get_quantile_freqs",
                 scores = as.integer(t(scores)),
                 n.genes = as.integer(nrow(scores)),
                 n.grp1 = as.integer(length(grp1)),
                 n.grp2 = as.integer(length(grp2)),
                 labels = as.integer(labels),
                 n.labels = as.integer(n.labels),
                 freqs = as.integer(rep(0,n.freqs*nrow(scores))),
                 PACKAGE="QuantCombine"
                 )

    df <- data.frame(matrix(output$freqs,ncol=n.freqs,byrow=TRUE),row.names = row.names(scores))
    names(df) <- rep(labels,2)
    df
}
