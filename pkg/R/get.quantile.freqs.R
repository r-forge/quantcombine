#######################################################################
## File: get.quantile.freqs.R
##
## Author: Peter S. Bazeley, University of Toledo
##
## Description: Code for
##
## The problem:
##
## Contains: get.quantile.freqs
#######################################################################



#######################################################################
## Function: get.quantile.freqs
##
## Description: For each gene
##
## Input:
##        x = Data frame for 1 dataset with genes as rows and sample
##            scores as columns
##        grp1 = vector indicating group 1 sample columns
##        grp2 = vector indicating group 1 sample columns
##
## Output:
##        data = data frame with genes as rows and quantile frequencies
##               (scores) as columns, 1 group after the other
##
#######################################################################


get.quantile.freqs <- function(data,grp1,grp2) {

    #num.quantiles
    #num.labels <- num.quantiles + 1
    #seq <- seq(length.out = num.labels) #1,2,3,...
    #labels <- seq - median(seq) #shift to -2,-1,0,...
    
    ## if num.labels is odd, will get half values from median subtraction, so add 0.5
    #if(round(labels[1]) != labels[1]) {
    #    labels <- labels + 0.5	
    #}
	
    #quantile.unit <- 1/num.labels
    #quantile.thresholds <- sapply(1:num.quantiles,function(x){quantile.unit*x})

    labels.range <- range(data[1,]) #each gene will always have a max and min value
    labels <- seq(from = labels.range[1], to = labels.range[2])

    
    # <- function(gene.index,grp,quantile.label) { sum(data[gene.index,grp] == quantile.label)}          
    ##sapply(row.names(data),freqs

    #for(gene.index in row.names(data)) {

     #   for(grp in 
    #}

    
    f1 <- function(gene.index) {
        #cat(paste("gene.index is ",gene.index,"\n",sep=""))

        
        f2 <- function(grp) {

            #print(grp)
            
            #cat(paste("grp is ",grp,"\n",sep=""))
            freqs <- function(quantile.label) { sum(data[gene.index,grp] == quantile.label)}          
            
            sapply(labels,freqs)
        }
        c(f2(grp1),f2(grp2))
    }

    df <- data.frame(t(sapply(row.names(data),f1)),row.names = row.names(data))
    names(df) <- rep(labels,2)
    df
}

#end of get.quantile.freqs.R
