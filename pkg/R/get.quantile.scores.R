#######################################################################
# File: get.quantile.scores.R
#
# Author: Peter S. Bazeley, University of Toledo
#
# Description: Code for
#
# The problem:
#
# Contains: get.quantile.scores,write.quantile.scores
#######################################################################



#######################################################################
# Function: get.quantile.scores
#
# Description: For each gene
#
# Input:
#        x <- Data frame for 1 dataset with genes as rows and samples
#            as columns
#        grp1 = vector indicating group 1 sample columns
#        grp2 = vector indicating group 1 sample columns
#
# Output:
#         object for 1 dataset containing:
#             data = data frame with genes as rows and quantile
#                    frequencies (scores) as columns, 1 group after
#                    the other
#             score.labels = vector indicating score labels
#             group.labels = vector indicating group labels
#
#######################################################################


#Setup quantiles and groups and get group frequencies

get.quantile.scores <- function(data,grp1,grp2,n.quantiles = 4) {
    
    n.labels <- n.quantiles + 1
    seq <- seq(length.out = n.labels) #1,2,3,...
    labels <- seq - median(seq) #shift to -2,-1,0,...
    
    ## if num.labels is odd, will get half values from median subtraction, so add 0.5
    if(round(labels[1]) != labels[1]) {
        labels <- labels + 0.5	
    }
	#print(labels)
    quantile.unit <- 1/n.labels
    quantile.thresholds <- sapply(1:n.quantiles,function(x){quantile.unit * x})
    
    f1 <- function(gene.index) {
        #cat(paste("gene.index is ",gene.index,"\n",sep=""))

        
        ##calculate quantiles for gene
        
        quantiles <- quantile(as.numeric(data[gene.index,]),quantile.thresholds)
	#print(quantiles)        

        f2 <- function(grp) {

            #print(grp)
            
            #cat(paste("grp is ",grp,"\n",sep=""))

            grp.exprs <- as.numeric(data[gene.index,grp])
            #print(grp.exprs) ##
            ##intialize matrix to store adjusted expression values
            #values <- matrix(0,ncol=1,nrow=length(grp))
            values <- matrix(0,ncol = length(grp),nrow = 1)

            #print(values)
            
            values[grp.exprs < quantiles[1]] <- labels[1]
            
            #print(values)

            index <- 2
	    for(index in c(2:n.quantiles)) {
                
                values[grp.exprs >= quantiles[index - 1] & grp.exprs < quantiles[index] ] <- labels[index]
		
                #print(values)

            }
            
            values[grp.exprs >= quantiles[n.quantiles]] <- labels[n.labels]
	    
            #print(values)
            
            values
        }
        c(f2(grp1),f2(grp2))
    }
    
    df <- data.frame(t(sapply(row.names(data),f1)),row.names = row.names(data))
    names(df) <- names(data)
    df
}


#######################################################################
# Function: write.quantile.scores
#
# Description: For each gene, 
#
# Input:
#        x = Data frame for 1 dataset with genes as rows and samples
#            as columns
#        grp1 = vector indicating group 1 sample columns
#        grp2 = vector indicating group 1 sample columns
#
# Output:
#         object for 1 dataset containing:
#             data = data frame with genes as rows and quantile
#                    frequencies (scores) as columns, 1 group after
#                    the other
#             score.labels = vector indicating score labels
#             group.labels = vector indicating group labels
#
#######################################################################









#end of get.quantile.scores.R
