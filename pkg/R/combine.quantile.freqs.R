#######################################################################
## File: combine.quantile.freqs.R
##
## Author: Peter S. Bazeley, University of Toledo
##
## Description: Code for
##
## The problem:
##
## Contains: combine.quantile.freqs
#######################################################################



#######################################################################
## Function: combine.quantile.freqs
##
## Description: For each gene
##
## Input:
##        d1.fr = Data frame for 1st dataset with genes as rows and score
##                frequencies as columns
##        d2.fr = Data frame for 2nd dataset with genes as rows and score
##                frequencies as columns
##        ID.map = Data frame with paired, mapped IDs for probes.
##
## Output:
##        data = data frame with genes as rows and quantile frequencies
##               (scores) as columns, 1 group after the other. Row 
##               names are, and any additional IDs are added as
##               additional columns on right.
##
#######################################################################


combine.quantile.freqs <- function(d1.fr,d2.fr,ID.map) {

    ##check that inputs are data frames; convert to if otherwise
    if(!is.data.frame(d1.fr)) {

        if(is.matrix(d1.fr)) {
            d1.fr.df = data.frame(d1.fr,row.names=dimnames(d1.fr)[[1]])
            names(d1.fr.df) = dimnames(d1.fr)[[2]]
            d1.fr = d1.fr.df
        }
        else {
            ##can it be a vector and have a name?
            stop("First frequency data table must be data frame or matrix.")
        }
    }
    if(!is.data.frame(d2.fr)) {
        
        if(is.matrix(d2.fr)) {
            d2.fr.df = data.frame(d2.fr,row.names=dimnames(d2.fr)[[1]])
            names(d2.fr.df) = dimnames(d2.fr)[[2]]
            d2.fr = d2.fr.df
        }
        else {
            ##can it be a vector and have a name?
            stop("Second frequency data table must be data frame or matrix.")
        }
    }
    if(!is.data.frame(ID.map)) {
        
        if(is.matrix(ID.map)) {
            ID.map.df = data.frame(ID.map,row.names=dimnames(ID.map)[[1]])
            names(ID.map.df) = dimnames(ID.map)[[2]]
            ID.map = ID.map.df
        }
        else {
            stop("ID map must be data frame or matrix.")
        }
    }

    ##check that inputs have gene names

    ##check that inputs have same number of columns and column labels

    
    ##check that d1 IDs match ID.map[,1]
    d1.rown = row.names(d1.fr) %in% row.names(ID.map)
    col1 = row.names(d1.fr) %in% ID.map[,1]
    col2 = row.names(d1.fr) %in% ID.map[,2]

    ##NOTE: i should make it order agnostic...?
    ##then make common column optional

    col1.sum = sum(col1==TRUE)
    rown.sum = sum(d1.rown==TRUE)
    
    if(col1.sum > 0) {

        if(rown.sum > 0) {
            stop("Both ID.map's row names and first column contain IDs that match the first frequency table's IDs. Please use a data frame or matrix that only contains one set of IDs found in the input frequency table.")
        }
    }
    else if (rown.sum > 0) {

    }
    else {
        stop("Neither ID.map's row names nor first column match first frequency table's row names.")
    }
    
    ##subset common IDs
    ##NOTE: verify that order in ID map prevails in subset
    d1.fr.common = d1.fr[as.character(ID.map.df[,1]),]
    d2.fr.common = d2.fr[as.character(ID.map.df[,2]),]

    ## if ID.map longer than data, will show up as NAs at bottom, so remove these
    d1.fr.common = d1.fr.common[!is.na.data.frame(d1.fr.common)[,1],]
    d2.fr.common = d2.fr.common[!is.na.data.frame(d2.fr.common)[,1],]
    
    ##check which columns have data and which have IDs
    f1 = function(x) {is.factor(d1.fr.common[1,x])}
    f2 = function(x) {is.factor(d2.fr.common[1,x])}

    d1.fr.common.IDs = sapply(1:ncol(d1.fr.common),f1)
    d2.fr.common.IDs = sapply(1:ncol(d2.fr.common),f2)

    
    ##store IDs from each data
    d1.fr.common.IDcols = d1.fr.common[,d1.fr.common.IDs]
    d2.fr.common.IDcols = d2.fr.common[,d2.fr.common.IDs]

    
    ##rename d2 rows to same as d1
    ##order should be aligned by subset
    ##row.names(d2.fr.common) = row.names(d1.fr.common)

    
    ##merge datasets
    ##d1's IDs will prevail
    d.fr.common = d1.fr.common[,!d1.fr.common.IDs] + d2.fr.common[,!d2.fr.common.IDs]

    
    ##add IDs to end of merged data
    d.fr.common = cbind(d.fr.common,d1.fr.common.IDcols)
    d.fr.common = cbind(d.fr.common,d2.fr.common.IDcols)

    ##d1's IDs prevailed in merge, so add d2's IDs to end
    d.fr.common = cbind(d.fr.common,ID.map.df[,2])
    

    d.fr.common
}


##end of combine.quantile.freqs.R
