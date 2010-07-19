##store list of GSM names, one per row, in a text file

collateGSMs <- function(gsms,sample.names=gsms) {

    if( (length(sample.names)>0) & (length(gsms) != length(sample.names)) ) {
        stop("gsms and sample.names arguments must be the same length.")
    }

    ##check that valid vector present
    if( (is.null(gsms)) | (length(gsms)==0) ) {
        stop("Please enter a vector of GSMs for the first argument.")
    }
    
    ### END OF INPUT CHECKS ###
    
    missed = character()
    dat.nrow = 0

    ##installed via Depends
    require(GEOquery)


    for (gsm.name in gsms) {

        if(length(grep("^gsm",gsm.name,ignore.case=TRUE))>0) {
            ##incase they have a column header
        }
        else {
            next;
        }

        ##remove .CEL/.cel at end
        if(length(grep(".CEL$",gsm.name,ignore.case=TRUE))>0) {
            gsm.name = strsplit(gsm.name,".")[[1]][1]
        }

        ##remove spaces at end
        if(length(grep(" ",gsm.name))>0) {
            gsm.name = gsub(" ","",gsm.name)
        }
        
        ##download file; will catch if a connection error (among other things?)
	gsm = try(getGEO(gsm.name),silent=TRUE)
	
        ##download may fail
	if( !(class(gsm) == "try-error") ) { 	

            ##check for Table1
            tab = Table(gsm)

            if(exists("tab")) {
                
                ##check for VALUE
                n = names(tab)

                if(length(grep("^VALUE",n))>0) {
                    
                    calls = tab[,c("ID_REF","VALUE")]

                    if(exists("calls")) {

                        calls.sorted = calls[order(calls$ID_REF),]
                                            
                        if(!exists("dat")) {
                            dat = data.frame(calls.sorted,row.names=1)
                            dat.nrow = nrow(dat)
                        }
                        else {
                            
                            if(nrow(calls.sorted) != dat.nrow) {
                                stop("Number of genes in samples differ, so cannot collate data.\n\n")
                            }
                            
                            df = data.frame(calls.sorted,row.names=1)
                            dat = cbind(dat,df)
                        }
                    }
                    else {
                        warning(paste("Error in extracting data for ",gsm.name,sep=""))
                    }
                }
                else {
                    warning(paste("MAS 5.0 expressions missing from ",gsm.name,sep=""))
                }
            }
            else {
                warning(paste("Data table missing for ",gsm.name,sep=""))
            }
	}
	else {
            missed = c(missed,gsm.name)	
	}
    }
    
    if( length(missed) != 0 ) {
        
        m.list = ""
        m.pos = vector()
        #pres = vector()


        
        for(i in 1:length(missed)) {
            m.list = paste(m.list," ",missed[i],sep="")
            m.pos = c(m.pos, grep(paste("^",missed[i],"$",sep=""),gsms))
        }

        pres = sample.names[-unique(m.pos)]
        
        warning(paste("The following files were not downloaded due to an error during the download: ",m.list,sep=""))
    } else {
        pres = sample.names
    }

    if(exists("dat")) {
        names(dat) = pres
        
        return(dat)
    } else {
        stop(paste("Unable to find data for GSMs: ",gsms,sep=""))
    }
}
