

CutAndPaste <- function(inputSequence,nmer,MiSl){
    
    results <- list()
    peptide = strsplit(inputSequence,"")[[1]]
    L = length(peptide)
    
    if((2*L)>nmer){
        
        # compute all PCP with length <=nmer
        cp = computeCPomplete(L,nmer)
        print("compute all PCP with length")
        
        # get all PCP with length == nmer
        index = which((cp[,2]-cp[,1]+1)==nmer)
        cpNmer = cp[index,]
        if(length(cpNmer)==2){cpNmer = matrix(cpNmer,1,2)}
        print("got all Nmer PCP")
        if(length(cpNmer)>0){
            CPseq = translateCP(cpNmer,peptide)
        }
        else{
            CPseq = numeric()
        }
        print("CP translated")
   
   FLAG = TRUE
   if(FLAG==TRUE){
        if(length(cp[,1])>1){
            
            sp = computeSPcomplete(cp,maxL=nmer,minL=nmer,MiSl=MiSl)
            print(paste("all SP:",length(sp[,1])))
            
            SPseq = translateSP(sp,peptide)
            print("All SPs translated")
            
            # x = removeDoubles(CPseq,cpNmer)
            CPseqClean = CPseq
            cpClean = cpNmer
            print(paste("CP without doubles:",length(CPseqClean)))
            
            # x = removeDoubles(SPseq,sp)
            SPseqClean = SPseq
            spClean = sp
            print(paste("SP without doubles:",length(SPseqClean)))
            
            x = removeCPfromSP(SPseqClean,spClean,CPseqClean)
            SPseqClean = x[[1]]
            spClean = x[[2]]
            print(paste("SP without CP:",length(SPseqClean)))
            
        }
        else{
            spClean = numeric()
            SPseqClean = numeric()
            cpClean = numeric()
            CPseqClean = numeric()
        }
        
        results[[1]] = spClean
        results[[2]] = SPseqClean
   }
        results[[3]] = cpClean
        results[[4]] = CPseqClean
        
    }
    else{
        results[[1]] = numeric()
        results[[2]] = numeric()
        results[[3]] = numeric()
        results[[4]] = numeric()
        
        if(L==nmer){
            results[[3]] = matrix(c(1,nmer),1,2)
            results[[4]] = inputSequence
        }
    }
    
    return(results)

    
}





# compute all PCP with length <=nmer
computeCPomplete <- function(L,nmer){
    
    maxL = nmer+1
    CP <- numeric()
    
    for(i in 1:L){
        temp2 = c(i:min(L,(i+maxL-2)))
        temp1 = rep(i,length(temp2))
        temp3 = cbind(temp1,temp2)
        CP = rbind(CP,temp3)
    }
    return(CP)
}


# translate PCP
translateCP <- function(CP,peptide){
    
    CPseq = rep(NA,dim(CP)[1])
    for(i in 1:dim(CP)[1]){
        CPseq[i] = paste(peptide[CP[i,1]:CP[i,2]],sep="",collapse="")
    }
    return(CPseq)
    
}


# compute all PSP with length == nmer
computeSPcomplete <- function(cp,maxL,minL,MiSl){
    
    SP <- numeric()
    N = dim(cp)[1]
    NN = 10**8
    
    SP = matrix(NA,NN,4)
    a = 1
    for(i in 1:N){
        
        temp1 = rep(cp[i,1],N)
        temp2 = rep(cp[i,2],N)
        temp3 = cp[,1]
        temp4 = cp[,2]
        
        temp = cbind(temp1,temp2,temp3,temp4)
        L = temp[,4]-temp[,3]+temp[,2]-temp[,1]+2
        
        # this line excludes trans PSP
        # ind = which(((temp[,3]-temp[,2])==1)|(L>maxL)|(L<minL)|((temp[,3]-temp[,2])>(MiSl+1))|((temp[,1]-temp[,4])>(MiSl+1))|((temp[,3]<=temp[,2])&(temp[,4]>=temp[,1])))

        # this line includes trans PSP
        ind = which(((temp[,3]-temp[,2])==1)|(L>maxL)|(L<minL)|((temp[,3]-temp[,2])>(MiSl+1))|((temp[,1]-temp[,4])>(MiSl+1)))
        
        
        if(length(ind)>0){
            temp = temp[-ind,]
        }
        
        if(length(temp)>4){
            if((a+dim(temp)[1]-1)>NN){
                SP = rbind(SP,matrix(NA,NN,4))
            }
        }
        if(length(temp)>4){
            if(dim(temp)[1]>0){
                SP[c(a:(a+dim(temp)[1]-1)),] = temp
            }
            
            a = a+dim(temp)[1]
            
        }
    }
    # remove all empty lines from SP
    if(a<NN){
        SP = SP[-c(a:NN),]
    }
    
    
 
    
    return(SP)
}




# translate PSP
translateSP <- function(SP,peptide){
    
    SPseq = rep(NA,dim(SP)[1])
    for(i in 1:dim(SP)[1]){
        SPseq[i] = paste(peptide[c(SP[i,1]:SP[i,2],SP[i,3]:SP[i,4])],sep="",collapse="")
    }
    
    return(SPseq)
}


# remove double sequences
removeDoubles <- function(x,y){
    
    result <- list()
    
    temp = unique(x)
    x2 = rep(NA,length(temp))
    y2 = matrix(NA,length(temp),dim(y)[2])
    for(i in 1:length(temp)){
        ind = which(x==temp[i])
        x2[i] = x[ind[1]]
        y2[i,] = y[ind[1],]
    }
    
    result[[1]] = x2
    result[[2]] = y2
    return(result)
    
}


# remove PCP from PSP list
removeCPfromSP <- function(x,y,z){
    
    result <- list()
    index = which(x%in%z)
    if(length(index)>0){
        result[[1]] = x[-index]
        result[[2]] = y[-index,]
    }
    else{
        result[[1]] = x
        result[[2]] = y
    }
    return(result)
    
}








##########################################

getNames <- function(x){
	
	Names = rep(NA,dim(x)[1])
	
    if(dim(x)[2]==2){
        for(i in 1:dim(x)[1]){
            Names[i] = paste("PCP_",paste(x[i,],sep="",collapse="_"),sep="")
        }
    }
    
    if(dim(x)[2]==4){
        for(i in 1:dim(x)[1]){
            Names[i] = paste("PSP_",paste(x[i,],sep="",collapse="_"),sep="")
        }
    }
    
	return(Names)
	
}


##################################

getNamesIC <- function(x,y,IC,prec,mutation){
	
	if(sum(dim(x)[1])>0){
		Names = rep(NA,dim(x)[1])
		for(i in 1:dim(x)[1]){
			Names[i] = paste(paste(x[i,],sep="",collapse=" ")," id =",y,prec[i]," mutation: ",mutation[i],"  IC=",IC[i],sep=" ")
		}
	}
	else{
		Names = rep(NA,1)
		for(i in 1:1){
			Names[i] = paste(paste(x,sep="",collapse=" ")," id =",y,prec[i]," mutation: ",mutation[i],"  IC=",IC[i],sep=" ")
		}	

	}
	
	
	return(Names)
	
}



