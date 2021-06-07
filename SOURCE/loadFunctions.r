library(plyr)
getPositions <- function(seq,substrate){
    
    

    #########################
    # PCP
    #########################

        
    l = nchar(seq)
    
    k = which((grepl(pattern=seq,x=substrate)==TRUE))
    if(length(k)>0){
    
        pcp = numeric()
        
        for(j in 1:length(k)){
            a = substrate
            x = strsplit(a,split=seq)[[1]]
            nn = nchar(x)
            n1 = rep(NA,(length(nn)-1))
            n2 = rep(NA,(length(nn)-1))
            for(r in 1:(length(x)-1)){
                n1[r] = sum(nn[1:r])+(r-1)*nchar(seq)+1
                n2[r] = n1[r]+nchar(seq)-1
            }
            pcp = rbind(pcp,cbind(n1,n2))
        }
        return(pcp)
    }

    
    
    #########################
    # PSP
    #########################

     
    
    
    ll = nchar(seq)


    pept = unlist(seq)
    N = nchar(seq)

    # split peptides to P matrix
    P = strsplit(pept,split="")[[1]]

    # get permutations of length N
    x = c(1:N)
    y = c(1:N)
    z = as.vector(outer(x,y,paste,sep="_"))
    q = matrix(NA,length(z),2)
    for(i in 1:length(z)){
        q[i,] = as.numeric(strsplit(z[i],split="_")[[1]])
    }

    qs = apply(q,1,sum)
    k = which(qs==N)
    q = q[k,]

    # loop over all peptides
    res2 <- list()
    res1 <- list()

    psp <- list()

    psp <- list()
    res1 <- list()
    res2 <- list()
    
    # generate all strings for searches
    S = matrix(NA,dim(q)[1],2)
    for(i in 1:dim(q)[1]){
        S[i,1] = paste(P[1:q[i,1]],sep="",collapse="")
        S[i,2] = paste(P[(q[i,1]+1):N],sep="",collapse="")
    }

    # search each entry in prot for the two corresponding fragments and extract acc and positions

    for(i in 1:dim(S)[1]){
        
        psp[[i]] <- list()
        res1[[i]] = which((grepl(pattern=S[i,1],x=substrate)==TRUE))
        res2[[i]] = which((grepl(pattern=S[i,2],x=substrate)==TRUE))
        
        kk = which(res1[[i]]%in%res2[[i]])
        k = res1[[i]][kk]
        if(length(k)>0){
           
            for(j in 1:length(k)){

                a = substrate

            
                x = strsplit(a,split=S[i,1])[[1]]
                nn = nchar(x)
                n1 = rep(NA,(length(nn)-1))
                n2 = rep(NA,(length(nn)-1))
                for(r in 1:(length(x)-1)){
                    n1[r] = sum(nn[1:r])+(r-1)*nchar(S[i,1])+1
                    n2[r] = n1[r]+nchar(S[i,1])-1
                }
                #check if substrate Cterm==S[i,1]
                len = nchar(S[i,1])
                y = paste(strsplit(a,split="")[[1]][(nchar(a)-len+1):nchar(a)],collapse="")
                if(S[i,1]==y){
                    n1 = c(n1,nchar(a)-len+1)
                    n2 = c(n2,nchar(a))
                }
                tmp = unique(apply(cbind(n1,n2),1,paste,collapse="_"))
                tmp2 = matrix(as.numeric(unlist(strsplit(tmp,split="_"))),length(tmp),2,byrow=TRUE)
                n1 = tmp2[,1]
                n2 = tmp2[,2]

                x = strsplit(a,split=S[i,2])[[1]]
                nn = nchar(x)
                n3 = rep(NA,(length(nn)-1))
                n4 = rep(NA,(length(nn)-1))
                for(r in 1:(length(x)-1)){
                    n3[r] = sum(nn[1:r])+(r-1)*nchar(S[i,2])+1
                    n4[r] = n3[r]+nchar(S[i,2])-1
                }
                #check if substrate Cterm==S[i,2]
                len = nchar(S[i,2])
                y = paste(strsplit(a,split="")[[1]][(nchar(a)-len+1):nchar(a)],collapse="")
                if(S[i,2]==y){
                    n3 = c(n3,nchar(a)-len+1)
                    n4 = c(n4,nchar(a))
                }
                tmp = unique(apply(cbind(n3,n4),1,paste,collapse="_"))
                tmp2 = matrix(as.numeric(unlist(strsplit(tmp,split="_"))),length(tmp),2,byrow=TRUE)
                n3 = tmp2[,1]
                n4 = tmp2[,2]

                # get all internal combinations and keep only those with intv<=25
                
                z = as.vector(outer(n2,n3,paste,sep="_"))
                y = matrix(NA,length(z),2)
                for(zz in 1:length(z)){
                    y[zz,] = as.numeric(strsplit(z[zz],split="_")[[1]])
                }
                intv = y[,2]-y[,1]-1
                x = which(intv<0)
                if(length(x)>0){ intv[x] = y[x,1]-y[x,2]+1-nchar(S[i,1])-nchar(S[i,2]) }
                x = which(intv<0)
            #    if(length(x)>0){ intv[x] = 1000 }

                select = which(intv<=5000)
                
                nnn = length(select)
                if(nnn>0){
                    psp[[i]][[j]] = matrix(NA,nnn,5)

                    for(j2 in 1:nnn){
                        
                        psp[[i]][[j]][j2,] = c(pept,y[select[j2],1]-nchar(S[i,1])+1,y[select[j2],1],y[select[j2],2],y[select[j2],2]+nchar(S[i,2])-1)
                    }
                }
                
            }
            
        }


    }
    
    # unlist results and return as unique matrix with all possible explanations as rows
    x = unlist(psp)
  #  psp = matrix(x,length(x)/5,5,byrow=FALSE)
    
    res = numeric()
    for(i in 1:length(psp)){
        if(length(psp[[i]])>0){
            for(j in 1:length(psp[[i]])){
                res = rbind(res,psp[[i]][[j]])
            }
        }
    }
    
    
   # print(res)
    
    return(res)

    
}
