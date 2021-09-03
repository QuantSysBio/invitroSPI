### invitroSPI ###
# description:  mapping of peptides to the substrate sequence accounting for potential multi-mappers
# input:        extracted PSMs of the final database
# output:       ProteasomeDB.csv
# authors:      JL, modified by HPR

library(seqinr)

print("-----------------------------")
print("4) MAPPING AND POSTPROCESSING")
print("-----------------------------")


### INPUT ###
load(snakemake@input[["PSMs"]])
source("SOURCE/loadFunctions.r")


### MAIN PART ##
d = extracted
pb = txtProgressBar(min = 0, max = dim(d)[1], style = 3)

for(i in 1:dim(d)[1]){
  setTxtProgressBar(pb, i)
  
  if(!(d$productType[i]=="CONT" | d$productType[i]=="CONT_synError")){
    s = gsub("I","L",as.vector(d$pepSeq[i]))
    substrate = gsub("I","L",as.vector(d$substrateSeq[i]))
    x = getPositions(s,substrate)
    
    
    #PCP
    if(dim(x)[2]==2){
      d$positions[i] = paste(apply(x,1,paste,collapse="_"),collapse=";")
    }
    
    #PSP
    if(dim(x)[2]>2){
      # print(i)
      if(dim(x)[2]==5 & dim(x)[1]>1){
        d$positions[i] = paste(apply(x[,-1],1,paste,collapse="_"),collapse=";")
      }
      if(dim(x)[2]==5 & dim(x)[1]==1){
        d$positions[i] = paste(x[,-1],collapse="_")
      }
      
      types = rep("cis",dim(x)[1])
      
      intv = as.numeric(x[,4])-as.numeric(x[,3])
      k = which(intv<=0)
      if(length(k)>0){
        types[k] = "revCis"
        
        k2 = which(as.numeric(x[k,2])<=as.numeric(x[k,5]) & as.numeric(x[k,3])>=as.numeric(x[k,4]))
        if(length(k2)>0){
          types[k[k2]] = "trans"
        }
        
      }
      
      
      d$spliceType[i] = paste(types,collapse=";")
      
    }
  }
  
}


### OUTPUT ###
write.csv(d, file = unlist(snakemake@output[["ProteasomeDB"]]),
          row.names = F)


