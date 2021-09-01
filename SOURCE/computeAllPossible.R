### invitroSPI ###
# description:  calculate all possible spliced and non-spliced peptides from substrate
# input:        sample list, config file
# output:       substrateID_allPSP.fasta
# authors:      JL, modified by HPR

library(dplyr)
library(seqinr)
source("SOURCE/SplicingFunctions_Reduced_trans.r")

print("-----------------------------------------------------")
print("COMPUTE ALL POSSIBLE SPLICED AND NON SPLICED PEPTIDES")
print("-----------------------------------------------------")


### INPUT ###
# load sample list
sample_list = read.csv(file = "INPUT/sample_list.csv",
                       sep = ";", header = T, stringsAsFactors = F)

# load config table
config = read.table("INPUT/config.yaml")
project_name = config[,2][which(config[,1] == "project_name:")] %>%
  as.character()

# filter sample list
sample_list = sample_list[sample_list$project_name == project_name, ]

substrates = sample_list[, c("substrateID", "substrateSeq")] %>%
  unique()

### MAIN PART ###
# create project subdirectory
if(!dir.exists(paste0("OUTPUT/", project_name, "/"))) {
  dir.create(paste0("OUTPUT/", project_name, "/"))
}

# calculate all posssible PCPs and PSPs
N = nrow(substrates)

for(pept in 1:N){
  
  if (grepl(pattern = ".fasta", substrates$substrateSeq[pept])) {
    
    cntSeq = read.fasta(file = paste0("INPUT/sequences/", substrates$substrateSeq[pept]),
                        seqtype = "AA", strip.desc = T) %>%
      unlist()
    
  } else {
    cntSeq = substrates$substrateSeq[pept] %>% strsplit("") %>% unlist()
  }
  
  
  Name = substrates$substrateID[pept]
  Names = Name
  write.fasta(paste(cntSeq,sep="",collapse=""),
              Names,
              file.out=paste0("OUTPUT/", project_name, "/", Name, "_allPSP.fasta"),
              open="w")
  
  print("################")
  print(substrates$substrateID[pept])
  print("###############")
  
  for(nmer in 4:(3*length(cntSeq))){
    #for(nmer in 6:50){
    
    print("-------------------------------")
    print(paste("Computing ",nmer,"mers",sep=""))
    print("-------------------------------")
    
    inputSequence1 = paste(cntSeq,sep="",collapse="")
    SP1s = CutAndPaste(inputSequence1,nmer=nmer,MiSl=500)
    
    
    if(length(SP1s[[4]])>0){
      Name = substrates$substrateID[pept]
      Names = getNames(SP1s[[3]])
      
      write.fasta(lapply(SP1s[[4]],s2c),
                  Names,
                  as.string=FALSE,
                  file.out=paste0("OUTPUT/", project_name, "/", Name, "_allPSP.fasta"),
                  open="a")
    }
    
    if(length(SP1s[[2]])>0){
      Name = substrates$substrateID[pept]
      Names = getNames(SP1s[[1]])
      
      write.fasta(lapply(SP1s[[2]],s2c),
                  Names,
                  as.string=FALSE,
                  file.out=paste0("OUTPUT/", project_name, "/", Name, "_allPSP.fasta"),
                  open="a")
    }
    
    
  }
  
  
}

print("-----------------------------------------------------")
print("COMPUTATION FINISHED SUCCESSFULLY")
print("-----------------------------------------------------")


