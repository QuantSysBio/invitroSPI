### invitroSPI ###
# description:  extract MSDB entries that are clearly identifiable as PCP/PSP
# input:        MSDB
# output:       filtered PSMs
# authors:      GS, JL, modified by HPR

library(stringr)
library(plyr)

print("------------------------------------------------")
print("2) EXTRACT HIGH-QUALITY PSMs FROM SEARCH RESULTS")
print("------------------------------------------------")

delta_score = snakemake@params[["delta_score"]]
ion_score = snakemake@params[["ion_score"]]
q_value = snakemake@params[["q_value"]]

print(paste0("Delta score: ", delta_score))
print(paste0("Mascot ion score: ", ion_score))
print(paste0("q-value: ", q_value))


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       sep = ";", header = T, stringsAsFactors = F)
load(snakemake@input[["MSDB"]])


### MAIN PART ###
columns = names(MSDB)
runIDs = unique(MSDB$runID)

allPSMs = list()

pb = txtProgressBar(min = 0, max = length(runIDs), style = 3)
for (i in 1:length(runIDs)) {
  
  setTxtProgressBar(pb, i)
  
  selected <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(selected) = columns
  
  runIDTable <- MSDB[MSDB$runID == runIDs[i],]
  scanNum <- as.character(runIDTable$scanNum)
  scanNum <- unique(scanNum)
  
  substrateX <- gsub("(I|L)", "X", runIDTable$substrateSeq[1])
  
  # Iterate through scanNum and sort by rank
  for (k in 1:length(scanNum)) {
    PSPcandidatesIndex = NULL
    PCPcandidatesIndex = NULL
    
    filteredScans = data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(filteredScans) <- columns
    
    #replace leucin(L) and isoleucin(I) with an X 
    scanNumIndex <- which(as.character(runIDTable$scanNum) %in% scanNum[k])
    scanNumTable <- runIDTable[scanNumIndex,]
    scanNumTable <- scanNumTable[order(scanNumTable$rank),]
    origPepSeqs <- data.frame("pepSeq" = scanNumTable$pepSeq)
    
    #scanNumTable$pepSeq <- gsub("(I|L)", "X", scanNumTable$pepSeq) 
    replaceSeq <- gsub("(I|L)", "X", scanNumTable$pepSeq)
    replacesubtrateX = rep(NA,length(replaceSeq))
    newPos <- matrix(NA,length(replaceSeq),2)
    newSeq = rep(NA,length(replaceSeq))
    
    for (sq in 1:length(replaceSeq))
    {
      replacesubtrateX[sq] <- grepl(replaceSeq[sq],substrateX)
      if (replacesubtrateX[sq] == 1){
        newPos[sq,] <- str_locate(substrateX, replaceSeq[sq])[1,]
        newSeq[sq] <-  paste(strsplit(scanNumTable$substrateSeq[1],split="")[[1]][newPos[sq,1]:newPos[sq,2]],sep="",collapse="")
      }
    }
    
    scanNumTable$productType[which(replacesubtrateX==1)] = 'PCP'
    scanNumTable$spliceType[which(replacesubtrateX==1)] = NA
    temp<-newPos[which(replacesubtrateX==1),]
    
    if (length(temp)==2){
      temp = matrix(temp,1,2)
    }
    
    scanNumTable$positions[which(replacesubtrateX==1)] = apply(temp,1,paste,sep="_",collapse="_")
    
    scanNumTable$pepSeq[which(replacesubtrateX==1)] = newSeq[which(replacesubtrateX==1)]
    
    # How many entries in filteredScans are rank 1?
    filteredScans = scanNumTable
    filteredScans <- filteredScans[order(filteredScans$rank),]
    #print(filteredScans)
    #browser()
    
    if (filteredScans[1, "ionScore"] != 0){
      topScore <- filteredScans[1, "ionScore"]
    }else{
      print("filteredScan ion score is 0")
      filteredScans <- filteredScans[0,]
      next
    }
    topRanked <- length(which(as.numeric(filteredScans$rank) == 1))
    if (topRanked == 0)
    {
      filteredScans <- filteredScans[0,]
      next
    }
    
    # Calculate delta scores
    
    # If more than 1 entry is rank 1, remove entire table because correct seq can not be determined
    if (topRanked > 1) {
      
      topIndices = which(scanNumTable$rank == 1)
      topSeqReplaced = replaceSeq[topIndices]
      lenTopRank = length(unique(topSeqReplaced))
      if (lenTopRank > 1 & length(scanNumTable$productType[which(scanNumTable$rank == 1)]=="PCP")>1){ #
        filteredScans <- filteredScans[0, ]
      }
      if (lenTopRank > 1 & length(scanNumTable$productType[which(scanNumTable$rank == 1)]=="PSP")>1){ #
        filteredScans <- filteredScans[0, ]
      }
      if (lenTopRank > 1 & length(scanNumTable$productType[which(scanNumTable$rank == 1)]=="PCP")==1){ #
        PCPindex = which(scanNumTable$rank == 1)
        filteredScans <- filteredScans[PCPindex, ]
      }
      
    } 
    
    if (topRanked == 1) {
      lenTopRank = 1
    }
    
    if (nrow(filteredScans) > 1 & lenTopRank == 1 & filteredScans[1, "productType"] == "PCP") {
      #lowerRanked <- 2:nrow(filteredScans)
      #filteredScans <- filteredScans[-lowerRanked,]
      filteredScans <- filteredScans[1,]
      # print("Top ranked is PCP. Deleting other entries")
      
      # If rank 1 entry is PSP check if a probable PCP is present in the lower ranks (PCPcandidates)
    }else if (nrow(filteredScans) > 1 & lenTopRank == 1 & filteredScans[1, "productType"] == "PSP") {
      
      PCPcandidatesIndex_temp <- which(as.numeric(filteredScans$rank) > 1 & as.character(filteredScans$productType) == "PCP")
      PSPcandidatesIndex_temp <- which(as.numeric(filteredScans$rank) > 1 & as.character(filteredScans$productType) == "PSP")
      
      if((length(PCPcandidatesIndex_temp)==0)&(length(PSPcandidatesIndex_temp)==0)){
        filteredScans <- filteredScans[0,]
        next
      }
      
      if (length(PSPcandidatesIndex_temp)>0)
      {
        PSPcandidatesIndex = which(as.numeric(filteredScans$rank)==min(as.numeric(filteredScans$rank)[PSPcandidatesIndex_temp]) & as.character(filteredScans$productType) == "PSP")[1]
      }
      else
      {
        PSPcandidatesIndex = NULL
      }
      if (length(PCPcandidatesIndex_temp)>0)
      {
        PCPcandidatesIndex = which(as.numeric(filteredScans$rank)==min(as.numeric(filteredScans$rank)[PCPcandidatesIndex_temp]) & as.character(filteredScans$productType) == "PCP")[1]
      }
      else
      {
        PCPcandidatesIndex = NULL
      }
      if (length(PCPcandidatesIndex) > 0) {
        PCPcandidatescore <- filteredScans[PCPcandidatesIndex, "ionScore"]
        # Keep top ranked PSP if deltascore >0.3 and remove rest
        if (1 - PCPcandidatescore / topScore > delta_score)
        {
          # Calculate delta scores
          #topScore <- filteredScans[1, "ionScore"]
          if (length(PSPcandidatesIndex)>0)
          {
            PSPcandidatescore <- filteredScans[PSPcandidatesIndex, "ionScore"]
            if (1 - PSPcandidatescore / topScore > delta_score)
            {
              filteredScans <- filteredScans [1,]
              # Take candidate PCP as new top ranked and remove rest
            }
            else if (1 - PSPcandidatescore / topScore <= delta_score)
            {
              filteredScans <- filteredScans[0,]
            }
            
          }
          if (length(PSPcandidatesIndex)==0){
            filteredScans <- filteredScans[1,]
          }
        }
        # Take candidate PCP as new top ranked and remove rest
        else if (1 - PCPcandidatescore / topScore <= delta_score)
        {
          filteredScans <- filteredScans[PCPcandidatesIndex,]
        }
        
      } # end if (length(PCPcandidatesIndex) > 0) 
      # if only PSPs
      if (length(PSPcandidatesIndex)>0 & length(PCPcandidatesIndex)==0){
        PSPcandidatescore <- filteredScans[PSPcandidatesIndex, "ionScore"]
        if (1 - PSPcandidatescore / topScore > delta_score)
        {
          filteredScans <- filteredScans [1,]
          
          # Take candidate PCP as new top ranked and remove rest
        }
        else if (1 - PSPcandidatescore / topScore <= delta_score)
        {
          filteredScans <- filteredScans[0,]
        }
        
      }
    }
    
    #filter out scans with low qvalue or scores
    if(nrow(filteredScans)> 0){
      if(as.double(filteredScans$qValue) >= q_value | as.double(filteredScans$ionScore <= ion_score)){
        
        filteredScans <- filteredScans[0,]
      }else{
        selected <- rbind(filteredScans, selected)
      }
    }
    
  }
  
  allPSMs[[i]] = selected
  
}

allPSMs = plyr::ldply(allPSMs)

### OUTPUT ###

save(allPSMs, file = unlist(snakemake@output[["allPSMs"]]))

