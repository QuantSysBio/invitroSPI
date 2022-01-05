### invitroSPI ###
# description:  extract MSDB entries that are clearly identifiable as PCP/PSP
# input:        MSDB
# output:       filtered PSMs
# authors:      JL, GS, HPR

library(stringr)
library(plyr)
library(dplyr)

print("------------------------------------------------")
print("2) EXTRACT HIGH-QUALITY PSMs FROM SEARCH RESULTS")
print("------------------------------------------------")

project_name = snakemake@params[["project_name"]]
delta_score = snakemake@params[["delta_score"]]
ion_score = snakemake@params[["ion_score"]]
q_value = snakemake@params[["q_value"]]

print(paste0("Delta score: ", delta_score))
print(paste0("Mascot ion score: ", ion_score))
print(paste0("q-value: ", q_value))


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       sep = ",", header = T, stringsAsFactors = F)
# filter sample list
sample_list = sample_list[sample_list$project_name == project_name, ]

load(snakemake@input[["MSDB"]])


### MAIN PART ###
columns = names(MSDB)
runIDs = unique(MSDB$runID)


allPSMs = list()
allPSMs_Delta = list()

pb = txtProgressBar(min = 0, max = length(runIDs), style = 3)
for (i in 1:length(runIDs)) {

  setTxtProgressBar(pb, i)
  
  selected <- data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(selected) = columns
  
  selectedDeltaRecord = data.frame(matrix(ncol = length(columns), nrow = 0))
  colnames(selectedDeltaRecord) <- columns
  
  runIDTable <- MSDB[MSDB$runID == runIDs[i],]
  scanNum <- as.character(runIDTable$scanNum)
  scanNum <- unique(scanNum)
  
  
  # Iterate through scanNum and sort by rank
  for (k in 1:length(scanNum)) {
    
    PSPcandidatesIndex = NULL
    PCPcandidatesIndex = NULL
    
    deltaRecord = data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(deltaRecord ) <- columns
    
    scanNumIndex <- which(as.character(runIDTable$scanNum) %in% scanNum[k])
    scanNumTable <- runIDTable[scanNumIndex,]
    scanNumTable <- scanNumTable[order(scanNumTable$rank),]
    
    # # How many entries in filteredScans are rank 1?
    # scanNumTable = unique(scanNumTable)  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    filteredScans = scanNumTable
    filteredScans = filteredScans[order(filteredScans$rank),]
    

    # extract the top score
    if (filteredScans[1, "ionScore"] != 0){
      topScore <- filteredScans[1, "ionScore"]
    } else{
      filteredScans <- filteredScans[0,]
    }
    
    # extract scans with rank 1
    allRanks = filteredScans %>%
      distinct() %>%
      dplyr::count(rank)
    
    NoTopRankedScans <- allRanks$n[allRanks$rank == 1]
    if (NoTopRankedScans == 0) {
      filteredScans <- filteredScans[0,]
    }
    
    # if more than 1 scan is rank 1, discard scam unless there is only a single PCP
    if (NoTopRankedScans > 1) {
      
      topIndices = which(filteredScans$rank == 1)
      NoTopRankedPeps = filteredScans$pepSeq[topIndices] %>%
        unique() %>%
        length()
      
      filteredScans = unique(filteredScans)
      PCPindex = which(filteredScans$productType[filteredScans$rank == 1]=="PCP")
      PSPindex = which(filteredScans$productType[filteredScans$rank == 1]=="PSP")
      
      # if more than one peptide is rank 1
      # keep scan only if there is just a single PCP, otherwise discard
      if (NoTopRankedPeps > 1) {
        
        if (length(PCPindex) == 1) {
          
          filteredScans = filteredScans[PCPindex,]
          
          } else {
            
            if (length(PSPindex) > 1) {
              index = which(filteredScans$productType[which(filteredScans$rank == 1)]=="PSP")
              PSPcandidatescores <- filteredScans[index, "ionScore"]
              ind = which(1-PSPcandidatescores/topScore <= delta_score)
              deltaRecord = rbind(deltaRecord,filteredScans[index[ind],])
            }
            
            filteredScans = filteredScans[0,]
            
        }
        
      }
      
    }
    # if there is only one scan on rank 1, there is also only one peptide
    
    # if there is only a single scan on rank 1
    # (nrow(filteredScans) will be larger than 1 since the previous if condition did not apply)
    # if the 1st rank is a PCP, assign it as such
    if (NoTopRankedScans == 1 & filteredScans[1, "productType"] == "PCP") {
      
      filteredScans = unique(filteredScans)
      filteredScans <- filteredScans[1,]
      
      # If rank 1 entry is PSP check if a likely PCP is present in the lower ranks (PCPcandidates)
      # also check if there are lower-ranked PSPs
    } else if (NoTopRankedScans == 1 & filteredScans[1, "productType"] == "PSP") {
      
      filteredScans = unique(filteredScans)
      
      # lower-ranked PCPs and PSPs
      PCPcandidatesIndex_temp <- which(as.numeric(filteredScans$rank) > 1 & as.character(filteredScans$productType) == "PCP")
      PSPcandidatesIndex_temp <- which(as.numeric(filteredScans$rank) > 1 & as.character(filteredScans$productType) == "PSP")
      
      # if there are no PSPs or PCPs with a lower rank, assign the PSP
      if(length(PCPcandidatesIndex_temp)==0 & length(PSPcandidatesIndex_temp)==0){
        filteredScans <- filteredScans[1,]
      }
      
      # if there are PSPs with a lower rank --> get highest-ranked PSM
      if (length(PSPcandidatesIndex_temp)>0) {
        #keeps only top PSP candidate
        PSPcandidatesIndex = which(as.numeric(filteredScans$rank)==min(as.numeric(filteredScans$rank)[PSPcandidatesIndex_temp]) & as.character(filteredScans$productType) == "PSP")[1]
        # keeps all PSP candidates
        PSPcandidatesIndexAll = PSPcandidatesIndex_temp
      } else {
        PSPcandidatesIndex = NULL
      }
      
      # if there are PCPs with a lower rank --> get highest-ranked PSM
      if (length(PCPcandidatesIndex_temp)>0) {
        PCPcandidatesIndex = which(as.numeric(filteredScans$rank)==min(as.numeric(filteredScans$rank)[PCPcandidatesIndex_temp]) & as.character(filteredScans$productType) == "PCP")[1]
      } else {
        PCPcandidatesIndex = NULL
      }
      
      # for lower-ranked PCPs --> calculate Delta score
      # Keep top ranked PSP if deltascore > 0.3 and remove rest
      if (length(PCPcandidatesIndex) > 0) {
        
        
        PCPcandidatescore <- filteredScans[PCPcandidatesIndex, "ionScore"]
        # if the Delta score between the top-ranked PSP and a lower-ranked PCP
        # is larger than the threshold
        if (1 - PCPcandidatescore / topScore > delta_score) {
          
          # if there are any lower-ranked PSPs
          # Calculate delta_score scores also for lower-ranked PSPs
          if (length(PSPcandidatesIndex)>0) {
            PSPcandidatescore <- filteredScans[PSPcandidatesIndex, "ionScore"]
            
            # lower-ranked PSP, but with very low ion score --> assign top-ranked PSP
            if (1 - PSPcandidatescore / topScore > delta_score) {
              filteredScans <- filteredScans [1,]
              
              # discard spectrum bc there are several PSP candidates with similar
              # ion score
            } else if (1 - PSPcandidatescore / topScore <= delta_score) {
              
              PSPcandidatescores <- filteredScans[PSPcandidatesIndexAll, "ionScore"]
              ind = which(1-PSPcandidatescores/topScore <= delta_score)
              deltaRecord = rbind(deltaRecord,filteredScans[PSPcandidatesIndexAll[ind],])
              
              filteredScans <- filteredScans[0,]
            }
            
            # if there are no lower-ranked PSPs --> assign the top-ranked PSP
          } else if (length(PSPcandidatesIndex)==0){
            filteredScans <- filteredScans[1,]
          }
          
          # if the Delta score is too small:
          # Take candidate PCP as new top ranked and remove rest
        } else if (1 - PCPcandidatescore / topScore <= delta_score) {
          filteredScans <- filteredScans[PCPcandidatesIndex,]
        }
        
        # end if there are both PCPs and PSPs in lower ranks
        # if there are lower ranked PSPs but no PCPs: continue
        
      } else if (length(PSPcandidatesIndex) > 0 & length(PCPcandidatesIndex)==0){
        PSPcandidatescore <- filteredScans[PSPcandidatesIndex, "ionScore"]
        
        # keep the candidate PSP if the Delta score is large enough
        if (1 - PSPcandidatescore / topScore > delta_score) {
          filteredScans <- filteredScans [1,]
          
          # if the Delta score is too small --> discard entire scan
        } else if (1 - PSPcandidatescore / topScore <= delta_score) {
          
          PSPcandidatescores <- filteredScans[PSPcandidatesIndexAll, "ionScore"]
          ind = which(1-PSPcandidatescores/topScore <= delta_score)
          deltaRecord = rbind(deltaRecord,filteredScans[PSPcandidatesIndexAll[ind],])
          
          filteredScans <- filteredScans[0,]
        }
        
      }
    }
    
    # apply filter for q-value and ion score
    if(nrow(filteredScans) > 0){
      if(as.double(filteredScans$qValue) >= q_value | as.double(filteredScans$ionScore <= ion_score)){
        filteredScans <- filteredScans[0,]
      }else{
        selected <- rbind(filteredScans, selected)
      }
    }
    
    if(nrow(deltaRecord)> 0){
      k = which(as.double(deltaRecord$qValue) < q_value & as.double(deltaRecord$ionScore > ion_score))
      deltaRecord <- deltaRecord[k,]
      
      if(nrow(deltaRecord)> 0){
        selectedDeltaRecord <- rbind(deltaRecord,selectedDeltaRecord)
      }
    }
    
  }
  
  allPSMs[[i]] = selected
  allPSMs_Delta[[i]] = selectedDeltaRecord

}

allPSMs = plyr::ldply(allPSMs)
allPSMs_Delta = plyr::ldply(allPSMs_Delta)

paste0("number of PSMs: ", nrow(allPSMs)) %>%
  print()


### OUTPUT ###

save(allPSMs, file = unlist(snakemake@output[["allPSMs"]]))
save(allPSMs_Delta, file = unlist(snakemake@output[["allPSMs_Delta"]]))


