### invitroSPI ###
# description:  removal of peptides that are present in the control data set from the database
# input:        filtered (all) PSMs
# output:       extracted PSMs of the final database
# authors:      GS, JL, HPR

library(stringr)
library(dplyr)

print("----------------------------")
print("3) IDENTIFY SYNTHESIS ERRORS")
print("----------------------------")


### INPUT ###
load(snakemake@input[["allPSMs"]])
keep_synErrors = snakemake@params[["keep_synErrors"]]

# gsub("I","L",subControls0$pepSeq) %>% unique() %>% length()
# gsub("I","L",ProteasomeDB$pepSeq[str_detect(ProteasomeDB$productType, "synError")]) %>% unique() %>% length()
# gsub("I","L",ProteasomeDB$pepSeq[str_detect(ProteasomeDB$productType, "_synErrorPrecursor")]) %>% unique() %>% length()
# gsub("I","L",ProteasomeDB$pepSeq[str_detect(ProteasomeDB$productType, "_synError$")]) %>% unique() %>% length()


### MAIN PART ###
subControls0 = allPSMs[allPSMs$digestTime == "CTRL", ]
extracted0 = allPSMs[allPSMs$digestTime != "CTRL", ]

unqTSNctrl<-unique(subControls0$substrateID)
unqTSNextr<-unique(extracted0$substrateID)
unqTSN = unqTSNextr

allExtracted = numeric()
for (i in 1:length(unqTSN))
{
  
  subIndex=which(subControls0$substrateID == unqTSN[i])
  extIndex=which(extracted0$substrateID == unqTSN[i])
  extracted<-extracted0[extIndex,]
  subControls<-subControls0[subIndex,]
  
  
  # Get all unique peptide seqs present in the substrate controls
  subControlSeq <- unique(gsub("I","L",subControls$pepSeq))
  subControlPCP = unique(gsub("I","L",subControls$pepSeq[subControls$productType == "PCP"]))
  
  # get SR2s of detected sequences
  extPos = str_split_fixed(extracted$positions,"[:punct:]",Inf)[,c(1:4)]
  extPos = apply(extPos,2,as.numeric)
  extracted = extracted %>%
    mutate(sr2s = gsub("I","L",str_sub(extracted$substrateSeq, extPos[,3], extPos[,4])),
           synErrSR2 = ifelse(sr2s %in% subControlPCP, "yes", "no"),
           synErrSR2 = ifelse(productType == "PCP", NA, synErrSR2)) %>%
    select(-sr2s)
  
  
  # Check which entries of extracted are also present in the substrate controls.
  errorIndex <- which(gsub("I","L",extracted$pepSeq) %in% subControlSeq)
  errorSeqs <- unique(gsub("I","L",extracted$pepSeq[errorIndex]))
  
  # We conservatively regard all matches as synthesis errors and delete them from extracted
  # extracted0 <- extracted0[-extIndex[errorIndex], ]
  
  # Match PSPs to errors
  PSPtoRemove <- c()
  for (b in 1:nrow(extracted)){
    candidateType <- extracted$productType[b]
    if (str_detect(candidateType,"PSP")){
      candidateSeq <- gsub("I","L",extracted$pepSeq[b])
      
      if (any(grepl(candidateSeq, subControlSeq))){
        PSPtoRemove <- c(PSPtoRemove, b)
      }
    }
  }
  
  mergeIndex=(unique(c(errorIndex, PSPtoRemove)))
  k = which(PSPtoRemove%in%errorIndex)
  if(length(k)>0){
    PSPtoRemove = PSPtoRemove[-k]
  }
  # Label PSPs that were found to match errors
  if (length(errorIndex)>0) {
    extracted$productType[errorIndex] = paste(extracted$productType[errorIndex],"_synError",sep="")
  }
  if (length(PSPtoRemove)>0) {
    extracted$productType[PSPtoRemove] = paste(extracted$productType[PSPtoRemove],"_synErrorPrecursor",sep="")
  }
  
  allExtracted = rbind(allExtracted,extracted)
  
} # the for loop

extracted = allExtracted
subControls = subControls0

# Also, check if timepoint 0 measurements or substrate controls are missing for any substrates
# controlSubstrates <- c(unique(subControls$substrateID), zeroSubstrates)
controlSubstrates <- unique(subControls$substrateID)
extractedSubstrates <- unique(extracted$substrateID)

missingSubIndex <- which(!(extractedSubstrates %in% controlSubstrates))

# remove synthesis errors
if (keep_synErrors == "no") {
     
    if(length(missingSubIndex) > 0){
      print("WARNING! The following substrates do not have any substrate control measurements:")
      
      print(extractedSubstrates[missingSubIndex])
      
      print("They are being removed from your table!")
      extracted <- extracted[-which(extracted$substrateID == extractedSubstrates[missingSubIndex[i]]),]
      
    }else{
      print("Success! All substrates in extracted are covered by substrate control measurements!")
    }

 
  print("Removing synthesis errors from the final dataset")
  
  kk = which(str_detect(extracted$productType, "_synError"))
  if (length(kk > 0)) {
    extracted = extracted[-kk, ]
    
    print("........................................................")
    paste0("number of synthesis errors that are being removed: ", length(kk)) %>%
      print()
    print("........................................................")
    
  }
  
}

### OUTPUT ###

save(extracted, file = unlist(snakemake@output[["PSMs"]]))

