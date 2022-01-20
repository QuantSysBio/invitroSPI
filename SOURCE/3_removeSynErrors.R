### invitroSPI ###
# description:  removal of peptides that are present in the control data set from the database
# input:        filtered (all) PSMs
# output:       extracted PSMs of the final database
# authors:      GS, JL, modified by HPR

library(stringr)

print("----------------------------")
print("3) IDENTIFY SYNTHESIS ERRORS")
print("----------------------------")


### INPUT ###
load(snakemake@input[["allPSMs"]])
keep_synErrors = snakemake@params[["keep_synErrors"]]


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
  subControlSeq <- gsub("I","L",unique(subControls$pepSeq))
  
  # Check which entries of extracted are also present in the substrate controls.
  errorIndex <- which(gsub("I","L",extracted$pepSeq) %in% gsub("I","L",subControlSeq))
  errorSeqs <- unique(gsub("I","L",extracted[errorIndex, "pepSeq"]))
  
  # We conservatively regard all matches as synthesis errors and delete them from extracted
  # extracted0 <- extracted0[-extIndex[errorIndex], ]
  
  # Match PSPs to errors
  PSPtoRemove <- c()
  for (b in 1:nrow(extracted)){
    candidateType <- extracted[b, "productType"]
    if (candidateType == "PSP"){
      candidateSeq <- gsub("I","L",extracted[b, "pepSeq"])
      
      if (nchar(candidateSeq) > 4 & any(grepl(candidateSeq, subControlSeq))){
        PSPtoRemove <- c(PSPtoRemove, b)
      }
    }
  }
  
  mergeIndex=(unique(c(errorIndex, PSPtoRemove)))
  k = which(PSPtoRemove%in%errorIndex)
  if(length(k)>0){
    PSPtoRemove = PSPtoRemove[-k]
  }
  # Remove PSPs that were found to match errors
  if (length(errorIndex)>0)
  {
    extracted$productType[errorIndex] = paste(extracted$productType[errorIndex],"_synError",sep="")
  }
  if (length(PSPtoRemove)>0)
  {
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

