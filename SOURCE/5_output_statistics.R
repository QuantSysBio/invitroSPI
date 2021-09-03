### invitroSPI ###
# description:  generate statistics of ProteasomeDB
# input:        ProteasomeDB
# output:       DB statistics
# authors:      HPR

source("SOURCE/invitroSPI_utils.R")

print("-------------------------------")
print("5) GENERATE DATABASE STATISTICS")
print("-------------------------------")


### INPUT ###
ProteasomeDB = read.csv(snakemake@input[["ProteasomeDB"]],
                        stringsAsFactors = F)
# ProteasomeDB = read.csv("../data/submission/ProteasomeDB.csv", stringsAsFactors = F)

### MAIN PART ###
##### stats #####
print("#################################")
print("assigned peptide-spectrum matches")
paste0("n = ", nrow(ProteasomeDB)) %>%
  print()

ProteasomeDB.psm = ProteasomeDB %>%
  ILredundancy() %>%
  filterPepLength() %>%
  removeSubstrateFromPCPs %>%
  disentangleMultimappers.Type() %>%
  DBcosmetics()
print("#################################")


print("###############")
print("unique peptides")
ProteasomeDB.unique = ProteasomeDB.psm %>%
  remSynthErrors() %>%
  uniquePeptides() %>%
  DBcosmetics()
print("###############")

###### plots #####
TypeFrequency = function(DB) {
  
  subSeqs = DB$substrateSeq %>% unique()
  freqs = list()
  
  for (s in 1:length(subSeqs)) {
    cntDB = DB[DB$substrateSeq == subSeqs[s], ]
    cntFreqs = cntDB$spliceType %>% table() %>% as.data.frame()
    cntFreqs$. = as.character(cntFreqs$.)
    
    # remove multimappers if present
    k = which(str_detect(cntFreqs$., "type_multi-mapper"))
    if (length(k) > 0) { cntFreqs = cntFreqs[-k, ] }
    
    cntFreqs$.[which(str_detect(cntFreqs$., "cis_"))] = "cis"
    cntFreqs$.[which(str_detect(cntFreqs$., "revCis_"))] = "revCis"
    if (any(str_detect(cntFreqs$., "trans_"))) { cntFreqs$.[which(str_detect(cntFreqs$., "trans_"))] = "trans" }
    
    names(cntFreqs)[1] = "type"
    cntFreqs$Freq = as.numeric(cntFreqs$Freq)
    
    cntFreqs = cntFreqs %>%
      group_by(type) %>%
      summarise(Freq = sum(Freq))
    cntFreqs$Freq = as.numeric(cntFreqs$Freq) / nrow(cntDB)
    freqs[[s]] = cntFreqs
    
  }
  
  freqs = plyr::ldply(freqs)
  names(freqs) = c("type", "value")
  
  freqs$type = factor(freqs$type,
                      levels = c("PCP", "cis", "revCis", "trans"))
  freqs$value = freqs$value * 100
  
  return(freqs)
}

freqs = TypeFrequency(DB = ProteasomeDB.unique)

### OUTPUT ###
pdf(file = unlist(snakemake@output[["DB_stats"]]),
    height = 6, width = 4)

boxplot(value~type, data = freqs,
        col = c(plottingCols["PCP"], plottingCols["cis"],
                plottingCols["revCis"], plottingCols["trans"]),
        main = "unique peptides",
        xlab = "product type",
        ylab = "frequency per substrate (%)")

dev.off()

