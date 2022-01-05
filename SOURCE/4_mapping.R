### invitroSPI ###
# description:  mapping of peptides to the substrate sequence accounting for potential multi-mappers
# input:        extracted PSMs of the final database
# output:       ProteasomeDB.csv
# authors:      JL, modified by HPR

library(seqinr)
source("SOURCE/loadFunctions.r")

print("-----------------------------")
print("4) MAPPING AND POSTPROCESSING")
print("-----------------------------")


### INPUT ###
load(snakemake@input[["PSMs"]])
load(snakemake@input[["allPSMs_Delta"]])



### MAIN PART ##
d = mapping(extracted)
d_delta = mapping(allPSMs_Delta)


### OUTPUT ###
write.csv(d, file = unlist(snakemake@output[["ProteasomeDB"]]),
          row.names = F)
write.csv(d_delta, file = unlist(snakemake@output[["DeltaPeptides"]]),
          row.names = F)

