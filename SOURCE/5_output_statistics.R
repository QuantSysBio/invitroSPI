### invitroSPI ###
# description:  generate statistics + plots of ProteasomeDB
# input:        ProteasomeDB
# output:       DB statistics, plots
# authors:      HPR

library(dplyr)
library(stringr)
library(grid)
library(ggplot2)
library(gridExtra)

theme_set(theme_bw())

source("SOURCE/invitroSPI_utils.R")
source("SOURCE/plottingFunctions.R")

print("-------------------------------")
print("5) GENERATE DATABASE STATISTICS")
print("-------------------------------")


### INPUT ###
ProteasomeDB = read.csv(snakemake@input[["ProteasomeDB"]],
                        stringsAsFactors = F)
# ProteasomeDB = read.csv("OUTPUT/test_data/ProteasomeDB.csv", stringsAsFactors = F)

S = ProteasomeDB$substrateID %>% unique()
tps = ProteasomeDB$digestTime %>% unique() %>% as.numeric() %>% sort()


### MAIN PART ###

# ----- 1) pre-processing -----

Peps = ProteasomeDB %>%
  remSynthErrors() %>%
  ILredundancy() %>%
  filterPepLength() %>%
  disentangleMultimappers.Type() %>%
  disentangleMultimappers.SRlen() %>%
  disentangleMultimappers.IVSeqLen() %>%
  disentangleMultimappers.AA() %>%
  DBcosmetics()

ProteasomeDB = ProteasomeDB %>%
  remSynthErrors()


# ----- 2) general stats + peptides over time -----

uniquePeps = Peps %>% uniquePeptides()

# plotNumberofPeptides(ProteasomeDB,
#                      outname = "OUTPUT/test_data/number_of_products.pdf")
plotNumberofPeptides(ProteasomeDB,
                     outname = unlist(snakemake@output[["number_of_products"]]))

# pdf("OUTPUT/test_data/DB_stats.pdf", height = 4, width = 6)
pdf(file = unlist(snakemake@output[["DB_stats"]]), height = 4, width = 6)
generalStats(uniquePeps, tp = "all") %>% grid.table()

for (t in 1:length(tps)) {
  cntDB = Peps[Peps$digestTime == tps[t], ] %>%
    uniquePeptides()
  grid.newpage()
  generalStats(DB = cntDB, tp = tps[t]) %>% grid.table()
  
}
dev.off()



# ----- 3) coverage maps -----

# pdf("OUTPUT/test_data/coverage_map.pdf", height = 12, width = 16)
# pdf(unlist(snakemake@output[["coverage_map"]]), height = 12, width = 16)
# for (s in 1:length(S)) {
#   
#   cntDB = uniquePeps[uniquePeps$substrateID == S[s], ]
#   out = plotCoverage(cntDB, name = S[s], tp = "all")
#   replayPlot(out)
#   
# }
# 
# dev.off()

# ----- 4) length distributions + product type frequencies -----

# pdf("OUTPUT/Delta03/length_distributions.pdf", height = 10, width = 14)
pdf(file = unlist(snakemake@output[["length_distributions"]]), height = 12, width = 16)

freq = TypeFrequency(uniquePeps)
pl = PepLength(uniquePeps, tp="all")
srlen = SRLength(uniquePeps, tp="all")
ivlen = IVSeqLength(uniquePeps, tp="all")

cntG = grid.arrange(freq,pl,srlen,ivlen,
                    nrow = 2, ncol = 2)

print(cntG)

for (t in 1:length(tps)) {
  
  cntDB = uniquePeps[uniquePeps$digestTime == tps[t], ]
  
  freq = TypeFrequency(cntDB)
  pl = PepLength(cntDB, tp=tps[t])
  srlen = SRLength(cntDB, tp=tps[t])
  ivlen = IVSeqLength(cntDB, tp=tps[t])
  
  cntG = grid.arrange(freq,pl,srlen,ivlen,
                      nrow = 2, ncol = 2)
  
  print(cntG)
  
}
dev.off()

