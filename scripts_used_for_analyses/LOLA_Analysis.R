---
#title: "Chemosensory gene LOLA analysis"
#original script modify by: "Sergio Olvera"
---

  source("http://bioconductor.org/biocLite.R")
biocLite("LOLA")
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")

devtools::install_github("nsheff/LOLA")
install.packages("path/to/LOLA", repos=NULL)

library(data.table)
library("LOLA")

setwd("Pathway")

regionDB = loadRegionDB("LOLA_database/Species_database")
names(regionDB)
ORregions = readBed("LOLA_database/Gene_regions_LOLA.bed")
universeSet= readBed("LOLA_database/Species_genome_LOLA.bed")                            
locResults = runLOLA(ORregions, universeSet, regionDB, cores=1)
colnames(locResults)
head(locResults)
locResults[order(support, decreasing=TRUE),]
locResults[order(maxRnk, decreasing=TRUE),]

writeCombinedEnrichment(locResults, outFolder= "Gene_file", includeSplits=TRUE)

#Extract the line with the enrinchment value
Collection <- locResults[1,]

#Extraction of the overlaping region
extract_overlap <- extractEnrichmentOverlaps(Collection, ORregions, regionDB)
extract_overlap

#
userSetsRedefined =	redefineUserSets(ORregions, universeSet)
resRedefined = runLOLA(userSetsRedefined, universeSet, regionDB, cores=1)

g = plotTopLOLAEnrichments(resRedefined)