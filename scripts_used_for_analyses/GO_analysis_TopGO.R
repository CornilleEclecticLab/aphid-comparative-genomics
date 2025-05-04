---
#title: "Gene Ontology analysis"
#original script modify by: "Sergio Olvera"
---
source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
biocLite("GO.db")
biocLite("biomaRt")
biocLite("Rgraphviz")
BiocManager::install("Rgraphviz")
# Load the required R packages
library(topGO)
library(GO.db)
library(biomaRtr)
library(Rgraphviz)

GOenrich=function(GO.file, genes, background, topGOFiles,
                  list.ont = c("BP", "MF", "CC"),
                  pval =  0.01, min.genes.tot=5, min.genes.signif=5, max.nb.genes=10){
  require(topGO)
  ZmB73_5a_xref.GO.topGO.FG <- read.delim(GO.file, header=FALSE, stringsAsFactors = F)
  geneID2GO <- readMappings(file = GO.file)
  
  ## Prepare background gene list
  geneNames <- ZmB73_5a_xref.GO.topGO.FG$V1
  geneNames <- geneNames[geneNames %in% background]
  
  ## Prepare gene list
  myInterestingGenes=unique(genes)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  
  ## Run GO enrichment analysis
  allRes=list() # List of GO enrichment analysis results
  GOtermsGenes=NULL # Table containing gene list for each enriched GO ID
  
  if(length(genes)>=min.genes.tot){ ## Run topGO only if the list of genes to be analyzed contain more than 5 genes
    
    for(ont in list.ont){ # Run topGo for each ontology
      GOdata <- new("topGOdata", ontology = ont, allGenes = geneList,
                    annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
      
      ## Compute stats using elim strategy
      test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
      resultElimFisher <- getSigGroups(GOdata, test.stat)
      allRes[[`ont`]] <- GenTable(GOdata,
                                  elim = resultElimFisher, orderBy = "elim",
                                  ranksOf = "elim", topNodes = length(resultElimFisher@score))
      allRes[[`ont`]]$ontology=rep(ont, nrow(allRes[[`ont`]]))
      
      #remove potential 0
      if (min(resultElimFisher@score)==0){
        resultElimFisher@score[which(resultElimFisher@score==0)]=min(resultElimFisher@score[-(which(resultElimFisher@score==0))])
      }
      
      ## Plot topGO output
      pdf(paste0(topGOFiles, "_", ont, ".pdf"), width=11.5, height=8)
      showSigOfNodes(GOdata, score(resultElimFisher), firstSigNodes = 10, useInfo ='all')
      dev.off()
      
      ## Write gene list for each enriched GO ID
      myterms <- allRes[[`ont`]]$GO.ID
      mygenes  <- genesInTerm(GOdata, myterms)
      
      for (j in 1:length(myterms)) {
        myterm <- myterms[j]
        mygenesforterm <- mygenes[myterm][[1]]
        myfactor <- mygenesforterm %in% myInterestingGenes # find the genes that are in the list of genes of interest
        mygenesforterm2 <- mygenesforterm[myfactor == TRUE]
        mygenesforterm2 <- paste(mygenesforterm2, collapse=',')
        GOtermsGenes <- rbind(GOtermsGenes,
                              data.frame("Ontology"=ont, "GOTerm"=myterm,"Genes"=mygenesforterm2))
      }
    }
    allResfinal=rbind(allRes$BP, allRes$MF, allRes$CC)
    ## write only significant results (3 genes minimum in the group + )
    write.table(allResfinal[(allResfinal$Significant >= min.genes.signif) & (allResfinal$elim <= pval),],
                file=paste0(topGOFiles, "_enrichment.txt"),
                quote=FALSE,row.names=FALSE, sep="\t")
    write.table(GOtermsGenes, file=paste0(topGOFiles, "_term_genes.txt"),
                quote=FALSE,row.names=FALSE, sep="\t")
    return(list(enrichment=allResfinal, genelist=GOtermsGenes))
  }
  else {return(list(enrichment="", genelist=""))}
}

setwd("pathway")
# I create a vector of the genes ids that are of interest to me.
genes_of_interest=c(read.table("pathway/file_genes_of_interest.txt")$V1)

# I create a vector contaning all the gene ids of my organfism
genes_all=unique(c(read.table("/pathway/file_gene_list.txt"))$V1)

#GOenrich(GO.file, genes, background, topGOFiles,list.ont = c("BP", "MF", "CC"),pval =  0.01, min.genes.tot=5, min.genes.signif=5, max.nb.genes=10)
#Parameters : 
# Go.file : file contaning as first column the genes ids and second column their Goterms
# genes : vector of our gene of interest
# background : vector of all gene gene ids of the organism
# topGOFiles : Prefix used to name all the output files

GOenrich("pathway/file_ID_GO.txt",genes_of_interest,genes_all,"Single_copy")
