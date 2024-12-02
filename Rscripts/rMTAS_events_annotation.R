library('biomaRt')
library('org.Mm.eg.db')
library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(ggVennDiagram)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library("VennDiagram")

get_enrichPathway <- function(symbols, pcutoff=0.05, ncat=10) {
  
  gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  gene = gene[!is.na(gene)]
  print(paste0(length(gene), ' mapped genes'))
  yy = enrichPathway(gene, pvalueCutoff=pcutoff, organism = "mouse", readable = TRUE)
  print(head(as.data.frame(yy)))
  if (length(yy) > 0) {
    print(dotplot(yy, showCategory=ncat))
    return(yy)
  }
  return(NULL)
}

get_enrichGO <- function(symbols, pcutoff=0.01, qcutoff=0.05, ncat=20, ont="BP") {
  
  gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  gene = gene[!is.na(gene)]
  print(paste0(length(gene), ' mapped genes'))
  
  ego <- enrichGO(gene          = gene,
                  OrgDb         = org.Mm.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = pcutoff,
                  qvalueCutoff  = qcutoff,
                  readable      = TRUE)
  head(ego)
  if (length(ego) > 0) {
    print(barplot(ego, showCategory=ncat))
    return(ego)
  }
  return(NULL)
}

annotate_genes <- function(mats, group) {
  
  if (group == 'DOI'){
    mats = mats[mats$IncLevelDifference < 0,] ## Treatment
  } else {
    mats = mats[mats$IncLevelDifference > 0,] ## CTRL
  }
  symbols = unique(mats$GeneID)
  get_enrichPathway(symbols)
  get_enrichGO(symbols)
  
}
