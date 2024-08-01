library(tidyverse)
library(hrbrthemes)
library(ggVennDiagram)
library("VennDiagram")
library(ggbio)
library(foreach)
library(doMC)

registerDoMC(cores = 4)


check_UCSC_seqnames <- function(obj) {
  
  GRCm39_table = read.delim("files/sequence_report.tsv", sep = "\t")
  
  if (grepl('NC_', levels(seqnames(obj))[1])) {
    obj <- obj[as.character(seqnames(obj)) %in% GRCm39_table$RefSeq.seq.accession, ]
    obj_df <- as.data.frame(obj)
    obj_df$seqnames <- plyr::mapvalues(x = obj_df$seqnames, from = GRCm39_table$RefSeq.seq.accession, to = GRCm39_table$UCSC.style.name)
    obj <- makeGRangesFromDataFrame(obj_df, keep.extra.columns = T)
  } 
  obj = keepStandardChromosomes(obj, pruning.mode="coarse")
}


build_coordinates_gr <- function(res, et, site='single') {
  
  if (et == "SE") {
    aux = res[, c("chr", "strand", "exonStart_0base", "exonEnd", "GeneID", "source", "IncLevelDifference", "status", "passed")]
  } else if (et == "A5SS") {
    aux = res[, c("chr", "strand", "shortES", "shortEE", "GeneID", "source", "IncLevelDifference","status","passed")]
  } else if (et == "A3SS") {
    aux = res[, c("chr", "strand", "shortES", "shortEE", "GeneID", "source", "IncLevelDifference","status","passed")]
  } else if (et == "MXE") {
    aux = res[, c("chr", "strand", "X2ndExonStart_0base", "X2ndExonEnd", "GeneID", "source", "IncLevelDifference","status","passed")]
  } else if (et == "RI") {
    aux = res[, c("chr", "strand", "riExonStart_0base", "riExonEnd", "GeneID", "source", "IncLevelDifference","status","passed")]
  }
  colnames(aux) <- c("chr", "strand", "start", "end", "name", "source", "IncLevelDifference","status","passed")
  aux$event_type = et
  aux_gr = makeGRangesFromDataFrame(aux, keep.extra.columns = T)
  aux_gr = check_UCSC_seqnames(aux_gr)
  
  return(aux_gr)
}

mark_filter_events <- function(res, qt=0.8, psi=0.1){
  #q1 = quantile(res$avgPSI_1, qt)
  #q2 = quantile(res$avgPSI_2, qt)
  res$PSI_quantile_passed = (res$avgPSI_1 >= qt & res$avgPSI_2 >= qt)
  res$IncLevelDifference_passed = abs(res$IncLevelDifference) >= psi
  
  res$passed = ifelse(res$PSI_quantile_passed & res$IncLevelDifference_passed, "yes", "no")
  return(res)
}

get_AS_associated_genes <- function(res) {
  data_ = res[, c('study','source', 'GeneID')]
  return(data_[!duplicated(data_),])
}

get_overlapping_genes <- function(x) {
  t=get.venn.partitions(x, keep.elements = T, force.unique = T)
  return(t$..values..[[1]])
}

jaccard_index <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

overlap_coef <- function(a, b){
  intersection = length(intersect(a, b))
  den = min(length(a), length(b))
  return (intersection/den)
}

count_overlapping_events <- function(gr) {
  ## baseline 
  gr_baseline = reduce(gr)
  #save_GR_bed_file(gr_baseline, dir = "files/", file_name = "baseline_events.bed", extra=FALSE)
  gr_list_source <- split(x = gr, f = gr$source)
  
  overlaps <- foreach(i=1:length(gr_list_source), .combine = cbind)%dopar%{
    gr_source <- gr_list_source[[i]]
    overlap <- countOverlaps(query = gr_baseline, subject = gr_source, type = 'within') > 0
    return(overlap)
  }
  colnames(overlaps) <- names(gr_list_source)
  gr_baseline_df = as.data.frame(gr[findOverlaps(gr_baseline, gr, select = "first"),])
  count = rowSums(overlaps)
  aux = cbind(gr_baseline_df[, -which(colnames(gr_baseline_df) == 'source')], overlaps)
  aux$counts = count
  gr_overlapps = makeGRangesFromDataFrame(aux, keep.extra.columns = T)
  return(gr_overlapps)
}







