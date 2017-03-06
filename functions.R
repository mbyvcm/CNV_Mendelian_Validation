# Christopher Medway
# March 2017
# Script for mendelian validation of CNVs using trios


penncnv_to_list <- function(x) {
  
  # read file for fam member
  x <-read.table(x, sep = "", stringsAsFactors = F)
  # split out chr:start-end
  l <- lapply(strsplit(x[,1], split = ':'), function(x) c(x[1],x[2]))
  # extract important values for making granges object
  chr   <- unlist(lapply(l, function(x) x[1]))
  start <- unlist(lapply(l, function(x) lapply(strsplit(x[2], split = '-'), function(x) x[1])))
  end   <- unlist(lapply(l, function(x) lapply(strsplit(x[2], split = '-'), function(x) x[2])))
  id    <- stringr::str_sub(string = basename(x$V5), start = 1, end = 5)
  # split into list by IID 
  df <- data.frame(id, chr, start, end, stringsAsFactors = F)
  return(split(df, f = df$id))
}


run_validation <- function(ALN, child, mother, father) {
  
  # extract df from list for single family
  child  <- child[[ALN]]
  mother <- mother[[ALN]]
  father <- father[[ALN]]
  
  # get vector of fractions 
  mat <- calculate_frac_overlap(child,mother)
  pat <- calculate_frac_overlap(child,father)
  
  # return childs df with mat and pat fractions added
  return(cbind(child,mat,pat))
}


calculate_frac_overlap <- function(child, parent) {
  
  par <- vector(mode = "numeric",length = dim(child)[1])
  
  # make query granges
  qgr <- makegranges(child)
  # make subject granges
  sgr <- makegranges(parent)
  
  # find overlaps
  ol   <- findOverlaps(qgr,sgr)
  # query row(s) containing hits 
  qrow <- queryHits(ol)
  # subject row(s) containing hits
  srow <- subjectHits(ol)
  
  # length of parentCNV + childCNV (inc flanking)
  total_length <- vector(mode = "numeric", length = length(ol))
  
  for (i in seq(length(ol))) {
    max <-  as.numeric(max(child[qrow[i],'end'], parent[srow[i],'end']))
    min <-  as.numeric(min(child[qrow[i],'start'], parent[srow[i],'start']))
    total_length[i] <- max - min 
  }
  
  # length of overlap only
  ollength <- width(intersect(qgr,sgr))
  
  # fraction of CNV that is overlapping (exact CNV = 1)
  par[qrow] <- ollength/total_length
  return(par)
}

makegranges <- function(x) {
  return(GRanges(seqnames = x$chr, ranges = IRanges(start = as.numeric(x$start), end = as.numeric(x$end) )))
}