
##  Script to retrieve DNA methylation or chromatin accessibility over specific genomic region ##

setwd("/mnt/hel/home/qinluo/data/scNMt")

library(argparse)
library(GenomicRanges)
library(doParallel)
library(dplyr)


# Initialize argument parser
p <- ArgumentParser(description='')

#CG was the methylation site for accessibility evaluation, GC is for methylation evaluation
p$add_argument('-c','--context', type="character",help='CG or GC')

#use mclapply for parallel computing, define number of cores according to machine availablity
p$add_argument('-n','--cores', type="integer",help='Number of cores')

#arguments define the region where to look at
p$add_argument('-c','--chromosome', type="character",help='chromosome with chr')
p$add_argument('-s','--start', type="integer",help='start position')
p$add_argument('-e','--end', type="integer",help='end position')

#define output
p$add_argument('-o','--output', type="character",help='output file')


# Read arguments
args <- p$parse_args(commandArgs(TRUE))



#define the target region as a genomic range object
ir_seq = IRanges(start=args$start, end=args$end)
gr_seq = GRanges(seqnames = args$chromosome,ranges=ir_seq)


#define the function to retrieve the values in the target region
retrieve_signal <- function(input) {
  
  rt = read.table(input,header=T)
  rt = cbind(rt,cell=as.character(input))
  
  ir = IRanges(start =rt[,2],end=rt[,2])
  gr = GRanges(seqnames = paste0("chr",rt[,1]),ranges=ir,mcols=as.data.frame(rt))
  
  hits1= findOverlaps(gr_seq,gr, ignore.strand=T)
  
  if (length(hits1)==0){
         table = NA
         }
  else {
   
        PE = gr[subjectHits(hits1)]
        PE = as.data.frame(mcols(PE))
 
        table = PE
        
        #output table contains the information of each detected site of each single cell
        colnames(table) = c("Chr","Pos","met_reads","nonmet_reads","rate","cell")
        }
  
   return(table)

}


#if context input was "CG", do analysis on the accessibility data
if (opts$context=="CG"){
    acc_files = list.files(pattern=".acc.tsv.gz",full.names = T)
    acc_name = list.files(pattern=".acc.tsv.gz",full.names = F)


    data_table = mclapply(acc_files,retrieve_signal,args$cores)
    na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
    data.table = na.omit.list(data_table)

    output =bind_rows(data.table)
    write.table(output,args$output,sep="\t",row.names = F,quote = F)
}



#if context input was "CG", do analysis on the accessibility data
if (opts$context=="GC"){
    acc_files = list.files(pattern=".meth.tsv.gz",full.names = T)
    acc_name = list.files(pattern=".meth.tsv.gz",full.names = F)

    data_table = mclapply(acc_files,retrieve_signal,args$cores)
    na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }
    data.table = na.omit.list(data_table)

    output =bind_rows(data.table)
    write.table(output,args$outputt,sep="\t",row.names = F,quote = F)
}
