
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biostrings", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("seqinr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

data_wrangling = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### READ top_file ----
  
  top_file<-as.data.frame(fread(file=opt$top_file, sep="\t", header=T), stringsAsFactors=F)

  cat("top_file_0\n")
  cat(str(top_file))
  cat("\n")
  
  #### READ collect_TF_motif_results ----
  
  collect_TF_motif_results<-as.data.frame(fread(file=opt$collect_TF_motif_results, sep="\t", header=T), stringsAsFactors=F)
  
  cat("collect_TF_motif_results_0\n")
  cat(str(collect_TF_motif_results))
  cat("\n")
  
  #### READ TF_search fasta ----

  fastaFile<-readDNAStringSet(file = opt$TF_search)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df_fasta <- data.frame(seq_name, sequence)
  
  cat("df_fasta_0\n")
  cat(str(df_fasta))
  cat("\n")
  
  df_fasta$query_region<-gsub("^>","",df_fasta$seq_name)
  
  cat("df_fasta_1\n")
  cat(str(df_fasta))
  cat("\n")
  
  #### 1 Select concordant motif to CHIPseq factor----
  
  collect_TF_motif_results_concordant<-collect_TF_motif_results[which(collect_TF_motif_results$Symbol == collect_TF_motif_results$CHIP_factor),]
  
  cat("collect_TF_motif_results_concordant_0\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  #### 2 merge by query region to obtain motif ----
  
  collect_TF_motif_results_concordant<-merge(collect_TF_motif_results_concordant,
                                             df_fasta,
                                             by="query_region")
  
  cat("collect_TF_motif_results_concordant_1\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  collect_TF_motif_results_concordant$motif_sequence<-substr(collect_TF_motif_results_concordant$sequence,collect_TF_motif_results_concordant$start,
                                                             collect_TF_motif_results_concordant$end-1)
  
  cat("collect_TF_motif_results_concordant_1.5\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  collect_TF_motif_results_concordant$motif_sequence_rev_comp<-NA
  
  FLAG_NA<-sum(is.na(collect_TF_motif_results_concordant$motif_sequence))
  
  cat("FLAG_NA_0\n")
  cat(str(FLAG_NA))
  cat("\n")
  
  
  
  
  #### Reverse complementary ----
  
  
  # http://www.biotechworld.it/bioinf/2016/02/04/dna-sequence-manipulation-in-r-how-do-i-get-the-reverse-complement-of-a-dna-sequence/
  
  baseMatch <- c("A", "T", "C", "G", "N")
  names(baseMatch) <- c("T", "A", "G", "C", "N")
  
  DEBUG<-0
  
  result<-data.frame()
  
  
  # Function to process each row
  process_row <- function(k, collect_TF_motif_results_concordant, DEBUG) {
    collect_TF_motif_results_concordant_sel <- collect_TF_motif_results_concordant[k,]
    
    if (DEBUG == 1) {
      cat("collect_TF_motif_results_concordant_sel_0\n")
      cat(str(collect_TF_motif_results_concordant_sel))
      cat("\n")
    }
    
    myDNA <- collect_TF_motif_results_concordant_sel$motif_sequence
    
    if (DEBUG == 1) {
      cat("myDNA_0\n")
      cat(sprintf(as.character(myDNA)))
      cat("\n")
    }
    
    # Reverse complement using Biostrings
    revDNA <- reverseComplement(DNAString(myDNA))
    
    if (DEBUG == 1) {
      cat("revDNA_0\n")
      cat(sprintf(as.character(revDNA)))
      cat("\n")
    }
    
    # Complement the reversed sequence (reverse complement is already done)
    complDNA <- as.character(revDNA)
    
    if (DEBUG == 1) {
      cat("complDNA_0\n")
      cat(sprintf(as.character(complDNA)))
      cat("\n")
    }
    
    collect_TF_motif_results_concordant_sel$motif_sequence_rev_comp <- complDNA
    
    return(collect_TF_motif_results_concordant_sel)  # Return the processed row
  }
  
  # Apply the function to each row using lapply
  result_list <- lapply(1:nrow(collect_TF_motif_results_concordant), 
                        process_row, 
                        collect_TF_motif_results_concordant = collect_TF_motif_results_concordant, 
                        DEBUG = 0)
  
  # Combine the list into a data frame
  result <- do.call(rbind, result_list)
  
  # for(k in 1:dim(collect_TF_motif_results_concordant)[1]){
  #   
  #   
  #   collect_TF_motif_results_concordant_sel<-collect_TF_motif_results_concordant[k,]
  #   
  #   if(DEBUG == 1){
  #     
  #     cat("collect_TF_motif_results_concordant_sel_0\n")
  #     cat(str(collect_TF_motif_results_concordant_sel))
  #     cat("\n")
  #   }
  #   
  #   myDNA <- collect_TF_motif_results_concordant_sel$motif_sequence
  #   
  #   if(DEBUG == 1){
  #     cat("myDNA_0\n")
  #     cat(sprintf(as.character(myDNA)))
  #     cat("\n")
  #   }
  #   
  #   
  #   revDNA <- sapply(1:nchar(myDNA), (function(i){
  #     pos <- nchar(myDNA) - i + 1
  #     substr(myDNA, pos, pos)
  #   }))
  #   
  #   revDNA <- paste(revDNA, collapse = "")
  #   
  #   if(DEBUG == 1){
  #     cat("revDNA_0\n")
  #     cat(sprintf(as.character(revDNA)))
  #     cat("\n")
  #   }
  #   
  #   complDNA <- paste(sapply(1:nchar(revDNA), (function(i){
  #     myBase <- substr(revDNA, i,i)
  #     baseMatch[myBase]
  #   })), collapse = "")
  #   
  #   if(DEBUG == 1){
  #     cat("complDNA_0\n")
  #     cat(sprintf(as.character(complDNA)))
  #     cat("\n")
  #   }
  #   
  #   collect_TF_motif_results_concordant_sel$motif_sequence_rev_comp<-complDNA
  #   
  #   result<-rbind(collect_TF_motif_results_concordant_sel,result)
  #   
  # 
  #   
  # }#k in 1:dim(collect_TF_motif_results_concordant
  # 
  # 
  # 
  # 
  # setwd(out)
  # 
  # 
  # write.table(collect_TF_motif_results_concordant,file="test.tsv", sep="\t", quote=F, row.names = F)
  
  cat("result_0\n")
  cat(str(result))
  cat("\n")
  
  
  #### subset prior to merge ----
  
  
  indx.keep<-which(colnames(result)%in%c("query_region","Symbol","source","start","end","Motif_ID","strand","score","cell_type","CHIP_factor","ensembl_gene_id","motif_sequence","motif_sequence_rev_comp"))
  
  
  cat("indx.keep_0\n")
  cat(str(indx.keep))
  cat("\n")
  
  
  result<-unique(result[,indx.keep])
  
  colnames(result)[which(colnames(result) == 'query_region')]<-'Peak_ID'
  colnames(result)[which(colnames(result) == 'start')]<-'motif_rel_start'
  colnames(result)[which(colnames(result) == 'end')]<-'motif_rel_end'
  colnames(result)[which(colnames(result) == 'score')]<-'motif_score'
  colnames(result)[which(colnames(result) == 'strand')]<-'motif_strand'
  
  cat("result_1\n")
  cat(str(result))
  cat("\n")
  
  #### merge with first top file ----
  
  
  collect_TF_motif_results_concordant<-merge(result,
                                             top_file,
                                             by=c("Peak_ID","cell_type","CHIP_factor","ensembl_gene_id"))
  
  cat("collect_TF_motif_results_concordant_4\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  collect_TF_motif_results_concordant$motif_start_38<-collect_TF_motif_results_concordant$Peak_start_38+collect_TF_motif_results_concordant$motif_rel_start-1
  
  collect_TF_motif_results_concordant$motif_end_38<-collect_TF_motif_results_concordant$Peak_start_38+collect_TF_motif_results_concordant$motif_rel_end -1
  
  cat("collect_TF_motif_results_concordant_4\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  collect_TF_motif_results_concordant<-collect_TF_motif_results_concordant[order(collect_TF_motif_results_concordant$cell_type,
                                                                                 collect_TF_motif_results_concordant$CHIP_factor,
                                                                                 collect_TF_motif_results_concordant$motif_score,
                                                                                 decreasing = T),]
  
  cat("collect_TF_motif_results_concordant_5\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  indx.reorder<-which(colnames(collect_TF_motif_results_concordant)%in%c("Peak_ID","cell_type","CHIP_factor",'Motif_ID','motif_score','motif_strand','source','chrom','motif_start_38','motif_end_38','motif_sequence','motif_sequence_rev_comp'))
  
  cat("indx.reorder_0\n")
  cat(str(indx.reorder))
  cat("\n")
  
  unselected_columns<-colnames(collect_TF_motif_results_concordant)[-which(colnames(collect_TF_motif_results_concordant)%in%c("Peak_ID","cell_type","CHIP_factor",'Motif_ID','motif_score','motif_strand','source','chrom','motif_start_38','motif_end_38','motif_sequence','motif_sequence_rev_comp'))]
  
  cat("unselected_columns_0\n")
  cat(str(unselected_columns))
  cat("\n")
  
  indx.not.reorder<-which(colnames(collect_TF_motif_results_concordant)%in%unselected_columns)
  
  cat("indx.not.reorder_0\n")
  cat(str(indx.not.reorder))
  cat("\n")
  
  collect_TF_motif_results_concordant<-collect_TF_motif_results_concordant[,c(indx.reorder,indx.not.reorder)]
  
  cat("collect_TF_motif_results_concordant_6\n")
  cat(str(collect_TF_motif_results_concordant))
  cat("\n")
  
  #### save ----
  
  setwd(out)
  
  write.table(collect_TF_motif_results_concordant,file="Final_table_CHIP_controls.tsv", sep="\t", quote=F, row.names = F)
  
  
}







printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--top_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--collect_TF_motif_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_search"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  data_wrangling(opt)


  
}


###########################################################################

system.time( main() )