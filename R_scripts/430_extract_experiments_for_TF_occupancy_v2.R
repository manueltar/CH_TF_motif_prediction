
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
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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
suppressMessages(library("liftOver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

collect = function(option_list)
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
  
  #### READ and transform table_sel ----
  
  table_sel = opt$table_sel
  
  cat("table_sel_\n")
  cat(sprintf(as.character(table_sel)))
  cat("\n")
  
 
  #### READ List_of_TFs_for_occupancy ----
  
  List_of_TFs_for_occupancy<-as.data.frame(fread(file=opt$List_of_TFs_for_occupancy, sep="\t", header=T), stringsAsFactors=F)

  cat("List_of_TFs_for_occupancy_0\n")
  cat(str(List_of_TFs_for_occupancy))
  cat("\n")
 
  
  #### READ fileList ----
  
  fileList<-as.data.frame(fread(file=opt$fileList, sep="\t", header=F), stringsAsFactors=F)
  
  colnames(fileList)<-c('filename','genome_assembly','Track_type_class','Track_type','Cell_type_class','Cell_type','Threshold','Experimental_IDs_included')

  cat("fileList_0\n")
  cat(str(fileList))
  cat("\n")
  
  
  #### LOOP TFs -----
  
  array_TFs<-unique(List_of_TFs_for_occupancy$ensembl_gene_id)
  
  cat("array_TFs_0\n")
  cat(str(array_TFs))
  cat("\n")
  
  DEBUG<-0
  
  Results<-data.frame()
  for(i in 1:length(array_TFs)){
    
    array_TFs_sel<-array_TFs[i]
    
    cat("---------------------------------->\t")
    cat(sprintf(as.character(array_TFs_sel)))
    cat("\n")
    
    List_of_TFs_for_occupancy_sel<-List_of_TFs_for_occupancy[which(List_of_TFs_for_occupancy$ensembl_gene_id == array_TFs_sel),]
    
    if(DEBUG == 1){
      
      cat("List_of_TFs_for_occupancy_sel_0\n")
      cat(str(List_of_TFs_for_occupancy_sel))
      cat("\n")
    }
    
    fileList_sel<-fileList[which(fileList$genome_assembly == 'hg19' &
                                   fileList$Cell_type_class == 'Blood' &
                                   fileList$Track_type%in%List_of_TFs_for_occupancy_sel$Symbol),]
    
    if(DEBUG == 1){
      
      cat("fileList_sel_0\n")
      cat(str(fileList_sel))
      cat("\n")
    }
    
    if(dim(fileList_sel)[1] >0){
      
      experiments<-unique(unlist(strsplit(fileList_sel$Experimental_IDs_included, split= ',')))
      
      if(DEBUG == 1){
        
        cat("experiments_0\n")
        cat(str(experiments))
        cat("\n")
      }
      
      
      a.df<-as.data.frame(cbind(rep(array_TFs_sel, length(experiments)),
                                rep(paste(List_of_TFs_for_occupancy_sel$Symbol, collapse=';'), length(experiments)),
                                experiments), stringsAsFactors=F)
      
      colnames(a.df)<-c('ensembl_gene_id','Symbol_string','experiments')
      
      if(DEBUG == 1){
        
        cat("a.df_0\n")
        cat(str(a.df))
        cat("\n")
      }
      
      Results<-unique(rbind(a.df,Results))
      
     
      
    }#dim(fileList_sel)[1] >0
  }#i in 1:length(array_TFs)
  
  
  cat("Results_0\n")
  cat(str(Results))
  cat("\n")
  
  
  ############# split and save ---------------------------
  
  
  n <- 1000
  
  nr <- dim(Results)[1]
  
  List_partitions<-split(Results, rep(1:ceiling(nr/n), each=n, length.out=nr))
  
  cat("List_partitions_0\n")
  cat(str(List_partitions))
  cat("\n")

  array_partitions<-names(List_partitions)  
  
  for(i in 1:length(array_partitions)){
    
    array_partitions_sel<-array_partitions[i]
    
    
    cat("---------------------------------->\t")
    cat(sprintf(as.character(array_partitions_sel)))
    cat("\n")
    
    
    df_partitioned<-List_partitions[[array_partitions_sel]]
    
    
    cat("df_partitioned_0\n")
    cat(str(df_partitioned))
    cat("\n")
    
    setwd(out)
    
    
    write.table(df_partitioned, file=paste('TF_experiments_search_file_',table_sel,"_",array_partitions_sel,".tsv",sep=''), sep="\t", quote=F, row.names = F, col.names = F)
    
    
  }#i in 1:length(array_partitions)
  
  


  
  
  
 
  
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
    make_option(c("--table_sel"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--motifs_collected"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--List_of_TFs_for_occupancy"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fileList"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--spanning_of_motif"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--allele"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  collect(opt)


  
}


###########################################################################

system.time( main() )