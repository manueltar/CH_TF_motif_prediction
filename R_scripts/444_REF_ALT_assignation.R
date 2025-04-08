
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
  
  #### READ and transform GWAS_file ----
  
  GWAS_file = as.data.frame(fread(file=opt$GWAS_file, sep="\t", header=T), stringsAsFactors=F)
  
  cat("GWAS_file_0\n")
  cat(str(GWAS_file))
  cat("\n")

  
  #### READ input_bed ----
  
  input_bed<-as.data.frame(fread(file=opt$input_bed, sep="\t", header=F), stringsAsFactors=F)
  colnames(input_bed)<-c('chr','start','end','name','score','strand')
 
  cat("input_bed_0\n")
  cat(str(input_bed))
  cat("\n")
 
  input_bed$snp<-gsub("__.+$","",input_bed$name)
  input_bed$tag<-gsub("^[^__]+__","",input_bed$name)
  input_bed$tag<-gsub("__.+$","",input_bed$tag)
  input_bed$rsid<-gsub("^[^__]+__","",input_bed$name)
  input_bed$rsid<-gsub(paste(unique(input_bed$tag), collapse="|"),"",input_bed$rsid)
  input_bed$rsid<-gsub("^__","",input_bed$rsid)
  
  cat("input_bed_1\n")
  cat(str(input_bed))
  cat("\n")
  
  #### READ fastaFile ----
  
  fastaFile<-readDNAStringSet(file = opt$input_fasta)
  seq_name = names(fastaFile)
  ref = paste(fastaFile)
  df_fasta <- data.frame(seq_name, ref)
  
  cat("df_fasta_0\n")
  cat(str(df_fasta))
  cat("\n")
  
  df_fasta$seq_name<-gsub("^>","",df_fasta$seq_name)
  df_fasta$chr<-gsub(":.+$","",df_fasta$seq_name)
  df_fasta$start<-gsub("^[^:]+:","",df_fasta$seq_name)
  df_fasta$start<-gsub("-.+$","",df_fasta$start)
  df_fasta$end<-gsub("^[^:]+:[0-9]+-","",df_fasta$seq_name)
  
  cat("df_fasta_1\n")
  cat(str(df_fasta))
  cat("\n")
  
  #### Merge fasta with bed ----
  
  input_bed<-unique(merge(input_bed,
                   df_fasta,
                   by=c('chr','start','end')))
  
  cat("input_bed_POST_merge_with_fasta_ref\n")
  cat(str(input_bed))
  cat("\n")
  
  #### Merge with GWAS_file ----
  
  input_bed<-unique(merge(input_bed,
                   GWAS_file,
                   by=c('chr','snp','rsid','tag')))
  
  cat("input_bed_POST_merge_with_GWAS_file\n")
  cat(str(input_bed))
  cat("\n")
  
  #### Assign ref and alt -----
  
  input_bed$alt<-NA
  
  indx.1<-which(input_bed$a0 == input_bed$ref)
  
  input_bed$alt[indx.1]<-input_bed$a1[indx.1]
  
  
  indx.2<-which(input_bed$a1 == input_bed$ref)
  
  input_bed$alt[indx.2]<-input_bed$a0[indx.2]
  

  cat("POST_ASSIGNATION\n")
  cat(str(input_bed))
  cat("\n")
  
  FLAG_unassigned<-sum(is.na(input_bed$alt))
  
  
  cat("FLAG_unassigned_0\n")
  cat(str(FLAG_unassigned))
  cat("\n")
  
  if(FLAG_unassigned >0)
  {
    table_NA<-input_bed[is.na(input_bed$alt),]
    
    cat("table_NA_0\n")
    cat(str(table_NA))
    cat("\n")
    
    
    input_bed_rest<-input_bed[!is.na(input_bed$alt),]
    
    cat("input_bed_rest_0\n")
    cat(str(input_bed_rest))
    cat("\n")
    
    
    table_NA_no_leftovers_for_alignment_check<-table_NA[which(table_NA$length > 1 & table_NA$tag != 'check_first'),]
    
    cat("table_NA_no_leftovers_for_alignment_check_0\n")
    cat(str(table_NA_no_leftovers_for_alignment_check))
    cat("\n")
    
    if(dim(table_NA_no_leftovers_for_alignment_check)[1] >0){
      
      setwd(out)
      
      write.table(input_bed, file='Unassigned_ref_and_alt.tsv', sep="\t",quote=F,row.names=F)
      
      stop("Unassigned alt variants persist\n")
      
      
      
    }else{
      
      setwd(out)

      write.table(input_bed_rest, file='config_file_assigned_ref_and_alt.tsv', sep="\t",quote=F,row.names=F)
      
      
    }# dim(table_NA_no_leftovers_for_alignment_check)[1] >0
    
    
  }else{
    
    setwd(out)
    
    write.table(input_bed, file='config_file_assigned_ref_and_alt.tsv', sep="\t",quote=F,row.names=F)
    
    
  }#FLAG_unassigned >0
  
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
    make_option(c("--input_bed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_fasta"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GWAS_file"), type="character", default=NULL, 
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
  
  data_wrangling(opt)
  
  
}


###########################################################################

system.time( main() )