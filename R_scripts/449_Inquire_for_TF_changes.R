
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
  
  #### READ and transform TF_REF ----
  
  TF_REF = opt$TF_REF
  
  cat("TF_REF_\n")
  cat(sprintf(as.character(TF_REF)))
  cat("\n")
  
  #### READ and transform TF_ALT ----
  
  TF_ALT = opt$TF_ALT
  
  cat("TF_ALT_\n")
  cat(sprintf(as.character(TF_ALT)))
  cat("\n")
  
  #### READ CHIP_variants_file ----
  
  CHIP_variants_file<-as.data.frame(fread(file=opt$CHIP_variants_file, sep="\t", header=T), stringsAsFactors=F)
  
  cat("CHIP_variants_file_0\n")
  cat(str(CHIP_variants_file))
  cat("\n")
  
 
  #### READ Final_table_TF_motif_prediction ----
  
  Final_table_TF_motif_prediction<-as.data.frame(fread(file=opt$Final_table_TF_motif_prediction, sep="\t", header=T), stringsAsFactors=F)

  cat("Final_table_TF_motif_prediction_0\n")
  cat(str(Final_table_TF_motif_prediction))
  cat("\n")
 
  
  Final_table_TF_motif_prediction_REF<-Final_table_TF_motif_prediction[which(Final_table_TF_motif_prediction$allele == "REF"),]
  
  cat("Final_table_TF_motif_prediction_REF_0\n")
  cat(str(Final_table_TF_motif_prediction_REF))
  cat("\n")
  
  Final_table_TF_motif_prediction_REF_sel<-Final_table_TF_motif_prediction_REF[which(Final_table_TF_motif_prediction_REF$Symbol == TF_REF & Final_table_TF_motif_prediction_REF$Intersect_SNP == "YES"),]
  
  cat("Final_table_TF_motif_prediction_REF_sel_0\n")
  cat(str(Final_table_TF_motif_prediction_REF_sel))
  cat("\n")
  
  
  Final_table_TF_motif_prediction_ALT<-Final_table_TF_motif_prediction[which(Final_table_TF_motif_prediction$allele == "ALT"),]
  
  cat("Final_table_TF_motif_prediction_ALT_0\n")
  cat(str(Final_table_TF_motif_prediction_ALT))
  cat("\n")
  
  Final_table_TF_motif_prediction_ALT_sel<-Final_table_TF_motif_prediction_ALT[which(Final_table_TF_motif_prediction_ALT$Symbol == TF_ALT & Final_table_TF_motif_prediction_ALT$Intersect_SNP == "YES"),]
  
  cat("Final_table_TF_motif_prediction_ALT_sel_0\n")
  cat(str(Final_table_TF_motif_prediction_ALT_sel))
  cat("\n")
  
  Final_table_TF_motif_prediction_ALT_unwanted<-Final_table_TF_motif_prediction_ALT[which(Final_table_TF_motif_prediction_ALT$Symbol == TF_REF & Final_table_TF_motif_prediction_ALT$Intersect_SNP == "YES"),]
  
  cat("Final_table_TF_motif_prediction_ALT_unwanted_0\n")
  cat(str(Final_table_TF_motif_prediction_ALT_unwanted))
  cat("\n")
  
  
  #### overlap 1 ----
  
  indx.1<-which(Final_table_TF_motif_prediction_REF_sel$VAR_38%in%Final_table_TF_motif_prediction_ALT_sel$VAR_38)
  
  VARS_sel_1<-unique(Final_table_TF_motif_prediction_REF_sel$VAR_38[indx.1])
  
  
  cat("VARS_sel_1_0\n")
  cat(str(VARS_sel_1))
  cat("\n")
  
  indx.2<-which(VARS_sel_1%in%Final_table_TF_motif_prediction_ALT_unwanted$VAR_38)
  
  
  VARS_sel_DEF<-VARS_sel_1[-indx.2]
  
  cat("VARS_sel_DEF_0\n")
  cat(str(VARS_sel_DEF))
  cat("\n")
  
  
  if(length(VARS_sel_DEF)>0){
    
    Final_table_TF_motif_prediction_SET<-Final_table_TF_motif_prediction[which(Final_table_TF_motif_prediction$VAR_38%in%VARS_sel_DEF),]
    
    cat("Final_table_TF_motif_prediction_SET_0\n")
    cat(str(Final_table_TF_motif_prediction_SET))
    cat("\n")
    
    Final_table_TF_motif_prediction_SET<-merge(CHIP_variants_file,
                                                Final_table_TF_motif_prediction_SET,
                                                by="snp")

    cat("Final_table_TF_motif_prediction_SET_1\n")
    cat(str(Final_table_TF_motif_prediction_SET))
    cat("\n")
    
    write.table(Final_table_TF_motif_prediction_SET, file=paste('TF_break_selected_',TF_REF,'_',TF_ALT,".tsv",sep=''), sep="\t", quote=F, row.names = F)
    
    
  }# length(VARS_sel_DEF)>0
  
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
    make_option(c("--Final_table_TF_motif_prediction"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_REF"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_ALT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--CHIP_variants_file"), type="character", default=NULL, 
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