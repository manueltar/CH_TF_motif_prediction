
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
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("gwasrapidd", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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

options(warn = 0)

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
  
  #### READ and transform GWAS_equivalence ----
  
  GWAS_equivalence = read.table(file=opt$GWAS_equivalence,sep="\t", header=T)
  
  cat("GWAS_equivalence_0\n")
  cat(str(GWAS_equivalence))
  cat("\n")
  
  #### READ and transform Threshold_PP ----
  
  Threshold_PP = opt$Threshold_PP
  
  cat("Threshold_PP_\n")
  cat(sprintf(as.character(Threshold_PP)))
  cat("\n")
  
 
  #### READ and transform indir ----
  
  indir = opt$indir
  
  cat("indir_0\n")
  cat(sprintf(as.character(indir)))
  cat("\n")
  
  
  setwd(indir)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("_susie_finemap\\.rds$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
  
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'LocusID'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  
  df_files$path<-paste(indir,files_df_sel,sep='')
  df_files$study_id<-gsub("^CHIP_GWAS_","",files_df_sel)
  df_files$study_id<-gsub("_.+$","",df_files$study_id)
  df_files$study_id<-paste("CHIP_GWAS_",df_files$study_id,sep='')
  
  
  df_files$LocusID<-gsub("^.+_locus_","",files_df_sel)
  df_files$LocusID<-gsub("_susie_finemap\\.rds$","",df_files$LocusID)
  
  df_files$LocusID<-paste(df_files$study_id,"locus",df_files$LocusID, sep='_')

  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
  
  ##### Reading loop -----
  
  LocusID_array<-df_files$LocusID
  
  cat("LocusID_array_0\n")
  cat(str(LocusID_array))
  cat("\n")
  

  DEBUG<-0
  
  list_results<-list()
  
  
  for(i in 1:length(LocusID_array))
  {
    LocusID_array_sel<-LocusID_array[i]
    
    cat("------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(LocusID_array_sel)))
    cat("\n")
    
    df_files_sel<-df_files[which(df_files$LocusID == LocusID_array_sel),]
    
    study_id_sel<-unique(df_files_sel$study_id)
    
    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(df_files_sel)[1] >0)
    {
      sel_file<-df_files_sel$path
      
      # cat("--->\t")
      # cat(sprintf(as.character(sel_file)))
      # cat("\n")
      
      SIZE_gate<-file.info(sel_file)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(sel_file))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          Result_DF<-readRDS(file=sel_file)
          
          
          FM_df<-Result_DF[['finemapping_lABFs']]
          
          if(DEBUG == 1)
          {
            cat("FM_df_0\n")
            cat(str(FM_df))
            cat("\n")
          }
          
          FM_df_sel<-FM_df[which(FM_df$is_cs == TRUE),]
          
          if(DEBUG == 1)
          {
            cat("FM_df_sel_0\n")
            cat(str(FM_df_sel))
            cat("\n")
          }
          
          if(dim(FM_df_sel)[1] >0){
            
            FM_df_sel$LocusID<-LocusID_array_sel
            FM_df_sel$study_id<-study_id_sel
            
            indx.col<-which(colnames(FM_df_sel) == 'SNP.PP')
            
            
            FM_df_sel_filtered<-FM_df_sel[which(FM_df_sel[,indx.col] > Threshold_PP),]
            
            if(dim(FM_df_sel_filtered)[1] >0){
              
              if(DEBUG == 1)
              {
                cat("FM_df_sel_filtered_0\n")
                cat(str(FM_df_sel_filtered))
                cat("\n")
              }
              
              list_results[[i]]<-FM_df_sel_filtered
              
            }# dim(FM_df_sel_filtered)[1] >0
            

            
          }#dim(FM_df_sel)[1] >0
          
         
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
    
  }# i in 1:length(LocusID_array)
  
  if(length(list_results) >0)
  {
    Results = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    cat("Results_0\n")
    cat(str(Results))
    cat("\n")
    
    Results<-merge(Results,
                   GWAS_equivalence[,c(which(colnames(GWAS_equivalence)%in%c('input','study_id')))],
                   by="study_id",
                   all.x=T)
    
    cat("Results_1\n")
    cat(str(Results))
    cat("\n")
    
    
    Results$GWAS_catalogID<-gsub("_buildGRCh38_formatted\\.tsv\\.gz","",Results$input)
    Results$GWAS_catalogID<-gsub("\\/ssu\\/bsssu\\/michelas_finemapping\\/input\\/","",Results$GWAS_catalogID)
    
    
    cat("Results_2\n")
    cat(str(Results))
    cat("\n")
    
    array_GWAS_catalogID<-unique(Results$GWAS_catalogID)
    
    cat("array_GWAS_catalogID_0\n")
    cat(str(array_GWAS_catalogID))
    cat("\n")
    
    
    df_traits<-data.frame()
    
    for(k in 1:length(array_GWAS_catalogID)){
      
      ID_sel<-array_GWAS_catalogID[k]
      
      cat("-------------------------------------------------------------------->\t")
      cat(sprintf(as.character(ID_sel)))
      cat("\n")
      
      study_retrieved<-get_studies(study_id = ID_sel)
      
      # cat("study_retrieved_0\n")
      # cat(str(study_retrieved))
      # cat("\n")
      
      tmp_df<-cbind(ID_sel,study_retrieved@studies$reported_trait)
      
      colnames(tmp_df)<-c("GWAS_catalogID","trait")
      
      cat("tmp_df_0\n")
      cat(str(tmp_df))
      cat("\n")
      
      df_traits<-rbind(tmp_df,df_traits)
      
    }#k in 1:length(array_GWAS_catalogID
    
    Results<-merge(Results,
                   df_traits,
                   by="GWAS_catalogID")
    
    cat("Results_3\n")
    cat(str(Results))
    cat("\n")
    

    Results<-Results[,-which(colnames(Results) == 'input')]
    
    cat("Results_4\n")
    cat(str(Results))
    cat("\n")
    
    ### save ----
    
    setwd(out)
    
    write.table(Results, file='CHIP_variants.tsv', sep="\t",quote=F,row.names=F)
    
  }# length(list_results) >0
  

}


export_config_file = function(option_list)
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
  
  setwd(out)
  
  Results<-read.table(file='CHIP_variants.tsv', sep="\t", header=T)
  
  cat("Results_0\n")
  cat(str(Results))
  cat("\n")
  
  Results$chr<-gsub(":.+$","",Results$snp)
  
  Results$pos<-gsub("^[^:]+:","",Results$snp)
  Results$pos<-as.integer(gsub(":.+$","",Results$pos))
  
  Results$a0<-gsub("^[^:]+:[^:]+:","",Results$snp)
  Results$a0<-gsub(":.+$","",Results$a0)
  
  Results$a1<-gsub("^[^:]+:[^:]+:[^:]+:","",Results$snp)
  
  
  cat("Results_1\n")
  cat(str(Results))
  cat("\n")
  
  config_df<-data.frame()
  
  DEBUG<-0
  
  for(i in 1:dim(Results)[1]){
    
    Results_sel<-Results[i,]
    
    if(DEBUG == 1){
      cat("Results_sel_1\n")
      cat(str(Results_sel))
      cat("\n")
    }
   
    
    a0_nchar<-nchar(Results_sel$a0)
    a1_nchar<-nchar(Results_sel$a1)
    max_char<-max(a0_nchar,a1_nchar)
    
    
    if(DEBUG == 1){
      cat("nchars\n")
      cat(sprintf(as.character(a0_nchar)))
      cat("\t")
      cat(sprintf(as.character(a1_nchar)))
      cat("\n")
      cat("max_char\n")
      cat(sprintf(as.character(max_char)))
      cat("\n")
    }
    
    
    if(max_char > 1){
      
      tmp<-as.data.frame(cbind(Results_sel$snp,Results_sel$rsid, Results_sel$chr, Results_sel$pos,Results_sel$a0,Results_sel$a1, 1, "for_possible_insertion"))
      
      colnames(tmp)<-c("snp",'rsid',"chr","pos","a0","a1","length","tag")
      
      config_df<-rbind(tmp, config_df)
      
      }#max_char > 1
    
    tmp<-as.data.frame(cbind(Results_sel$snp,Results_sel$rsid, Results_sel$chr, Results_sel$pos,Results_sel$a0,Results_sel$a1, max_char, "check_first"))
    
    colnames(tmp)<-c("snp",'rsid',"chr","pos","a0","a1","length","tag")
    
    config_df<-rbind(tmp, config_df)
    
    
  }#i in 1:dim(Results)[1]
  
  config_df$pos<-as.integer(config_df$pos)
  config_df$length<-as.integer(config_df$length)
  
  config_df$length_adjusted<-config_df$length - 1 #### It goes one base away from the position if I don't do this
  
  
  cat("config_df_0\n")
  cat(str(config_df))
  cat("\n")
  
  config_df<-unique(config_df)
  
  cat("config_df_1\n")
  cat(str(config_df))
  cat("\n")
  
  #### GR object ----
  
  gr_VARs <- GRanges(
    seqnames = as.character(config_df$chr),
    strand='+',
    ranges=IRanges(
      start=config_df$pos,
      end=config_df$pos+config_df$length_adjusted,
      names=paste(config_df$snp, config_df$tag, config_df$rsid, sep = "__")))
  
  cat("gr_VARs_0\n")
  cat(str(gr_VARs))
  cat("\n")
  
  
  ##### SAVE ==================
  
  setwd(out)
  
  write.table(config_df, file='config_file.tsv', sep="\t",quote=F,row.names=F)
  
  export.bed(gr_VARs,con=paste('VARS','.bed',sep=''))
  
  
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
    make_option(c("--indir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GWAS_equivalence"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_PP"), type="numeric", default=NULL, 
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
  export_config_file(opt)
  
  
}


###########################################################################

system.time( main() )