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
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("reshape2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicFeatures", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("OrganismDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Homo.sapiens", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("gwascat", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("R.utils", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))



opt = NULL

filter_AbC = function(option_list)
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
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform AbC_cell_types ----
  
  AbC_cell_types = unlist(strsplit(opt$AbC_cell_types, split=","))
  
  cat("AbC_cell_types_\n")
  cat(sprintf(as.character(AbC_cell_types)))
  cat("\n")
  
  #### READ and transform Threshold_AbC ----
  
  Threshold_AbC = opt$Threshold_AbC
  
  cat("Threshold_AbC\n")
  cat(sprintf(as.character(Threshold_AbC)))
  cat("\n")
  
 
  
  #### Read AbC_scores ----
  
  
  AbC_scores<-as.data.frame(fread(file=opt$AbC_scores,sep="\t",header=T, stringsAsFactors = F))
  
  cat("AbC_scores_0\n")
  str(AbC_scores)
  cat("\n")
  
  
  AbC_scores_selected_cell_types<-AbC_scores[which(AbC_scores$CellType%in%AbC_cell_types),]
  
  cat("AbC_scores_selected_cell_types_0\n")
  str(AbC_scores_selected_cell_types)
  cat("\n")
  
  
  AbC_scores_selected_cell_types_Thresholded<-AbC_scores_selected_cell_types[which(AbC_scores_selected_cell_types$ABC.Score >= Threshold_AbC),]
  
  cat("AbC_scores_selected_cell_types_Thresholded_0\n")
  str(AbC_scores_selected_cell_types_Thresholded)
  cat("\n")
  
  
  ##### save -----
  
  
  setwd(out)
  
  saveRDS(AbC_scores_selected_cell_types_Thresholded, file='AbC_selected_thresholded.rds')
  
  write.table(AbC_scores_selected_cell_types_Thresholded, file='AbC_selected_thresholded.tsv', sep="\t",quote=F,row.names = F)
 
}

intersect_filtered_AbC = function(option_list)
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
  
  cat("OUT\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
 
  #### Read the Table_of_variants file -----
  
  Table_of_variants<-as.data.frame(fread(file=opt$Table_of_variants, sep="\t", header = T), stringsAsFactors = F) # /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv
  
  
  cat("Table_of_variants_0\n")
  cat(str(Table_of_variants))
  cat("\n")
  
  Table_of_variants$VAR_38<-paste(Table_of_variants$chr,Table_of_variants$pos,Table_of_variants$ref,Table_of_variants$alt,sep='_')
  
  cat("Table_of_variants_1\n")
  cat(str(Table_of_variants))
  cat("\n")
  
  gr_Table_of_variants <- GRanges(
    seqnames = as.character(Table_of_variants$chr),
    ranges=IRanges(
      start=as.numeric(Table_of_variants$pos),
      end=as.numeric(Table_of_variants$pos),
      name=Table_of_variants$VAR_38))
  
  #### Read AbC_selected_thresholded ----
  
  setwd(out)
  
  AbC_selected_thresholded<-readRDS(file='AbC_selected_thresholded.rds')
  
  cat("AbC_selected_thresholded_0\n")
  str(AbC_selected_thresholded)
  cat("\n")
  
  gr_AbC_selected_thresholded <-GRanges(
    seqnames = as.character(AbC_selected_thresholded$chr),
    name2=AbC_selected_thresholded$TargetGene,
    name3=AbC_selected_thresholded$CellType,
    name4=AbC_selected_thresholded$ABC.Score,
    ranges=IRanges(
      start=as.numeric(AbC_selected_thresholded$start),
      end=as.numeric(AbC_selected_thresholded$end),
      names = AbC_selected_thresholded$name))
  
  
  cat("gr_AbC_selected_thresholded_0\n")
  str(gr_AbC_selected_thresholded)
  cat("\n")
  
  
  
  
  #### LiftOver 38 -> 37 VARs ----
  
  
  VAR_in_38_df<-data.frame(chr=as.character(seqnames(gr_Table_of_variants)),
                           pos=end(gr_Table_of_variants),
                           VAR_38=names(gr_Table_of_variants),
                           stringsAsFactors = F)
  
  cat("VAR_in_38_df_0\n")
  str(VAR_in_38_df)
  cat("\n")
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/home/manuel.tardaguila/reference_files/hg38ToHg19.over.chain")
  
  seqlevelsStyle(gr_Table_of_variants) = "UCSC"  # necessary
  gr_Table_of_variants37 = liftOver(gr_Table_of_variants, ch)
  gr_Table_of_variants37 = unlist(gr_Table_of_variants37)
  genome(gr_Table_of_variants37) = "hg19"
  
  if(length(gr_Table_of_variants37) >0)
  {
    
    chr_37<-as.character(seqnames(gr_Table_of_variants37))
    names_37<-as.character(names(gr_Table_of_variants37))
    
    
    
    
    
    
    VAR_in_37_df<-data.frame(chr=as.character(seqnames(gr_Table_of_variants37)),
                             pos_37=end(gr_Table_of_variants37),
                             VAR_38=names(gr_Table_of_variants37),
                             stringsAsFactors = F)
    
    
    
    cat("VAR_in_37_df_0\n")
    str(VAR_in_37_df)
    cat("\n")
    
    
    VAR_DEF<-unique(merge(VAR_in_38_df,
                          VAR_in_37_df,
                          by=c("chr","VAR_38"),
                          all=TRUE))
    
    
    
    cat("VAR_DEF_2\n")
    str(VAR_DEF)
    cat("\n")
    # 
    
    
    
    VAR_DEF<-merge(VAR_DEF,
                   Table_of_variants,
                   by=c("chr","pos","VAR_38"),
                   all=T)
    
    cat("VAR_DEF_3\n")
    str(VAR_DEF)
    cat("\n")
    
    
    
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_Table_of_variants37) >0
  
  #### gr_VARS_37 ----
  
  
  
  
  gr_VARS_37 <- GRanges(
    seqnames = VAR_DEF$chr,
    ranges=IRanges(
      start=as.numeric(VAR_DEF$pos_37),
      end=as.numeric(VAR_DEF$pos_37),
      name=VAR_DEF$VAR_38))
  
  cat("gr_VARS_37_0\n")
  cat(str(gr_VARS_37))
  cat("\n")
  
  
  

  
  #### Intersect TSS to genes with Regulatory features ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_VARS_37,gr_AbC_selected_thresholded,
                    ignore.strand = TRUE)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_AbC_selected_thresholded<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_AbC_selected_thresholded\n")
    cat(str(subjectHits_AbC_selected_thresholded))
    cat("\n")
  }
  
  queryHits_VARS<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_VARS\n")
    cat(str(queryHits_VARS))
    cat("\n")
  }
  
  VARS_df <- data.frame(chr=as.character(seqnames(gr_VARS_37)),
                                                  pos_37=as.integer(end(gr_VARS_37)),
                                                  VAR_38=names(gr_VARS_37), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("VARS_df_0\n")
    cat(str(VARS_df))
    cat("\n")
  }
  
  VARS_df_hits<-VARS_df[queryHits_VARS,]
  
  if(DEBUG == 1)
  {
    cat("VARS_df_hits_0\n")
    cat(str(VARS_df_hits))
    cat("\n")
  }
 
  
  AbC_scores_df_VARS <- data.frame(chr=as.character(seqnames(gr_AbC_selected_thresholded)),
                                                      AbC_start_37=as.integer(start(gr_AbC_selected_thresholded)),
                                                      AbC_end_37=as.integer(end(gr_AbC_selected_thresholded)),
                                                      Symbol=as.character(gr_AbC_selected_thresholded$name2),
                                                      Cell_Type=as.character(gr_AbC_selected_thresholded$name3),
                                                      ABC_score=as.numeric(gr_AbC_selected_thresholded$name4),
                                                      ABC_name=names(gr_AbC_selected_thresholded),
                                                      stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("AbC_scores_df_VARS_0\n")
    cat(str(AbC_scores_df_VARS))
    cat("\n")
  }
  
  AbC_scores_df_VARS_hits<-AbC_scores_df_VARS[subjectHits_AbC_selected_thresholded,]
  
  if(dim(AbC_scores_df_VARS_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("AbC_scores_df_VARS_hits_0\n")
      cat(str(AbC_scores_df_VARS_hits))
      cat("\n")
    }
    
    AbC_scores_df_VARS_hits<-cbind(AbC_scores_df_VARS_hits,VARS_df_hits[,-1])
    
    if(DEBUG == 1)
    {
      cat("AbC_scores_df_VARS_hits_1\n")
      cat(str(AbC_scores_df_VARS_hits))
      cat("\n")
      cat(str(unique(AbC_scores_df_VARS_hits$VAR_38)))
      cat("\n")
    }
   
    AbC_scores_df_VARS_hits<-merge(unique(Table_of_variants[,which(colnames(Table_of_variants)%in%c('VAR_38','snp','rsid'))]),
                                   AbC_scores_df_VARS_hits,
                                   by=c("VAR_38"))
    
    if(DEBUG == 1)
    {
      cat("AbC_scores_df_VARS_hits_2\n")
      cat(str(AbC_scores_df_VARS_hits))
      cat("\n")
      cat(str(unique(AbC_scores_df_VARS_hits$VAR_38)))
      cat("\n")
    }
   
    
    # #### SAVE ----
    
    setwd(out)
    
    write.table(AbC_scores_df_VARS_hits, file="AbC_hits.tsv", sep="\t", quote = F, row.names = F)
    

    
    
    
  }#dim(AbC_scores_df_VARS_hits)[1] >0
  
  
 
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
  traceback()
  options(show.error.locations = TRUE)
  
  option_list <- list(
    make_option(c("--AbC_scores"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_of_variants"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--AbC_cell_types"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_AbC"), type="numeric", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
       make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
        make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
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
  
 
 filter_AbC(opt)
 intersect_filtered_AbC(opt)
  
}




###########################################################################

system.time( main() )
