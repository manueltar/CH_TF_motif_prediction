
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("reshape2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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
suppressMessages(library("vroom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
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
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("DESeq2", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("R.utils", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("liftOver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rCNV", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)


data_wrangling_and_intersect = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ CHIP_variants_file ----
  
  CHIP_variants_file<-as.data.frame(fread(file=opt$CHIP_variants_file, sep="\t", header=T), stringsAsFactors=F)
  
  cat("CHIP_variants_file_0\n")
  cat(str(CHIP_variants_file))
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
  
  
  gr_VARS_37 <- GRanges(
    seqnames = VAR_DEF$chr,
    ranges=IRanges(
      start=as.numeric(VAR_DEF$pos_37),
      end=as.numeric(VAR_DEF$pos_37),
      name=VAR_DEF$VAR_38))
  
  cat("gr_VARS_37_0\n")
  cat(str(gr_VARS_37))
  cat("\n")
  
 
  
  #### Read the ENCODE file -----
  
  ENCODE_file<-as.data.frame(fread(file=opt$ENCODE_file, sep="\t", header = T))
  
  
  cat("ENCODE_file_0\n")
  cat(str(ENCODE_file))
  cat("\n")
  
  #### Select the Dnase peaks -----
  
  ENCODE_DNase<-ENCODE_file[grep('Dnase', ENCODE_file$file),]
  
  cat("ENCODE_DNase_0\n")
  cat(str(ENCODE_DNase))
  cat("\n")
  
  gr_ENCODE_DNase <- GRanges(
    seqnames = ENCODE_DNase$chrom,
    name2=as.character(ENCODE_DNase$cell_type),
    name3=as.character(ENCODE_DNase$score),
    ranges=IRanges(
      start=ENCODE_DNase$chromStart,
      end=ENCODE_DNase$chromEnd))
  
  
  #### find overlap with SNP ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_VARS_37,gr_ENCODE_DNase)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_ENCODE_DNase<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_ENCODE_DNase\n")
    cat(str(subjectHits_ENCODE_DNase))
    cat("\n")
  }
  
  queryHits_VAR<-queryHits(m)
  
  if(DEBUG == 1)
  {
    cat("queryHits_VAR\n")
    cat(str(queryHits_VAR))
    cat("\n")
  }
  
  VAR_df <- data.frame(chr=as.character(seqnames(gr_VARS_37)),
                       pos_37=as.integer(end(gr_VARS_37)),
                       VAR_38=names(gr_VARS_37), stringsAsFactors = F)
  
  if(DEBUG == 1)
  {
    cat("VAR_df_0\n")
    cat(str(VAR_df))
    cat("\n")
  }
  
  VAR_df_hits<-VAR_df[queryHits_VAR,]
  
  if(DEBUG == 1)
  {
    cat("VAR_df_hits_0\n")
    cat(str(VAR_df_hits))
    cat("\n")
  }
  
  ENCODE_DNase_df <- data.frame(chr=as.character(seqnames(gr_ENCODE_DNase)),
                            Peak_start=as.integer(start(gr_ENCODE_DNase)),
                            Peak_end=as.integer(end(gr_ENCODE_DNase)),
                            cell_type=as.character(gr_ENCODE_DNase$name2),
                            score=as.numeric(gr_ENCODE_DNase$name3),
                            stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("ENCODE_DNase_df_0\n")
    cat(str(ENCODE_DNase_df))
    cat("\n")
  }
  
  
  
  ENCODE_DNase_df_hits<-ENCODE_DNase_df[subjectHits_ENCODE_DNase,]
  
  if(dim(ENCODE_DNase_df_hits)[1] >0)
  {
    
    
    if(DEBUG == 1)
    {
      cat("ENCODE_DNase_df_hits_0\n")
      cat(str(ENCODE_DNase_df_hits))
      cat("\n")
    }
    
    Final_df<-cbind(VAR_df_hits,
                    ENCODE_DNase_df_hits)
    
    colnames(Final_df)[which(colnames(Final_df) == 'score')]<-'DNase score'
    
    
    if(DEBUG == 1)
    {
      cat("Final_df_0\n")
      cat(str(Final_df))
      cat("\n")
    }
    
    Final_df<-Final_df[,-1]
    
    
    if(DEBUG == 1)
    {
      cat("Final_df_1\n")
      cat(str(Final_df))
      cat("\n")
    }
    
   
    
    Final_df<-merge(Final_df,
                      Table_of_variants,
                      by=c("chr","VAR_38"))
    
    if(DEBUG == 1)
    {
      cat("Final_df_1\n")
      cat(str(Final_df))
      cat("\n")
    }
    
    Final_df_subset<-Final_df[,which(colnames(Final_df)%in%c('snp','VAR_38','Peak_start','Peak_end',"cell_type","DNase score"))]
    
    if(DEBUG == 1)
    {
      cat("Final_df_subset_0\n")
      cat(str(Final_df_subset))
      cat("\n")
    }
    
    Final_df_subset<-merge(Final_df_subset,
                             CHIP_variants_file,
                             by="snp")
    if(DEBUG == 1)
    {
      cat("Final_df_subset_1\n")
      cat(str(Final_df_subset))
      cat("\n")
    }
    
    setwd(out)  
    
    
    write.table(Final_df_subset, file=paste('DNase_overlap',".tsv",sep=''), sep="\t", quote=F, row.names = F)
    
    
  }#dim(ENCODE_DNase_df_hits)[1] >0
  
  
  
  
  
  
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
    make_option(c("--ENCODE_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_of_variants"), type="character", default=NULL, 
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
  
  data_wrangling_and_intersect(opt)

  
}

###########################################################################

system.time( main() )