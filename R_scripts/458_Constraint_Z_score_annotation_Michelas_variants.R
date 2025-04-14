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

intersect_Constraint_Z = function(option_list)
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
  
  #### Read Constraint_Z ----
  
  
  Constraint_Z<-as.data.frame(fread(file=opt$Constraint_Z,sep="\t",header=T, stringsAsFactors = F))
  
  cat("Constraint_Z\n")
  str(Constraint_Z)
  cat("\n")
  
  gr_Constraint_Z <-unique(GRanges(
    seqnames = as.character(Constraint_Z$chrom),
    possible= as.numeric(Constraint_Z$possible),
    expected= as.numeric(Constraint_Z$expected),
    observed= as.numeric(Constraint_Z$observed),
    oe= as.numeric(Constraint_Z$oe),
    z= as.numeric(Constraint_Z$z),
    ranges=IRanges(
      start=as.numeric(Constraint_Z$start),
      end=as.numeric(Constraint_Z$end),
      names = Constraint_Z$element_id)))
  
 
  #### Read the Table_of_variants file -----
  
  Table_of_variants<-as.data.frame(fread(file=opt$Table_of_variants, sep="\t", header = T), stringsAsFactors = F) # /group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv
  
  
  cat("Table_of_variants_0\n")
  cat(str(Table_of_variants))
  cat("\n")
  
  Table_of_variants$VAR_38<-paste(Table_of_variants$chr,Table_of_variants$pos,Table_of_variants$ref,Table_of_variants$alt,sep='_')
  
  cat("Table_of_variants_1\n")
  cat(str(Table_of_variants))
  cat("\n")
  
  gr_VARS <- GRanges(
    seqnames = as.character(Table_of_variants$chr),
    ranges=IRanges(
      start=as.numeric(Table_of_variants$pos),
      end=as.numeric(Table_of_variants$pos),
      name=Table_of_variants$VAR_38))
  

 
 
  #### Intersect TSS to genes with Regulatory features ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_VARS,gr_Constraint_Z,
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
  
  VARS_df <- data.frame(chr=as.character(seqnames(gr_VARS)),
                                                  pos=as.integer(end(gr_VARS)),
                                                  VAR_38=names(gr_VARS), stringsAsFactors = F)
  
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
  
  
  Constraint_Z_df_VARS <- data.frame(chr=as.character(seqnames(gr_Constraint_Z)),
                                                      region_start=as.integer(start(gr_Constraint_Z)),
                                                      region_end=as.integer(end(gr_Constraint_Z)),
                                                       possible=as.numeric(gr_Constraint_Z$possible),
                                                       expected=as.numeric(gr_Constraint_Z$expected),
                                                       observed=as.numeric(gr_Constraint_Z$observed),
                                                       oe=as.numeric(gr_Constraint_Z$oe),
                                                       z=as.numeric(gr_Constraint_Z$z),
                                                       element_id=names(gr_Constraint_Z),
                                                      stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("Constraint_Z_df_VARS_0\n")
    cat(str(Constraint_Z_df_VARS))
    cat("\n")
  }
  
  Constraint_Z_df_VARS_hits<-Constraint_Z_df_VARS[subjectHits_AbC_selected_thresholded,]
  
  if(dim(Constraint_Z_df_VARS_hits)[1] >0)
  {
    if(DEBUG == 1)
    {
      cat("Constraint_Z_df_VARS_hits_0\n")
      cat(str(Constraint_Z_df_VARS_hits))
      cat("\n")
    }
    
    Constraint_Z_df_VARS_hits<-cbind(Constraint_Z_df_VARS_hits,VARS_df_hits[,-1])
    
    if(DEBUG == 1)
    {
      cat("Constraint_Z_df_VARS_hits_1\n")
      cat(str(Constraint_Z_df_VARS_hits))
      cat("\n")
      cat(str(unique(Constraint_Z_df_VARS_hits$VAR_38)))
      cat("\n")
    }
   
    Constraint_Z_df_VARS_hits<-merge(unique(Table_of_variants[,which(colnames(Table_of_variants)%in%c('VAR_38','snp','rsid'))]),
                                   Constraint_Z_df_VARS_hits,
                                   by=c("VAR_38"))
    
    if(DEBUG == 1)
    {
      cat("Constraint_Z_df_VARS_hits_2\n")
      cat(str(Constraint_Z_df_VARS_hits))
      cat("\n")
      cat(str(unique(Constraint_Z_df_VARS_hits$VAR_38)))
      cat("\n")
    }
   
    
    # #### SAVE ----
    
    setwd(out)
    
    write.table(Constraint_Z_df_VARS_hits, file="Constraint_Z_hits.tsv", sep="\t", quote = F, row.names = F)
    

    
    
    
  }#dim(Constraint_Z_df_VARS_hits)[1] >0
  
  
 
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
    make_option(c("--Constraint_Z"), type="character", default=NULL,
                metavar="FILE.txt",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_of_variants"), type="character", default=NULL,
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
  
 
 intersect_Constraint_Z(opt)
  
}




###########################################################################

system.time( main() )
