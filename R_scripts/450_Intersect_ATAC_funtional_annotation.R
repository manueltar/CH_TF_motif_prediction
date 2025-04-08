
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


opt = NULL

options(warn = 1)

intersect_with_ALL = function(option_list)
{
  suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  library("R.oo", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/")
  library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/")
  suppressMessages(library("liftOver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("rCNV", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  
  
  
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
  
  #### READ ATAC_RAW ----
  
  
  
  
  gr_VARS_37 <- GRanges(
    seqnames = VAR_DEF$chr,
    ranges=IRanges(
      start=as.numeric(VAR_DEF$pos_37),
      end=as.numeric(VAR_DEF$pos_37),
      name=VAR_DEF$VAR_38))
  
  cat("gr_VARS_37_0\n")
  cat(str(gr_VARS_37))
  cat("\n")
  
  

  
  #### READ ATAC_peaks ----
  
  
  ATAC_peaks<-as.data.frame(fread(file=opt$ATAC_peaks, sep="\t", header =F), stringAsFactors=F)
  
  colnames(ATAC_peaks)<-c('chr','start','end')
  
  
  cat("ATAC_peaks_0\n")
  cat(str(ATAC_peaks))
  cat("\n")
  
  ATAC_counts<-as.data.frame(fread(file=opt$ATAC_counts, sep="\t", header =T), stringAsFactors=F)
  
  colnames(ATAC_counts)[which(colnames(ATAC_counts) == 'GMP-A')]<-'GMP_A'
  colnames(ATAC_counts)[which(colnames(ATAC_counts) == 'GMP-B')]<-'GMP_B'
  colnames(ATAC_counts)[which(colnames(ATAC_counts) == 'GMP-C')]<-'GMP_C'
  
  
  cat("ATAC_counts_0\n")
  cat(str(ATAC_counts))
  cat("\n")

  
  ATAC_DEF<-cbind(ATAC_peaks,ATAC_counts)
  
  cat("ATAC_DEF_0\n")
  cat(str(ATAC_DEF))
  cat("\n")
  
           
  
  gr_ATAC_DEF <- GRanges(
    seqnames = ATAC_DEF$chr,
    name2=as.character(ATAC_DEF$B),
    name3=as.character(ATAC_DEF$CD4),
    name4=as.character(ATAC_DEF$CD8),
    name5=as.character(ATAC_DEF$CLP),
    name6=as.character(ATAC_DEF$CMP),
    name7=as.character(ATAC_DEF$Ery),
    name8=as.character(ATAC_DEF$GMP_A),
    name9=as.character(ATAC_DEF$GMP_B),
    name10=as.character(ATAC_DEF$GMP_C),
    name11=as.character(ATAC_DEF$HSC),
    name12=as.character(ATAC_DEF$LMPP),
    name13=as.character(ATAC_DEF$mDC),
    name14=as.character(ATAC_DEF$Mega),
    name15=as.character(ATAC_DEF$MEP),
    name16=as.character(ATAC_DEF$Mono),
    name17=as.character(ATAC_DEF$MPP),
    name18=as.character(ATAC_DEF$NK),
    name19=as.character(ATAC_DEF$pDC),
    ranges=IRanges(
      start=ATAC_DEF$start,
      end=ATAC_DEF$end))
  
  
  #### find overlap with SNP ----
  
  DEBUG <- 1
  
  m <- findOverlaps(gr_VARS_37,gr_ATAC_DEF)
  
  if(DEBUG == 1)
  {
    cat("m\n")
    cat(str(m))
    cat("\n")
  }
  
  subjectHits_ATAC_DEF<-subjectHits(m)
  
  if(DEBUG == 1)
  {
    cat("subjectHits_ATAC_DEF\n")
    cat(str(subjectHits_ATAC_DEF))
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
  
  ATAC_DEF_df <- data.frame(chr=as.character(seqnames(gr_ATAC_DEF)),
                            Peak_start=as.integer(start(gr_ATAC_DEF)),
                            Peak_end=as.integer(end(gr_ATAC_DEF)),
                            B=as.numeric(gr_ATAC_DEF$name2),
                            CD4=as.numeric(gr_ATAC_DEF$name3),
                            CD8=as.numeric(gr_ATAC_DEF$name4),
                            CLP=as.numeric(gr_ATAC_DEF$name5),
                            CMP=as.numeric(gr_ATAC_DEF$name6),
                            Ery=as.numeric(gr_ATAC_DEF$name7),
                            GMP_A=as.numeric(gr_ATAC_DEF$name8),
                            GMP_B=as.numeric(gr_ATAC_DEF$name9),
                            GMP_C=as.numeric(gr_ATAC_DEF$name10),
                            HSC=as.numeric(gr_ATAC_DEF$name11),
                            LMPP=as.numeric(gr_ATAC_DEF$name12),
                            mDC=as.numeric(gr_ATAC_DEF$name13),
                            Mega=as.numeric(gr_ATAC_DEF$name14),
                            MEP=as.numeric(gr_ATAC_DEF$name15),
                            Mono=as.numeric(gr_ATAC_DEF$name16),
                            MPP=as.numeric(gr_ATAC_DEF$name17),
                            NK=as.numeric(gr_ATAC_DEF$name18),
                            pDC=as.numeric(gr_ATAC_DEF$name19),
                            stringsAsFactors = F)
  
  
  if(DEBUG == 1)
  {
    cat("ATAC_DEF_df_0\n")
    cat(str(ATAC_DEF_df))
    cat("\n")
  }
  
  
  
  ATAC_DEF_df_hits<-ATAC_DEF_df[subjectHits_ATAC_DEF,]
  
  if(dim(ATAC_DEF_df_hits)[1] >0)
  {
    

    if(DEBUG == 1)
    {
      cat("ATAC_DEF_df_hits_0\n")
      cat(str(ATAC_DEF_df_hits))
      cat("\n")
    }
    
    Final_df<-cbind(VAR_df_hits,
                    ATAC_DEF_df_hits)
    
    
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
    
    Final_df.m<-melt(Final_df, id.vars=c("pos_37","VAR_38","chr","Peak_start","Peak_end"), variable.name= "cell_type", value.name="ATACseq counts")
    
    
    if(DEBUG == 1)
    {
      cat("Final_df.m_0\n")
      cat(str(Final_df.m))
      cat("\n")
    }
    
    Final_df.m<-merge(Final_df.m,
                      Table_of_variants,
                      by="VAR_38")
    
    if(DEBUG == 1)
    {
      cat("Final_df.m_1\n")
      cat(str(Final_df.m))
      cat("\n")
    }
  
    setwd(out)  
      
      
    write.table(Final_df.m, file=paste('ATAC_overlap',".tsv",sep=''), sep="\t", quote=F, row.names = F)
  
 
  }#dim(ATAC_DEF_df_hits)[1] >0
  

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
    make_option(c("--ATAC_peaks"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_counts"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_of_variants"), type="character", default=NULL, 
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
  
  intersect_with_ALL(opt)

  
}


###########################################################################

system.time( main() )