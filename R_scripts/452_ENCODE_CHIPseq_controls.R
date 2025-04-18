
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
  
  #### READ and transform TF_to_search ----
  
  TF_to_search = unlist(strsplit(opt$TF_to_search, split=","))
  
  cat("TF_to_search_\n")
  cat(sprintf(as.character(TF_to_search)))
  cat("\n")
  
  #### READ and transform selected_cell_types ----
  
  selected_cell_types = unlist(strsplit(opt$selected_cell_types, split=","))
  
  cat("selected_cell_types_\n")
  cat(sprintf(as.character(selected_cell_types)))
  cat("\n")
  
  #### Read the ENCODE file -----
  
  ENCODE_file<-as.data.frame(fread(file=opt$ENCODE_file, sep="\t", header = T))
  
  
  cat("ENCODE_file_0\n")
  cat(str(ENCODE_file))
  cat("\n")
  
  ENCODE_file$Peak_ID<-paste(ENCODE_file$cell_type,paste(paste(ENCODE_file$chrom,
                                                                                     ENCODE_file$chromStart,
                                                                                     ENCODE_file$chromEnd,
                                                                                     sep="_"),ENCODE_file$CHIP_factor,sep=";"),sep="|")
  
  cat("ENCODE_file_1\n")
  cat(str(ENCODE_file))
  cat("\n")
  
  
  ENCODE_file_sel<-ENCODE_file[which(ENCODE_file$CHIP_factor%in%TF_to_search),]
  
  cat("ENCODE_file_sel_0\n")
  cat(str(ENCODE_file_sel))
  cat("\n")
  
  if(dim(ENCODE_file_sel)[1] >0){
    
    ENCODE_file_sel_CT_sel<-ENCODE_file_sel[which(ENCODE_file_sel$cell_type%in%selected_cell_types),]
    
   
    
    if(dim(ENCODE_file_sel_CT_sel)[1] >0){
      
      
      cat("ENCODE_file_sel_CT_sel_0\n")
      cat(str(ENCODE_file_sel_CT_sel))
      cat("\n")
      
      
      
      # cat("Found CHIPseq evidence lifting it over to GRCh38\n")
      
      gr_ENCODE <- GRanges(
        seqnames = as.character(ENCODE_file_sel_CT_sel$chrom),
        ranges=IRanges(
          start=as.numeric(ENCODE_file_sel_CT_sel$chromStart),
          end=as.numeric(ENCODE_file_sel_CT_sel$chromEnd),
          name=as.character(ENCODE_file_sel_CT_sel$Peak_ID)))
      
      
      #### LiftOver ENCODE data to GRCh38----
      
     
      
      ENCODE_file_sel_CT_sel_df<-data.frame(
                          chr=as.character(seqnames(gr_ENCODE)),
                          Peak_ID=names(gr_ENCODE),
                          Peak_start_37=start(gr_ENCODE),
                          Peak_end_37=end(gr_ENCODE),
                          stringsAsFactors = F)
      
      cat("ENCODE_file_sel_CT_sel_df_0\n")
      str(ENCODE_file_sel_CT_sel_df)
      cat("\n")
      
      #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
      ch = import.chain("/home/manuel.tardaguila/reference_files/hg19ToHg38.over.chain")
      
      seqlevelsStyle(gr_ENCODE) = "UCSC"  # necessary
      gr_ENCODE38 = liftOver(gr_ENCODE, ch)
      gr_ENCODE38 = unlist(gr_ENCODE38)
      genome(gr_ENCODE38) = "hg38"
      
      if(length(gr_ENCODE38) >0)
      {
        
        ENCODE_file_sel_CT_sel_38_df<-data.frame(
                               chr=as.character(seqnames(gr_ENCODE38)),
                               Peak_start_38=start(gr_ENCODE38),
                               Peak_end_38=end(gr_ENCODE38),
                               Peak_ID=names(gr_ENCODE38),
                               stringsAsFactors = F)
        
        
        cat("ENCODE_file_sel_CT_sel_38_df_0\n")
        str(ENCODE_file_sel_CT_sel_38_df)
        cat("\n")
        
        ENCODE_file_sel_CT_sel_df<-merge(ENCODE_file_sel_CT_sel_38_df,
                                         ENCODE_file_sel_CT_sel_df,
                                         by=c("chr","Peak_ID"))
        
        cat("ENCODE_file_sel_CT_sel_df_1\n")
        str(ENCODE_file_sel_CT_sel_df)
        cat("\n")
        
        
        ENCODE_file_sel_CT_sel_df<-merge(ENCODE_file_sel_CT_sel_df,
                                         ENCODE_file,
                                         by=c("Peak_ID"))
        
        cat("ENCODE_file_sel_CT_sel_df_2\n")
        str(ENCODE_file_sel_CT_sel_df)
        cat("\n")
        
        
        test<- ENCODE_file_sel_CT_sel_df %>%
          group_by(cell_type,CHIP_factor) %>%
          top_n(n = 10, wt = score)
        
        

        cat("test_0\n")
        str(test)
        cat("\n")
        
        
        
        
        ENCODE_file_sel_CT_sel_df_top10 <- as.data.frame(test, stringsAsFactors=F)
        
        cat("ENCODE_file_sel_CT_sel_df_top10_0\n")
        str(ENCODE_file_sel_CT_sel_df_top10)
        cat("\n")      
        
        
        # ENCODE_file_sel_CT_sel_df <- ENCODE_file_sel_CT_sel_df[order(ENCODE_file_sel_CT_sel_df$score, decreasing = TRUE), ]  # Top N highest values by group
        # 
        # cat("ENCODE_file_sel_CT_sel_df_3\n")
        # str(ENCODE_file_sel_CT_sel_df)
        # cat("\n")
        # 
        # ENCODE_file_sel_CT_sel_df.dt <- data.table(ENCODE_file_sel_CT_sel_df, key = c("cell_type","CHIP_factor"))
        # 
        # ENCODE_file_sel_CT_sel_df_top10 <- as.data.frame(ENCODE_file_sel_CT_sel_df.dt[ , head(.SD, 10), by = key(ENCODE_file_sel_CT_sel_df.dt)], stringsAsFactors=F)
        # 
        # cat("ENCODE_file_sel_CT_sel_df_top10_0\n")
        # str(ENCODE_file_sel_CT_sel_df_top10)
        # cat("\n")        
       
        
        setwd(out)
        
        write.table(ENCODE_file_sel_CT_sel_df_top10, 
                    file=paste("ENCODE_top10_controls_",paste(TF_to_search,collapse="_" ),".tsv",sep=""), 
                    row.names = F, quote=F, sep="\t")
        
        #### GR object ----
        
        gr_controls <- GRanges(
          seqnames = as.character(michelas_table$chr),
          strand='+',
          ranges=IRanges(
            start=michelas_table$Peak_start_38,
            end=michelas_table$Peak_end_38,
            names=michelas_table$Peak_ID))
        
        cat("gr_controls_0\n")
        cat(str(gr_controls))
        cat("\n")
        
        #### export bed ----
        
        setwd(out)
        
        export.bed(gr_controls,con=paste('controls','.bed',sep=''))
        
      }else{
        
        
        stop("NO_LIFT_OVER\n")
        
      }# length(gr_ENCODE38) >0
      
      
    }#dim(ENCODE_file_sel_CT_sel)[1] >0
    
  }#dim(ENCODE_file_sel)[1] >0
 
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
    make_option(c("--selected_cell_types"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ENCODE_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_to_search"), type="character", default=NULL, 
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