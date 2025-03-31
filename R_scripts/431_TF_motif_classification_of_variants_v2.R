
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
  
  #### READ and transform table_sel ----
  
  table_sel = opt$table_sel
  
  cat("table_sel_\n")
  cat(sprintf(as.character(table_sel)))
  cat("\n")
  
  TF_motif_prediction<-as.data.frame(fread(file=opt$TF_motif_prediction, sep="\t", header = F, fill=TRUE), stringsAsFactors = F)
  
  colnames(TF_motif_prediction)<-c('chr','start','end','ensembl_gene_id')
  
  cat("TF_motif_prediction_0\n")
  cat(str(TF_motif_prediction))
  cat("\n")
  
  
  TF_motif_prediction_NO_NA<-TF_motif_prediction[!is.na(TF_motif_prediction$ensembl_gene_id),]
  
  cat("TF_motif_prediction_NO_NA_0\n")
  cat(str(TF_motif_prediction_NO_NA))
  cat("\n")
  
  
  
  
  TF_motif_prediction_NO_NA$Peak_ID<-paste(TF_motif_prediction_NO_NA$ensembl_gene_id,paste(TF_motif_prediction_NO_NA$chr,TF_motif_prediction_NO_NA$start,TF_motif_prediction_NO_NA$end, sep="_"), sep='__')
  
  
  cat("TF_motif_prediction_NO_NA_1\n")
  cat(str(TF_motif_prediction_NO_NA))
  cat("\n")
  
  
  gr_TF_motif_prediction_NO_NA <- GRanges(
    seqnames = as.character(TF_motif_prediction_NO_NA$chr),
    ranges=IRanges(
      start=as.numeric(TF_motif_prediction_NO_NA$start),
      end=as.numeric(TF_motif_prediction_NO_NA$end),
      name=TF_motif_prediction_NO_NA$Peak_ID))
  
  
  #### LiftOver 38 -> 37 ----
  

  motif_in_38_df<-data.frame(chr=as.character(seqnames(gr_TF_motif_prediction_NO_NA)),
                     start=start(gr_TF_motif_prediction_NO_NA),
                     end=end(gr_TF_motif_prediction_NO_NA),
                     Peak_ID=names(gr_TF_motif_prediction_NO_NA),
                     stringsAsFactors = F)
  
  cat("motif_in_38_df_0\n")
  str(motif_in_38_df)
  cat("\n")
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/home/manuel.tardaguila/reference_files/hg38ToHg19.over.chain")
  
  seqlevelsStyle(gr_TF_motif_prediction_NO_NA) = "UCSC"  # necessary
  gr_TF_motif_prediction_NO_NA37 = liftOver(gr_TF_motif_prediction_NO_NA, ch)
  gr_TF_motif_prediction_NO_NA37 = unlist(gr_TF_motif_prediction_NO_NA37)
  genome(gr_TF_motif_prediction_NO_NA37) = "hg19"
  
  if(length(gr_TF_motif_prediction_NO_NA37) >0)
  {
    
    chr_37<-as.character(seqnames(gr_TF_motif_prediction_NO_NA37))
    names_37<-as.character(names(gr_TF_motif_prediction_NO_NA37))
    
    
    
    
    
    
    motif_in_37_df<-data.frame(chr=as.character(seqnames(gr_TF_motif_prediction_NO_NA37)),
                          start_37=start(gr_TF_motif_prediction_NO_NA37),
                          end_37=end(gr_TF_motif_prediction_NO_NA37),
                          Peak_ID=names(gr_TF_motif_prediction_NO_NA37),
                          stringsAsFactors = F)
    

    
    cat("motif_in_37_df_0\n")
    str(motif_in_37_df)
    cat("\n")
    
    
    motif_DEF<-unique(merge(motif_in_38_df,
                             motif_in_37_df,
                             by=c("chr","Peak_ID"),
                            all=TRUE))
    
    
    
    cat("motif_DEF_2\n")
    str(motif_DEF)
    cat("\n")
    # 
    
    
    
    motif_DEF<-merge(motif_DEF,
              TF_motif_prediction_NO_NA,
              by=c("chr","start","end","Peak_ID"),
              all=T)
    
    cat("motif_DEF_3\n")
    str(motif_DEF)
    cat("\n")
    
   
    
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_TF_motif_prediction_NO_NA37) >0
  
  
  
  
  #### READ and transform indir_chip_atlas_search_results ----
  
  indir_chip_atlas_search_results = opt$indir_chip_atlas_search_results
  
  cat("indir_chip_atlas_search_results_\n")
  cat(sprintf(as.character(indir_chip_atlas_search_results)))
  cat("\n")
  
  #### Read in all the results of the meta-analisis ----
  
  file_list <- list.files(path=indir_chip_atlas_search_results, include.dirs = FALSE)
  
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  
  indexes_sel <- grep("\\.out\\.gz",file_list)
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  colnames(file_list_sel)<-"file"
  
 
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$ensembl_gene_id<-gsub("_[^\\.]+\\.out\\.gz","",file_list_sel$file)
  
  cat("file_list_sel_0.25\n")
  cat(str(file_list_sel))
  cat("\n")
  
  
  ############# LOOP --------------------------
  
  List_RESULTS<-list()
 
  
  DEBUG<-0
  
  Result<-data.frame()
  
  for(i in 1:dim(file_list_sel)[1]){
    
    read_file_sel<-file_list_sel$file[i]
    ensembl_gene_id_sel<-file_list_sel$ensembl_gene_id[i]

    
    cat("-------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(read_file_sel)))
    cat("\t")
    cat(sprintf(as.character(ensembl_gene_id_sel)))
    cat("\n")
   
    
    setwd(indir_chip_atlas_search_results)
    
    
    
    SIZE_gate<-file.info(read_file_sel)$size
    
    if(DEBUG ==1)
    {
      cat("SIZE_gate\n")
      cat(str(SIZE_gate))
      cat("\n")
    }
    
    
    
    if(SIZE_gate> 0)
    {
      
      
      LINE_gate<-length(readLines(read_file_sel))
      
      if(DEBUG ==1)
      {
        cat("LINE_gate\n")
        cat(str(LINE_gate))
        cat("\n")
      }
      
      if(LINE_gate> 0)
      {
        
        
        
        df<-as.data.frame(fread(file=read_file_sel, sep="\t", header = F, fill=TRUE), stringsAsFactors = F)
        colnames(df)<-c('chr','start','end','experiment','score')
        
        
        if(DEBUG ==1)
        {
          cat("df_0\n")
          cat(str(df))
          cat("\n")
        }
        
        df$ensembl_gene_id<-ensembl_gene_id_sel
        
        if(DEBUG ==1)
        {
          cat("df_1\n")
          cat(str(df))
          cat("\n")
        }
        
        gr_df <- GRanges(
          seqnames = as.character(df$chr),
          name2=as.character(df$experiment),
          ranges=IRanges(
            start=as.numeric(df$start),
            end=as.numeric(df$end),
            name=df$ensembl_gene_id))
        
        
        
        
        
        motif_DEF_sel<-motif_DEF[which(motif_DEF$ensembl_gene_id == ensembl_gene_id_sel),]
        
        if(DEBUG ==1)
        {
          cat("motif_DEF_sel_0\n")
          cat(str(motif_DEF_sel))
          cat("\n")
        }
        
        gr_motif_DEF_sel <- GRanges(
          seqnames = as.character(motif_DEF_sel$chr),
          start_38= as.integer(motif_DEF_sel$start),
          end_38= as.integer(motif_DEF_sel$end),
          ranges=IRanges(
            start=as.numeric(motif_DEF_sel$start_37),
            end=as.numeric(motif_DEF_sel$end_37),
            name=motif_DEF_sel$ensembl_gene_id))
        
        m <- findOverlaps(gr_motif_DEF_sel,gr_df)
        
        if(DEBUG == 1)
        {
          cat("m\n")
          cat(str(m))
          cat("\n")
        }
        
        subjectHits_df<-subjectHits(m)
        
        if(DEBUG == 1)
        {
          cat("subjectHits_df\n")
          cat(str(subjectHits_df))
          cat("\n")
        }
        
        queryHits_motif_DEF_sel<-queryHits(m)
        
        if(DEBUG == 1)
        {
          cat("queryHits_motif_DEF_sel\n")
          cat(str(queryHits_motif_DEF_sel))
          cat("\n")
        }
        
        if(length(queryHits_motif_DEF_sel) >0){
          
          df_df <- data.frame(chr=as.character(seqnames(gr_df)),
                                                        start=as.integer(start(gr_df)),
                                                        experiment=as.character(gr_df$name2),
                                                        end=as.integer(end(gr_df)),
                                                        ensembl_gene_id=names(gr_df), stringsAsFactors = F)
          
          if(DEBUG == 1)
          {
            cat("df_df_0\n")
            cat(str(df_df))
            cat("\n")
          }
          
          df_df_hits<-df_df[subjectHits_df,]
          
          if(DEBUG == 1)
          {
            cat("df_df_hits_0\n")
            cat(str(df_df_hits))
            cat("\n")
          }
          
          
          
          motif_DEF_sel_df <- data.frame(chr=as.character(seqnames(gr_motif_DEF_sel)),
                                         start_37=as.integer(start(gr_motif_DEF_sel)),
                                         end_37=as.integer(end(gr_motif_DEF_sel)),
                                         start=as.integer(gr_motif_DEF_sel$start_38),
                                         end=as.integer(gr_motif_DEF_sel$end_38),
                                                      ensembl_gene_id=names(gr_motif_DEF_sel), stringsAsFactors = F)
          
          if(DEBUG == 1)
          {
            cat("motif_DEF_sel_df_0\n")
            cat(str(motif_DEF_sel_df))
            cat("\n")
          }
          
          motif_DEF_sel_df_hits<-motif_DEF_sel_df[queryHits_motif_DEF_sel,]
          
          
          if(DEBUG == 1)
          {
            cat("motif_DEF_sel_df_hits_0\n")
            cat(str(motif_DEF_sel_df_hits))
            cat("\n")
          }
          
          motif_DEF_sel_df_hits$chip_seq_experiment<-df_df_hits$experiment
          
          if(DEBUG == 1)
          {
            cat("motif_DEF_sel_df_hits_1\n")
            cat(str(motif_DEF_sel_df_hits))
            cat("\n")
          }
          
          motif_DEF_sel_df_hits.dt<-data.table(motif_DEF_sel_df_hits,
                                                              key=c("chr","start","end","start_37","end_37","ensembl_gene_id"))
          
          motif_DEF_sel_df_hits_collapsed<-as.data.frame(motif_DEF_sel_df_hits.dt[,.(chip_seq_experiment_string=paste(unique(chip_seq_experiment),collapse=';')),
                                                                                                  by=key(motif_DEF_sel_df_hits.dt)], stringsAsFactors=F)
          
          
          if(DEBUG == 1)
          {
            cat("motif_DEF_sel_df_hits_collapsed_0\n")
            cat(str(motif_DEF_sel_df_hits_collapsed))
            cat("\n")
          }
          
          Result<-rbind(motif_DEF_sel_df_hits_collapsed,Result)
          
          
          
          
        }#length(queryHits_motif_DEF_sel) >0
        
      }#LINE_gate
      else{
        
        cat(sprintf("empty_file\n"))
        cat(sprintf(as.character(read_file_sel)))
        cat("\n")
        
        
      }
    }#SIZE_gates
  }#i in 1:dim(file_list_sel)[1]
  
  

  cat("Result_0\n")
  cat(str(Result))
  cat("\n")
  
  
  ####################### SAVE ###############################################################
 
  
  setwd(out)
  
  write.table(Result, 
              file=paste("TF_motifs_supported_by_chipseq_",table_sel,".tsv",sep=""), 
              row.names = F, quote=F, sep="\t")
  
 
}

merge_with_final_table = function(option_list)
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
  
  #### READ and transform out2 ----
  
  out2 = opt$out2
  
  cat("out2_\n")
  cat(sprintf(as.character(out2)))
  cat("\n")
  
  #### READ and transform table_sel ----
  
  table_sel = opt$table_sel
  
  cat("table_sel_\n")
  cat(sprintf(as.character(table_sel)))
  cat("\n")
  
  
  #### Read TF_motifs_with_an_experimental_HIT ----
  
  TF_motifs_with_an_experimental_HIT<-as.data.frame(fread(file=opt$TF_motifs_with_an_experimental_HIT, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  cat("TF_motifs_with_an_experimental_HIT_0\n")
  cat(str(TF_motifs_with_an_experimental_HIT))
  cat("\n")
  
  colnames(TF_motifs_with_an_experimental_HIT)[which(colnames(TF_motifs_with_an_experimental_HIT) == 'start')]<-'motif_start'
  colnames(TF_motifs_with_an_experimental_HIT)[which(colnames(TF_motifs_with_an_experimental_HIT) == 'end')]<-'motif_end'
  colnames(TF_motifs_with_an_experimental_HIT)[which(colnames(TF_motifs_with_an_experimental_HIT) == 'start_37')]<-'motif_start_37'
  colnames(TF_motifs_with_an_experimental_HIT)[which(colnames(TF_motifs_with_an_experimental_HIT) == 'end_37')]<-'motif_end_37'
  
  
  cat("TF_motifs_with_an_experimental_HIT_1\n")
  cat(str(TF_motifs_with_an_experimental_HIT))
  cat("\n")
  
  #### Read Final_table_TF_motif_prediction ----
  
  Final_table_TF_motif_prediction<-as.data.frame(fread(file=opt$Final_table_TF_motif_prediction, sep="\t", header = T, fill=TRUE), stringsAsFactors = F)
  
  cat("Final_table_TF_motif_prediction_0\n")
  cat(str(Final_table_TF_motif_prediction))
  cat("\n")
  
  #### Merge and do the final table----
  
  Final_table_TF_motif_prediction<-merge(Final_table_TF_motif_prediction,
                                              TF_motifs_with_an_experimental_HIT,
                                              by=c('chr','motif_start','motif_end','ensembl_gene_id'),
                                              all.x=T)
  
  cat("Final_table_TF_motif_prediction_1\n")
  cat(str(Final_table_TF_motif_prediction))
  cat("\n")
  
  ############################################################################### SAVE   ###############################################################################
  
  setwd(out2)
  
  write.table(Final_table_TF_motif_prediction, 
              file=paste("Final_table_TF_motif_prediction_with_chipseq_support_",table_sel,".tsv",sep=""), 
              row.names = F, quote=F, sep="\t")
  
  
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
    make_option(c("--TF_motif_prediction"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--indir_chip_atlas_search_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Final_table_TF_motif_prediction"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_motifs_with_an_experimental_HIT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out2"), type="character", default=NULL, 
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
  merge_with_final_table(opt)
    
  
}

###########################################################################

system.time( main() )