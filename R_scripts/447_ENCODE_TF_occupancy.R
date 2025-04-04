
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
  
  #### Read the Final_table_TF_motif_prediction file -----
  
  Final_table_TF_motif_prediction<-as.data.frame(fread(file=opt$Final_table_TF_motif_prediction, sep="\t", header = T), stringsAsFactors = F)
  

  cat("Final_table_TF_motif_prediction_0\n")
  cat(str(Final_table_TF_motif_prediction))
  cat("\n")
  
  #### Read the ENCODE file -----
  
  ENCODE_file<-as.data.frame(fread(file=opt$ENCODE_file, sep="\t", header = T))
  
  
  cat("ENCODE_file_0\n")
  cat(str(ENCODE_file))
  cat("\n")
  
  #### Discard genes without ENSG ------------------------------------------
  
  
  Final_table_TF_motif_prediction_NO_NA<-Final_table_TF_motif_prediction[!is.na(Final_table_TF_motif_prediction$ensembl_gene_id),]
  
  cat("Final_table_TF_motif_prediction_NO_NA_0\n")
  cat(str(Final_table_TF_motif_prediction_NO_NA))
  cat("\n")
  
  
  Final_table_TF_motif_prediction_NO_NA_sub<-Final_table_TF_motif_prediction_NO_NA[,-which(colnames(Final_table_TF_motif_prediction_NO_NA)%in%c('motif_start_37','motif_end_37'))]
  
  
  cat("Final_table_TF_motif_prediction_NO_NA_sub_0\n")
  cat(str(Final_table_TF_motif_prediction_NO_NA_sub))
  cat("\n")
  
  
  
  Final_table_TF_motif_prediction_NO_NA_sub$Peak_ID<-paste(Final_table_TF_motif_prediction_NO_NA_sub$ensembl_gene_id,
                                                       paste(Final_table_TF_motif_prediction_NO_NA_sub$chr,
                                                             Final_table_TF_motif_prediction_NO_NA_sub$motif_start,
                                                             Final_table_TF_motif_prediction_NO_NA_sub$motif_end, sep="_"), sep='__')
  
  
  cat("Final_table_TF_motif_prediction_NO_NA_sub_1\n")
  cat(str(Final_table_TF_motif_prediction_NO_NA_sub))
  cat("\n")
  
  
  gr_Final_table_TF_motif_prediction_NO_NA_sub <- GRanges(
    seqnames = as.character(Final_table_TF_motif_prediction_NO_NA_sub$chr),
    ranges=IRanges(
      start=as.numeric(Final_table_TF_motif_prediction_NO_NA_sub$motif_start),
      end=as.numeric(Final_table_TF_motif_prediction_NO_NA_sub$motif_end),
      name=Final_table_TF_motif_prediction_NO_NA_sub$Peak_ID))
  
  
  #### LiftOver 38 -> 37 motifs ----
  
  
  motif_in_38_df<-data.frame(chr=as.character(seqnames(gr_Final_table_TF_motif_prediction_NO_NA_sub)),
                             motif_start=start(gr_Final_table_TF_motif_prediction_NO_NA_sub),
                             motif_end=end(gr_Final_table_TF_motif_prediction_NO_NA_sub),
                             Peak_ID=names(gr_Final_table_TF_motif_prediction_NO_NA_sub),
                             stringsAsFactors = F)
  
  cat("motif_in_38_df_0\n")
  str(motif_in_38_df)
  cat("\n")
  
  #path = system.file(package="liftOver", "extdata", "Hg19Tohg38.over.chain")
  ch = import.chain("/home/manuel.tardaguila/reference_files/hg38ToHg19.over.chain")
  
  seqlevelsStyle(gr_Final_table_TF_motif_prediction_NO_NA_sub) = "UCSC"  # necessary
  gr_Final_table_TF_motif_prediction_NO_NA_sub37 = liftOver(gr_Final_table_TF_motif_prediction_NO_NA_sub, ch)
  gr_Final_table_TF_motif_prediction_NO_NA_sub37 = unlist(gr_Final_table_TF_motif_prediction_NO_NA_sub37)
  genome(gr_Final_table_TF_motif_prediction_NO_NA_sub37) = "hg19"
  
  if(length(gr_Final_table_TF_motif_prediction_NO_NA_sub37) >0)
  {
    
    chr_37<-as.character(seqnames(gr_Final_table_TF_motif_prediction_NO_NA_sub37))
    names_37<-as.character(names(gr_Final_table_TF_motif_prediction_NO_NA_sub37))
    
    
    
    
    
    
    motif_in_37_df<-data.frame(chr=as.character(seqnames(gr_Final_table_TF_motif_prediction_NO_NA_sub37)),
                               motif_start_37=start(gr_Final_table_TF_motif_prediction_NO_NA_sub37),
                               motif_end_37=end(gr_Final_table_TF_motif_prediction_NO_NA_sub37),
                               Peak_ID=names(gr_Final_table_TF_motif_prediction_NO_NA_sub37),
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
                     Final_table_TF_motif_prediction_NO_NA_sub,
                     by=c("chr","motif_start","motif_end","Peak_ID"),
                     all=T)
    
    cat("motif_DEF_3\n")
    str(motif_DEF)
    cat("\n")
    
    
    
    
  }else{
    
    
    stop("NO_LIFT_OVER\n")
    
  }# length(gr_Final_table_TF_motif_prediction_NO_NA_sub37) >0
  
  

  
  
  ############# LOOP --------------------------
  
  List_RESULTS<-list()
 
  
  DEBUG<-0
  
  compressed_results<-data.frame()
  
  detailed_results<-data.frame()
  
  array_ENSG<-unique(motif_DEF$ensembl_gene_id)
  
  START<-1
  
  for(i in START:length(array_ENSG)){
    
    ensembl_gene_id_sel<-array_ENSG[i]

    
    cat("-------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(ensembl_gene_id_sel)))
    cat("\n")
    
    ENCODE_file_sel<-ENCODE_file[which(ENCODE_file$ensembl_gene_id == ensembl_gene_id_sel),]
    
    if(dim(ENCODE_file_sel)[1] >0){
      
      if(DEBUG ==1)
      {
        cat("ENCODE_file_sel\n")
        cat(str(ENCODE_file_sel))
        cat("\n")
      }
      
          gr_ENCODE <- GRanges(
            seqnames = as.character(ENCODE_file_sel$chrom),
            cell_type=as.character(ENCODE_file_sel$cell_type),
            experiment=as.character(ENCODE_file_sel$file),
            score=as.numeric(ENCODE_file_sel$score),
            ranges=IRanges(
              start=as.numeric(ENCODE_file_sel$chromStart),
              end=as.numeric(ENCODE_file_sel$chromEnd),
              name=ENCODE_file_sel$ensembl_gene_id))
      
      motif_DEF_sel<-motif_DEF[which(motif_DEF$ensembl_gene_id == ensembl_gene_id_sel),]
      
      
        
        if(DEBUG ==1)
        {
          cat("motif_DEF_sel\n")
          cat(str(motif_DEF_sel))
          cat("\n")
        }
      
          gr_motif_DEF_sel <- GRanges(
            seqnames = as.character(motif_DEF_sel$chr),
            motif_start= as.integer(motif_DEF_sel$motif_start),
            motif_end= as.integer(motif_DEF_sel$motif_end),
            ensembl_gene_id= as.character(motif_DEF_sel$ensembl_gene_id),
            query_region= as.character(motif_DEF_sel$query_region),
            Symbol= as.character(motif_DEF_sel$Symbol),
            source= as.character(motif_DEF_sel$source),
            Motif_ID= as.character(motif_DEF_sel$Motif_ID),
            Intersect_SNP= as.character(motif_DEF_sel$Intersect_SNP),
            ranges=IRanges(
              start=as.numeric(motif_DEF_sel$motif_start_37),
              end=as.numeric(motif_DEF_sel$motif_end_37),
              strand=as.character(motif_DEF_sel$strand),
              name=motif_DEF_sel$Peak_ID))
          
          if(DEBUG == 1)
          {
            cat("gr_motif_DEF_sel_0\n")
            cat(str(gr_motif_DEF_sel))
            cat("\n")
          }
      
      
              m <- findOverlaps(gr_motif_DEF_sel,gr_ENCODE, ignore.strand = TRUE)

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

                ENCODE_HIT_df <- data.frame(chr=as.character(seqnames(gr_ENCODE)),
                                    Peak_start=as.integer(start(gr_ENCODE)),
                                    Peak_end=as.integer(end(gr_ENCODE)),
                                    cell_type=as.character(gr_ENCODE$cell_type),
                                    experiment=as.character(gr_ENCODE$experiment),
                                    score=as.numeric(gr_ENCODE$score), stringsAsFactors = F)

                if(DEBUG == 1)
                {
                  cat("ENCODE_HIT_df_0\n")
                  cat(str(ENCODE_HIT_df))
                  cat("\n")
                }

                ENCODE_HIT_df_hits<-ENCODE_HIT_df[subjectHits_df,]

                if(DEBUG == 1)
                {
                  cat("ENCODE_HIT_df_hits_0\n")
                  cat(str(ENCODE_HIT_df_hits))
                  cat("\n")
                }
                
                
                motif_HIT_df <- data.frame(chr=as.character(seqnames(gr_motif_DEF_sel)),
                                           motif_start_37=as.integer(start(gr_motif_DEF_sel)),
                                           motif_end_37=as.integer(end(gr_motif_DEF_sel)),
                                           strand=strand(gr_motif_DEF_sel),
                                           Peak_ID=names(gr_motif_DEF_sel),
                                           motif_start=as.integer(gr_motif_DEF_sel$motif_start),
                                           motif_end=as.integer(gr_motif_DEF_sel$motif_end),
                                           ensembl_gene_id=as.character(gr_motif_DEF_sel$ensembl_gene_id),
                                           query_region=as.character(gr_motif_DEF_sel$query_region),
                                           Symbol=as.character(gr_motif_DEF_sel$Symbol),
                                           source=as.character(gr_motif_DEF_sel$source),
                                           Motif_ID=as.character(gr_motif_DEF_sel$Motif_ID),
                                           Intersect_SNP=as.character(gr_motif_DEF_sel$Intersect_SNP),
                                             stringsAsFactors = F)
                
                if(DEBUG == 1)
                {
                  cat("motif_HIT_df_SIMPLER\n")
                  cat(str(motif_HIT_df))
                  cat("\n")
                }
               

                motif_HIT_df_hits<-motif_HIT_df[queryHits_motif_DEF_sel,]


                if(DEBUG == 1)
                {
                  cat("motif_HIT_df_hits_0\n")
                  cat(str(motif_HIT_df_hits))
                  cat("\n")
                }
                
                tmp_detailed<-unique(cbind(motif_HIT_df_hits,ENCODE_HIT_df_hits))
                
                
                
                if(DEBUG == 1)
                {
                  cat("tmp_detailed_1\n")
                  cat(str(tmp_detailed))
                  cat("\n")
                }
                
                
                detailed_results<-rbind(tmp_detailed,detailed_results)

                
                
                motif_HIT_df_hits$ENCODE_chip_seq_experiment<-ENCODE_HIT_df_hits$experiment
                motif_HIT_df_hits$ENCODE_cell_type<-ENCODE_HIT_df_hits$cell_type
                
                
                if(DEBUG == 1)
                {
                  cat("motif_HIT_df_hits_1\n")
                  cat(str(motif_HIT_df_hits))
                  cat("\n")
                }
                
                
               
                
                motif_HIT_df_hits.dt<-data.table(motif_HIT_df_hits,
                                                 key=c("chr","motif_start_37","motif_end_37","strand","motif_start","motif_end","ensembl_gene_id"))
                
                motif_HIT_df_hits_collapsed<-as.data.frame(motif_HIT_df_hits.dt[,.(ENCODE_chip_seq_experiment_string=paste(unique(paste(ENCODE_cell_type,ENCODE_chip_seq_experiment, sep="|")),collapse=';')),
                                                                                by=key(motif_HIT_df_hits.dt)], stringsAsFactors=F)
                
                
                if(DEBUG == 1)
                {
                  cat("motif_HIT_df_hits_collapsed_0\n")
                  cat(str(motif_HIT_df_hits_collapsed))
                  cat("\n")
                }
                
                compressed_results<-rbind(motif_HIT_df_hits_collapsed,compressed_results)
                
                

              }#length(queryHits_motif_DEF_sel) >0
      
    }#dim(ENCODE_file_sel)[1] >0
  }#i in START:length(array_ENSG)
    
    
  
  
  cat("detailed_results_0\n")
  cat(str(detailed_results))
  cat("\n")

  cat("compressed_results_0\n")
  cat(str(compressed_results))
  cat("\n")
  
  if(dim(detailed_results)[1] >0){
    
    motif_DEF<-merge(motif_DEF,
                     compressed_results,
                     by=c("chr","motif_start_37","motif_end_37","strand","motif_start","motif_end","ensembl_gene_id"),
                     all.x=T)
    
    cat("motif_DEF_0\n")
    cat(str(motif_DEF))
    cat("\n")
    
    
    ####################### SAVE ###############################################################
    
    
    setwd(out)
    
    write.table(detailed_results, 
                file=paste("TF_motifs_supported_by_ENCODE_chipseq_",table_sel,".tsv",sep=""), 
                row.names = F, quote=F, sep="\t")
    
    
    setwd(out2)
    
    write.table(motif_DEF, 
                file=paste("Final_table_Final_table_TF_motif_prediction_with_chipseq_support_plus_ENCODE_",table_sel,".tsv",sep=""), 
                row.names = F, quote=F, sep="\t")
  }else{
    
    
    setwd(out2)
    
    write.table(motif_DEF, 
                file=paste("Final_table_Final_table_TF_motif_prediction_with_chipseq_support_plus_ENCODE_",table_sel,".tsv",sep=""), 
                row.names = F, quote=F, sep="\t")
    
  }# dim(detailed_results)[1] >0
 
  
 
  
 
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
    make_option(c("--Final_table_TF_motif_prediction"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ENCODE_file"), type="character", default=NULL, 
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

  
}

###########################################################################

system.time( main() )