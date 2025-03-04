
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
 
  #### READ input_REF ----
  
  input_REF<-as.data.frame(fread(file=opt$input_REF, sep="\t", header=T), stringsAsFactors=F)

  cat("input_REF_0\n")
  cat(str(input_REF))
  cat("\n")
 
  
  #### READ input_ALT ----
  
  input_ALT<-as.data.frame(fread(file=opt$input_ALT, sep="\t", header=T), stringsAsFactors=F)

  cat("input_ALT_0\n")
  cat(str(input_ALT))
  cat("\n")
 
  
  
  #### Merge all ----
  
  ALL_df<-rbind(input_REF,input_ALT)
  
  ALL_df<-ALL_df[order(ALL_df$query_region,ALL_df$score),]
  
  cat("ALL_df_0\n")
  cat(str(ALL_df))
  cat("\n")
  
 
  
  
  ##### ALL_df subset -----
  
  ALL_df_subset<-unique(ALL_df[,which(colnames(ALL_df)%in%c('chr','motif_start','motif_end','ensembl_gene_id'))])
  
  
  cat("ALL_df_subset_0\n")
  cat(str(ALL_df_subset))
  cat("\n")
  
  
  
  #### bed file ,I DON'T TRUST export.bed ----
  
  
  Bed_df<-ALL_df_subset
  colnames(Bed_df)[which(colnames(Bed_df) == 'ensembl_gene_id')]<-'name'
  
  indx.reorder<-c(which(colnames(Bed_df) == 'chr'),
                  which(colnames(Bed_df) == 'motif_start'),
                  which(colnames(Bed_df) == 'motif_end'),
                  which(colnames(Bed_df) == 'name'))
  
  cat("indx.reorder_0\n")
  cat(str(indx.reorder))
  cat("\n")
  cat(sprintf(as.character(indx.reorder)))
  cat("\n")
  
  Bed_df<-Bed_df[,indx.reorder]
  
  cat("Bed_df_0\n")
  cat(str(Bed_df))
  cat("\n")
  
  Bed_df$chr<-factor(Bed_df$chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                                 "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                 "chr22","chr23","chrX","chrY"), ordered=T)
  
  
  Bed_df<-Bed_df[order(Bed_df$chr,Bed_df$motif_start, decreasing = F),]
  
  
  cat("Bed_df_1\n")
  cat(str(Bed_df))
  cat("\n")
  
  #### order the bedfile
  
  
  ##### select all the possible TFs in my file -----
  
  names_subset<-as.data.frame(unique(ALL_df[!is.na(ALL_df$ensembl_gene_id),c(which(colnames(ALL_df) == 'ensembl_gene_id'))]), stringsAsFactors=F)
  colnames(names_subset)<-'ensembl_gene_id'
  
  cat("names_subset_0\n")
  cat(str(names_subset))
  cat("\n")
  
 
  
  multiVals <- function(x) paste(x,collapse=";")
  
  
  names_subset$Symbol <- mapIds(org.Hs.eg.db, keys=names_subset$ensembl_gene_id, keytype="ENSEMBL",
                            column="SYMBOL", multiVals=multiVals)
  
  cat("names_subset_1\n")
  cat(str(names_subset))
  cat("\n")
  
  names_subset_long<-unique(as.data.frame(cSplit(names_subset,sep = ';', direction = "long",splitCols = "Symbol"),stringsAsFactors=F))
  
  cat("names_subset_long_0\n")
  cat(str(names_subset_long))
  cat("\n")
  
  
  ############# SAVE ---------------------------
  
  setwd(out)
  
  # export.bed(gr_ALL_df,con='TF_motif_prediction_.bed') #### I DON'T TRUST export.bed
  
  write.table(Bed_df, file=paste('TF_motif_prediction',".bed",sep=''), sep="\t", quote=F, row.names = F, col.names = F)
  write.table(ALL_df, file=paste('Final_table_TF_motif_prediction',".tsv",sep=''), sep="\t", quote=F, row.names = F)
  write.table(names_subset_long, file=paste('List_of_TFs_for_occupancy',".tsv",sep=''), sep="\t", quote=F, row.names = F)
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
    make_option(c("--motifs_collected"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_REF"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_ALT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--spanning_of_motif"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--allele"), type="numeric", default=NULL, 
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