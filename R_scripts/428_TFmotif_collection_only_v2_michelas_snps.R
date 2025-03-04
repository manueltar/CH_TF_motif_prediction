
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
  
  #### READ and transform allele ----
  
  allele = opt$allele
  
  cat("allele_\n")
  cat(sprintf(as.character(allele)))
  cat("\n")
  
  #### READ and transform spanning_of_motif ----
  
  spanning_of_motif = opt$spanning_of_motif
  
  cat("spanning_of_motif_0\n")
  cat(str(spanning_of_motif))
  cat("\n")
  
  #### READ input_file_HOMER ----
  
  input_file_HOMER<-as.data.frame(fread(file=opt$input_file_HOMER, sep="\t", header=F, skip=5), stringsAsFactors=F)
  colnames(input_file_HOMER)<-c('query_region','start','end','Motif_ID','score','strand')
  
  cat("input_file_HOMER_0\n")
  cat(str(input_file_HOMER))
  cat("\n")
  cat(str(unique(input_file_HOMER$Motif_ID)))
  cat("\n")
  
  input_file_HOMER$Symbol<-input_file_HOMER$Motif_ID
  
  input_file_HOMER$Symbol<-gsub("/Homer$","",input_file_HOMER$Symbol)
  
  
  input_file_HOMER_long<-unique(as.data.frame(cSplit(input_file_HOMER,sep = '/', direction = "long",splitCols = "Symbol"),stringsAsFactors=F))
  
  
  
  input_file_HOMER_long$Symbol<-gsub("^[^_]+_","",input_file_HOMER_long$Symbol)
  input_file_HOMER_long$Symbol<-gsub("\\(.+$","",input_file_HOMER_long$Symbol)
  input_file_HOMER_long$Symbol<-gsub("-halfsite","",input_file_HOMER_long$Symbol)
  
  input_file_HOMER_long$source<-'HOMER'
  
  
  cat("input_file_HOMER_long_3\n")
  cat(str(input_file_HOMER_long))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_HOMER_long$Symbol))))
  cat("\n")
  
  input_file_HOMER_long_CHIP_seq<-input_file_HOMER_long[grep("ChIP-Seq",input_file_HOMER_long$Symbol),]
  
  
  input_file_HOMER_long_CHIP_seq$Symbol<-gsub("^[^-]+-","",input_file_HOMER_long_CHIP_seq$Symbol)
  input_file_HOMER_long_CHIP_seq$Symbol<-gsub("-.+$","",input_file_HOMER_long_CHIP_seq$Symbol)

  cat("input_file_HOMER_long_CHIP_seq_0\n")
  cat(str(input_file_HOMER_long_CHIP_seq))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_HOMER_long_CHIP_seq$Symbol))))
  cat("\n")
  
  
  
  
  df_HOMER<-unique(input_file_HOMER_long_CHIP_seq[,c(which(colnames(input_file_HOMER_long_CHIP_seq) == 'query_region'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'Symbol'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'source'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'start'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'end'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'Motif_ID'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'strand'),
                                which(colnames(input_file_HOMER_long_CHIP_seq) == 'score'))])
  
  cat("df_HOMER_0\n")
  cat(str(df_HOMER))
  cat("\n")
  cat(str(unique(df_HOMER$Symbol)))
  cat("\n")
  
  #### READ input_file_JASPAR ----
  
  input_file_JASPAR<-as.data.frame(fread(file=opt$input_file_JASPAR, sep="\t", header=F, skip=5), stringsAsFactors=F)
  colnames(input_file_JASPAR)<-c('query_region','start','end','Motif_ID','score','strand')
  
  cat("input_file_JASPAR_0\n")
  cat(str(input_file_JASPAR))
  cat("\n")
  cat(str(unique(input_file_JASPAR$Motif_ID)))
  cat("\n")
  
  input_file_JASPAR$Symbol<-input_file_JASPAR$Motif_ID
  
  input_file_JASPAR$Symbol<-gsub("^[^_]+_","",input_file_JASPAR$Symbol)
  input_file_JASPAR$Symbol<-gsub("\\(.+$","",input_file_JASPAR$Symbol)
  
    
  # cat("input_file_JASPAR_2\n")
  # cat(str(input_file_JASPAR))
  # cat("\n")
  # cat(sprintf(as.character(unique(input_file_JASPAR$Symbol))))
  # cat("\n")
  
  input_file_JASPAR.dt<-data.table(input_file_JASPAR, key=c('query_region','start','end','Motif_ID','score','strand'))
  
  input_file_JASPAR_DEF<-as.data.frame(input_file_JASPAR.dt[,.(Symbol=unlist(strsplit(Symbol, split='::'))), by = key(input_file_JASPAR.dt)], stringsAsFactors=F)
  input_file_JASPAR_DEF$source<-'JASPAR'
  
  
  cat("input_file_JASPAR_DEF_2\n")
  cat(str(input_file_JASPAR_DEF))
  cat("\n")
  cat(sprintf(as.character(unique(input_file_JASPAR_DEF$Symbol))))
  cat("\n")
  
  df_JASPAR<-unique(input_file_JASPAR_DEF[,c(which(colnames(input_file_JASPAR_DEF) == 'query_region'),
                                             which(colnames(input_file_JASPAR_DEF) == 'Symbol'),
                                             which(colnames(input_file_JASPAR_DEF) == 'source'),
                                             which(colnames(input_file_JASPAR_DEF) == 'start'),
                                             which(colnames(input_file_JASPAR_DEF) == 'end'),
                                             which(colnames(input_file_JASPAR_DEF) == 'Motif_ID'),
                                             which(colnames(input_file_JASPAR_DEF) == 'strand'),
                                             which(colnames(input_file_JASPAR_DEF) == 'score'))])
  
  cat("df_JASPAR_0\n")
  cat(str(df_JASPAR))
  cat("\n")
  cat(str(unique(df_JASPAR$Symbol)))
  cat("\n")
  
  
  
  #### Merge all ----
  
  ALL_df<-rbind(df_HOMER,df_JASPAR)
  
  cat("ALL_df_0\n")
  cat(str(ALL_df))
  cat("\n")
  cat(str(unique(ALL_df$Symbol)))
  cat("\n")
  
  ALL_df_broken1<-as.data.frame(cSplit(ALL_df, splitCols= "query_region", sep = "|", direction = "wide", fixed = TRUE, drop=F), stringsAsFactors=F)
  
  cat("ALL_df_broken1_0\n")
  cat(str(ALL_df_broken1))
  cat("\n")
  cat(str(unique(ALL_df_broken1$Symbol)))
  cat("\n")
  
  ALL_df_broken2<-as.data.frame(cSplit(ALL_df_broken1, splitCols= "query_region_2", sep = "_", direction = "wide", fixed = TRUE, drop=F), stringsAsFactors=F)
  
  cat("ALL_df_broken2_0\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  cat(str(unique(ALL_df_broken2$Symbol)))
  cat("\n")
  
  
  colnames(ALL_df_broken2)[which(colnames(ALL_df_broken2) == 'query_region_1')]<-'rsid'
  colnames(ALL_df_broken2)[which(colnames(ALL_df_broken2) == 'query_region_2')]<-'VAR_38'
  colnames(ALL_df_broken2)[which(colnames(ALL_df_broken2) == 'query_region_2_1')]<-'chr'
  colnames(ALL_df_broken2)[which(colnames(ALL_df_broken2) == 'query_region_2_2')]<-'pos'
  colnames(ALL_df_broken2)[which(colnames(ALL_df_broken2) == 'query_region_2_3')]<-'ref'
  colnames(ALL_df_broken2)[which(colnames(ALL_df_broken2) == 'query_region_2_4')]<-'alt'
  
  
  cat("ALL_df_broken2_1\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  cat(str(unique(ALL_df_broken2$Symbol)))
  cat("\n")
  
  
  ALL_df_broken2$ensembl_gene_id <- mapIds(org.Hs.eg.db, keys=ALL_df_broken2$Symbol, keytype="SYMBOL",
                           column="ENSEMBL", multiVals="first")
  
  cat("ALL_df_broken2_1\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  
  ### harcoded recovery of PU.1 and others due to alias problems with the motif dBs ----
  
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'PU.1')]<-'ENSG00000066336'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'PU.1-IRF')]<-'ENSG00000066336'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'PU.1:IRF8')]<-'ENSG00000066336'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'E2A')]<-'ENSG00000071564'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'HEB')]<-'ENSG00000140262'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'RARa')]<-'ENSG00000131759'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'TCFL2')]<-'ENSG00000099949'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'THRb')]<-'ENSG00000151090'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'BMYB')]<-'ENSG00000101057'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'HIF-1a')]<-'ENSG00000100644'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'HIF-1b')]<-'ENSG00000143437'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'EAR2')]<-'ENSG00000160113'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'NFkB-p50,p52')]<-'ENSG00000109320,ENSG00000077150'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'RXR')]<-'ENSG00000186350'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'ZNF143|STAF')]<-'ENSG00000166478'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'AP-1')]<-'ENSG00000175592'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'TEAD')]<-'ENSG00000187079'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'EKLF')]<-'ENSG00000105610'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'AMYB')]<-'ENSG00000185697'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'HNF6')]<-'ENSG00000169856'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'ETS')]<-'ENSG00000134954'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol%in%c('NFkB-p65-Rel','NFkB-p65'))]<-'ENSG00000173039'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'HIF2a')]<-'ENSG00000116016'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'NPAS')]<-'ENSG00000130751'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'bHLHE40')]<-'ENSG00000134107'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'bHLHE41')]<-'ENSG00000123095'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'BMAL1')]<-'ENSG00000133794'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'RORa')]<-'ENSG00000069667'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'RORg')]<-'ENSG00000143365'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'NFAT')]<-'ENSG00000131196'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'CEBP')]<-'ENSG00000245848'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'RUNX')]<-'ENSG00000159216'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'LXRE')]<-'ENSG00000025434'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'NFY')]<-'ENSG00000001167'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'STAT5')]<-'ENSG00000126561'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'E2F')]<-'ENSG00000101412'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'RARg')]<-'ENSG00000172819'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'EBF')]<-'ENSG00000164330'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'CEBP:CEBP')]<-'ENSG00000245848'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'BORIS')]<-'ENSG00000124092'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'CENPBD1')]<-'ENSG00000177946'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'ZSCAN5')]<-'ENSG00000131848'
  ALL_df_broken2$ensembl_gene_id[which(ALL_df_broken2$Symbol == 'HKR1')]<-'ENSG00000181666'
  
  
  cat("ALL_df_broken2_2\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  
  # #### RMV motifs without TF id ----
  # 
  # ALL_df_broken2_NO_NA<-ALL_df_broken2[!is.na(ALL_df_broken2$ensembl_gene_id),]
  # 
  # cat("ALL_df_broken2_NO_NA_0\n")
  # cat(str(ALL_df_broken2_NO_NA))
  # cat("\n")
  # cat(str(unique(ALL_df_broken2_NO_NA$Symbol)))
  # cat("\n")
  
  #################### Add 1 to start + end because ALL_df_broken2 is 0 based  ----------------------
  
  ALL_df_broken2$start<-ALL_df_broken2$start+1
  ALL_df_broken2$end<-ALL_df_broken2$end+1
  
  cat("ALL_df_broken2_1\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  cat(str(unique(ALL_df_broken2$Symbol)))
  cat("\n")
  
  
  ############### Intersect the SNP position or not ---------------
  
  ALL_df_broken2$Intersect_SNP<-NA
  
  ALL_df_broken2$Intersect_SNP[which(ALL_df_broken2$start <= spanning_of_motif & ALL_df_broken2$end >= spanning_of_motif)]<-'YES'
  ALL_df_broken2$Intersect_SNP[-which(ALL_df_broken2$start <= spanning_of_motif & ALL_df_broken2$end >= spanning_of_motif)]<-'NO'
  
  cat("ALL_df_broken2_2\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  cat(str(unique(ALL_df_broken2$Symbol)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ALL_df_broken2$Intersect_SNP))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ALL_df_broken2$Intersect_SNP)))))
  cat("\n")
  
  ############### relative shifts to calculate absolute coordinates of the motif ---------------
  
  ALL_df_broken2$rel_shift_dw<-ALL_df_broken2$start - spanning_of_motif
  ALL_df_broken2$rel_shift_up<-ALL_df_broken2$end - spanning_of_motif
  

  cat("ALL_df_broken2_3\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  
  
  ALL_df_broken2$motif_start<-ALL_df_broken2$pos+ALL_df_broken2$rel_shift_dw
  ALL_df_broken2$motif_end<-ALL_df_broken2$pos+ALL_df_broken2$rel_shift_up
  
  ALL_df_broken2$allele<-allele
  

  cat("ALL_df_broken2_4\n")
  cat(str(ALL_df_broken2))
  cat("\n")
  
  ##### save -----
  
  setwd(out)
  write.table(ALL_df_broken2, file=paste(type,".tsv",sep=''), sep="\t", quote=F, row.names = F)
  
  
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
    make_option(c("--input_file_HOMER"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--input_file_JASPAR"), type="character", default=NULL, 
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