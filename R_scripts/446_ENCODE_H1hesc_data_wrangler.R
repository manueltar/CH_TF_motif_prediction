
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
suppressMessages(library("R.utils", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

options(warn = 1)

multiVals <- function(x) paste(x,collapse=";")

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
  
  #### READ and transform cell_type_array ----
  
  cell_type_array = unlist(strsplit(opt$cell_type_array, split=","))
  
  cat("cell_type_array_0\n")
  cat(sprintf(as.character(cell_type_array)))
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
  
  files_df_sel<-files_df[grep("UniPk\\.txt\\.gz$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
  
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'filename'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  df_files$path<-paste(indir,files_df_sel,sep='')
  df_files$filename<-files_df_sel
  df_files$filename<-gsub("\\.txt\\.gz$","",df_files$filename)
  
  
  
  
  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
  
  starting_point<-data.frame()
  
  for(i in 1:length(cell_type_array)){
    
    cell_type_array_sel<-cell_type_array[i]
    
    cat("--------------------------------------->\t")
    cat(sprintf(as.character(cell_type_array_sel)))
    cat("\n")
    
    master_indx<-grep(cell_type_array_sel,df_files$filename)
    
    df_ct_sel<-df_files[master_indx,]
    
    cat("df_ct_sel_1\n")
    cat(str(df_ct_sel))
    cat("\n")
    
    if(dim(df_ct_sel)[1] >0){
      
      df_ct_sel$CHIP_factor<-gsub(paste('.+',cell_type_array_sel,sep=''),"",df_ct_sel$filename)
      
      
      df_ct_sel$CHIP_factor<-gsub("a[0-9]{3,}a","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("Pcr[0-9]+.+$","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("V[0-9]+UniPk$","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("sc[0-9]{3,}","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("UniPk$","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("Iggrab.*$","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("Ucd.*$","",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("[0-9]{5,}","",df_ct_sel$CHIP_factor)
      
      
      manually_searched_patterns<-paste(c("a301218a","39875","ab26049","a300109a","sc14700","a301218a","200401194","nb6001263"), collapse="|")
      
      df_ct_sel$CHIP_factor<-gsub(manually_searched_patterns,"",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("Cmyc","myc",df_ct_sel$CHIP_factor)
      df_ct_sel$CHIP_factor<-gsub("Ifn[ag][0-9]+","",df_ct_sel$CHIP_factor)
      
      
      
      
      
      df_ct_sel$CHIP_factor<-toupper(df_ct_sel$CHIP_factor)
      
      df_ct_sel$CHIP_factor[which(df_ct_sel$CHIP_factor == 'PU1')]<-'SP1'
      
      df_ct_sel$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=df_ct_sel$CHIP_factor, keytype="SYMBOL",column="ENSEMBL", multiVals=multiVals)
      
      
      if(cell_type_array_sel == 'H1hesc'){
        
        df_ct_sel$cell_type<-'hESC_H1'
        
      }
      
      if(cell_type_array_sel == 'K562'){
        
        df_ct_sel$cell_type<-'K562'
        
      }
      
      if(cell_type_array_sel == 'Hl60'){
        
        df_ct_sel$cell_type<-'HL60'
        
      }
      
      starting_point<-rbind(df_ct_sel,starting_point)
      
    }#dim(df_ct_sel)[1] >0
    
    
    
  }#i in 1:length(cell_type_array)
 
  
  cat("starting_point_0\n")
  cat(str(starting_point))
  cat("\n")
  
 
  
  
  
  indx<-which(starting_point$CHIP_factor == "")
  
  starting_point$CHIP_factor[indx]<-'DNASE'
  starting_point$ensembl_gene_id[indx]<-NA
  
  
  cat("starting_point_3\n")
  cat(str(starting_point))
  cat("\n")
  
 
  setwd(out)
  
  write.table(starting_point, file='test.tsv', sep="\t",quote=F,row.names=F)
  
 
  ##### Reading loop -----
  
  filename_array<-starting_point$filename
  
  cat("filename_array_0\n")
  cat(str(filename_array))
  cat("\n")
  

  DEBUG<-0
  
  list_results<-list()
  
  
  for(i in 1:length(filename_array))
  {
    filename_array_sel<-filename_array[i]
    
    cat("------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(filename_array_sel)))
    cat("\n")
    
    starting_point_sel<-starting_point[which(starting_point$filename == filename_array_sel),]
    

    if(DEBUG == 1)
    {
      cat("starting_point_sel_0\n")
      cat(str(starting_point_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(starting_point_sel)[1] >0)
    {
      sel_file<-starting_point_sel$path
      
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
          CHIP_result<-as.data.frame(fread(file=sel_file, sep="\t", header=F), stringsAsFactors=F)
          
          colnames(CHIP_result)<-c("bin","chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak") # See  https://genome.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeAwgTfbsUniform&hgta_table=wgEncodeAwgTfbsHaibH1hescBcl11aPcr1xUniPk&hgta_doSchema=describe+table+schema
          

          if(DEBUG == 1)
          {
            cat("CHIP_result_0\n")
            cat(str(CHIP_result))
            cat("\n")
          }
          
          CHIP_result$CHIP_factor<-starting_point_sel$CHIP_factor
          CHIP_result$ensembl_gene_id<-starting_point_sel$ensembl_gene_id
          CHIP_result$cell_type<-starting_point_sel$cell_type
          CHIP_result$file<-starting_point_sel$filename
          
          indx.reorder<-c(which(colnames(CHIP_result)%in%c('CHIP_factor','ensembl_gene_id','cell_type','file')),
                          which(colnames(CHIP_result)%in%c("bin","chrom","chromStart","chromEnd","name","score","strand",
                                                           "signalValue","pValue","qValue","peak")))
         
          CHIP_result<-CHIP_result[,indx.reorder]
         
          if(DEBUG == 1)
          {
            cat("CHIP_result_1\n")
            cat(str(CHIP_result))
            cat("\n")
          }
          
         list_results[[i]]<-CHIP_result
         
        
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(starting_point_sel)[1] >0
    
  }# i in 1:length(filename_array)
  
  if(length(list_results) >0)
  {
    Results = unique(as.data.frame(data.table::rbindlist(list_results, fill = T)))
    
    cat("Results_0\n")
    cat(str(Results))
    cat("\n")
    
  
    
    ### save ----
    
    setwd(out)
    
    write.table(Results, file='ENCODE_CHIP_data.tsv', sep="\t",quote=F,row.names=F)
    
    gzip('ENCODE_CHIP_data.tsv',overwrite = TRUE)
    
  }# length(list_results) >0
  

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
    make_option(c("--cell_type_array"), type="character", default=NULL, 
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

  
}


###########################################################################

system.time( main() )