#!/bin/bash
 
eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2
List_of_TFs_for_occupancy=$3

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Master_path_analysis=$(echo "$output_dir")

 
Log_files=$(echo "$Master_path_analysis""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files

Search_result=$(echo "$Master_path_analysis""/""Search_result/")
Search_result_Log_files=$(echo "$Master_path_analysis""/""Search_result/""Log_files""/")

rm -rf $Search_result
mkdir -p $Search_result
rm -rf $Search_result_Log_files
mkdir -p $Search_result_Log_files


# ####################################################################################################################################################################### search_for_experiments_TF_occupancy

# type=$(echo "search_for_experiments_TF_occupancy")
# outfile_search_for_experiments_TF_occupancy=$(echo "$Log_files""outfile_1_""$type"".log")
# touch $outfile_search_for_experiments_TF_occupancy
# echo -n "" > $outfile_search_for_experiments_TF_occupancy
# name_search_for_experiments_TF_occupancy=$(echo "$type""_job")
# seff_name=$(echo "seff""_""$type")


# Rscript_search_for_experiments_TF_occupancy=$(echo "$Rscripts_path""430_extract_experiments_for_TF_occupancy.R")
 
 
 
# fileList=$(echo "/group/soranzo/manuel.tardaguila/MPRA_explore_results/TF_motif_analysis/reference_files/fileList.tab")


# myjobid_search_for_experiments_TF_occupancy=$(sbatch --output=$outfile_search_for_experiments_TF_occupancy --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024 --parsable --job-name $name_search_for_experiments_TF_occupancy --wrap="Rscript $Rscript_search_for_experiments_TF_occupancy --List_of_TFs_for_occupancy $List_of_TFs_for_occupancy --fileList $fileList --type $type --out $output_dir")
# myjobid_seff_search_for_experiments_TF_occupancy=$(sbatch --dependency=afterany:$myjobid_search_for_experiments_TF_occupancy --open-mode=append --output=$outfile_search_for_experiments_TF_occupancy --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_search_for_experiments_TF_occupancy >> $outfile_search_for_experiments_TF_occupancy")

# exit

config_file=$(echo "$output_dir""TF_experiments_search_file_1.tsv") #limit of 1000 jobs submitted per batch submission 3 config files to explore all the experiments
#config_file=$(echo "$output_dir""TF_experiments_search_file_2.tsv")
#config_file=$(echo "$output_dir""TF_experiments_search_file_3.tsv")



echo "$config_file"

declare -a arr

end_of_file=0
counter=0
array=()
array_Merge=()





while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?


  ((counter++))

#    echo "counter:$counter  LINE:$line"

  if [ $counter != 0 ]; then

      if [[ "$line" =~ [^[:space:]] ]]; then

        a=($(echo "$line" | tr "\t" '\n'))

        ENSG=${a[0]}
        experiment=${a[2]}
        symbol_string=${a[1]}



      echo "ENSG:$ENSG  experiment:$experiment  symbol_string:$symbol_string"

      tag=$(echo "$ENSG""__""$experiment")

      echo "------------------------------------------------------------------------------------->$tag"


      type=$(echo "$tag""_search_ChIP_atlas")
      outfile_search_ChIP_atlas=$(echo "$Search_result_Log_files""outfile_""$type"".log")
      touch $outfile_search_ChIP_atlas
      echo -n "" > $outfile_search_ChIP_atlas
      name_search_ChIP_atlas=$(echo "$type""_job")
      seff_name=$(echo "seff""_""$type")

      ChIP_atlas=$(echo "/group/soranzo/manuel.tardaguila/CH/TF_motif_analysis/TF_occupancy_filter/reference/allPeaks_light.hg19.10.bed.gz")
      outfile=$(echo "$Search_result""$tag"".tsv")


      myjobid_search_ChIP_atlas=$(sbatch --output=$outfile_search_ChIP_atlas --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128 --parsable --job-name $name_search_ChIP_atlas --wrap="zcat $ChIP_atlas|grep $experiment > $outfile")
      myjobid_seff_search_ChIP_atlas=$(sbatch --dependency=afterany:$myjobid_search_ChIP_atlas --open-mode=append --output=$outfile_search_ChIP_atlas --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_search_ChIP_atlas >> $outfile_search_ChIP_atlas
")

      exit



      fi #no spaces
    fi #counter > 1
done < "$config_file"

exit


####################################################################################################################

#$ cut -f1 TF_experiments_search_file.tsv|sort|uniq -c|sort -rg|awk -F" " '{print $2}' > ENSG_file.txt


ENSG_file=$(echo "$output_dir""ENSG_file.txt")


echo "$ENSG_file"

declare -a arr

end_of_file=0
counter=0
array=()
array_Merge=()

Search_result=$(echo "$Master_path_analysis""/""Search_result/")

Search_result_Log_files=$(echo "$Master_path_analysis""/""Search_result/""Log_files""/") 

while [[ $end_of_file == 0 ]]
do
  read -r line
  end_of_file=$?


  ((counter++))

#    echo "counter:$counter  LINE:$line"

  if [ $counter != 0 ]; then

      if [[ "$line" =~ [^[:space:]] ]]; then

        a=($(echo "$line" | tr "\t" '\n'))

        ENSG=${a[0]}
	outfile=$(echo "$ENSG""__""FINAL_file.out.gz")

	echo "--------------------->ENSG:$ENSG"
	search_param=$(echo "$ENSG""__""*"'.tsv')
	echo "$search_param"

	outfile=$(echo "$Search_result""$ENSG"".out.gz")

	echo "$outfile"

	type=$(echo "compress""_""$ENSG")
	outfile_compress_results_ChIP_atlas=$(echo "$Search_result_Log_files""outfile_""$type"".log")
	touch $outfile_compress_results_ChIP_atlas
	echo -n "" > $outfile_compress_results_ChIP_atlas
	name_compress_results_ChIP_atlas=$(echo "$type""_job")
	seff_name=$(echo "seff""_""$type")




	cd $Search_result

	myjobid_compress_results_ChIP_atlas=$(sbatch --output=$outfile_compress_results_ChIP_atlas --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128 --parsable --job-name $name_compress_results_ChIP_atlas --wrap="cat $search_param |gzip -9 > $outfile")
       myjobid_seff_compress_results_ChIP_atlas=$(sbatch --dependency=afterany:$myjobid_compress_results_ChIP_atlas --open-mode=append --output=$outfile_compress_results_ChIP_atlas --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_compress_results_ChIP_atlas >> $outfile_compress_results_ChIP_atlas")



      fi #no spaces
    fi #counter > 1
done < "$ENSG_file"

exit

####################################################################################################################################################################### Match_TF_occupancy

type=$(echo "Match_TF_occupancy")
outfile_Match_TF_occupancy=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_Match_TF_occupancy
echo -n "" > $outfile_Match_TF_occupancy
name_Match_TF_occupancy=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


Rscript_Match_TF_occupancy=$(echo "$Rscripts_path""342_MPRA_bc_cloned_collapse.R")
 
 
TF_motif_prediction_MPRA=$(echo "/group/soranzo/manuel.tardaguila/MPRA_explore_results/TF_motif_analysis/""TF_motif_prediction_MPRA.bed")
indir_chip_atlas_search_results=$(echo "/group/soranzo/manuel.tardaguila/MPRA_explore_results/TF_motif_analysis/TF_occupancy_filter/Search_result/")
TF_motifs_with_an_experimental_HIT=$(echo "/group/soranzo/manuel.tardaguila/MPRA_explore_results/TF_motif_analysis/TF_occupancy_filter/TF_motifs_supported_by_chipseq.tsv")
Final_table_TF_motif_prediction=$(echo "/group/soranzo/manuel.tardaguila/MPRA_explore_results/TF_motif_analysis/Final_table_TF_motif_prediction_MPRA.tsv")
out2=$(echo "/group/soranzo/manuel.tardaguila/MPRA_explore_results/TF_motif_analysis/")

myjobid_Match_TF_occupancy=$(sbatch --output=$outfile_Match_TF_occupancy --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024 --parsable --job-name $name_Match_TF_occupancy --wrap="Rscript $Rscript_Match_TF_occupancy --TF_motif_prediction_MPRA $TF_motif_prediction_MPRA --indir_chip_atlas_search_results $indir_chip_atlas_search_results --TF_motifs_with_an_experimental_HIT $TF_motifs_with_an_experimental_HIT --Final_table_TF_motif_prediction $Final_table_TF_motif_prediction --type $type --out $output_dir --out2 $out2")
myjobid_seff_Match_TF_occupancy=$(sbatch --dependency=afterany:$myjobid_Match_TF_occupancy --open-mode=append --output=$outfile_Match_TF_occupancy --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Match_TF_occupancy >> $outfile_Match_TF_occupancy")
