#!/bin/bash
 
eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Master_path_analysis=$(echo "$output_dir")

 
Log_files=$(echo "$Master_path_analysis""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

####################### Peaks are in hg19 liftOver
####################### Peaks are in hg19 liftOver
####################### Peaks are in hg19 liftOver

type=$(echo "ENCODE_Match_TF_occupancy")
outfile_ENCODE_Match_TF_occupancy=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_ENCODE_Match_TF_occupancy
echo -n "" > $outfile_ENCODE_Match_TF_occupancy
name_ENCODE_Match_TF_occupancy=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")


Rscript_ENCODE_Match_TF_occupancy=$(echo "$Rscripts_path""452_ENCODE_CHIPseq_controls.R")

ENCODE_file=$(echo "/group/soranzo/manuel.tardaguila/CH/reference_files/ENCODE/CHIPseq/ENCODE_CHIP_data.tsv.gz")
selected_cell_types=$(echo "hESC_H1,K562")
TF_to_search=$(echo "BCL11A,GATA6")



myjobid_ENCODE_Match_TF_occupancy=$(sbatch --output=$outfile_ENCODE_Match_TF_occupancy --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4096 --parsable --job-name $name_ENCODE_Match_TF_occupancy --wrap="Rscript $Rscript_ENCODE_Match_TF_occupancy --ENCODE_file $ENCODE_file --selected_cell_types $selected_cell_types --TF_to_search $TF_to_search --type $type --out $output_dir")
myjobid_seff_ENCODE_Match_TF_occupancy=$(sbatch --dependency=afterany:$myjobid_ENCODE_Match_TF_occupancy --open-mode=append --output=$outfile_ENCODE_Match_TF_occupancy --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_ENCODE_Match_TF_occupancy >> $outfile_ENCODE_Match_TF_occupancy")

