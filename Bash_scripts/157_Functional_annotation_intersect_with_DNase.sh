#!/bin/bash

MASTER_ROUTE=$1
analysis=$2

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0

eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

#### Intersect_DNase ####


type=$(echo "Intersect_DNase")
outfile_Intersect_DNase=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Intersect_DNase
echo -n "" > $outfile_Intersect_DNase
name_Intersect_DNase=$(echo "$type""_job")


Rscript_Intersect_DNase=$(echo "$Rscripts_path""451_ENCODE_DNase.R")

 
ENCODE_file=$(echo "/group/soranzo/manuel.tardaguila/CH/reference_files/ENCODE/CHIPseq/ENCODE_CHIP_data.tsv.gz")
Table_of_variants=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv")
CHIP_variants_file=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/CHIP_variants.tsv")

myjobid_Intersect_DNase=$(sbatch  --job-name=$name_Intersect_DNase --output=$outfile_Intersect_DNase --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=3 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Intersect_DNase --ENCODE_file $ENCODE_file --Table_of_variants $Table_of_variants --CHIP_variants_file $CHIP_variants_file  --type $type --out $output_dir")
myjobid_seff_Intersect_DNase=$(sbatch --dependency=afterany:$myjobid_Intersect_DNase --open-mode=append --output=$outfile_Intersect_DNase --job-name=$(echo "seff""_""$name_Intersect_DNase") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect_DNase >> $outfile_Intersect_DNase")

