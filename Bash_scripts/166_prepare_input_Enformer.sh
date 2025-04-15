#!/bin/bash>

MASTER_ROUTE=$1
output_dir=$MASTER_ROUTE

MASTER_ROUTE=$1
analysis=$2


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""/""$analysis""/")

Log_files=$(echo "$output_dir""Log_files""/")

rm -rf $Log_files
mkdir -p $Log_files
 
#######################################################################prepare_input_Enformer ###################################


type=$(echo "prepare_input_Enformer""_""$analysis")
outfile_prepare_input_Enformer=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_prepare_input_Enformer
echo -n "" > $outfile_prepare_input_Enformer
name_prepare_input_Enformer=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_prepare_input_Enformer=$(echo "$Rscripts_path""460_Prepare_Enformer_input.R")

Table_of_variants=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv")

myjobid_prepare_input_Enformer=$(sbatch --job-name=$name_prepare_input_Enformer --output=$outfile_prepare_input_Enformer --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_prepare_input_Enformer --Table_of_variants $Table_of_variants --type $type --out $output_dir")
myjobid_seff_prepare_input_Enformer=$(sbatch --dependency=afterany:$myjobid_prepare_input_Enformer --open-mode=append --output=$outfile_prepare_input_Enformer --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_prepare_input_Enformer >> $outfile_prepare_input_Enformer")
 

