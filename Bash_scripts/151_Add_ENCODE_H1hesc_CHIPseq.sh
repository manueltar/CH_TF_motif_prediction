#!/bin/bash
 
MASTER_ROUTE=$1
analysis=$2


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

#### Gather_data ###################################

type=$(echo "Gather_data""_""$analysis")
outfile_Gather_data=$(echo "$Log_files""outfile_1_""$type"".out")
touch $outfile_Gather_data
echo -n "" > $outfile_Gather_data
name_Gather_data=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Gather_data=$(echo "$Rscripts_path""446_ENCODE_H1hesc_data_wrangler.R")
indir=$(echo "/group/soranzo/manuel.tardaguila/CH/reference_files/ENCODE/H1hesc/")
cell_type_array=$(echo 'H1hesc,K562,Hl60')

myjobid_Gather_data=$(sbatch --job-name=$name_Gather_data --output=$outfile_Gather_data --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_Gather_data --indir $indir --cell_type_array $cell_type_array --cell_type_array $cell_type_array --type $type --out $output_dir")
myjobid_seff_Gather_data=$(sbatch --dependency=afterany:$myjobid_Gather_data --open-mode=append --output=$outfile_Gather_data --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Gather_data >> $outfile_Gather_data")


