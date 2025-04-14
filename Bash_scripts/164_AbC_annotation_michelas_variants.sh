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
 
#######################################################################binder_of_scores_AbC ###################################


type=$(echo "binder_of_scores_AbC""_""$analysis")
outfile_binder_of_scores_AbC=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_binder_of_scores_AbC
echo -n "" > $outfile_binder_of_scores_AbC
name_binder_of_scores_AbC=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_binder_of_scores_AbC=$(echo "$Rscripts_path""457_AbC_annotation_Michelas_variants.R")



Table_of_variants=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/config_file_assigned_ref_and_alt.tsv")
AbC_scores=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/AbC_Engreitz/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz")
AbC_cell_types=$(echo 'K562-Roadmap,H1-hESC-Roadmap,erythroblast-Corces2016,CD8-positive_alpha-beta_T_cell-Corces2016,CD4-positive_helper_T_cell-Corces2016,CD14-positive_monocyte_treated_with_LPS_d6-Novakovic2016,CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016,CD14-positive_monocyte_treated_with_RPMI_4h-Novakovic2016,CD14-positive_monocyte_treated_with_LPS_1h-Novakovic2016,CD14-positive_monocyte_treated_with_BG_4h-Novakovic2016,CD14-positive_monocyte_treated_with_RPMI_d1-Novakovic2016,CD14-positive_monocyte_treated_with_BG_d1-Novakovic2016,natural_killer_cell-Corces2016,CD4-positive_helper_T_cell-ENCODE,CD56-positive_natural_killer_cells-Roadmap,CD14-positive_monocyte_treated_with_BG_1h-Novakovic2016,CD19-positive_B_cell-Roadmap,CD3-positive_T_cell-Roadmap,CD14-positive_monocyte-Novakovic2016,CD14-positive_monocyte_treated_with_RPMI_1h-Novakovic2016,B_cell-ENCODE,CD14-positive_monocyte-ENCODE,CD8-positive_alpha-beta_T_cell-ENCODE,T-cell-ENCODE,CD14-positive_monocyte_treated_with_LPS_d1-Novakovic2016,dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_6_hour-Garber2017,CD14-positive_monocytes-Roadmap,dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_4_hour-Garber2017,CD34-positive_mobilized-Roadmap,dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_30_minute-Garber2017,dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_1_hour-Garber2017,dendritic_cell_treated_with_Lipopolysaccharide_0_ng-mL_for_0_hour-Garber2017,megakaryocyte-erythroid_progenitor-Corces2016,dendritic_cell_treated_with_Lipopolysaccharide_100_ng-mL_for_2_hour-Garber2017,CD14-positive_monocyte_treated_with_RPMI_d6-Novakovic2016,CD14-positive_monocyte_treated_with_BG_d6-Novakovic2016')
Threshold_AbC=$(echo "0.01")


myjobid_binder_of_scores_AbC=$(sbatch --job-name=$name_binder_of_scores_AbC --output=$outfile_binder_of_scores_AbC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_binder_of_scores_AbC --Table_of_variants $Table_of_variants --AbC_scores $AbC_scores --AbC_cell_types $AbC_cell_types --Threshold_AbC $Threshold_AbC --type $type --out $output_dir")
myjobid_seff_binder_of_scores_AbC=$(sbatch --dependency=afterany:$myjobid_binder_of_scores_AbC --open-mode=append --output=$outfile_binder_of_scores_AbC --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_binder_of_scores_AbC >> $outfile_binder_of_scores_AbC")
 

