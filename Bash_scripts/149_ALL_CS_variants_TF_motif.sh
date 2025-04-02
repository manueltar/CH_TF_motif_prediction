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

#### Gather_and_filter ###################################

type=$(echo "Gather_and_filter""_""$analysis")
outfile_Gather_and_filter=$(echo "$Log_files""outfile_1_""$type"".out")
touch $outfile_Gather_and_filter
echo -n "" > $outfile_Gather_and_filter
name_Gather_and_filter=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Gather_and_filter=$(echo "$Rscripts_path""443_TF_motif_CS_gather_and_filter_variants.R")
indir=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/")
Threshold_PP=$(echo '0.001')
GWAS_equivalence=$(echo "/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/hcoloc_submission_table_CHIP.tsv")

myjobid_Gather_and_filter=$(sbatch --job-name=$name_Gather_and_filter --output=$outfile_Gather_and_filter --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Gather_and_filter --indir $indir --Threshold_PP $Threshold_PP --GWAS_equivalence $GWAS_equivalence --type $type --out $output_dir")
myjobid_seff_Gather_and_filter=$(sbatch --dependency=afterany:$myjobid_Gather_and_filter --open-mode=append --output=$outfile_Gather_and_filter --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Gather_and_filter >> $outfile_Gather_and_filter")


#### bedtools_getfasta #############################


type=$(echo "bedtools_getfasta""_""$analysis")
outfile_bedtools_getfasta=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_bedtools_getfasta
echo -n "" > $outfile_bedtools_getfasta
name_bedtools_getfasta=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

module load bedtools2/2.31.0

reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2022.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")
input_bed=$(echo "$output_dir""VARS"".bed")
output_fasta=$(echo "$output_dir""intervals"".fasta")

myjobid_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_seff_Gather_and_filter --job-name=$name_bedtools_getfasta --output=$outfile_bedtools_getfasta --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=512M --parsable --wrap="bedtools getfasta -fi $reference_genome -bed $input_bed > $output_fasta")
myjobid_seff_bedtools_getfasta=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --open-mode=append --output=$outfile_bedtools_getfasta --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bedtools_getfasta >> $outfile_bedtools_getfasta")

#### printer_fasta_files #############################


type=$(echo "printer_fasta_files""_""$analysis")
outfile_printer_fasta_files=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_printer_fasta_files
echo -n "" > $outfile_printer_fasta_files
name_printer_fasta_files=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_printer_fasta_files=$(echo "$Rscripts_path""444_REF_ALT_assignation.R")

input_bed=$(echo "$output_dir""VARS"".bed")
input_fasta=$(echo "$output_dir""intervals"".fasta")
GWAS_file=$(echo "$output_dir""config_file.tsv")

# --dependency=afterany:$myjobid_bedtools_getfasta

myjobid_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_bedtools_getfasta --job-name=$name_printer_fasta_files --output=$outfile_printer_fasta_files --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_printer_fasta_files --input_bed $input_bed --input_fasta $input_fasta --GWAS_file $GWAS_file --type $type --out $output_dir")
myjobid_seff_printer_fasta_files=$(sbatch --dependency=afterany:$myjobid_printer_fasta_files --open-mode=append --output=$outfile_printer_fasta_files --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_printer_fasta_files >> $outfile_printer_fasta_files")

