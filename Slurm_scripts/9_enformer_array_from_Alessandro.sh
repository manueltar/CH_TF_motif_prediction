#!/bin/bash
#SBATCH --job-name=enformer_run
#SBATCH --partition=gpuq
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --output=%x_%A_%a.log
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --array=1-1529%50

echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

echo "init conda"
eval "$(conda shell.bash hook)"

module load cuda/11.5.1 cudnn/8.3.1.22-11.5-cuda-11.5.1

conda activate enformer_20231026

END=$SLURM_ARRAY_TASK_ID
# # END=600
# START=$(( $END - 199 ))

# input and output folder
INPUT_FOLDER="/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/Enformer/input"
OUTPUT_FOLDER="/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/Enformer/output"
SCRIPT_FOLDER="/ssu/gassu/GAU_tools/enformer_variant_score"
REF_FOLDER="/processing_data/reference_datasets/iGenomes/2023.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta"
LOG_FOLDER="/group/soranzo/manuel.tardaguila/CH/ALL_variants_in_CS/Enformer/Log_files/"

cd ${LOG_FOLDER}

# input variant is a file with one line per variant in GRCh38 in the chr_pos_ref_alt format, the array works for 200 variants per file

input_variant=$( ls ${INPUT_FOLDER}/input_list* | tail -n+${END} | head -1 )
echo $input_variant
# create the output folder
mkdir ${OUTPUT_FOLDER}/out_enformer_${END}

# launch the script
python ${SCRIPT_FOLDER}/enformer_predict.py --variants ${input_variant} --output ${OUTPUT_FOLDER}/out_enformer_${END} --ref_genome ${REF_FOLDER}/genome.fa | tee ${LOG_FOLDER}/enformer_${END}.log

echo "========================"
echo "Completed: $(date)"
