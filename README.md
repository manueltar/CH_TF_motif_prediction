# To get the TF predictions run

$ bash ~/Scripts/Wraper_scripts/140_TF_motif_prediction_michela_v3.sh /group/soranzo/manuel.tardaguila/CH/ TF_motif_analysis 15 15

# To run the TF occupancy finder

1. Download the peaks file from the ChiP_ATLAS: https://github.com/inutano/chip-atlas/wiki#downloads_doc

$ wget https://chip-atlas.dbcls.jp/data/hg19/allPeaks_light/allPeaks_light.hg19.10.bed.gz

2. Run

$ nohup bash ~/Scripts/Wraper_scripts/141_TF_occupancy.sh /group/soranzo/manuel.tardaguila/CH/TF_motif_analysis/ TF_occupancy_filter &

