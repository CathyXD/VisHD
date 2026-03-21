#!/bin/bash --login
#SBATCH --job-name=test
#SBATCH --output=test.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=0-02:00:00

export PATH=/home/sweng/SpaceRanger/spaceranger-4.0.1:$PATH
spaceranger count --id="VisiumHD_FFPE_LUT_245_20" \
          --transcriptome=/home/sweng/SpaceRanger/refdata-gex-GRCh38-2024-A \
          --probe-set=/home/sweng/SpaceRanger/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv \
          --fastqs=/pipeline/Runs/NextSeq/251120_VH01624_389_222CWN5NX/ProjectFolders/Project_Sirui-Weng/Sample_LUT-245-20 \
          --image=/home/sweng/VisHD/scan/20_40x.ome.tif \
          --cytaimage=/home/sweng/VisHD/alignment/assay_CAVG10520_2025-11-13_12-25-15_H1-CWVHCK2_1762998856_CytAssist/CAVG10520_2025-11-13_12-54-16_2025-11-13_12-25-15_H1-CWVHCK2_A1_LUT-245-20.tif \
          --slide=H1-CWVHCK2 \
          --area=B1 \
          --reorient-images=true \
          --localcores=16 \
          --localmem=128