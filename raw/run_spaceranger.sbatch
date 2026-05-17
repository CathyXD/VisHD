#!/bin/bash --login
#SBATCH --job-name=test
#SBATCH --output=test.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=32
#SBATCH --mem=50G
#SBATCH --time=0-04:00:00

export PATH=/home/sweng/SpaceRanger/spaceranger-4.1.0:$PATH
          
          
spaceranger count --id="LUT-245-17" \
          --transcriptome=/home/sweng/SpaceRanger/refdata-gex-GRCh38-2024-A \
          --probe-set=/home/sweng/SpaceRanger/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv \
          --fastqs=/pipeline/Runs/NextSeq/260217_VH01624_417_222HNKLNX/ProjectFolders/Project_Sirui-Weng/Sample_LUT-245-17 \
          --image=/home/sweng/VisHD/raw/images/LUT-245-17_40x.ome.tif \
          --cytaimage=/home/sweng/VisHD/raw/scans/CAVG10520_2026-02-12_15-02-31_2026-02-12_14-38-00_H1-JFJRPH4_D1_LUT-245-17.tif \
          --slide=H1-JFJRPH4 \
          --area=D1 \
          --create-bam=false \
          --loupe-alignment=/home/sweng/VisHD/raw/alignment/H1-JFJRPH4-D1-fiducials-image-registration.json \
          --localcores=32 \
          --localmem=50
          
spaceranger count --id="LUT-245-10" \
          --transcriptome=/home/sweng/SpaceRanger/refdata-gex-GRCh38-2024-A \
          --probe-set=/home/sweng/SpaceRanger/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv \
          --fastqs=/pipeline/Runs/NextSeq/260217_VH01624_417_222HNKLNX/ProjectFolders/Project_Sirui-Weng/Sample_LUT-245-10 \
          --image=/home/sweng/VisHD/raw/images/LUT-245-10_40x.ome.tif \
          --cytaimage=/home/sweng/VisHD/raw/scans/CAVG10520_2026-02-12_15-02-31_2026-02-12_14-38-00_H1-JFJRPH4_A1_LUT-245-10.tif \
          --slide=H1-JFJRPH4 \
          --area=A1 \
          --create-bam=false \
          --loupe-alignment=/home/sweng/VisHD/raw/alignment/H1-JFJRPH4-A1-fiducials-image-registration.json \
          --localcores=32 \
          --localmem=50
          