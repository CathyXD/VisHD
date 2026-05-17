#!/bin/bash --login
#SBATCH --job-name=split
#SBATCH --output="slurmout/3.split-%j-%a.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=3:00:00
#SBATCH --array=1-8

module load ood_r/4.5.0
Rscript ~/VisHD/3.1.tumour_normal_split.R "$SLURM_ARRAY_TASK_ID"