#!/bin/bash --login
#SBATCH --job-name=ataclone
#SBATCH --output="slurmout/2.ATAClone-%j-%a.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=03:00:00
#SBATCH --array=8

module load ood_r/4.5.0
Rscript ~/VisHD/2.ATAClone.R "$SLURM_ARRAY_TASK_ID"