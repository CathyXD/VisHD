#!/bin/bash --login
#SBATCH --job-name=vishd
#SBATCH --output=slurmout/0.generate_%a_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --array=2-3,6-7
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=08:00:00

module load ood_r/4.5.0
Rscript ~/VisHD/0.generate_srt.R