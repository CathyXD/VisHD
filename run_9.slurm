#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=5:00:00
#SBATCH --job-name=scimilarity_anno
#SBATCH --output="slurmout/9.normalcell_annotation-%j.out"
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --partition=work
#SBATCH --account=pawsey1172

conda_env="scimilarity"
source ${MYSOFTWARE}/miniconda3/bin/activate ${conda_env}

python "${MYSCRATCH}/VisHD/9.1.normalcell_annotation.py"
