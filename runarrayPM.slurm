#!/bin/bash --login
#SBATCH --job-name=ataclone
#SBATCH --output="slurmout/2.ATAClone-%j-%a.out"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=10:20:00
#SBATCH --array=1-8

module load apptainer
dir="/home/sweng"
image="docker://swweng/spatial:4.5v6"
cd ${dir}
imagename=${image##*/}
imagename=${imagename/:/_}.sif
  

# Create a user-specific tmp directory to avoid clashes between users
tmp_dir="/tmp/tmp_$USER"
mkdir -p $tmp_dir

srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
apptainer exec -c \
-B ${tmp_dir}:/tmp \
-B ${dir}:$HOME \
-B ${tmp_dir}:/var \
${imagename} \
Rscript ~/VisHD/2.infercnv_ATAClone.R "$SLURM_ARRAY_TASK_ID"