#!/bin/bash -l
# Allocate slurm resources, edit as necessary
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=150G
#SBATCH --time=10:00:00
#SBATCH --job-name=infercnv
#SBATCH --array=2,8
#SBATCH --output="slurmout/3.2.infercnv_subclone-%j-%a.out"
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --partition=work

dir="${MYSCRATCH}"
# Set the image and tag we want to use
image="docker://swweng/spatial:4.5v6"

# Load Singularity. The version may change over time.
module load singularity/4.1.0-slurm

# You should not need to edit the lines below
# Prepare the working directory
mkdir -p $dir
cd ${dir}

# Get the image filename
imagename=${image##*/}
  imagename=${imagename/:/_}.sif
  
  # Create a user-specific tmp directory to avoid clashes between users
  tmp_dir="/tmp/tmp_$USER"
  mkdir -p $tmp_dir
  
  srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
  singularity exec -c \
  -B ${tmp_dir}:/tmp \
  -B ${dir}:$HOME \
  -B ${tmp_dir}:/var \
  ${imagename} \
  Rscript ~/VisHD/3.2.infercnv_subclone.R "$SLURM_ARRAY_TASK_ID"