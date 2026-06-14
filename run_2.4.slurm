#!/bin/bash -l
# Allocate slurm resources, edit as necessary
#SBATCH --account=pawsey1172
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --job-name=second_round_check
#SBATCH --array=1-8
#SBATCH --partition=work
#SBATCH --output="slurmout/2.4.second_round_check-%j-%a.out"

dir="${MYSCRATCH}"
libpath="/scratch/pawsey1172/sweng/R_Library/4.5"
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
  Rscript ~/VisHD/2.4.second_round_check.R "$SLURM_ARRAY_TASK_ID"
