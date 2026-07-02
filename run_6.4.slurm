#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=03:00:00
#SBATCH --job-name=signature_analysis
#SBATCH --output="slurmout/6.4.signature_analysis-%j.out"
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --account=pawsey1172
#SBATCH --partition=work

dir="${MYSCRATCH}"
image="docker://swweng/spatial:4.5v6"

module load singularity/4.1.0-slurm

mkdir -p $dir
cd ${dir}

imagename=${image##*/}
imagename=${imagename/:/_}.sif

tmp_dir="/tmp/tmp_$USER"
mkdir -p $tmp_dir

srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
  singularity exec -c \
  -B ${tmp_dir}:/tmp \
  -B ${dir}:$HOME \
  -B ${tmp_dir}:/var \
  ${imagename} \
  Rscript ~/VisHD/6.4.signature_analysis.R
