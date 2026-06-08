#!/bin/bash -l
# Setonix Singularity submission for 5.2.integrate_DESeq2_deg.R
# Cross-sample pseudobulk DESeq2 (DT vs CB) + cohort B vs A (DT) analysis.
# Memory peak is the 8-way Seurat merge + JoinLayers; reruns hit the cache.
#SBATCH --account=pawsey1172
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH --job-name=integrate_DESeq2
#SBATCH --partition=work
#SBATCH --output="slurmout/5.2.integrate_DESeq2_deg-%j.out"

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
  Rscript ~/VisHD/5.2.integrate_DESeq2_deg.R
