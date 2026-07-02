#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --job-name=normal_scvi
#SBATCH --output="slurmout/8.2.normal_scvi-%j.out"
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

# Extra bind: keep /scratch/pawsey1172/sweng visible inside the container at
# its host path so scVIIntegration can locate
# /scratch/pawsey1172/sweng/conda/envs/scvi-tools/ as written in the script.
srun -N $SLURM_JOB_NUM_NODES -n $SLURM_NTASKS -c $SLURM_CPUS_PER_TASK \
  singularity exec -c \
  -B ${tmp_dir}:/tmp \
  -B ${dir}:$HOME \
  -B ${tmp_dir}:/var \
  -B /scratch/pawsey1172/sweng:/scratch/pawsey1172/sweng \
  ${imagename} \
  Rscript ~/VisHD/7.3.normalcell_integration_scvi.r
