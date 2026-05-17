#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=00:20:00
#SBATCH --job-name=spatial_archetype
#SBATCH --output="slurmout/spatial_archetype-%j-%a.out"
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --partition=work
#SBATCH --account=pawsey1172
#SBATCH --array=1-8

conda_env="partipy"
source ${MYSOFTWARE}/miniconda3/bin/activate ${conda_env}

samples=(LUT-245-07 LUT-245-09 LUT-245-10 LUT-245-11 LUT-245-15 LUT-245-16 LUT-245-17 LUT-245-20)
sample_id=${samples[$SLURM_ARRAY_TASK_ID - 1]}

outdir="${MYSCRATCH}/VisHD/${sample_id}/tumour/archetype_result"
mkdir -p "$outdir"

export SAMPLE_IDX=$SLURM_ARRAY_TASK_ID

jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.timeout=3600 \
  --output-dir "$outdir" \
  "${MYSCRATCH}/VisHD/6.1archetypal_analysis_tumour.ipynb"
