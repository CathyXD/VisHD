#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=1:00:00
#SBATCH --job-name=archetype_tumour
#SBATCH --array=1,3,6
#SBATCH --output="slurmout/6.1.archetypal_analysis_tumour-%j-%a.out"
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --partition=work
#SBATCH --account=pawsey1172

conda_env="partipy"
source ${MYSOFTWARE}/miniconda3/bin/activate ${conda_env}

# Per-sample selection is driven by the SAMPLE_IDX env var read inside the notebook
export SAMPLE_IDX=${SLURM_ARRAY_TASK_ID}

# Executed-notebook copies go to a per-task name so array tasks don't collide
exec_dir="${MYSCRATCH}/VisHD/notebooks_executed"
mkdir -p "$exec_dir"

jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.timeout=14400 \
  --output "6.1.archetypal_analysis_tumour_task${SLURM_ARRAY_TASK_ID}.ipynb" \
  --output-dir "$exec_dir" \
  "${MYSCRATCH}/VisHD/6.1.archetypal_analysis_tumour.ipynb"
