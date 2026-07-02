#!/bin/bash -l
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=5:00:00
#SBATCH --job-name=integrate_archetype
#SBATCH --output="slurmout/integrate_archetype-%j.out"
#SBATCH --mail-user=Sirui.Weng@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --partition=work
#SBATCH --account=pawsey1172

conda_env="partipy"
source ${MYSOFTWARE}/miniconda3/bin/activate ${conda_env}

outdir="${MYSCRATCH}/VisHD/6.3.integrate_tumour_archetype"
mkdir -p "$outdir"

jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.timeout=14400 \
  --output-dir "$outdir" \
  "${MYSCRATCH}/VisHD/6.5.integrate_tumour_archetype.ipynb"
