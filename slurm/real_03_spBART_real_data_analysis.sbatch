#!/bin/bash
#SBATCH --job-name=spBART_realdata
#SBATCH --output=logs/real_data_analysis/%j.out
#SBATCH --error=logs/real_data_analysis/%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
################################################################################
# real_03_spBART_real_data_analysis.sbatch  ->  real_data_analysis/03_real_data_analysis.R
# Submit from slurm/:   cd slurm && sbatch real_03_spBART_real_data_analysis.sbatch
# Requires the modified BART package (run 00_install_BART.sbatch once first).
################################################################################
set -eo pipefail   # no -u: Lmod "module" can reference unset variables
if [ "$(basename "${SLURM_SUBMIT_DIR}")" != "slurm" ]; then
  echo "ERROR: submit from slurm/ (cd slurm && sbatch real_03_spBART_real_data_analysis.sbatch)" >&2; exit 1
fi
SPBART_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
echo "Job ${SLURM_JOB_ID} task ${SLURM_ARRAY_TASK_ID:-NA} on $(hostname) at $(date)"
module load R 2>/dev/null || echo "No R module found, using system R"   # adapt to your cluster
Rscript "${SPBART_ROOT}/real_data_analysis/03_real_data_analysis.R"
echo "Finished at $(date)"
