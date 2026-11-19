#!/bin/bash
#SBATCH --job-name=BART_realdata
#SBATCH --output=logs/BART_application/%j.out
#SBATCH --error=logs/BART_application/%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
################################################################################
# real_05_BART_application.sbatch  ->  real_data_analysis/05_BART_application.R
# Submit from slurm/:   cd slurm && sbatch real_05_BART_application.sbatch
# Requires the modified BART package (run 00_install_BART.sbatch once first).
################################################################################
set -eo pipefail   # no -u: Lmod "module" can reference unset variables
if [ "$(basename "${SLURM_SUBMIT_DIR}")" != "slurm" ]; then
  echo "ERROR: submit from slurm/ (cd slurm && sbatch real_05_BART_application.sbatch)" >&2; exit 1
fi
SPBART_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
echo "Job ${SLURM_JOB_ID} task ${SLURM_ARRAY_TASK_ID:-NA} on $(hostname) at $(date)"
module load R 2>/dev/null || echo "No R module found, using system R"   # adapt to your cluster
Rscript "${SPBART_ROOT}/real_data_analysis/05_BART_application.R"
echo "Finished at $(date)"
