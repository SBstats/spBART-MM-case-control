#!/bin/bash
#SBATCH --job-name=S08_dgm_ext
#SBATCH --output=logs/08_dgm_extensions/task_%a.out
#SBATCH --error=logs/08_dgm_extensions/task_%a.err
#SBATCH --array=1-800
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
# SPBART_UNITS=800 SPBART_UNIT_HOURS=8
################################################################################
# sim_08_study_dgm_extensions.sbatch  ->  simulation/08_study_dgm_extensions.R
#
# 800 units of work (replicates / settings); a unit needs at most ~8 h.
# Submit through submit_all.sh, which reads the QOS limit on queued jobs and
# sets --array, --time and REPS_PER_TASK so everything fits (no throttle).
# Each array task runs REPS_PER_TASK units one after another; the R script sees
# each unit as SLURM_ARRAY_TASK_ID. Finished units are skipped, so resubmitting
# resumes where a task stopped. Direct "sbatch sim_08_study_dgm_extensions.sbatch" runs 1 unit per task.
################################################################################
set -o pipefail    # no -e/-u: one failed unit must not stop the others; Lmod "module" can reference unset variables
if [ "$(basename "${SLURM_SUBMIT_DIR}")" != "slurm" ]; then
  echo "ERROR: submit from slurm/" >&2; exit 1
fi
SPBART_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
N_UNITS=800
REPS_PER_TASK=${REPS_PER_TASK:-1}
echo "Job ${SLURM_JOB_ID} array task ${SLURM_ARRAY_TASK_ID} (${REPS_PER_TASK} units) on $(hostname) at $(date)"
module load R 2>/dev/null || echo "No R module found, using system R"   # adapt to your cluster
failed=0
for k in $(seq 0 $((REPS_PER_TASK - 1))); do
  UNIT=$(( (SLURM_ARRAY_TASK_ID - 1) * REPS_PER_TASK + k + 1 ))
  [ "${UNIT}" -le "${N_UNITS}" ] || continue
  echo "=== unit ${UNIT} of ${N_UNITS} (started $(date)) ==="
  SLURM_ARRAY_TASK_ID=${UNIT} Rscript "${SPBART_ROOT}/simulation/08_study_dgm_extensions.R" || { echo "unit ${UNIT} FAILED"; failed=1; }
done
echo "Finished at $(date)"
exit ${failed}
