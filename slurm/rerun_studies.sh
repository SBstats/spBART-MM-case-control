#!/bin/bash
################################################################################
# rerun_studies.sh -- rerun selected simulation studies.
#
#   bash rerun_studies.sh [--fresh] [--dry-run] job1.sbatch [job2.sbatch ...]
#
#   (default)  fill gaps only: finished results are kept, and every array task
#              skips replicates whose rep_###.RDS already exists
#   --fresh    rerun ALL replicates: the studies' existing results are first
#              MOVED (not deleted) to simulation/results_archive/<timestamp>/
#   --dry-run  show what would be moved and the submission plan; change nothing
#
# Submission goes through submit_all.sh, which sizes the arrays to the QOS
# limit on queued jobs (no throttle) and needs no install/prepare dependency.
# Find what is missing first with:  Rscript check_missing.R
################################################################################
set -o pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
FRESH=0; DRY=0; JOBS=()
for a in "$@"; do
  case "$a" in --fresh) FRESH=1 ;; --dry-run) DRY=1 ;; *) JOBS+=("$a") ;; esac
done
[ ${#JOBS[@]} -gt 0 ] || { echo "Usage: bash rerun_studies.sh [--fresh] [--dry-run] job1.sbatch ..." >&2; exit 1; }

RES="../simulation/results"
results_of() {  # result folders (relative to simulation/results) written by a job
  case "$1" in
    sim_01_TrueModel1_regular_BART.sbatch)        ls -d TrueModel1/*/*/regular_BART ;;
    sim_02_TrueModel1_spBART.sbatch)              ls -d TrueModel1/*/*/spBART ;;
    sim_03_TrueModel2_regular_BART.sbatch)        ls -d TrueModel2/*/*/regular_BART ;;
    sim_04_TrueModel2_spBART.sbatch)              ls -d TrueModel2/*/*/spBART ;;
    sim_06_TrueModel2_logistic_regression.sbatch) ls -d TrueModel2/*/*/logistic_reg TrueModel2/*/*/probit_reg ;;
    sim_07_study_sigma_beta.sbatch)               ls -d 07_sigma_beta ;;
    sim_08_study_dgm_extensions.sbatch)           ls -d 08_dgm_extensions ;;
    sim_09_study_semisynthetic.sbatch)            ls -d 09_semisynthetic ;;
    sim_10_study_full_pipeline.sbatch)            ls -d 10_full_pipeline ;;
    sim_11_study_multichain.sbatch)               ls -d 11_multichain ;;
    sim_12_study_bmi_grid.sbatch)                 ls -d 12_bmi_grid ;;
    sim_14_study_confounder_adjustment.sbatch)    ls -d 14_confounder_adjustment ;;
    sim_15_study_ultrahigh_dimension.sbatch)      ls -d 15_ultrahigh_dimension ;;
    *) return 1 ;;
  esac
}

known() { case "$1" in sim_0[1-4]_*|sim_06_*|sim_0[7-9]_*|sim_1[0-24-5]_*) return 0 ;; *) return 1 ;; esac; }
for j in "${JOBS[@]}"; do
  [ -f "$j" ] || { echo "ERROR: no job file $j" >&2; exit 1; }
  if [ "$FRESH" -eq 1 ] && ! known "$j"; then echo "ERROR: --fresh does not know the result folders of $j" >&2; exit 1; fi
done

#-- 1. --fresh: move the studies' existing results into an archive ---------------
if [ "$FRESH" -eq 1 ]; then
  STAMP=$(date +%Y%m%d_%H%M%S)
  ARCH="../simulation/results_archive/${STAMP}"
  echo "Fresh rerun: existing results of these studies go to simulation/results_archive/${STAMP}/"
  for j in "${JOBS[@]}"; do
    dirs=$(cd "$RES" 2>/dev/null && results_of "$j" 2>/dev/null)
    if [ -z "$dirs" ]; then echo "  $j: no existing results"; continue; fi
    for d in $dirs; do
      n=$(find "$RES/$d" -name 'rep_*.RDS' | wc -l | tr -d ' ')
      if [ "$DRY" -eq 1 ]; then
        echo "  would move ${d} (${n} result files)"
      else
        mkdir -p "$ARCH/$(dirname "$d")" && mv "$RES/$d" "$ARCH/$d" || { echo "ERROR: could not move $d" >&2; exit 1; }
        echo "  moved ${d} (${n} result files)"
      fi
    done
  done
fi

#-- 2. Submit through submit_all.sh (QOS-aware array sizing, no dependencies) ----
if [ "$DRY" -eq 1 ]; then
  bash submit_all.sh --dry-run "${JOBS[@]}"
else
  bash submit_all.sh "${JOBS[@]}"
fi
