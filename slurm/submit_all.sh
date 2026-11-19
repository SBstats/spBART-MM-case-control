#!/bin/bash
################################################################################
# submit_all.sh -- submit EVERY job at once, with dependencies:
#
#   00_install_BART                       (BART 2.9.9.1 + all R packages, self-test)
#     ├── 00b_prepare_analysis_base       (real X / Z for semi-synthetic + BART application)
#     │     ├── sim_07 sigma_beta grid    (includes a semi-synthetic scenario)
#     │     ├── sim_09 semi-synthetic
#     │     ├── sim_11 multi-chain        (includes a semi-synthetic scenario)
#     │     └── real_05 BART application
#     ├── sim_01-04, sim_06               (DGM 1 and DGM 2 studies)
#     ├── sim_08 further DGMs
#     ├── sim_10 full pipeline
#     ├── sim_12 BMI-effect grid
#     ├── sim_14 confounder adjustment    (8 covariates + age-proxy null genes)
#     ├── sim_15 ultra-high dimension     (p = 5000, 10000, no screening)
#     ├── real_03 spBART real data        (Section 4)
#     ├── real_04 elastic net vs spBART
#     ├── real_07 multi-chain stability   (Section 4.2)
#     └── real_08 BC-only analysis        (Sections 4.1-4.2)
#
# QOS limit on queued jobs, handled automatically (no array throttle):
#   every array task counts as one queued job. The script reads the QOS limit
#   (sacctmgr MaxSubmitPU) and the jobs you already have queued, then picks the
#   smallest REPS_PER_TASK (units of work per array task) that fits everything,
#   and sets --array and --time (unit hours x REPS_PER_TASK, capped at the QOS
#   MaxWall) for each array. Nothing is submitted if it cannot fit.
#
# Usage:   cd slurm && bash submit_all.sh
#          bash submit_all.sh --dry-run      # show the plan, submit nothing
#          bash submit_all.sh [--dry-run] job1.sbatch job2.sbatch ...
#              submit only these jobs, without install/prepare dependencies
#              (for reruns once the install and data preparation have run)
# Override the detected limit with SPBART_MAX_JOBS=<n> if sacctmgr is not readable.
# Account / QOS: the job files do not name an account. Provide your own cluster
# credentials for every submission:
#   SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh ...
# SPBART_QOS is optional; without it, also set SPBART_MAX_JOBS=<your queued-job limit>.
################################################################################
set -o pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
DRY=0; ONLY=()
for a in "$@"; do if [ "$a" = "--dry-run" ]; then DRY=1; else ONLY+=("$a"); fi; done
MARGIN=25                                   # queued-job slots kept free

FIRST=(00_install_BART.sbatch)
AFTER_INSTALL=(00b_prepare_analysis_base.sbatch
  sim_01_TrueModel1_regular_BART.sbatch sim_02_TrueModel1_spBART.sbatch
  sim_03_TrueModel2_regular_BART.sbatch sim_04_TrueModel2_spBART.sbatch
  sim_06_TrueModel2_logistic_regression.sbatch sim_08_study_dgm_extensions.sbatch
  sim_10_study_full_pipeline.sbatch sim_12_study_bmi_grid.sbatch
  sim_14_study_confounder_adjustment.sbatch sim_15_study_ultrahigh_dimension.sbatch
  real_03_spBART_real_data_analysis.sbatch real_04_comparative_ElasticNet_vs_spBART.sbatch
  real_07_multichain_stability.sbatch real_08_BConly_analysis.sbatch)
AFTER_PREP=(sim_07_study_sigma_beta.sbatch sim_09_study_semisynthetic.sbatch
  sim_11_study_multichain.sbatch real_05_BART_application.sbatch)
ALL=("${FIRST[@]}" "${AFTER_INSTALL[@]}" "${AFTER_PREP[@]}")
[ ${#ONLY[@]} -gt 0 ] && ALL=("${ONLY[@]}")

units_of() { grep -h '^# SPBART_UNITS=' "$1" | sed 's/.*SPBART_UNITS=\([0-9]*\).*/\1/'; }
hours_of() { grep -h '^# SPBART_UNITS=' "$1" | sed 's/.*SPBART_UNIT_HOURS=\([0-9]*\).*/\1/'; }
wall_to_hours() {  # D-HH:MM:SS | HH:MM:SS | minutes -> whole hours (empty if unlimited/unknown)
  local w=$1 d=0 h m s
  [ -z "$w" ] && return
  if [[ "$w" == *-* ]]; then d=${w%%-*}; w=${w#*-}; fi
  IFS=: read -r h m s <<< "$w"
  if [ -z "$m" ]; then echo $(( 10#$h / 60 )); else echo $(( 10#$d * 24 + 10#$h )); fi
}
array_tasks() { echo $(( ($(units_of "$1") + REPS - 1) / REPS )); }
array_hours() {
  local h=$(( $(hours_of "$1") * REPS ))
  [ -n "$MAXWALL_H" ] && [ "$h" -gt "$MAXWALL_H" ] && h=$MAXWALL_H
  echo "$h"
}

#-- 1. Work to schedule ---------------------------------------------------------
N_SINGLE=0; TOTAL_UNITS=0
for f in "${ALL[@]}"; do
  [ -f "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
  u=$(units_of "$f")
  if [ -z "$u" ]; then N_SINGLE=$((N_SINGLE + 1)); else TOTAL_UNITS=$((TOTAL_UNITS + u)); fi
done

#-- 2. Limits -------------------------------------------------------------------
ACCOUNT=${SPBART_ACCOUNT:-}
QOS=${SPBART_QOS:-}
if [ -z "$ACCOUNT" ]; then
  echo "ERROR: no SLURM account given. Provide your own cluster credentials, e.g." >&2
  echo "       SPBART_ACCOUNT=<your_account> SPBART_QOS=<your_qos> bash submit_all.sh" >&2
  echo "       (SPBART_QOS is optional; without it, also set SPBART_MAX_JOBS=<queued-job limit>)" >&2
  exit 1
fi
if [ -n "$QOS" ]; then
  LIMIT=${SPBART_MAX_JOBS:-$(sacctmgr -n -P show qos "$QOS" format=MaxSubmitPU 2>/dev/null | head -1)}
  MAXWALL_H=$(wall_to_hours "$(sacctmgr -n -P show qos "$QOS" format=MaxWall 2>/dev/null | head -1)")
  QUEUED=$(squeue -u "$USER" -h -r -q "$QOS" 2>/dev/null | wc -l | tr -d ' ')   # your jobs in this QOS
else
  LIMIT=${SPBART_MAX_JOBS:-}
  MAXWALL_H=""
  QUEUED=$(squeue -u "$USER" -h -r 2>/dev/null | wc -l | tr -d ' ')              # all your queued jobs
fi
if ! [[ "$LIMIT" =~ ^[0-9]+$ ]]; then
  echo "ERROR: could not read your queued-job limit${QOS:+ for QOS ${QOS}}. Run again with SPBART_MAX_JOBS=<limit>." >&2
  [ -n "$QOS" ] && echo "       (check with: sacctmgr show qos ${QOS} format=Name,MaxSubmitPU,MaxWall)" >&2
  exit 1
fi
AVAIL=$((LIMIT - QUEUED - N_SINGLE - MARGIN))
echo "Account ${ACCOUNT}, QOS ${QOS:-none}: limit ${LIMIT} queued jobs per user, max wall ${MAXWALL_H:-unlimited} h; you have ${QUEUED} queued in it"
echo "Work: ${TOTAL_UNITS} array units + ${N_SINGLE} single jobs; slots available for arrays: ${AVAIL}"
[ "$AVAIL" -gt 0 ] || { echo "ERROR: no free slots. Cancel queued jobs first." >&2; exit 1; }

#-- 3. Smallest REPS_PER_TASK that fits ------------------------------------------
REPS=1
while :; do
  N=0
  for f in "${ALL[@]}"; do [ -n "$(units_of "$f")" ] && N=$((N + $(array_tasks "$f"))); done
  [ "$N" -le "$AVAIL" ] && break
  REPS=$((REPS + 1))
done
echo "Plan: REPS_PER_TASK=${REPS} -> ${N} array tasks + ${N_SINGLE} single jobs = $((N + N_SINGLE)) queued jobs"
for f in "${ALL[@]}"; do
  [ -n "$(units_of "$f")" ] || continue
  if [ -n "$MAXWALL_H" ] && [ "$(hours_of "$f")" -gt "$MAXWALL_H" ]; then
    echo "ERROR: $f needs $(hours_of "$f") h per unit, more than the QOS max wall (${MAXWALL_H} h)." >&2; exit 1
  fi
  note=""; [ $(( $(hours_of "$f") * REPS )) -gt "$(array_hours "$f")" ] && note="(capped at QOS max wall; unfinished units resume on resubmission)"
  printf '  %-46s array 1-%-4s time %3s h %s\n' "$f" "$(array_tasks "$f")" "$(array_hours "$f")" "$note"
done
[ "$DRY" -eq 1 ] && { echo "Dry run: nothing submitted."; exit 0; }

#-- 4. Submit; stop at the first failure and say how to undo ---------------------
SUBMITTED=()
sub() {  # sub <dependency id or ""> <sbatch file>
  local dep=$1 f=$2 opts=(--account="${ACCOUNT}") id
  [ -n "$QOS" ] && opts+=(--qos="${QOS}")
  # SLURM does not create missing log folders: make the ones this job writes to
  grep -h '^#SBATCH --\(output\|error\)=' "$f" | sed 's/^#SBATCH --[a-z]*=//' | while read -r p; do mkdir -p "$(dirname "$p")"; done
  [ -n "$dep" ] && opts+=(--dependency=afterok:"$dep")
  if [ -n "$(units_of "$f")" ]; then
    opts+=(--array=1-"$(array_tasks "$f")" --time="$(array_hours "$f")":00:00 --export=ALL,REPS_PER_TASK="${REPS}")
  fi
  id=$(sbatch --parsable "${opts[@]}" "$f")
  if [ $? -ne 0 ] || ! [[ "${id%%;*}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: submitting $f failed." >&2
    [ ${#SUBMITTED[@]} -gt 0 ] && echo "Already submitted; to undo:  scancel ${SUBMITTED[*]}" >&2
    exit 1
  fi
  id=${id%%;*}; SUBMITTED+=("$id"); printf '%-48s %s\n' "$f" "$id"
}
if [ ${#ONLY[@]} -gt 0 ]; then
  for f in "${ONLY[@]}"; do sub "" "$f"; done
else
  sub "" "${FIRST[0]}"; INSTALL=${SUBMITTED[0]}
  for f in "${AFTER_INSTALL[@]}"; do sub "${INSTALL}" "$f"; done
  PREP=${SUBMITTED[1]}                                   # 00b_prepare_analysis_base
  for f in "${AFTER_PREP[@]}"; do sub "${PREP}" "$f"; done
fi
echo "All ${#SUBMITTED[@]} submissions accepted. Monitor with: squeue -u \$USER"
echo "To cancel everything from this submission:  scancel ${SUBMITTED[*]}"
