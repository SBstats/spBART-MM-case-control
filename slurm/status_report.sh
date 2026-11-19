#!/bin/bash
################################################################################
# status_report.sh -- one report of what ran, what failed and what is stuck.
# Writes status_report.txt (useful when asking for help). Changes nothing.
#
# Usage:  cd slurm && bash status_report.sh [start date, default: 3 days ago]
################################################################################
cd "$(dirname "${BASH_SOURCE[0]}")"
SINCE=${1:-$(date -d '3 days ago' +%Y-%m-%d 2>/dev/null || date -v-3d +%Y-%m-%d)}
OUT=status_report.txt
NAMES="install_BART,prep_base,TM1_regBART,TM1_spBART,TM2_regBART,TM2_spBART,TM2_logistic,S07_sigma2,S08_dgm_ext,S09_semisyn,S10_pipeline,S11_multichain,S12_bmi_grid,S14_confound,S15_ultrahigh,spBART_realdata,spBART_EN_comp,BART_realdata,spBART_multichain,spBART_BConly"

newest() { ls -t "$1"/*."$2" 2>/dev/null | head -1; }
show_tail() { [ -n "$1" ] && [ -s "$1" ] && { echo "--- $1 (last ${2:-40} lines)"; tail -n "${2:-40}" "$1"; }; }

{
echo "=== spBART status report: $(date) on $(hostname); jobs since ${SINCE} ==="

echo; echo "### 1. Queue now (job name / state / reason: number of jobs)"
squeue -u "$USER" -h -r -o "%j|%T|%r" | sort | uniq -c | sort -k2

echo; echo "### 2. Finished or running since ${SINCE} (job name / state: number of jobs or array tasks)"
sacct -u "$USER" -S "$SINCE" -X -n -P --name="$NAMES" --format=JobName,State | \
  sed 's/CANCELLED by [0-9]*/CANCELLED/' | sort | uniq -c | sort -k2

echo; echo "### 3. Install and data-preparation jobs"
sacct -u "$USER" -S "$SINCE" -X --name=install_BART,prep_base \
  --format=JobID%14,JobName%14,State%22,ExitCode,Elapsed,MaxRSS,Start,End
for d in install prepare_base; do
  show_tail "$(newest logs/$d out)" 40
  show_tail "$(newest logs/$d err)" 40
done

echo; echo "### 4. First failed array tasks (up to 2 per job name)"
sacct -u "$USER" -S "$SINCE" -X -n -P --name="$NAMES" --format=JobID,JobName,State,ExitCode | \
  awk -F'|' '$3 ~ /FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL/' | \
  awk -F'|' '{ if (++n[$2] <= 2) print }' | while IFS='|' read -r id name state code; do
    echo "--- ${name} ${id}: ${state} (exit ${code})"
    task=${id##*_}
    for f in $(grep -l "job-name=${name}$" *.sbatch 2>/dev/null); do
      dir=$(grep -h '^#SBATCH --output=' "$f" | sed 's#.*--output=logs/\([^/]*\)/.*#\1#')
      show_tail "logs/${dir}/task_${task}.err" 25
    done
  done

echo; echo "### 5. Result files written so far"
for d in ../simulation/results/*/ ; do
  [ -d "$d" ] && printf '%-60s %s RDS files\n' "$d" "$(find "$d" -name 'rep_*.RDS' | wc -l)"
done
ls -d ../real_data_analysis/results_* 2>/dev/null

echo; echo "### 6. Jobs stuck on a dependency that can never be met (cancel after fixing)"
STUCK=$(squeue -u "$USER" -h -o "%A %r" | awk '$2 ~ /DependencyNeverSatisfied/ {print $1}' | sort -u | tr '\n' ' ')
echo "${STUCK:-none}"
[ -n "$STUCK" ] && echo "To cancel them:  scancel ${STUCK}"
} > "$OUT" 2>&1

cat "$OUT"
echo; echo "Saved to $(pwd)/${OUT}"
