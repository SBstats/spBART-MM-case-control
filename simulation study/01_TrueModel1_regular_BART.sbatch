#!/bin/bash
#SBATCH --job-name=TM1_regBART
#SBATCH --output=logs/TrueModel1_regular_BART/rep_%a.out
#SBATCH --error=logs/TrueModel1_regular_BART/rep_%a.err
#SBATCH --array=1-500
#SBATCH --time=05:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=
#SBATCH --qos=


################################################################################
# 01_TrueModel1_regular_BART.sbatch
#
# SLURM batch script for running Regular BART simulation on True Model 1
# Each array task runs one replicate (1-500) for BOTH sample sizes (n=1000, n=1500)
#
# USAGE:
#   sbatch 01_TrueModel1_regular_BART.sbatch
#
# MONITORING:
#   squeue -u $USER                    # Check job status
#   sacct -j <JOB_ID> --format=...     # Check completed jobs
#
# ESTIMATED TIME:
#   - Per replicate: ~1-2 hours (5 folds x 2 sample sizes)
#   - Total wall time: 4 hours (with buffer)
#
# MEMORY:
#   - 8GB per task (sufficient for probit BART with ~100-200 screened genes)
#
################################################################################

# Print job info
echo "=============================================="
echo "SLURM Job ID: $SLURM_JOB_ID"
echo "SLURM Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "SLURM Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Running on host: $(hostname)"
echo "Working directory: $(pwd)"
echo "Start time: $(date)"
echo "=============================================="

# Load R module (adjust based on your cluster's module system)
module load R/4.3.1 2>/dev/null || module load R 2>/dev/null || echo "No R module found, using system R"

# Use SLURM_SUBMIT_DIR - the directory from which sbatch was called
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"

# Create logs subdirectory if it doesn't exist
mkdir -p "${SCRIPT_DIR}/logs/TrueModel1_regular_BART"

# Change to script directory
cd "${SCRIPT_DIR}"

# Run the R script
echo "Starting R script for replicate $SLURM_ARRAY_TASK_ID..."
echo "Script directory: ${SCRIPT_DIR}"
Rscript "${SCRIPT_DIR}/01_TrueModel1_regular_BART.R"

# Print completion info
echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
