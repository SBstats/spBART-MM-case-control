#!/bin/bash
#SBATCH --job-name=TM2_logistic
#SBATCH --output=logs/TrueModel2_logistic_reg/rep_%a.out
#SBATCH --error=logs/TrueModel2_logistic_reg/rep_%a.err
#SBATCH --array=1-500
#SBATCH --time=05:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=
#SBATCH --qos=

################################################################################
# 06_TrueModel2_logistic_regression.sbatch
#
# SLURM batch script for running Logistic Regression on True Model 2
# Each array task runs one replicate (1-500) for ALL n,p combinations
#
# USAGE:
#   sbatch 06_TrueModel2_logistic_regression.sbatch
#
# MONITORING:
#   squeue -u $USER                    # Check job status
#   sacct -j <JOB_ID> --format=...     # Check completed jobs
#
# ESTIMATED TIME:
#   - Per replicate: ~15-30 minutes (3 p-values x 2 sample sizes)
#   - Total wall time: 2 hours (with buffer)
#
# MEMORY:
#   - 16GB per task (sufficient for logistic regression with up to 5000 genes)
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
mkdir -p "${SCRIPT_DIR}/logs/TrueModel2_logistic_reg"

# Change to script directory
cd "${SCRIPT_DIR}"

# Run the R script
echo "Starting R script for replicate $SLURM_ARRAY_TASK_ID..."
echo "Script directory: ${SCRIPT_DIR}"
Rscript "${SCRIPT_DIR}/06_TrueModel2_logistic_regression.R"

# Print completion info
echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
