#!/bin/bash
#SBATCH --job-name=TM1_spBART
#SBATCH --output=logs/TrueModel1_spBART/rep_%a.out
#SBATCH --error=logs/TrueModel1_spBART/rep_%a.err
#SBATCH --array=1-500
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=
#SBATCH --qos=


################################################################################
# 02_TrueModel1_spBART.sbatch
#
# SLURM batch script for running Semi-parametric BART simulation on True Model 1
# Each array task runs one replicate (1-500) for BOTH sample sizes (n=1000, n=1500)
#
# NOTE: spBART is slower than regular BART due to custom Gibbs sampler
#       (n_burn=1000, n_iter=1000 per fold x 5 folds x 2 sample sizes)
#
# USAGE:
#   sbatch 02_TrueModel1_spBART.sbatch
#
# MONITORING:
#   squeue -u $USER                    # Check job status
#   sacct -j <JOB_ID> --format=...     # Check completed jobs
#
# ESTIMATED TIME:
#   - Per replicate: ~2-4 hours (5 folds x 2 sample sizes, 2000 iterations each)
#   - Total wall time: 6 hours (with buffer)
#
# MEMORY:
#   - 16GB per task (sufficient for spBART with ~100-200 screened genes)
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
mkdir -p "${SCRIPT_DIR}/logs/TrueModel1_spBART"

# Change to script directory
cd "${SCRIPT_DIR}"

# Run the R script
echo "Starting R script for replicate $SLURM_ARRAY_TASK_ID..."
echo "Script directory: ${SCRIPT_DIR}"
Rscript "${SCRIPT_DIR}/02_TrueModel1_spBART.R"

# Print completion info
echo "=============================================="
echo "End time: $(date)"
echo "Exit status: $?"
echo "=============================================="
