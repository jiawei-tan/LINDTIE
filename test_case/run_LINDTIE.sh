#!/bin/bash

#SBATCH --job-name=run_LINDTIE
#SBATCH --partition=regular
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your-email@example.com
#SBATCH --output=script_output/%x_%J.out
#SBATCH --error=script_output/%x_%J.err

# Exit on the first failure. Without this, a broken environment still exits 0 and
# SLURM reports the job as successful despite having done nothing.
set -euo pipefail

# `module` is a shell function, and sbatch runs a non-login shell, so the module
# system must be initialised explicitly. Adjust for your site if needed
# (Lmod clusters use /usr/share/lmod/lmod/init/bash).
source /etc/profile.d/modules.sh

module load nextflow/25.04.2 singularity/4.1.5

# ----------------------------------------------------------------------------
# Where container images are downloaded
#
# LINDTIE runs each step inside a container, so the first run downloads 13
# container images (roughly 3 GB in total). Two separate caches are involved:
#
#   1. The finished images. On this cluster the `nextflow` module already
#      places these on scratch for you (via NXF_SINGULARITY_CACHEDIR), so
#      there is nothing to do. On other systems, set that variable yourself.
#
#   2. A temporary staging area that Singularity uses *while* downloading.
#      This one defaults to your home directory (~/.singularity/cache), which
#      usually has a small disk quota. Once it fills up, downloads fail
#      part-way through with "disk quota exceeded".
#
# The lines below move that staging area onto scratch, where there is room.
# Change the path to somewhere with plenty of free space on your own system.
# It is only needed while downloading, so it is safe to delete afterwards.
# ----------------------------------------------------------------------------
export SINGULARITY_CACHEDIR=/vast/scratch/users/$USER/nextflow/singularity_stage
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"  # for clusters that use Apptainer
mkdir -p "$SINGULARITY_CACHEDIR"

# modify the path to the LINDTIE base directory
LINDTIE_dir=/path/to/your/LINDTIE

nextflow run "$LINDTIE_dir/main.nf" -params-file "$LINDTIE_dir/params.yaml" -profile singularity
