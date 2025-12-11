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

module load nextflow/25.04.2 singularity/4.1.5

# modify the path to the LINDTIE base directory
LINDTIE_dir=/path/to/your/LINDTIE

nextflow run $LINDTIE_dir/main.nf -params-file $LINDTIE_dir/params.yaml -profile singularity 