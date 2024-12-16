#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=02:00:00
#SBATCH --job-name=smoothingLength_gal0
#SBATCH --output=smoothingLength_gal0.out
#SBATCH --error=smoothingLength_gal0.err

module purge 
ml python/3.11.5

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/hybridInterpolator_smoothingLength

export PYTHONUNBUFFERED=1 # this is the write .out file while the code is running 
python hybridInterpolator_usingFline.py gal0 firebox 0.0 40