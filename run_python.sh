#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=z0_linearNDInterpolator_gal10
#SBATCH --output=z0_linearNDInterpolator_gal10.out
#SBATCH --error=z0_linearNDInterpolator_gal10.err

module purge 
ml python/3.11.5

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius

export PYTHONUNBUFFERED=1 # this is the write .out file while the code is running 
python linearNDInterpolator_usingIline.py gal10 firebox 0.0 40