#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=hybridInterp_line_m12i
#SBATCH --output=hybridInterp_line_m12i.out
#SBATCH --error=hybridInterp_line_m12i.err

module purge 
ml python/3.11.5

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius

export PYTHONUNBUFFERED=1 # this is the write .out file while the code is running 
python hybridInterpolator_usingFline.py m12i_res7100_md zoom_in 0.0 40