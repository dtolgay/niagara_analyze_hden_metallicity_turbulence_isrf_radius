#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=metal
#SBATCH --output=metal.out
#SBATCH --error=metal.err

module purge 
ml python/3.11.5

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/hybridInterpolator_debugging_z3


python hybridInterpolator_usingFline_metal_mult_10^1.48.py m12i_res7100_md zoom_in 3.0 40