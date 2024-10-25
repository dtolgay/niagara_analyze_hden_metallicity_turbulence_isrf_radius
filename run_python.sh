#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=z0_linearNDInterpolator
#SBATCH --output=z0_linearNDInterpolator.out
#SBATCH --error=z0_linearNDInterpolator.err

module purge 
ml python/3.11.5

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius

python linearNDInterpolator_usingIline.py m12i_res7100_md zoom_in 0.0