#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=z3_otherProperties
#SBATCH --output=z3_otherProperties.out
#SBATCH --error=z3_otherProperties.err

module purge 
ml python/3.11.5

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius

python nearestCloudyRun_otherProperties.py m12i_res7100_md zoom_in 3.0