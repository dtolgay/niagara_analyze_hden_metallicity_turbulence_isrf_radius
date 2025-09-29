#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=20:00:00
#SBATCH --job-name=test_particles
#SBATCH --output=test_particles.out
#SBATCH --error=test_particles.err


module purge 
# Load the venv 
source ~/.venv_all/bin/activate

# Go to the working directory
cd /scratch/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/interpolate

# Run Python 

python interpolating_for_gas_particles.py temperature &
python interpolating_for_gas_particles.py hden &

wait