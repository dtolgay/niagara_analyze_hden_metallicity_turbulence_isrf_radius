#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=1:00:00
#SBATCH --job-name=one_job_line_emissions_RBFInterpolator
#SBATCH --output=one_job_line_emissions_RBFInterpolator.out
#SBATCH --error=one_job_line_emissions_RBFInterpolator.err


module purge 
module load python/3.11.5 

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/interpolate

number_of_background_galaxies=40
redshift=0.0
target="line_emissions" # temperature, line_emissions, abundance

# Function to wait for all background processes to finish
wait_for_jobs() {
    for job in $(jobs -p)
    do
        wait $job
    done
}

time python interpolating_for_gas_particles.py gal4 firebox 0.0 line_emissions &
# time python interpolating_for_gas_particles.py gal2 firebox 0.0 line_emissions &
wait_for_jobs