#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=23:00:00
#SBATCH --job-name=one_job_RBFInterpolator
#SBATCH --output=one_job_RBFInterpolator.out
#SBATCH --error=one_job_RBFInterpolator.err


module purge 
module load python/3.11.5 

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/interpolate

number_of_background_galaxies=40
redshift=0.0

# Start the python job and write the output real time to output file 
stdbuf -oL -eL python interpolating_for_gas_particles_delete.py m12i_res7100_md zoom_in $redshift temperature > temperature_output.log 2>&1 &
stdbuf -oL -eL python interpolating_for_gas_particles_delete.py m12i_res7100_md zoom_in $redshift abundance > abundance_output.log 2>&1 &
stdbuf -oL -eL python interpolating_for_gas_particles_delete.py m12i_res7100_md zoom_in $redshift line_emissions > line_emissions_output.log 2>&1 &

wait