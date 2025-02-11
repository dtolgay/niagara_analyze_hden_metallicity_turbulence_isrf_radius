#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=20:00:00
#SBATCH --job-name=line_emissions_regularGridInterpolator
#SBATCH --output=line_emissions_regularGridInterpolator.out
#SBATCH --error=line_emissions_regularGridInterpolator.err


module purge 
module load python/3.11.5 

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/regularGridInterpolator

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

# ####### firebox
# # Counter for every 10 galaxies
# counter=0

# for i in {0..1000}; do
# # for i in {0..51}; do
#     python interpolating_for_gas_particles.py gal$i firebox $redshift $target &

#     # Increment counter
#     ((counter++))

#     if [ $counter -ge $number_of_background_galaxies ]; then
#         wait_for_jobs
#         counter=0
#     fi
# done

# # Wait for the last set of jobs to finish
# wait_for_jobs

# ####### zoom_in
# # Counter for every 10 galaxies
# counter=0

# # List of galaxy names
# galaxy_names=(
#     "m12b_res7100_md" 
#     "m12c_res7100_md"
#     "m12f_res7100_md"
#     "m12i_res7100_md"
#     "m12m_res7100_md"
#     "m12q_res7100_md"
#     "m12r_res7100_md"
#     "m12w_res7100_md"
#     "m11d_r7100"
#     "m11e_r7100"
#     "m11h_r7100"
#     "m11i_r7100"
#     "m11q_r7100"    
# )


# for galaxy in "${galaxy_names[@]}"; do
#     python interpolating_for_gas_particles.py $galaxy zoom_in $redshift $target &

#     # Increment counter
#     ((counter++))

#     if [ $counter -ge $number_of_background_galaxies ]; then
#         wait_for_jobs
#         counter=0
#     fi
# done

# # Wait for the last set of jobs to finish
# wait_for_jobs


####### particle_split
galaxy_names=(
    "m12i_r880_md"    
)


for galaxy in "${galaxy_names[@]}"; do
    python interpolating_for_gas_particles.py $galaxy particle_split $redshift $target &

    # Increment counter
    ((counter++))

    if [ $counter -ge $number_of_background_galaxies ]; then
        wait_for_jobs
        counter=0
    fi
done

# Wait for the last set of jobs to finish
wait_for_jobs
