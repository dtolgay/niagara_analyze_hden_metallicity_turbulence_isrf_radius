#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=particleSplit_nearestCloudyRun_usingIline
#SBATCH --output=particleSplit_nearestCloudyRun_usingIline.out
#SBATCH --error=particleSplit_nearestCloudyRun_usingIline.err


module purge 
module load python/3.11.5 

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius

number_of_background_galaxies=1
redshift=0.0

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

# for i in {600..999}; do
#     python nearestCloudyRun_usingIline.py gal$i firebox $redshift &

#     # Increment counter
#     ((counter++))

#     # Every 10th galaxy, wait for all background jobs to finish
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
#     python nearestCloudyRun_usingIline.py $galaxy zoom_in $redshift &

#     # Increment counter
#     ((counter++))

#     # Every 10th galaxy, wait for all background jobs to finish
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
    python nearestCloudyRun_usingIline.py $galaxy particle_split $redshift &

    # Increment counter
    ((counter++))

    # Every 10th galaxy, wait for all background jobs to finish
    if [ $counter -ge $number_of_background_galaxies ]; then
        wait_for_jobs
        counter=0
    fi
done

# Wait for the last set of jobs to finish
wait_for_jobs
