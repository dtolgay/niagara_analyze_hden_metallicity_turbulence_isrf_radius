#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=12:00:00
#SBATCH --job-name=smoothingLength
#SBATCH --output=smoothingLength.out
#SBATCH --error=smoothingLength.err


module purge 
module load python/3.11.5 

cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/hybridInterpolator_smoothingLength

number_of_background_galaxies=1
redshift=0.0
number_of_processors_per_galaxy=40

# Function to wait for all background processes to finish
wait_for_jobs() {
    for job in $(jobs -p)
    do
        wait $job
    done
}


## Single galaxy 
# python hybridInterpolator_usingFline.py m12i_res7100_md particle_split 0.0 40

# ####### firebox
# # Counter for every 10 galaxies
# counter=0

# # for i in {600..999}; do
# for i in {30..35}; do
#     python hybridInterpolator_usingFline.py gal$i firebox $redshift $number_of_processors_per_galaxy &

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
#     python hybridInterpolator_usingFline.py $galaxy zoom_in $redshift $number_of_processors_per_galaxy &

#     # Increment counter
#     ((counter++))

#     if [ $counter -ge $number_of_background_galaxies ]; then
#         wait_for_jobs
#         counter=0
#     fi
# done

# # Wait for the last set of jobs to finish
# wait_for_jobs


# ####### particle_split
# galaxy_names=(
#     "m12i_r880_md"    
# )


# for galaxy in "${galaxy_names[@]}"; do
#     python hybridInterpolator_usingFline.py $galaxy particle_split $redshift $number_of_processors_per_galaxy &

#     # Increment counter
#     ((counter++))

#     if [ $counter -ge $number_of_background_galaxies ]; then
#         wait_for_jobs
#         counter=0
#     fi
# done

# # Wait for the last set of jobs to finish
# wait_for_jobs