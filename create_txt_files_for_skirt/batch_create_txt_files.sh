#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=23:00:00
#SBATCH --job-name=z2.0_create_txt_files
#SBATCH --output=z2.0_create_txt_files.out
#SBATCH --error=z2.0_create_txt_files.err

# rrg-rbond-ac 
# rrg-murray-ac

# Set NUMEXPR_MAX_THREADS to match the number of processors per node requested
export NUMEXPR_MAX_THREADS=80


# Load the venv 
module purge 
source /scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/.venv/bin/activate


cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/create_txt_files_for_skirt

number_of_background_galaxies=40

# Function to wait for all background processes to finish
wait_for_jobs() {
    for job in $(jobs -p)
    do
        wait $job
    done
}


redshift=2.0

# ######################################################## zoom_in_dtolgay
# counter=0

# # List of galaxy names
# galaxy_names=(
#     "m11d_r7100"
#     "m11e_r7100"
#     "m11h_r7100"
#     "m11i_r7100"
#     "m11q_r7100"
# )


# for galaxy in "${galaxy_names[@]}"; do
#     python create_txt_files_for_skirt.py $galaxy zoom_in_tolgay 99 $redshift &

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


# ######################################################## zoom_in
# # Counter for every $number_of_background_galaxies galaxies
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
# )


# for galaxy in "${galaxy_names[@]}"; do
#     python create_txt_files_for_skirt.py $galaxy zoom_in 99 $redshift &

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



######################################################## firebox
# Counter for every $number_of_background_galaxies galaxies
counter=0

for i in {1..1000}
do
    # python create_txt_files_for_skirt_average_sobolev_smoothing_length.py dummy firebox $i &
    python create_txt_files_for_skirt.py dummy firebox $i $redshift &
    

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


# ######################################################## particle_split
# # Counter for every $number_of_background_galaxies galaxies
# counter=0

# # List of galaxy names
# galaxy_names=(
#     "m12i_r880_md" 
# )


# for galaxy in "${galaxy_names[@]}"; do
#     # python create_txt_files_for_skirt_average_sobolev_smoothing_length.py $galaxy particle_split 99 &
#     python create_txt_files_for_skirt.py $galaxy particle_split 99 $redshift &

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