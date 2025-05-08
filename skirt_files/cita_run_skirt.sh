#!/bin/bash -l
#PBS -l nodes=1:ppn=128
#PBS -l walltime=23:50:00
#PBS -N z1_zoom_in
#PBS -o z1_zoom_in.out
#PBS -e z1_zoom_in.err
#PBS -q starq 
#PBS -r n
#PBS -j oe


redshift=1.0
directory_name=voronoi_1e5


# Function to wait for the number of running jobs to drop below a limit
wait_for_jobs() {
    local max_jobs=$1
    while : ; do
        # Count the number of running jobs
        local current_jobs=$(jobs -rp | wc -l)
        # If the number of running jobs is less than the max, break the loop
        if [[ $current_jobs -lt $max_jobs ]]; then
            break
        fi
        sleep 5
    done
}

# ######################################################################################################
# # firebox

# cd $post_processing_fire_outputs
# cd skirt/runs_hden_radius/firebox/z${redshift}

# # Max number of concurrent jobs
# MAX_JOBS=10

# # Loop through directories
# for i in {900..1000}; do
#     dir="gal$i/$directory_name"
#     if [[ -d $dir ]]; then
#         cd $dir
        
#         # Check for file existence with '_grid_radiationField_wavelengths' in the name
#         if ! ls *_grid_radiationField_wavelengths* 1> /dev/null 2>&1; then
#             skirt "gal$i.ski" &  # Run skirt in the background if no such file exists
#         fi
        
#         cd - > /dev/null    # Return to the previous directory without printing the path

#         # Wait if the number of running jobs reaches MAX_JOBS
#         wait_for_jobs $MAX_JOBS
#     fi
# done

# # Wait for all background jobs to finish before exiting the script
# wait


# ######################################################################################################
# zoom_in

cd $post_processing_fire_outputs
cd skirt/runs_hden_radius/zoom_in/z$redshift

# Max number of concurrent jobs
MAX_JOBS=7

# # List of galaxy names
galaxy_names=(
    "m12b_res7100_md" 
    "m12c_res7100_md"
    "m12f_res7100_md"
    "m12i_res7100_md"
    "m12m_res7100_md"
    "m12q_res7100_md"
    "m12r_res7100_md"
    "m12w_res7100_md"
    "m11d_r7100"
    "m11e_r7100"
    "m11h_r7100"
    "m11i_r7100"
    "m11q_r7100"    
)


# Loop through galaxy names
for galaxy in "${galaxy_names[@]}"; do

    dir="$galaxy/$directory_name"
    # dir="gal$i/average_sobolev_smoothingLength" 
    if [[ -d $dir ]]; then
        cd $dir
        
        # Check for file existence with '_grid_radiationField_wavelengths' in the name
        if ! ls *_grid_radiationField_wavelengths* 1> /dev/null 2>&1; then
            skirt "$galaxy.ski" &  # Run skirt in the background
        fi
        
        cd - > /dev/null    # Return to the previous directory without printing the path

        # Wait if the number of running jobs reaches MAX_JOBS
        wait_for_jobs $MAX_JOBS
    fi
done


# Wait for all background jobs to finish before exiting the script
wait
