#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=10:00:00
#SBATCH --job-name=run_skirt
#SBATCH --output=run_skirt.out
#SBATCH --error=run_skirt.err

# modules 
module purge



redshift=2.0
directory_name=voronoi_1e5
MAX_JOBS=5


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

######################################################################################################
# firebox

cd $dongwoo_directory
cd runs_hden_radius/firebox/z${redshift}

# Max number of concurrent jobs

# Loop through directories
for i in {900..1000}; do
    dir="gal$i/$directory_name"
    if [[ -d $dir ]]; then
        cd $dir
        
        # Check for file existence with '_grid_radiationField_wavelengths' in the name
        if ! ls *_grid_radiationField_wavelengths* 1> /dev/null 2>&1; then
            skirt "gal$i.ski" &  # Run skirt in the background if no such file exists
        fi
        
        cd - > /dev/null    # Return to the previous directory without printing the path

        # Wait if the number of running jobs reaches MAX_JOBS
        wait_for_jobs $MAX_JOBS
    fi
done

# Wait for all background jobs to finish before exiting the script
wait