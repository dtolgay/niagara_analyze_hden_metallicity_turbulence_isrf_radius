#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=hybridInterpolator_otherProperties
#SBATCH --output=hybridInterpolator_otherProperties.out
#SBATCH --error=hybridInterpolator_otherProperties.err


module purge 
module load python/3.11.5 

cd ~
cd /home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius

# Read arguments (e.g., start and end indices of galaxies)
start_idx=$1
end_idx=$2

number_of_background_galaxies=1
redshift=0.0
number_of_processors_per_galaxy=40


wait_for_jobs() {
    for job in $(jobs -p)
    do
        wait $job
    done
}


####### firebox
# Counter for every 10 galaxies
counter=0

# Run the loop for the specified range of galaxy indices
for i in $(seq $start_idx $end_idx); do
    python hybridInterpolator_otherProperties.py gal$i firebox $redshift $number_of_processors_per_galaxy &
    # Wait for jobs every 10 galaxies
    ((counter++))
    if [ $counter -ge $number_of_background_galaxies ]; then
        wait_for_jobs
        counter=0
    fi
done

# Wait for the last set of jobs to finish
wait_for_jobs
