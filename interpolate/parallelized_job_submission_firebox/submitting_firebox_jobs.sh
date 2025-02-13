#!/bin/bash

# Number of galaxies processed per job
number_of_galaxies_to_process=50
number_of_submitted_jobs=20 # 1000/50


# Loop to create and submit jobs
for ((i=0; i<number_of_submitted_jobs; i++)); do
    start=$((i * number_of_galaxies_to_process))
    end=$((start + number_of_galaxies_to_process))
    echo "Submitting job for galaxies $start to $end"
    
    # Submit the job
    sbatch batch_linearNDInterpolator_argumentNeeded.sh $start $end
done
