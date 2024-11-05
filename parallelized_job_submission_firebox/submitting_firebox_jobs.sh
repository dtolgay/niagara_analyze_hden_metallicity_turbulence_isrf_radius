#!/bin/bash

# # Loop to create and submit 50 jobs, each processing 20 galaxies
# for i in {0..49}; do
#     start=$((i * 20))
#     end=$((start + 19))
#     echo "Submitting job for galaxies $start to $end"
#     sbatch batch_linearNDInterpolator_argumentNeeded.sh $start $end
# done


# Loop to create and submit 50 jobs, each processing 20 galaxies
for i in {0..49}; do
    start=$((i * 2))
    end=$((start + 1))
    echo "Submitting job for galaxies $start to $end"
    sbatch batch_linearNDInterpolator_argumentNeeded.sh $start $end
done


# # Loop to create and submit 50 jobs, each processing 2 galaxies
# for i in {0..49}; do
#     start=$(( i * 2 + 1 ))
#     end=$(( start + 1 ))
#     echo "Submitting job for galaxies $start to $end"
#     sbatch batch_linearNDInterpolator_argumentNeeded.sh $start $end
# done