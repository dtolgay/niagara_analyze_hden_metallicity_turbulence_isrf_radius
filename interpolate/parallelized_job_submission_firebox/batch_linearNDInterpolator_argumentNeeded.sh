#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --time=23:00:00
#SBATCH --job-name=z2_temperature
#SBATCH --output=z2_temperature.out
#SBATCH --error=z2_temperature.err



module purge 
module load python/3.11.5 

cd /scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/interpolate

# Read arguments (e.g., start and end indices of galaxies)
start_idx=$1
end_idx=$2


number_of_background_galaxies=20
redshift=2.0
target="temperature" # temperature, line_emissions, abundance

wait_for_jobs() {
    for job in $(jobs -p)
    do
        wait $job
    done
}


####### firebox
# Counter for every 10 galaxies
counter=0

for i in $(seq $start_idx $end_idx); do
    python interpolating_for_gas_particles.py gal$i firebox $redshift $target &

    # Increment counter
    ((counter++))

    if [ $counter -ge $number_of_background_galaxies ]; then
        wait_for_jobs
        counter=0
    fi
done

# Wait for the last set of jobs to finish
wait_for_jobs
