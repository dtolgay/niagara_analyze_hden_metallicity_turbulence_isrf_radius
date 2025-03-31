#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=05:00:00
#SBATCH --job-name=seperated_firebox_galaxies
#SBATCH --output=seperated_firebox_galaxies.out
#SBATCH --error=seperated_firebox_galaxies.err

# cd to working directory 
cd $SLURM_SUBMIT_DIR

module purge 
ml python/3.11.5

# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 line_emissions 2 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 temperature 2 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 abundance 2 &

# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 line_emissions 5 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 temperature 5 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 abundance 5 &

# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 line_emissions 10 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 temperature 10 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 abundance 10 &

# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 line_emissions 0.5 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 temperature 0.5 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 abundance 0.5 &

# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 line_emissions 0.1 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 temperature 0.1 &
# python -u debug_interpolating_for_gas_particles.py m12i_res7100_md zoom_in 0.0 abundance 0.1 &


python -u interpolating_for_gas_particles.py gal0 firebox 0.0 line_emissions &
python -u interpolating_for_gas_particles.py gal0 firebox 0.0 temperature &
python -u interpolating_for_gas_particles.py gal0 firebox 0.0 abundance &

wait