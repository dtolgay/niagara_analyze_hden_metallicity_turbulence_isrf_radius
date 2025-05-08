#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=01:00:00
#SBATCH --job-name=run_skirt
#SBATCH --output=run_skirt.out
#SBATCH --error=run_skirt.err

# Load modules 
module purge
ml gcc/8.3.0 openmpi/4.1.4

cd /gpfs/fs0/scratch/r/rbond/dongwooc/scratch_rwa/doga/runs_hden_radius/firebox/z2.0/gal0/voronoi_1e5

mpirun -np 40 skirt gal0.ski