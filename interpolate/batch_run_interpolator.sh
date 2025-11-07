#!/bin/bash
#SBATCH --account=rrg-rbond-ac
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --time=10:00:00
#SBATCH --job-name=z3_temperature_interpolate
#SBATCH --output=%x.out
#SBATCH --error=%x.err

# Go to working directory
cd /scratch/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/interpolate


## Calculation for parallel njob_parallel 
nnodes=1
ncores=192
njob_parallel=32 # 2, 4, 8, 16, 32, 64, 192 -- 64
cores_per_job=$(echo "$ncores / $njob_parallel" | bc)   # integer truncation
echo "Each job will use $cores_per_job cores."

# Set the threads for libraries 
export OMP_NUM_THREADS=$cores_per_job
export MKL_NUM_THREADS=$cores_per_job
export OPENBLAS_NUM_THREADS=$cores_per_job
export NUMEXPR_NUM_THREADS=$cores_per_job
echo "Set OMP_NUM_THREADS, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS, NUMEXPR_NUM_THREADS to $cores_per_job."


# Load the modules
module purge
source ~/.venv_all/bin/activate

# Clean up the old error messages and folders 
source clean.sh

### Set up the parameters 
redshift=3.0
target="temperature" #luminosity_from_luminosity_per_mass_by_dividing_to_4pi, abundance, temperature

outdir="timings_z$redshift/$target"
# Check if directory exists. If exists delete it and create a new one
[ -d "$outdir" ] && rm -rf "$outdir"
# Create directory
mkdir -p "$outdir"

zoom_in_galaxies=(
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


# Build commands.tsv: two columns: galaxy suite
{
  # Firebox 
  for i in {0..1000}; do
    echo "gal$i firebox $redshift $target" 
  done
  # Zoom-in
  for g in "${zoom_in_galaxies[@]}" ; do
    echo "$g zoom_in $redshift $target"
  done
  # High res zoom-in
  for g in m12i_r880_md ; do 
    echo "$g particle_split $redshift $target" 
  done    
} > $outdir/commands.tsv


parallel --jobs "$njob_parallel" --colsep ' ' --joblog "$outdir/parallel.log" \
  '/usr/bin/time -v python interpolating_for_gas_particles.py {1} {2} {3} {4} > '"$outdir"'/{1}.out 2> '"$outdir"'/{1}.time' \
  :::: "$outdir/commands.tsv"

echo "All tasks completed."