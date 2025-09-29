# Go to working directory
cd /scratch/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/interpolate


## Calculation for parallel njob_parallel 
nnodes=1
ncores=192
njob_parallel=64 # 2, 4, 8, 16, 32
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
redshift=0.0
target="luminosity_from_luminosity_per_mass_by_dividing_to_4pi"

outdir="timings_z$redshift"
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
  for i in {0..61}; do
    echo "gal$i firebox $redshift $target" 
  done
  # # Zoom-in
  # for g in "${zoom_in_galaxies[@]}" ; do
  #   echo "$g zoom_in $redshift $target"
  # done
  # # High res zoom-in
  # for g in m12i_r880_md ; do 
  #   echo "$g particle_split $redshift $target" 
  # done    
} > $outdir/commands.tsv

# Use xargs to run up to $njob_parallel concurrent srun tasks
# Each line has two fields: {galaxy} {suite}
xargs -a "$outdir/commands.tsv" -n4 -P "$njob_parallel" -I{} bash -lc '
  set -- {};
  galaxy="$1"; suite="$2"; redshift_in_file="$3"; target_="$4";
  srun --exclusive -N1 -n1 -c '"$cores_per_job"' /usr/bin/time -v \
    python interpolating_for_gas_particles.py "$galaxy" "$suite" "$redshift_in_file" "$target_" \
    > '"$outdir"'/"$galaxy".out 2> '"$outdir"'/"$galaxy".time
'
