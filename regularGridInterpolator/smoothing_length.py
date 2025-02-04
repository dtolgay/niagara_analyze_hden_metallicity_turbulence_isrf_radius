import os
import time
import pandas as pd # type: ignore
import numpy as np
import sys 
import joblib # type: ignore

def main(files_info, interpolators, interpolator_type):
    # Check if file exits. If exists, stop running with a message indicating that the file already exists.
    if os.path.exists(files_info["write_file_name"]):
        print(f"{files_info['write_file_name']} already exists. Stopping the run.")
        return None
    
    # Read the centers
    gas_particles = read_cloudy_gas_particles(files_info)

    # Read the interpolator
    interpolator = joblib.load(interpolators[interpolator_type]["interpolator_path"])

    # Prepare the gas particles for interpolation
    interpolation_columns_for_gas_particles = [
        "metallicity",
        "hden",
        "turbulence",
        "isrf",
        "smoothing_length",
    ]

    gas_particles = take_log_the_interpolation_centers_for_gas_particles(gas_particles, interpolation_columns_for_gas_particles)

    # measure the time it takes to interpolate in minutes
    time_start = time.time()
    # Interpolate the gas particles
    gas_particles = interpolate(
        interpolator = interpolator, 
        gas_particles = gas_particles, 
        interpolation_columns_for_gas_particles = interpolation_columns_for_gas_particles, 
        target_columns = interpolators[interpolator_type]["target_columns"]
    )
    time_end = time.time()
    print(f"Interpolation took {(time_end-time_start)/60:.2f} minutes.")

    
    print(gas_particles.head())

    return None 


def interpolate(interpolator, gas_particles, interpolation_columns_for_gas_particles, target_columns):
    
    print(gas_particles[interpolation_columns_for_gas_particles].values)

    # Interpolate the gas particles
    interpolated_values = interpolator(
        gas_particles[interpolation_columns_for_gas_particles].values
    )

    # Put the results into dataframe 
    gas_particles[target_columns] = interpolated_values

    return gas_particles


def read_cloudy_gas_particles(files_info):
    """
    Read the cloudy gas particles
    """

    # Read the cloudy gas particles
    gas_column_names = [
        "x",                                    # pc 
        "y",                                    # pc 
        "z",                                    # pc 
        "smoothing_length",                     # pc 
        "mass",                                 # Msolar
        "metallicity",                          # Zsolar
        "temperature",                          # Kelvin
        "vx",                                   # km/s
        "vy",                                   # km/s
        "vz",                                   # km/s
        "hden",                                 # 1/cm3
        "radius",                               # pc 
        "sfr",                                  # Msolar / year
        "turbulence",                           # km/s
        "density",                              # gr/cm3
        "mu_theoretical",                       # 1
        "average_sobolev_smoothingLength",      # pc 
        "index",                                # 1
        "isrf",                                 # G0
    ]

    gas_particles = pd.DataFrame(
        np.loadtxt(f'{files_info["gas_particles_path"]}/cloudy_gas_particles.txt'),
        columns=gas_column_names,
    )    

    return gas_particles


def take_log_the_interpolation_centers_for_gas_particles(gas_particles, interpolation_columns_for_gas_particles):
    
    for column in interpolation_columns_for_gas_particles:
        gas_particles[f"log_{column}"] = np.log10(gas_particles[column])

    return gas_particles

if __name__ == "__main__":

    interpolators_base_fdir = "/scratch/m/murray/dtolgay/cloudy_runs/interpolators"

    interpolators = {
        "temperature": {
            "file_name": "temperatures.txt",
            "target_columns": ["Th2", "Tco", "T"],
            "interpolator_path": f"{interpolators_base_fdir}/RegularGridInterpolator_linear_temperature.joblib",
            "write_file_name": "temperature_smoothingLength_regularGridInterpolator_linear",
        },
        "line_emissions": {
            "file_name": "I_line_values_without_reversing.txt",
            "target_columns": [
                "I_ly_alpha",
                "I_h_alpha",
                "I_h_beta",
                "I_co_10",
                "I_co_21",
                "I_co_32",
                "I_co_43",
                "I_co_54",
                "I_co_65",
                "I_co_76",
                "I_co_87",
                "I_13co",
                "I_c2",
                "I_o3_88",
                "I_o3_5006",
                "I_o3_4958",
            ],
            "interpolator_path": f"{interpolators_base_fdir}/RegularGridInterpolator_linear_flux.joblib",
            "write_file_name": "line_emissions_smoothingLength_regularGridInterpolator_linear",
        },
        "abundance": {
            "file_name": "other_properties.csv",
            "target_columns": [
                "fh2",
                "fCO",
            ],
            "interpolator_path": f"{interpolators_base_fdir}/RegularGridInterpolator_linear_abundance.joblib",
            "write_file_name": "abundance_smoothingLength_regularGridInterpolator_linear",
        }
    }

    base_fdir =  "/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/runs_hden_radius",
    galaxy_name = sys.argv[1]
    galaxy_type = sys.argv[2]
    redshift = sys.argv[3]
    interpolator_type = sys.argv[4]

    galaxy_info = {
        "base_fdir": base_fdir,
        "galaxy_name": galaxy_name,
        "galaxy_type": galaxy_type,
        "redshift": redshift,
        "direcotory": "voronoi_1e6",
    }

    files_info = {
        "write_file_name": f"{base_fdir}/{interpolator_type}_regularGridInterpolator_linear_smoothingLength.txt",
        # "gas_particles_path": f"{base_fdir}/{galaxy_info['galaxy_type']}/z{galaxy_info['redshift']}/{galaxy_info['galaxy_name']}/{galaxy_info['directory']}",
        "gas_particles_path": f"/scratch/m/murray/dtolgay/cloudy_runs/z_0/m12i_res7100_md_test",
    
    }


    main(files_info, interpolators, interpolator_type)