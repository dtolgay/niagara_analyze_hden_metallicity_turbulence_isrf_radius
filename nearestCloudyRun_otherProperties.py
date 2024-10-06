import sys
sys.path.append("/home/m/murray/dtolgay/scratch")

import numpy as np 
import pandas as pd 
import os
from scipy.spatial import KDTree
from time import time 

from tools import constants

# Global variables
epsilon = 1e-30


def main(galaxy_name, galaxy_type, redshift):

    # directory_name = "lichen_voronoi_1e6"
    directory_name = "voronoi_1e6"
    # directory_name = "voronoi_1e6_improved_wavelength_bin"
    # directory_name = "trial1"

    print(
        f"------------------------------------------ {galaxy_name} ------------------------------------------"
    )

    ## Check if file exits. If it exists do not continue running the code, if not run the code.
    cloudy_gas_particles_file_directory = f"/home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"
    # cloudy_gas_particles_file_directory = f"/home/m/murray/dtolgay/scratch/cloudy_runs/z_3/m12f_res7100_md_test"

    write_file_path = f"{cloudy_gas_particles_file_directory}/otherProperties_nearestCloudyRun.txt"

    print("\n")
    if os.path.isfile(write_file_path):
        print(f"{write_file_path} exits, the code is stopped.")
        return 0
    else: 
        print(f"{write_file_path} doen't exist. Continuing....")


    # Read gas particles 
    gas_particles_df, gas_column_names = read_cloudy_gas_particles(cloudy_gas_particles_file_directory)

    ################ Read training data particles 
    properties_column_names = [
        "fh2",
        "fCO",
    ]
    
    # 1st set of run
    train_data_base_file_dir_1 = "/scratch/m/murray/dtolgay/cloudy_runs/z_0"
    train_data_main_directory_1 = "cr_1_CO87_CII_H_O3/cr_1_CO87_CII_H_O3_metallicity_above_minus_2" 

    train_data_df_1, properties_column_names_with_log = read_training_data(
        base_file_dir = train_data_base_file_dir_1, 
        main_directory = train_data_main_directory_1, 
        file_name = "other_properties.csv", 
        properties_column_names = properties_column_names,
    )    

    # 2nd set of run
    train_data_base_file_dir_2 = "/scratch/m/murray/dtolgay/cloudy_runs/z_0"
    train_data_main_directory_2 = "cr_1_CO87_CII_H_O3/cr_1_CO87_CII_H_O3_metallicity_minus2_minus3point5" 

    train_data_df_2, properties_column_names_with_log = read_training_data(
        base_file_dir = train_data_base_file_dir_2, 
        main_directory = train_data_main_directory_2, 
        file_name = "other_properties.csv", 
        properties_column_names = properties_column_names,
    )    


    # Concattanete two dataframes 
    train_data_df = pd.concat([train_data_df_2, train_data_df_1])
    train_data_file_paths = [f"{train_data_base_file_dir_1}/{train_data_main_directory_1}", f"{train_data_base_file_dir_2}/{train_data_main_directory_2}"]

    # train_data_file_paths = [f"{train_data_base_file_dir_1}/{train_data_main_directory_1}"]
    # train_data_df = train_data_df_1


    ########
    # Interpolate
    gas_indices_interpolatedValues = interpolate_otherProperties(
        gas_particles_df=gas_particles_df, 
        train_data_df=train_data_df, 
        properties_column_names_with_log=properties_column_names_with_log
        )

    column_names = ['index'] + properties_column_names # Now this is not log because I took the exponential when I am interpolating 
    gas_indices_luminosities_df = pd.DataFrame(gas_indices_interpolatedValues, columns=column_names)

    ### 
    # Merge two dataframes
    if len(gas_indices_luminosities_df) == len(gas_particles_df):
        print("Lengths of luminosities and gas particles are the same. Merging can be done.")
        merged_df = gas_particles_df.merge(gas_indices_luminosities_df, how='left', on='index', validate='one_to_one') # Check if it is one to one 
    else:
        print("Lengths of luminosities and gas particles are NOT same. Exiting with code 3...")
        exit(3) 

    # Write to a file
    write_to_a_file(
        write_file_path = write_file_path, 
        train_data_file_paths = train_data_file_paths,
        gas_column_names = gas_column_names, 
        properties_column_names = properties_column_names, 
        merged_df = merged_df
        )            

    return 0


def read_cloudy_gas_particles(cloudy_gas_particles_file_directory):
    # Define the column names based on your description
    gas_column_names = [
        "x",
        "y",
        "z",
        "smoothing_length",
        "mass",
        "metallicity",
        "temperature",
        "vx",
        "vy",
        "vz",
        "hden",
        "radius",
        "sfr",
        "turbulence",
        "density",
        "mu_theoretical",
        "average_sobolev_smoothingLength",
        "index",
        "isrf",
    ]

    gas_particles_df = pd.read_csv(
        f"{cloudy_gas_particles_file_directory}/cloudy_gas_particles.txt",
        delim_whitespace=True,
        comment="#",
        names=gas_column_names,
    )

    gas_particles_df['dummy_radius'] = gas_particles_df['smoothing_length'] / 2 # TODO: Delete 

    # isrf of the gas particles can be zero, therefore set them equal to a very small number
    gas_particles_df.loc[gas_particles_df["isrf"] == 0, "isrf"] = 1e-30

    # Extend the dataframe by adding the log of the parameters
    gas_particles_df[
        [
            "log_metallicity",
            "log_density",
            "log_turbulence",
            "log_isrf",
            "log_hden",
            "log_radius",
            "log_smoothing_length",
            "log_average_sobolev_smoothingLength",
            "log_dummy_radius", # TODO: Delete
        ]
    ] = np.log10(
        gas_particles_df[
            [
                "metallicity",
                "density",
                "turbulence",
                "isrf",
                "hden",
                "radius",
                "smoothing_length",
                "average_sobolev_smoothingLength",
                "dummy_radius", # TODO: Delete
            ]
        ]
    )  # Take the log of the gas properties and interpolate using these logarithmic values.  

    print(f"{cloudy_gas_particles_file_directory}/cloudy_gas_particles.txt read and dataframe is created!")      
    
    return gas_particles_df, gas_column_names 

def read_training_data(base_file_dir, main_directory, file_name, properties_column_names):

    #################################################
    # Get the trained data
    print("Training data is started to be read.")

    # Read file
    path2TrainingData = f"{base_file_dir}/{main_directory}/{file_name}"
    unprocessed_train_data = pd.read_csv(path2TrainingData) 

    ############## Process the cloudy data 
    # Take the log of the properties 
    properties_column_names_with_log = []
    for property in properties_column_names:
        unprocessed_train_data[property] += epsilon # Add a very small number 
        unprocessed_train_data[f"log_{property}"] = np.log10(unprocessed_train_data[property])
        properties_column_names_with_log.append(f"log_{property}")

    # Discard all nan values 
    print("Dropping NaN containing lines")
    processed_train_data = unprocessed_train_data.dropna()        
    train_data_df = processed_train_data.drop(properties_column_names, axis=1) # Drop the columns which log is not taken.

    # # Double check if there is any NaN
    # if (np.isnan(train_data_df.values).any()):
    #     print("Still there are NaN values. Exiting with code 1...")
    #     exit(1)
    # elif (np.isinf(train_data_df.values).any()):
    #     print("Still there are inf values. Exiting with code 2...")
    #     exit(2)

    ######
    # Add the column density data to interpolate that too 
    train_data_df['log_column_density'] = np.log10(
        (10**train_data_df['log_hden'] / constants.cm2pc**3) * (10**train_data_df['log_radius']) * (constants.mu_h * constants.proton_mass * constants.kg2Msolar)
    ) # Msolar / pc^2

    print(f"{path2TrainingData} is read.")


    return train_data_df, properties_column_names_with_log

def interpolate_otherProperties(gas_particles_df, train_data_df, properties_column_names_with_log):

    print("I am in the interpolate_otherProperties")

    train_data_column_names = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",    
    ]    

    tree = KDTree(
        train_data_df[train_data_column_names].to_numpy(),
    ) # Create a tree for the training data

    scale_length = [
        "log_average_sobolev_smoothingLength"
    ]

    gas_data_column_names = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",      
    ] + scale_length

    gas_indices_luminosities = []
    
    intial_index = gas_particles_df.iloc[0]['index']
    for index, gas in gas_particles_df.iterrows():
        if intial_index == 0:
            if (gas['index'] % int(1e5) == 1):
                print(f"{gas['index']} finished. Left {len(gas_particles_df) - gas['index']}")

        # Find the closest train_data 
        distances, cloesest_run_indices = tree.query(gas[gas_data_column_names].to_numpy(), k=1) 

        # Interpolated value -> Actually this is not an interpolation but this is nearestND type interpolation
        log_interpolated_value = train_data_df.iloc[cloesest_run_indices][properties_column_names_with_log].to_numpy()
        interpolated_values = 10**log_interpolated_value

        gas_indices_luminosities.append(
            np.concatenate(([gas['index']], interpolated_values))
        )

    return gas_indices_luminosities


def write_to_a_file(write_file_path, train_data_file_paths, gas_column_names, properties_column_names, merged_df):

    train_data_file_paths_str = "\n".join(train_data_file_paths)

    header = f"""
    Estimated according to:
    ---------------------
    log_metallicity
    log_hden
    log_turbulence
    log_isrf
    log_average_sobolev_smoothingLength
    ---------------------

    Used training centers:
    ---------------------
    {train_data_file_paths_str}
    ---------------------

    Column 0: x-coordinate (pc)
    Column 1: y-coordinate (pc)
    Column 2: z-coordinate (pc)
    Column 3: smoothing length (pc)
    Column 4: mass (Msolar)
    Column 5: metallicity (Zsolar)
    Column 6: temperature (K)
    Column 7: vx (km/s)
    Column 8: vy (km/s)
    Column 9: vz (km/s)
    Column 10: hydrogen density (cm^-3)
    Column 11: radius (pc)
    Column 12: sfr (Msolar/yr)
    Column 13: turbulence (km/s)
    Column 14: density (gr/cm^-3)
    Column 15: mu_theoretical (1)
    Column 16: average_sobolev_smoothingLength (pc)
    Column 17: index [1]
    Column 18: isrf [G0]
    Column 19: fh2 [1]
    Column 20: fCO [1] Σco / ΣH2
    """

    write_df = merged_df[gas_column_names + properties_column_names]

    np.savetxt(fname=write_file_path, X=write_df, fmt="%.5e", header=header)

    print(f"File saved to: {write_file_path}")

    return 0

if __name__ == "__main__":
    galaxy_name = sys.argv[1]
    galaxy_type = sys.argv[2]
    redshift = sys.argv[3]

    main(galaxy_name, galaxy_type, redshift)