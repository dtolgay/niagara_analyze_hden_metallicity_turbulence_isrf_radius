from itertools import chain
import sys
sys.path.append("/home/m/murray/dtolgay/scratch")
from tools import constants

import numpy as np 
import pandas as pd 
import os

from scipy.spatial import KDTree
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

from time import time 
from concurrent.futures import ProcessPoolExecutor


# Global variables
epsilon = 1e-30
# mu = 1.2
mu = 1.38 # Krumholz and Gnedin 



def main(galaxy_name, galaxy_type, redshift, max_workers, write_interpolator_info=False):

    start = time()

    directory_name = "voronoi_1e6"
    # directory_name = "voronoi_1e5"

    print(
        f"------------------------------------------ {galaxy_name} ------------------------------------------"
    )

    ## Check if file exits. If it exists do not continue running the code, if not run the code.
    cloudy_gas_particles_file_directory = f"/home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"
    # cloudy_gas_particles_file_directory = f"/home/m/murray/dtolgay/scratch/cloudy_runs/z_3/m12f_res7100_md_test"

    write_file_path = f"{cloudy_gas_particles_file_directory}/L_line_averageSobolevH_hybridInterpolator_flux2Luminosity.txt"

    print("\n")
    if os.path.isfile(write_file_path):
        print(f"{write_file_path} exits, the code is stopped.")
        return 0
    else: 
        print(f"{write_file_path} doen't exist. Continuing....")


    # Read gas particles 
    gas_particles_df, gas_column_names = read_cloudy_gas_particles(cloudy_gas_particles_file_directory)

    base_line_names = [
        "ly_alpha",
        "h_alpha",
        "h_beta",
        "co_10",
        "co_21",
        "co_32",
        "co_43",
        "co_54",
        "co_65",
        "co_76",
        "co_87",
        "13co",
        "c2",
        "o3_88",
        "o3_5006",
        "o3_4958",        
    ]

    ################ Read training data particles 

    # 1st set of run
    train_data_base_file_dir_1 = "/scratch/m/murray/dtolgay/cloudy_runs/z_0"
    train_data_main_directory_1 = "cr_1_CO87_CII_H_O3/cr_1_CO87_CII_H_O3_metallicity_above_minus_2" 

    train_data_df_1, line_names_with_log = read_training_data(
        base_file_dir = train_data_base_file_dir_1, 
        main_directory = train_data_main_directory_1, 
        file_name = "I_line_values_without_reversing.txt", 
        base_line_names = base_line_names
    )    

    # 2nd set of run
    train_data_base_file_dir_2 = "/scratch/m/murray/dtolgay/cloudy_runs/z_0"
    train_data_main_directory_2 = "cr_1_CO87_CII_H_O3/cr_1_CO87_CII_H_O3_metallicity_minus2_minus3point5" 

    train_data_df_2, line_names_with_log = read_training_data(
        base_file_dir = train_data_base_file_dir_2, 
        main_directory = train_data_main_directory_2, 
        file_name = "I_line_values_without_reversing.txt", 
        base_line_names = base_line_names
    )    

    train_data_file_paths = [f"{train_data_base_file_dir_1}/{train_data_main_directory_1}", f"{train_data_base_file_dir_2}/{train_data_main_directory_2}"]

    # Concattanete two dataframes 
    train_data_df = pd.concat([train_data_df_2, train_data_df_1])

    ########

    # train_data_base_file_dir_test = "/home/m/murray/dtolgay/scratch/cloudy_runs/z_3/"
    # train_data_main_director_test = "m12f_res7100_md_test"
    # train_data_df, line_names_with_log = read_training_data(
    #     base_file_dir = train_data_base_file_dir_test, 
    #     main_directory = train_data_main_director_test, 
    #     file_name = "I_line_values_without_reversing.txt", 
    #     base_line_names = base_line_names
    # )   

    # train_data_file_paths = [f"{train_data_base_file_dir_test}/{train_data_main_director_test}"]

    # Split dataframe into several dataframes to run the parallely. 
    gas_particles_df_chunks = split_dataframe(
            df=gas_particles_df,
            max_workers=max_workers, 
        )
    
    ################

    # # Calculate the line luminosities from intensity data 
    # gas_indices_luminosities = []
    # for gas_particles_df_chunk in gas_particles_df_chunks:
    #     gas_indices_luminosities.append(
    #         calculate_Lline(
    #             gas_particles_df_chunk, 
    #             train_data_df, 
    #             line_names_with_log
    #         )
    #     )

    # Calculate the line luminosities from intensity data in parallel
    # used_interpolator_info_chunks = []
    # gas_indices_luminosities_chunks = []
    gas_indices_luminosities = []
    used_interpolator_info = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(calculate_Lline, gas_particles_df_chunk, train_data_df, line_names_with_log)
            for gas_particles_df_chunk in gas_particles_df_chunks
        ]
        for future in futures:
            result = future.result()
            gas_indices_luminosities.extend(result[0])
            used_interpolator_info.extend(result[1])    
        
            
    # Create df of the retuned luminosities 
    log_line_names = []
    for line_name in base_line_names:
        log_line_names.append(f"L_{line_name}")

    column_names = ['index'] + log_line_names

    gas_indices_luminosities_df = pd.DataFrame(gas_indices_luminosities, columns=column_names).sort_values(by="index", ascending=True)
    used_interpolator_info_df = pd.DataFrame(used_interpolator_info, columns=['index', 'used_interpolator'])


    # Change the unit of the calculated CO luminosities 
    CO_lines_and_wavelengths = {
        "L_co_10": 2600.05e-6,  # meter
        "L_co_21": 1300.05e-6,
        "L_co_32": 866.727e-6,
        "L_co_43": 650.074e-6,
        "L_co_54": 520.08e-6,
        "L_co_65": 433.438e-6,
        "L_co_76": 371.549e-6,
        "L_co_87": 325.137e-6,
        "L_13co": 2719.67e-6,
    }
    
    gas_indices_luminosities_df = change_unit_of_CO_emission(
        gas_indices_luminosities_df=gas_indices_luminosities_df,
        lines_and_wavelengths=CO_lines_and_wavelengths,
    )

    # Merge two dataframes
    if len(gas_indices_luminosities_df) == len(gas_particles_df):
        print("Lengths of luminosities and gas particles are the same. Merging can be done.")
        merged_df = gas_particles_df.merge(gas_indices_luminosities_df, how='inner', on='index', validate='one_to_one') # Check if it is one to one 
    else:
        print("Lengths of luminosities and gas particles are NOT same. Exiting with code 3...")
        exit(3)

    # Write to a file
    write_to_a_file(
        write_file_path = write_file_path, 
        train_data_file_paths = train_data_file_paths,
        gas_column_names = gas_column_names, 
        base_line_names = base_line_names, 
        merged_df = merged_df
        )
    

    # Write used interpolator to another file
    if write_interpolator_info:
        fpath_write_interpolator_info = f"{cloudy_gas_particles_file_directory}/hybridInterpolator_info.csv"
        used_interpolator_info_df.to_csv(
            fpath_write_interpolator_info
        )

    end = time()

    print(f"Code took {np.round((end - start) / 60, 3)} minutes")

    return 0

######## Functions

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

def read_training_data(base_file_dir, main_directory, file_name, base_line_names):

    #################################################
    # Get the trained data

    print("Training data is started to be read.")

    line_names = []
    for line_name in base_line_names:
        line_names.append(f"I_{line_name}")

    column_names = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",
    ]  + line_names

    # Read file
    path2TrainingData = f"{base_file_dir}/{main_directory}/{file_name}"
    unprocessed_train_data = pd.DataFrame(
        np.loadtxt(fname=path2TrainingData),
        columns=column_names,
    )

    ############## Process the cloudy data 
    # Discard all nan values 
    print("Dropping NaN containing lines")
    unprocessed_train_data = unprocessed_train_data.dropna()

    # Check if all intensities are positive and set 0 values to epsilon
    print(f"Check if all intensities are positive. Then set 0 values to {epsilon}")
    all_positive_columns = (unprocessed_train_data[line_names] >= 0).all().all()
    if all_positive_columns:
        print(f"All of the intensity values are non-negative. Continuing...")
    else:
        # Set values smaller or equal to zero to epsilon in specified columns
        for col in line_names:
            unprocessed_train_data[col] = unprocessed_train_data[col].map(lambda x: epsilon if x <= 0 else x)
        print(f"Not all intensities are are non-negative. Setting them to epsilon")


    line_names_with_log = []
    for column in line_names:
        unprocessed_train_data[f"log_{column}"] = np.log10(unprocessed_train_data[column])
        line_names_with_log.append(f"log_{column}") # Store the new line names


    train_data_df = unprocessed_train_data[[
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",
        ] + line_names_with_log]  # Only use the log of the line luminosities    

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
        (10**train_data_df['log_hden'] / constants.cm2pc**3) * (10**train_data_df['log_radius']) * (mu * constants.proton_mass * constants.kg2Msolar)
    ) # Msolar / pc^2

    print(f"{path2TrainingData} is read.")


    return train_data_df, line_names_with_log

def split_dataframe(df, max_workers):
    # Create different chunks of of dataframe to run them parallely
    n = len(df)
    chunk_size = -(
        -n // max_workers
    )  # Ceiling division to ensure all rows are included

    # Split the dataframe into chunks and store in an array
    return [df[i : i + chunk_size] for i in range(0, n, chunk_size)]

def prepare_interpolator(k, gas, gas_data_column_names, tree, train_data_df, train_data_column_names, line_names_with_log, interpolator="LinearNDInterpolator"):
    # Query the tree for neighbors
    distances, indices = tree.query(gas[gas_data_column_names].to_numpy(), k=k)
    
    # Set up linearNDInterpolator
    if interpolator == "LinearNDInterpolator":  
        interpolator = LinearNDInterpolator(
            points=train_data_df.iloc[indices][train_data_column_names].to_numpy(),
            values=train_data_df.iloc[indices][line_names_with_log].to_numpy()
        )
    elif interpolator == "NearestNDInterpolator":
        interpolator = NearestNDInterpolator(
            train_data_df.iloc[indices][train_data_column_names].to_numpy(),
            train_data_df.iloc[indices][line_names_with_log].to_numpy()
        )
    else:
        return None
    
    return interpolator

def calculate_Lline(gas_particles_df, train_data_df, line_names_with_log):

    print("I am in the calculate_Lline")
    
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
    used_interpolator_info = []

    intial_index = gas_particles_df.iloc[0]['index']
    for index, gas in gas_particles_df.iterrows():
        if intial_index == 0:
            if (gas['index'] % int(1e5) == 1):
                print(f"{gas['index']} finished. Left {len(gas_particles_df) - gas['index']}")

        # List of k values to try in order
        k_values = [50, 100, 500, 1000, 2000, 3000, 5000, int(1e4)]

        for k in k_values:
            try:
                # Try to query and interpolate with the current k value
                interpolator = prepare_interpolator(
                        k, gas, gas_data_column_names, tree, train_data_df, train_data_column_names, line_names_with_log, interpolator="LinearNDInterpolator"
                    )
                
                # Interpolate to check if there are NaN values 
                interpolated_intensities = 10**interpolator(gas[gas_data_column_names])[0] # It returns an array of arrays. That's why [0] is done.

                # If there exist any NaN change iterate to the next k value:
                if np.isnan(interpolated_intensities).any(): 
                    if k < 300:
                        continue
                    else:
                        interpolator = prepare_interpolator(
                                k, 
                                gas, 
                                gas_data_column_names, 
                                tree, 
                                train_data_df, 
                                train_data_column_names, 
                                line_names_with_log,
                                interpolator="NearestNDInterpolator"
                            )
                        interpolated_intensities = 10**interpolator(gas[gas_data_column_names])[0] # It returns an array of arrays. That's why [0] is done.
                        used_interpolator = "NearestND"
                        break
                else: 
                    used_interpolator = "LinearND"
                    break  # Break out of the loop if and there exist no NaN values 
            except Exception as e:
                # If it fails with the current k, continue to the next one
                continue
                    
        if interpolator == None:
            print(f"Error: interpolator is None for index: {gas['index']}")
            exit(99)

        # Calculate area 
        area = gas['mass'] * constants.M_sun2gr / (gas['density'] * (10**gas[scale_length[0]] * constants.pc2cm)) 

        # Calculate luminosity
        luminosities = interpolated_intensities * area

        gas_indices_luminosities.append(
            np.concatenate(([gas['index']], luminosities))
        )
        
        used_interpolator_info.append([gas['index'], used_interpolator])

    return gas_indices_luminosities, used_interpolator_info 

def meters_to_Ghz_calculator(wavelength_in_meters):
    c = 299792458  # m/s
    frequency_in_Ghz = c / wavelength_in_meters * 1e-9
    return frequency_in_Ghz

def return_ergs_per_second2radio_units(rest_frequency):
    ergs_per_second2solar_luminosity = constants.ergs2Lsolar # TODO: Everyone is citing different values. Ask the value of conversion from erg/s -> Lsolar
    solar_luminosity2radio_units = (3e-11 * (rest_frequency**3)) ** (-1)  # Rest frequency should be in Ghz
    ergs_per_second2radio_units = (ergs_per_second2solar_luminosity * solar_luminosity2radio_units)

    return ergs_per_second2radio_units

def change_unit_of_CO_emission(gas_indices_luminosities_df, lines_and_wavelengths):

    for line in lines_and_wavelengths: 
        conversion_factor = return_ergs_per_second2radio_units(
            rest_frequency=meters_to_Ghz_calculator(lines_and_wavelengths[line])
        )
        gas_indices_luminosities_df[line] *= conversion_factor

    return gas_indices_luminosities_df

def write_to_a_file(write_file_path, train_data_file_paths, gas_column_names, base_line_names, merged_df):

    train_data_file_paths_str = "\n".join(train_data_file_paths)

    header = f"""
    Gas particles for {galaxy_name} galaxy

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
    Column 19: L_ly_alpha [erg s^-1]
    Column 20: L_h_alpha [erg s^-1]
    Column 21: L_h_beta [erg s^-1]
    Column 22: L_co_10 [K km s^-1 pc^2]
    Column 23: L_co_21 [K km s^-1 pc^2]
    Column 24: L_co_32 [K km s^-1 pc^2]
    Column 25: L_co_43 [K km s^-1 pc^2]
    Column 26: L_co_54 [K km s^-1 pc^2]
    Column 27: L_co_65 [K km s^-1 pc^2]
    Column 28: L_co_76 [K km s^-1 pc^2]
    Column 29: L_co_87 [K km s^-1 pc^2]
    Column 30: L_13co [K km s^-1 pc^2]
    Column 31: L_c2 [erg s^-1]
    Column 32: L_o3_88 [erg s^-1]
    Column 33: L_o3_5006 [erg s^-1]
    Column 34: L_o3_4958 [erg s^-1]
    """

    line_names = []
    for line_name in base_line_names:
        line_names.append(f"L_{line_name}")

    write_df = merged_df[gas_column_names + line_names]

    np.savetxt(fname=write_file_path, X=write_df, fmt="%.8e", header=header)

    print(f"File saved to: {write_file_path}")

    return 0

if __name__ == "__main__":
    galaxy_name = sys.argv[1]
    galaxy_type = sys.argv[2]
    redshift = sys.argv[3]
    max_workers = int(sys.argv[4]) 

    main(galaxy_name, galaxy_type, redshift, max_workers, write_interpolator_info=True)