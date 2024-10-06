# Imported modules
import sys
sys.path.append("/home/m/murray/dtolgay/scratch")

import numpy as np
import pandas as pd
from time import time
from concurrent.futures import ProcessPoolExecutor
import os # Check if file exits. 

from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import KDTree

import sys  # To run the batch script

from tools import constants


###########################################################################################################################################################################################################################

# Global variables
epsilon = 1e-30
# mu = 1.2
mu = 1.38 # Krumholz and Gnedin 

###########################################################################################################################################################################################################################


def main(galaxy_name, galaxy_type, redshift, max_workers):
    print(
        f"------------------------------------------ {galaxy_name} ------------------------------------------"
    )

    ## Check if file exits. If it exists do not continue running the code, if not run the code.
    cloudy_gas_particles_file_directory = f"/home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{redshift}/{galaxy_name}/voronoi_1e6"

    write_file_name = f"{cloudy_gas_particles_file_directory}/L_line_averageSobolevSmoothing_cloudyRun2_usingIline.txt"
    # write_file_name = f"{cloudy_gas_particles_file_directory}/L_line_dummyRadius_cloudyRun2_usingIline.txt"

    print("\n")
    if os.path.isfile(write_file_name):
        print(f"{write_file_name} exits, the code is stopped.")
        return 0
    else: 
        print(f"{write_file_name} doen't exist. Continuing....")


    start_code = time()

    #################################################
    # Get the gas particles

    print("Gas particles is started to read")

    # Define the column names based on your description
    gas_column_names = [
        "x-coordinate",
        "y-coordinate",
        "z-coordinate",
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

    # TODO: Delete 
    gas_particles_df['dummy_radius'] = gas_particles_df['smoothing_length'] / 2

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

    print(
        f"{cloudy_gas_particles_file_directory}/cloudy_gas_particles.txt read and dataframe is created!"
    )

    #################################################
    # Get the trained data

    print("Training data is started to be read.")

    ######

    # cr1_train_path = "/fs/lustre/scratch/dtolgay/cloudy_runs/z_0/cr_1_train/L_line_values.txt"
    # cr1_train_low_isrf_path = "/fs/lustre/scratch/dtolgay/cloudy_runs/z_0/cr_1_train_low_isrf/L_line_values.txt"

    # train_data_normal_isrf_df = read_training_data(file_path = cr1_train_path)
    # train_data_low_isrf_df = read_training_data(file_path = cr1_train_low_isrf_path)

    # # Append two dataframes together
    # train_data_df = pd.concat([train_data_normal_isrf_df, train_data_low_isrf_df], ignore_index = True)

    ######

    base_file_dir = "/home/m/murray/dtolgay/scratch/cloudy_runs/z_0"
    main_direcory = "cr_1_CO87_CII_H_O3_metallicity_above_minus_2"

    line_names = [
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
    ]

    column_names = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",
    ]  + line_names

    # Read file
    unprocessed_train_data = pd.DataFrame(
        np.loadtxt(fname=f"{base_file_dir}/{main_direcory}/I_line_values_without_reversing.txt", skiprows=1),
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
        unprocessed_train_data[line_names] = unprocessed_train_data[line_names].applymap(lambda x: epsilon if x <= 0 else x)
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

    # Double check if there is any NaN
    if (np.isnan(train_data_df.values).any()):
        print("Still there are NaN values. Exiting with code 1...")
        exit(1)
    elif (np.isinf(train_data_df.values).any()):
        print("Still there are inf values. Exiting with code 2...")
        exit(2)

    ######
    # Add the column density data to interpolate that too 
    train_data_df['log_column_density'] = np.log10(
        (10**train_data_df['log_hden'] / constants.cm2pc**3) * (10**train_data_df['log_radius']) * (mu * constants.proton_mass * constants.kg2Msolar)
    ) # Msolar / pc^2


    print(f"Training data read and data frame is created.")


    #################################################

    # Check the maxima values for the parameters. If the minima is higher than the gas particle value, cast the value of gas particle to that minima value.
    # If the maxima is lower than the gas particle value cast gas particle value to maxima. Don't extrapolate
    parameter_names = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",
    ]

    for parameter in parameter_names:
        min_train = min(train_data_df[parameter])
        max_train = max(train_data_df[parameter])

        if min(gas_particles_df[parameter]) < min_train:
            print(f"\nmin boundary problem in {parameter}")
            print(
                "Setting the gas particles having values smaller than the min_train equal to the min_train."
            )
            gas_particles_df.loc[
                gas_particles_df[parameter] < min_train, parameter
            ] = min_train * 9 / 10
        if max(gas_particles_df[parameter] > max_train):
            print(f"\nmax boundary problem in {parameter}")
            print(
                "Setting the gas particles having values bigger than the max_train equal to the max_train."
            )
            gas_particles_df.loc[
                gas_particles_df[parameter] > max_train, parameter
            ] = max_train * 9 / 10

    # Do it manually for smoothing length
    min_train = min(train_data_df["log_radius"])
    max_train = max(train_data_df["log_radius"])
    parameter = "log_smoothing_length"

    if min(gas_particles_df[parameter]) < min_train:
        print(f"\nmin boundary problem in {parameter}")
        print(
            "Setting the gas particles having values smaller than the min_train equal to the min_train."
        )
        gas_particles_df.loc[
            gas_particles_df[parameter] < min_train, parameter
        ] = min_train * 9 / 10
    if max(gas_particles_df[parameter] > max_train):
        print(f"\nmax boundary problem in {parameter}")
        print(
            "Setting the gas particles having values bigger than the max_train equal to the max_train."
        )
        gas_particles_df.loc[
            gas_particles_df[parameter] > max_train, parameter
        ] = max_train * 9 / 10

    # Do it manually for log_average_sobolev_smoothingLength
    min_train = min(train_data_df["log_radius"])
    max_train = max(train_data_df["log_radius"])
    parameter = "log_average_sobolev_smoothingLength"

    if min(gas_particles_df[parameter]) < min_train:
        print(f"\nmin boundary problem in {parameter}")
        print(
            "Setting the gas particles having values smaller than the min_train equal to the min_train."
        )
        gas_particles_df.loc[
            gas_particles_df[parameter] < min_train, parameter
        ] = min_train * 9 / 10
    if max(gas_particles_df[parameter] > max_train):
        print(f"\nmax boundary problem in {parameter}")
        print(
            "Setting the gas particles having values bigger than the max_train equal to the max_train."
        )
        gas_particles_df.loc[
            gas_particles_df[parameter] > max_train, parameter
        ] = max_train * 9 / 10

    ##############################################################################

    # Assuming gas_particles_df is your DataFrame and max_workers is defined
    gas_particles_df_array = split_dataframe(gas_particles_df, max_workers)

    ##############################################################################
    print(f"Line luminosity calculation started! It will take some time")

    timer_luminosity_prediction_start = time()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit each chunk of the DataFrame to the executor
        futures = [
            executor.submit(
                process_chunk,
                chunk,
                train_data_df,
                line_names_with_log,
            )
            for chunk in gas_particles_df_array
        ]

        # Collect the results
        results = [future.result() for future in futures]

    timer_luminosity_prediction_end = time()
    print("Line luminosity finished.")
    print(f"It took {round(timer_luminosity_prediction_end - timer_luminosity_prediction_start, 2)} seconds")

    print("Now results array is altered by flattening it and constructed the new array with the desired structure.")
    

    # Flattening the results list to a single dimension
    flattened_results = [item for sublist in results for item in sublist]    
    print("Results array are flattened")

    # #################################################
    print("Merging gas particles and calculated line luminosities")
    index_Lline_df = pd.DataFrame(
        flattened_results,
        columns=["index"] + line_names,
    )

    merged_df = gas_particles_df.merge(index_Lline_df, on="index", how="left")
    print("Merged!")


    ################ 
    # Interpolating other parameters
    other_columns_to_interpolate = [
        'log_column_density'
    ]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit each chunk of the DataFrame to the executor
        futures = [
            executor.submit(
                interpolate_other_galactic_properties,
                chunk,
                train_data_df,
                other_columns_to_interpolate,
            )
            for chunk in gas_particles_df_array
        ]

        # Collect the results
        results = [future.result() for future in futures]    

    # Flattening the results list to a single dimension
    flattened_results = [item for sublist in results for item in sublist]    
    print("Results array are flattened")        

    #################################################
    print("Merging gas particles and interpolated other parameters")
    index_other_parameters = pd.DataFrame(
        flattened_results,
        columns=[
            "index",
            'log_column_density'
        ],
    )

    merged_df = merged_df.merge(index_other_parameters, on="index", how="left")
    print("Merged!")


    ################
    # Change the units for the CO emission. Although there is I, in the merged_df, the units of lines are the luminosity units. I have multiplied them with the area.
    CO_lines_and_wavelengths = {
        "I_co_10": 2600.05e-6,  # meter
        "I_co_21": 1300.05e-6,
        "I_co_32": 866.727e-6,
        "I_co_43": 650.074e-6,
        "I_co_54": 520.08e-6,
        "I_co_65": 433.438e-6,
        "I_co_76": 371.549e-6,
        "I_co_87": 325.137e-6,
        "I_13co": 2719.67e-6,
    }

    for line in CO_lines_and_wavelengths: 
        conversion_factor = return_ergs_per_second2radio_units(
            rest_frequency=meters_to_Ghz_calculator(CO_lines_and_wavelengths[line])
        )
        merged_df[line] *= conversion_factor

    #################################################
    # Write to a file

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
    {base_file_dir}/{main_direcory}
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
    Column 35: log_column_density_cloudy [log(Msolar pc^-2)]
    """

    write_df = merged_df[gas_column_names + line_names + other_columns_to_interpolate]


    np.savetxt(fname=write_file_name, X=write_df, fmt="%.5e", header=header)

    print(f"File saved to: {write_file_name}")

    end_code = time()
    print(f"The code took {round((end_code - start_code) / 60, 2)} minutes to run")

    return 0


###########################################################################################################################################################################################################################

# Create different chunks of of dataframe to run them parallely
def split_dataframe(df, max_workers):
    n = len(df)
    chunk_size = -(
        -n // max_workers
    )  # Ceiling division to ensure all rows are included

    # Split the dataframe into chunks and store in an array
    return [df[i : i + chunk_size] for i in range(0, n, chunk_size)]


def process_chunk(gas_particles_df_chunk, train_data, line_luminosities_column_names):

    # prediction_center_column_names_for_train_data is for the training data
    prediction_center_column_names_for_train_data = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",
    ]

    # prediction_center_column_names_for_gas is for the gas particles
    scale_length = [
        # "log_radius",
        "log_average_sobolev_smoothingLength",
        # "log_smoothing_length",
        # "log_dummy_radius",
        ]

    prediction_center_column_names_for_gas = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
    ] + scale_length  

    tree = KDTree(
        train_data[prediction_center_column_names_for_train_data].to_numpy(),
    ) # Create a tree for the training data        

    try:
        index_Lline = []

        i = 0
        for index, gas_df in gas_particles_df_chunk.iterrows(): 
            if i == 0:
                initial_index = index

            predict_center = gas_df[prediction_center_column_names_for_gas]     

            # Get the indices of the closest training datapoints
            neighbour_indices = tree.query(predict_center, k=200, workers=1)[1] # Do not return distance with [1]
            # Get the closest training datapoints
            neighbours_train_data = train_data.iloc[neighbour_indices]

            interpolator = LinearNDInterpolator(
                neighbours_train_data[prediction_center_column_names_for_train_data].to_numpy(), 
                neighbours_train_data[line_luminosities_column_names].to_numpy(), 
            )                   

            interpolated_values = interpolator(predict_center)[0]  # log(erg s^-1 cm^-2)
            interpolated_intensity = 10**interpolated_values # erg s^-1 cm^-2

            area = gas_df['mass'] * constants.M_sun2gr / (gas_df[scale_length[0]] * constants.pc2cm * gas_df['density']) # cm^2
            L_line_from_interpolated_values = interpolated_intensity * area # erg s^-1 

            index_Lline.append(
                np.append(index, [L_line_from_interpolated_values[i] for i in range(len(L_line_from_interpolated_values))])
            )

            # Increase the counter
            i += 1 

            # Only print for the initial chunk. Don't need to print everything for all of the concurrent processes.
            if (initial_index == 0) and (i % 1e3 == 1):
                print(f"\nPrinting ony done for the first chunk.")
                print(f"{i} done. Left {len(gas_particles_df_chunk) - i}\n")

        return index_Lline

    except Exception as e:
        print("An exception occurred:", e)

        index_Lline.append(
            np.append(index, [np.nan for i in range(len(line_luminosities_column_names))])
        )

    return index_Lline

def interpolate_other_galactic_properties(gas_particles_df_chunk, train_data, columns_to_interpolate):

    # prediction_center_column_names_for_train_data is for the training data
    prediction_center_column_names_for_train_data = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
        "log_radius",
    ]

    # prediction_center_column_names_for_gas is for the gas particles
    scale_length = [
        # "log_radius",
        "log_average_sobolev_smoothingLength",
        # "log_smoothing_length",
        # "log_dummy_radius",
        ]

    prediction_center_column_names_for_gas = [
        "log_metallicity",
        "log_hden",
        "log_turbulence",
        "log_isrf",
    ] + scale_length  

    tree = KDTree(
        train_data[prediction_center_column_names_for_train_data].to_numpy(),
    ) # Create a tree for the training data        

    try:
        index_interpolated_values = []

        i = 0
        for index, gas_df in gas_particles_df_chunk.iterrows(): 
            if i == 0:
                initial_index = index

            predict_center = gas_df[prediction_center_column_names_for_gas]     

            # Get the indices of the closest training datapoints
            neighbour_indices = tree.query(predict_center, k=200, workers=1)[1] # Do not return distance with [1]
            # Get the closest training datapoints
            neighbours_train_data = train_data.iloc[neighbour_indices]

            interpolator = LinearNDInterpolator(
                neighbours_train_data[prediction_center_column_names_for_train_data].to_numpy(), 
                neighbours_train_data[columns_to_interpolate].to_numpy(), 
            )                   

            interpolated_values = interpolator(predict_center)[0]  # log(....)

            index_interpolated_values.append(
                np.append(index, [interpolated_values[i] for i in range(len(interpolated_values))])
            )

            # Increase the counter
            i += 1 

            # Only print for the initial chunk. Don't need to print everything for all of the concurrent processes.
            if (initial_index == 0) and (i % 1e3 == 1):
                print(f"\nPrinting ony done for the first chunk.")
                print(f"{i} done. Left {len(gas_particles_df_chunk) - i}\n")

        return index_interpolated_values

    except Exception as e:
        print("An exception occurred:", e)

        index_interpolated_values.append(
            np.append(index, [np.nan for i in range(len(columns_to_interpolate))])
        )

    return index_interpolated_values


def meters_to_Ghz_calculator(wavelength_in_meters):
    c = 299792458  # m/s
    frequency_in_Ghz = c / wavelength_in_meters * 1e-9
    return frequency_in_Ghz

def return_ergs_per_second2radio_units(rest_frequency):
    ergs_per_second2solar_luminosity = (3.826e33) ** (-1)  ## TODO: Probably there is a typo in the number. 3.846??
    solar_luminosity2radio_units = (3e-11 * (rest_frequency**3)) ** (-1)  # Rest frequency should be in Ghz
    ergs_per_second2radio_units = (ergs_per_second2solar_luminosity * solar_luminosity2radio_units)

    return ergs_per_second2radio_units

###########################################################################################################################################################################################################################

if __name__ == "__main__":
    galaxy_name = sys.argv[1]
    galaxy_type = sys.argv[2]
    redshift = sys.argv[3]
    max_workers = int(sys.argv[4])

    main(galaxy_name, galaxy_type, redshift, max_workers)

    ################################################################################################################################################################################################

    # # #######################
    # # zoom_in and zoom_in_tolgay
    # galaxy_names = [
    #     "m12b_res7100_md",
    #     "m12c_res7100_md",
    #     "m12f_res7100_md",
    #     "m12i_res7100_md",
    #     "m12m_res7100_md",
    #     "m12q_res7100_md",
    #     "m12r_res7100_md",
    #     "m12w_res7100_md",
    #     "m11d_r7100",
    #     "m11e_r7100",
    #     "m11h_r7100",
    #     "m11i_r7100",
    #     "m11q_r7100",
    # ]

    # for galaxy_name in galaxy_names:
    #     galaxy_type = "zoom_in"
    #     redshift = "3.0"
    #     max_workers = 20

    #     main(galaxy_name, galaxy_type, redshift, max_workers)

    # #######################
    # # Firebox
    # # for file_number in range(29):
    # for file_number in [29]:

    #     # FIREBox Galaxy
    #     galaxy_name = f"gal{file_number}"
    #     galaxy_type = "firebox"
    #     redshift = "3.0"
    #     max_workers = 30

    #     main(galaxy_name, galaxy_type, redshift, max_workers)

    ################################################################################################################################################################################################
