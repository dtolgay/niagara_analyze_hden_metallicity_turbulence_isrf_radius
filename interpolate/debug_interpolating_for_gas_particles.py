from math import ceil
import sys 
sys.path.append("/scratch/m/murray/dtolgay")
from tools import constants # type: ignore
import os
import pandas as pd # type: ignore
import numpy as np # type: ignore
import joblib # type: ignore

def main(files_info, interpolators, interpolator_name, interpolator_target_type, isrf_mult_factor):

    # print(f" --------------------------------------- {files_info['galaxy_name']} --------------------------------------- ")

    # Check if file exits. If exists, stop running with a message indicating that the file already exists.
    if os.path.exists(files_info["write_file_name"]):
        print(f"{files_info['write_file_name']} already exists. Stopping the run.")
        return None
    
    ### Read the centers
    gas_particles = read_cloudy_gas_particles(files_info)

    # TODO: Delete 
    gas_particles['isrf'] = gas_particles['isrf'] * isrf_mult_factor 

    # Note the gas particles columns before the interpolation
    gas_particles_base_columns = gas_particles.columns # This is the columns before the interpolation
    gas_particles_columns_write_to_file = list(gas_particles_base_columns) + list(interpolators[interpolator_target_type]["target_columns"]) # This is the full list of columns after the interpolation that will be written to the file 

    # Read the interpolator
    interpolator = joblib.load(interpolators[interpolator_target_type]["interpolator_path"])

    # Prepare the gas particles for interpolation
    base_interpolation_columns_for_gas_particles = [
        "metallicity",
        "hden",
        "turbulence",
        "isrf",
        "smoothing_length",
    ]
    used_shielding_length_name = base_interpolation_columns_for_gas_particles[-1]

    gas_particles, interpolation_columns_for_gas_particles = take_log_the_interpolation_centers_for_gas_particles(gas_particles, base_interpolation_columns_for_gas_particles)

    ########## Set the boundaries of gas particles within the limits of the interpolator
    if interpolator_name == "RegularGridInterpolator_linear":
        # If the radiation field of gas particles is smaller than the boundary of the interpolator,
        # set it to the minimum value of the interpolator's training data
        x_min = [grid.min() for grid in interpolator.grid]
        x_max = [grid.max() for grid in interpolator.grid]

    elif interpolator_name in ["GaussianProcessRegressor_rbf", "GaussianProcessRegressor_rational_quadratic"]:
        # If the radiation field of gas particles is smaller than the boundary of the interpolator,
        # set it to the minimum value of the interpolator's training data
        x_min = interpolator.X_train_.min(axis=0)
        x_max = interpolator.X_train_.max(axis=0)


    elif interpolator_name == "RBFInterpolator":
        # TODO: Don't know how to extract the information from the interpolator so set xmin=-inf xmax=inf
        x_min = [-2., -5.,  -3., -5.,  0.]
        x_max = [1.,  4.,  3.,  5.,  3.5,]

    # Put the boundaries of the interpolator into a dictionary
    print("\n\n Interpolator boundaries:")
    print(f"centers: {interpolation_columns_for_gas_particles}")
    print(f"x_min: {x_min}") 
    print(f"x_max: {x_max}")
    print("\n\n")    
    interpolator_boundaries = dict(zip(interpolation_columns_for_gas_particles, zip(x_min, x_max)))

    # If the properties of the gas particles are outside of the boundaries of the interpolator, 
    # set them to the boundaries of the interpolator
    for i, column in enumerate(interpolation_columns_for_gas_particles):
        gas_particles[column] = gas_particles[column].clip(
            lower=interpolator_boundaries[column][0], 
            upper=interpolator_boundaries[column][1]
        )

    # measure the time it takes to interpolate in minutes
    # Interpolate the gas particles
    gas_particles = interpolate(
        interpolator_name = interpolator_name,
        interpolator = interpolator, 
        gas_particles = gas_particles, 
        interpolation_columns_for_gas_particles = interpolation_columns_for_gas_particles, 
        target_columns = interpolators[interpolator_target_type]["target_columns"]
    )

    ##########

    # Interpolate the remaining NaN values with the NearestNDInterpolator
    # If gas particles with target columns have NaN indices interpolate them with the NearestNDInterpolator
    if gas_particles[interpolators[interpolator_target_type]["target_columns"]].isna().any().any():
        print("There are NaN values in the target columns. Interpolating with NearestNDInterpolator.")
        gas_particles = interpolate_for_nan_indices(
            data = gas_particles, 
            x_columns = interpolation_columns_for_gas_particles, 
            target_columns = interpolators[interpolator_target_type]["target_columns"], 
            interpolator_file_path = files_info["NearestNDInterpolator_file_path"]
            )

    # Calculate the line luminosities for the interpolated values if the interpolator type is line_emissions
    # Change the unit of the calculated CO luminosities   
    if interpolator_target_type == "line_emissions":  

        gas_particles = calculate_line_luminosities_for_interpolated_values(
            gas_particles = gas_particles, 
            target_columns = interpolators[interpolator_target_type]["target_columns"],
            used_shielding_length_name = used_shielding_length_name,
        )

        # Change the unit of the CO emissions from erg/s to K km/s pc2
        # Define CO wavelengths 
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

        gas_particles = change_unit_of_CO_emission(
            gas_particles = gas_particles, 
            lines_and_wavelengths = CO_lines_and_wavelengths
        )
    
    # Write the results to a file
    write_to_a_file(
        gas_particles=gas_particles,
        interpolators=interpolators, 
        interpolator_target_type=interpolator_target_type, 
        files_info=files_info,
        gas_particles_columns_write_to_file = gas_particles_columns_write_to_file,
        used_shielding_length_name = used_shielding_length_name,
    )

    return None 

def interpolate_for_nan_indices(data, x_columns, target_columns, interpolator_file_path):
    """
    Interpolates NaN values in the specified target columns of a DataFrame using a pre-trained interpolator.
    Parameters:
    data (pd.DataFrame): The input DataFrame containing NaN values to be interpolated.
    x_columns (list of str): List of column names to be used as features for interpolation.
    target_columns (list of str): List of column names where NaN values need to be interpolated.
    interpolator_file_path (str): Path to the file containing the pre-trained interpolator.
    Returns:
    pd.DataFrame: A DataFrame with NaN values in the target columns interpolated.
    """

    
    # Get the NaN rows 
    nan_rows = data[data.isna().any(axis=1)].copy()

    # Get the centers to interpolate 
    centers = nan_rows[x_columns].copy() 

    # Get the interpolator
    interpolator = joblib.load(interpolator_file_path)

    # interpolate 
    interpolated_values = interpolator(centers)

    # Create dataframe to store the values 
    interpolated_data = pd.DataFrame(centers, columns=x_columns)
    interpolated_data[target_columns] = interpolated_values
    
    
    # Merge the dataframes such that initially NaN values will be the interpolated values 
    data_merged = data.combine_first(interpolated_data)

    return data_merged

def interpolate(interpolator_name, interpolator, gas_particles, interpolation_columns_for_gas_particles, target_columns):

    if interpolator_name in ["RegularGridInterpolator_linear", "NearestNDInterpolator", "RBFInterpolator"]:
        # Interpolate the gas particles
        log_interpolated_values = interpolator(
            gas_particles[interpolation_columns_for_gas_particles].values
        )

    elif interpolator_name in ["GaussianProcessRegressor_rbf", "GaussianProcessRegressor_rational_quadratic"]:
        # Interpolate the gas particles
        log_interpolated_values = interpolator.predict(
            gas_particles[interpolation_columns_for_gas_particles].values # If I don't use values, it gives an error because of the column names where I use log_radius, but I interpolate according to log_smoothing_length
        )

    # Take the exponent of the interpolated values
    interpolated_values = 10 ** log_interpolated_values

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
    
    columns_names = []
    for column in interpolation_columns_for_gas_particles:
        new_column_name = f"log_{column}"
        columns_names.append(new_column_name)
        gas_particles[new_column_name] = np.log10(gas_particles[column])

    return gas_particles, columns_names

def calculate_line_luminosities_for_interpolated_values(gas_particles, target_columns, used_shielding_length_name):

    # Multiply with the gas particles area to find the luminosity of the line
    gas_particles['area'] = gas_particles['mass'] * constants.M_sun2gr / (gas_particles['density'] * (gas_particles[used_shielding_length_name] * constants.pc2cm))  # cm2

    # Multiply the interpolated fluxes with the area to find the luminosity of the line
    gas_particles[target_columns] = gas_particles[target_columns].multiply(gas_particles['area'], axis=0)

    return gas_particles

def meters_to_Ghz_calculator(wavelength_in_meters):
    c = 299792458  # m/s
    frequency_in_Ghz = c / wavelength_in_meters * 1e-9
    return frequency_in_Ghz

def return_ergs_per_second2radio_units(rest_frequency):
    ergs_per_second2solar_luminosity = constants.ergs2Lsolar # TODO: Everyone is citing different values. Ask the value of conversion from erg/s -> Lsolar
    solar_luminosity2radio_units = (3e-11 * (rest_frequency**3)) ** (-1)  # Rest frequency should be in Ghz
    ergs_per_second2radio_units = (ergs_per_second2solar_luminosity * solar_luminosity2radio_units)

    return ergs_per_second2radio_units

def change_unit_of_CO_emission(gas_particles, lines_and_wavelengths):

    for line in lines_and_wavelengths: 

        conversion_factor = return_ergs_per_second2radio_units(
            rest_frequency=meters_to_Ghz_calculator(lines_and_wavelengths[line])
        )
        gas_particles[line] *= conversion_factor

    return gas_particles

def write_to_a_file(gas_particles, interpolators, interpolator_target_type, files_info, gas_particles_columns_write_to_file, used_shielding_length_name):

    # Add the base header
    header = f"""
    Estimated according to:
    ---------------------
    log_metallicity
    log_hden
    log_turbulence
    log_isrf
    {used_shielding_length_name}
    ---------------------

    Used training centers:
    ---------------------
    {interpolators[interpolator_target_type]["interpolator_path"]}
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
        Column 18: isrf [G0]"""
    
    # Add the interpolator specific header
    if interpolator_target_type == "line_emissions":
        header_interpolater_specific = f"""
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
    elif interpolator_target_type == "abundance":
        header_interpolater_specific = f"""
        Column 19: fh2 [1]
        Column 20: fCO [1] Σco / ΣH2
        """
    elif interpolator_target_type == "temperature":
        header_interpolater_specific = f"""
        Column 19: Th2 [K]
        Column 20: Tco [K]
        Column 21: T [K]
        """

    header += header_interpolater_specific


    write_df = gas_particles[gas_particles_columns_write_to_file].copy()
    write_file_path = files_info["write_file_name"]
    number_of_significant_digits = ceil(np.log10(len(write_df)))  
 
    np.savetxt(fname=write_file_path, X=write_df, fmt=f"%.{number_of_significant_digits}e", header=header)

    print(f"File saved to: {write_file_path}")
    
    return None

if __name__ == "__main__":

    interpolators_base_fdir = "/scratch/m/murray/dtolgay/cloudy_runs/interpolators"

    interpolator_names = [
        # "GaussianProcessRegressor_rbf",
        # "GaussianProcessRegressor_rational_quadratic",
        # "LinearNDInterpolator",
        # "NearestNDInterpolator",
        "RBFInterpolator",
        # "RegularGridInterpolator_linear",
    ]
    interpolator_name = interpolator_names[0] 

    interpolators = {
        "temperature": {
            "file_name": "temperatures.txt",
            "target_columns": ["Th2", "Tco", "T"],
            "interpolator_path": f"{interpolators_base_fdir}/{interpolator_name}_temperature.joblib",
            "write_file_name": f"temperature_smoothingLength_{interpolator_name}",
            "interpolator_identifier_name": "temperature",
        },
        "line_emissions": {
            "file_name": "I_line_values_without_reversing.txt",
            "target_columns": [
                "L_ly_alpha",
                "L_h_alpha",
                "L_h_beta",
                "L_co_10",
                "L_co_21",
                "L_co_32",
                "L_co_43",
                "L_co_54",
                "L_co_65",
                "L_co_76",
                "L_co_87",
                "L_13co",
                "L_c2",
                "L_o3_88",
                "L_o3_5006",
                "L_o3_4958",
            ],
            "interpolator_path": f"{interpolators_base_fdir}/{interpolator_name}_flux.joblib",
            "write_file_name": f"line_emissions_smoothingLength_{interpolator_name}",
            "interpolator_identifier_name": "flux",
        },
        "abundance": {
            "file_name": "other_properties.csv",
            "target_columns": [
                "fh2",
                "fCO",
            ],
            "interpolator_path": f"{interpolators_base_fdir}/{interpolator_name}_abundance.joblib",
            "write_file_name": f"abundance_smoothingLength_{interpolator_name}",
            "interpolator_identifier_name": "abundance",
        }
    }

    base_fdir =  "/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/runs_hden_radius"
    galaxy_name = sys.argv[1]
    galaxy_type = sys.argv[2]
    redshift = sys.argv[3]
    interpolator_target_type = sys.argv[4] # temperature, line_emissions, abundance
    isrf_mult_factor = float(sys.argv[5]) # TODO:


    galaxy_info = {
        "base_fdir": base_fdir,
        "galaxy_name": galaxy_name,
        "galaxy_type": galaxy_type,
        "redshift": redshift,
        "directory": "voronoi_1e6",
    }

    gas_particles_path = f"{base_fdir}/{galaxy_info['galaxy_type']}/z{galaxy_info['redshift']}/{galaxy_info['galaxy_name']}/{galaxy_info['directory']}"
    # gas_particles_path = f"/scratch/m/murray/dtolgay/cloudy_runs/z_0/m12i_res7100_md_test"

    files_info = {
        "write_file_name": f"{gas_particles_path}/{interpolator_target_type}_{interpolator_name}_smoothingLength_isrfMult{isrf_mult_factor}.txt", # TODO: Change the name of the file
        "gas_particles_path": gas_particles_path,
        "NearestNDInterpolator_file_path": f"/scratch/m/murray/dtolgay/cloudy_runs/interpolators/NearestNDInterpolator_{interpolators[interpolator_target_type]['interpolator_identifier_name']}.joblib", 
        "galaxy_name": galaxy_name,
    }

    main(files_info, interpolators, interpolator_name, interpolator_target_type, isrf_mult_factor)