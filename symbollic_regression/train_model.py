import numpy as np
import pandas as pd
from pysr import PySRRegressor
import os


def main(which_redshifts:list, save_dir:str):

    ############################## Read the data ##############################

    base_dir = "/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files"

    # Define the file path
    file_path = f"{base_dir}/analyze_hden_metallicity_turbulence_isrf_radius/data"
    # file_name = "galactic_properties_averageSobolevH_hybridInterpolator_z0_usingIvalues.csv"

    file_name = "galactic_properties_smoothingLength_RBFInterpolator_z0_usingFvalues_voronoi_1e5.csv"
    z0 = read_galactic_properties(file_path=file_path, file_name=file_name)

    file_name = "galactic_properties_smoothingLength_RBFInterpolator_z1_usingFvalues_voronoi_1e5.csv"
    z1 = read_galactic_properties(file_path=file_path, file_name=file_name)

    file_name = "galactic_properties_smoothingLength_RBFInterpolator_z2_usingFvalues_voronoi_1e5.csv"
    z2 = read_galactic_properties(file_path=file_path, file_name=file_name)

    file_name = "galactic_properties_smoothingLength_RBFInterpolator_z3_usingFvalues_voronoi_1e5.csv"
    z3 = read_galactic_properties(file_path=file_path, file_name=file_name)


    # Merge the three redshifts 
    merged_df = pd.concat([z0, z1, z2, z3], ignore_index=True)

    # use only the redshifts that are in the list
    merged_df = merged_df[merged_df['redshift'].isin(which_redshifts)].copy()

    # Take logarithm of the data 
    data_df = take_log_of_the_data(data=merged_df)
    data_df_original = data_df.copy()

    ############################## Define the X columns and the target ##############################

    # Dictionary mapping feature column names to their LaTeX-formatted symbols
    x_vars = {
        "log_metallicity_gas_half_light_visible_z500pc": {
            "symbol": r"$\log_{10}(Z_{\mathrm{gas}}(R_{50}))$"
        },
        "log_metalicity_star_half_light_visible_z500pc": {
            "symbol": r"$\log_{10}(Z_{\star}(R_{50}))$"
        },
        "log_sfr_10Myr": {
            "symbol": r"$\log_{10}(\mathrm{SFR}_{10\,\mathrm{Myr}})$"
        },
        "log_star_mass": {
            "symbol": r"$\log_{10}(M_\star)$"
        },
        "log_gas_mass": {
            "symbol": r"$\log_{10}(M_{\mathrm{gas}})$"
        }, 
        "log_Pgas_half_light_visible_z500pc": {
            "symbol": r"$\log_{10}(P_{\mathrm{gas}}(R_{50}))$"
        },
        "log_Pstar_half_light_visible_z500pc": {
            "symbol": r"$\log_{10}(P_{\mathrm{star}}(R_{50}))$"
        },
        "log_Ptotal_half_light_visible_z500pc": {
            "symbol": r"$\log_{10}(P_{\mathrm{tot}}(R_{50}))$"
        },
        "log_halo_mass": {
            "symbol": r"$\log_{10}(M_{\mathrm{halo}})$"
        },
        "log_h2_weighted_temperature_mass_average_half_light_visible_z500pc":{
            "symbol": "Th2"
        },
    }

    # List of relevant columns for analysis
    x_columns = list(x_vars.keys())

    target_column = "log_L_co_10"

    X = data_df[x_columns].values
    y = data_df[target_column].values



    ############################## Save the model ##############################
    #  
    # Define the path to save the file 
    save_base_fdir = "/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/notebooks/symbollic_regression/outputs"
    save_dir = "model1"
    save_fdir = f"{save_base_fdir}/{save_dir}"

    # Check if file exits. If not create. If exists, delete the folder and create a new one
    if not os.path.exists(save_fdir):
        os.makedirs(save_fdir)
    else:
        import shutil
        shutil.rmtree(save_fdir)
        os.makedirs(save_fdir)

    # Define the settings for PySRRegressor
    settings_for_pysrRegressor= {
        "niterations": 40,
        "binary_operators": ["+", "-", "*", "/"],
        "unary_operators": ["pow_10(x) = 10^x"],  # Julia definition
        "extra_sympy_mappings": {
            "pow_10": lambda x: 10**x  # Python syntax
        },
        "model_selection": "best",      # Choose the best model by loss
        "elementwise_loss": "loss(prediction, target) = (prediction - target)^2",  # Mean squared error
        "verbosity": 0,
        "warm_start": True,
        "output_directory": save_fdir
    }


    # Define the symbolic regressor
    model = PySRRegressor(
        niterations             = settings_for_pysrRegressor['niterations'],                # Number of evolutionary iterations
        binary_operators        = settings_for_pysrRegressor['binary_operators'],           # Binary operators to use
        unary_operators         = settings_for_pysrRegressor['unary_operators'],            # Unary operators to use
        extra_sympy_mappings    = settings_for_pysrRegressor['extra_sympy_mappings'],       # Extra sympy mappings for custom functions
        model_selection         = settings_for_pysrRegressor['model_selection'],            # Model selection strategy
        elementwise_loss        = settings_for_pysrRegressor['elementwise_loss'],           # Loss function
        verbosity               = settings_for_pysrRegressor['verbosity'],                  # Verbosity level
        warm_start              = settings_for_pysrRegressor['warm_start'],                 # Whether to warm start the model
        output_directory        = settings_for_pysrRegressor['output_directory'],           # Directory to save the output
    )

    # Fit the model
    model.fit(X,y)

    # Write the options to a txt file. 
    with open(f"{save_fdir}/options.txt", "w") as f:
        for key, value in settings_for_pysrRegressor.items():
            f.write(f"{key}: {value}\n")

    # Append the used x and y columns to a txt file
    with open(f"{save_fdir}/options.txt", "a") as f:
        f.write(f"X columns: {x_columns}\n")
        f.write(f"Target column: {target_column}\n")

    # Write the redshifts used to a txt file
    redshifts_used = np.unique(data_df['redshift'].values)
    with open(f"{save_fdir}/options.txt", "a") as f:
        f.write(f"Redshifts used: {redshifts_used}\n")

    return None


def read_galactic_properties(file_path, file_name):

    print("I am in the read_galactic_properties function")

    # Read the DataFrame, skipping the header lines
    galaxies = pd.read_csv(
        f"{file_path}/{file_name}", 
        sep=',', 
    )


    galaxies['alpha_co_cloudy'] = galaxies['h2_mass_cloudy'] / galaxies['L_co_10'] 
    galaxies['X_co_cloudy'] = galaxies['alpha_co_cloudy'] * 6.3e19 


    galaxies['alpha_co_semi_analytical'] = galaxies['h2_mass_semi_analytical'] / galaxies['L_co_10']
    galaxies['X_co_semi_analytical'] = galaxies['alpha_co_semi_analytical'] * 6.3e19 

    return galaxies 


def take_log_of_the_data(data):

    print("I am in the take_log_of_the_data function")
    
    # Taking the logarithm of the data
    log_data = pd.DataFrame()

    not_log_taken_columns = [
        'name', 'galaxy_type', 'redshift', 'number_of_NaN_indices', 'alpha_co_cloudy', 'X_co_cloudy', 'alpha_co_semi_analytical', 'X_co_semi_analytical'
        ]
    
    columns = [col for col in data.columns if col not in not_log_taken_columns]
    
    for key in data.keys():
        if key in columns:
            log_data[f"log_{key}"] = np.log10(data[key])
        else: 
            log_data[key] = data[key]


    # Delete Nan and inf
    log_data = log_data[~log_data.isin([np.inf, -np.inf, np.nan]).any(axis=1)] 
    
    return log_data


if __name__ == "__main__":

    use_redshifts = [0, 1, 2, 3]
    save_dir = "z0_z1_z2_z3"

    main(which_redshifts=use_redshifts, save_dir=save_dir)