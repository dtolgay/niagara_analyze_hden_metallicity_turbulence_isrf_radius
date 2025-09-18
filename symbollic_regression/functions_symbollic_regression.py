import pysr
import numpy as np 
import pandas as pd 
from matplotlib import pyplot as plt

def calculate_std_for_every_model(model: pysr.sr.PySRRegressor, X, y_expected):

    for i, model_number in enumerate(model.equations_.index):        
        # Predict the expected y values 
        y_pred = model.predict(X = X, index = model_number)

        # Calculate the residuals 
        residuals = y_expected - y_pred

        # Calculate the std 
        std = np.std(residuals)

        # Store the std in the model
        model.equations_.loc[model_number, "std"] = std

    return model

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
    numeric_cols = log_data.select_dtypes(include=[np.number]).columns
    log_data = log_data[np.isfinite(log_data[numeric_cols]).all(axis=1)]
    
    return log_data

def plot_psyr_results(
        model: pysr.sr.PySRRegressor, 
        model_number:int, 
        actual_y:np.ndarray, 
        X:np.ndarray, 
        ax_main:plt.axis, 
        ax_resid:plt.axis, 
        x_columns:list, 
        x_vars:dict,
        convert_to_log_scale: bool = False,
    ):

    """
    Plots the predicted results from a PySRRegressor model against the actual values and visualizes the residuals.
    Parameters:
        model (pysr.sr.PySRRegressor): The symbolic regression model used for predictions.
        model_number (int): The index of the model to be used for predictions.
        actual_y (np.ndarray): The actual target values to compare against predictions.
        X (np.ndarray): The input features used for making predictions.
        ax_main (plt.axis): The axis object for the main prediction plot.
        ax_resid (plt.axis): The axis object for the residuals plot.
        x_columns (list): A list of column names corresponding to the input features.
        x_vars (dict): A dictionary mapping column names to their respective LaTeX symbols.
        convert_to_log_scale (bool, optional): Calculates the residuals by taking the log10. Also converts linear values to log values. Defaults to False.
    Returns:
        None: This function does not return any value; it only generates plots.
    """

    y_pred = model.predict(X, model_number)
    equation_raw = str(model.equations_.iloc[model_number]['sympy_format'])
    if convert_to_log_scale:
        residuals = np.log10(actual_y) - np.log10(y_pred)
        # Delete the nan's 
        non_nan_indices = np.isfinite(residuals)        
        residuals = residuals[non_nan_indices]
        actual_y = np.log10(actual_y[non_nan_indices])
        y_pred = np.log10(y_pred[non_nan_indices])
    else:
        residuals = actual_y - y_pred


    # Find the std of the residuals 
    std = np.std(residuals)
    
    # x = y 
    x_dummy = np.linspace(min(actual_y), max(actual_y), num=100)
    y_dummy = x_dummy
    
    # Replace x0, x1, ... with LaTeX labels
    for i, col in enumerate(x_columns):
        label = x_vars[col]["symbol"]
        equation_raw = equation_raw.replace(f'x{i}', label)

    # Plot actual vs predicted
    ax_main.scatter(actual_y, y_pred, label=equation_raw, s=10)
    ax_main.plot(x_dummy, y_dummy, label="x=y", color = "red")
    ax_main.set_ylabel("Prediction")
    ax_main.set_title(f"Model {model_number} σ = {std:.2f}")
    ax_main.legend(fontsize='x-small', loc='upper left')
    ax_main.tick_params(labelbottom=False)
    # Plot residuals
    ax_resid.scatter(actual_y, residuals, color='gray', s=10)
    ax_resid.axhline(0, color='black', linestyle='--', linewidth=1)
    ax_resid.set_xlabel("FIRE Calculation")
    ax_resid.set_ylabel("Residuals")

    # if log_scale:
    #     ax_main.set_xscale('log')
    #     ax_main.set_yscale('log')
    #     ax_resid.set_xscale('log')

        
    ax_resid.set_ylim([-2, 2])

    return None