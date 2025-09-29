import sys
sys.path.append("/scratch/m/murray/dtolgay/")

import numpy as np
import pandas as pd 
from time import time

from tools import functions, constants  # type: ignore 
from tools import functions_readfiles as readfiles # type: ignore
from tools.functions_create_galactic_properties import sfr_calculator, halo_mass_calculator, mass_average # type: ignore
from tools.functions_create_galactic_properties import get_halo_mass_from_previously_calculated_file # type: ignore

import astropy.units as units
import galactic_props_functions as galactic_prop_funcs

def main(redshift, running_cluster:str):

    # base file name to start searching files
    directory_name = "voronoi_1e5"
    galaxies_list = []

    ####################################### For particle split 
    galaxy_names = [
        "m12i_r880_md"
    ]

    for galaxy_name in galaxy_names:
        try: 
            galaxies_list.append(get_galactic_properties_for_a_galaxy(
                galaxy_name=galaxy_name, 
                galaxy_type="particle_split", 
                redshift=redshift, 
                directory_name=directory_name
            ))
        except Exception as e:
            print(f"Exception occured: \n{e}")    

    ####################################### For firebox 

    for i in range(int(1e3)):
    # for i in range(int(10)):
    # for i in [0]:
        try:
            galaxies_list.append(get_galactic_properties_for_a_galaxy(
                galaxy_name=f"gal{i}", 
                galaxy_type="firebox", 
                redshift=redshift, 
                directory_name=directory_name
            ))
            
        except Exception as e:
            print(f"Exception occured: \n{e}") 

    ####################################### For zoom_in
    galaxy_names = [
        "m12b_res7100_md", 
        "m12c_res7100_md",
        "m12f_res7100_md",
        "m12i_res7100_md",
        "m12m_res7100_md",
        "m12q_res7100_md",
        "m12r_res7100_md",
        "m12w_res7100_md",
        "m11d_r7100",
        "m11e_r7100",
        "m11h_r7100",
        "m11i_r7100",
        "m11q_r7100",            
    ]

    for galaxy_name in galaxy_names:
        try:
            galaxies_list.append(get_galactic_properties_for_a_galaxy(
                galaxy_name=galaxy_name, 
                galaxy_type="zoom_in", 
                redshift=redshift, 
                directory_name=directory_name
            ))
            
        except Exception as e:
            print(f"Exception occured: \n{e}")   

    ####################################### Turn list to a pandas dataframe
    galaxies = pd.DataFrame(galaxies_list)


    ####################################### Write to a file

    if running_cluster == "cita":
        save_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/data"    
    elif running_cluster == "niagara":
        save_dir = "/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/data"

    file_path = f"{save_dir}/galactic_properties_smoothingLength_RBFInterpolator_z{int(float(redshift))}_usingFvalues_{directory_name}.csv"

    # Append the DataFrame to the file
    galaxies.to_csv(file_path, mode='w', sep=',', index=False, float_format='%.5e')

    print(f"File saved to {file_path}")


    return 0

#### Functions defined below 


def get_galactic_properties_for_a_galaxy(galaxy_name:str, galaxy_type:str, redshift:str, directory_name:str, running_cluster:str):
    
    print(f"--------------  {galaxy_name}  --------------")
    
    start = time()

    if running_cluster == "cita":
        base_fdir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
        read_dir = f"{base_fdir}/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"
        skirt_files_dir = read_dir
    elif running_cluster == "niagara":
        base_fdir = "/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/runs_hden_radius"
        read_dir = f"{base_fdir}/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"
        skirt_files_dir = "/gpfs/fs0/scratch/r/rbond/dongwooc/scratch_rwa/doga/runs_hden_radius"


    file_names = {
        "interpolated_lines": "line_emissions_RBFInterpolator_smoothingLength.txt",
        "abundance": "abundance_RBFInterpolator_smoothingLength.txt",
        "temperature": "temperature_RBFInterpolator_smoothingLength.txt",
        "semi_analytical": "semi_analytical_smoothingLength_cf_1.txt",
        "i0_total_fits": f"{skirt_files_dir}/{galaxy_name}_i0_total.fits",
        "skirt_wavelengths": f"{skirt_files_dir}/{galaxy_name}_grid_radiationField_wavelengths.dat",
    }

    
    # Read gas particles
    base_dir="/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
    path_without_file_name = f'{base_dir}/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}'

    # gas, lines = readfiles.read_interpolated_Lline(
    #     galaxy_name = galaxy_name, 
    #     galaxy_type = galaxy_type, 
    #     redshift = redshift, 
    #     directory_name = directory_name, 
    #     file_name = file_names['interpolated_lines']
    # )

    gas, lines = readfiles.read_interpolated_files_usingFilePath2(
        path = f"{path_without_file_name}/{file_names['interpolated_lines']}",
        interpolation_type = "line_emissions",
    )

    # gas_other_properties = readfiles.read_otherProperties(
    #         galaxy_name = galaxy_name, 
    #         galaxy_type = galaxy_type, 
    #         redshift = redshift, 
    #         directory_name = directory_name, 
    #         file_name = file_names['otherProperties'],
    #     )   
    gas_abundances, file_specific_columns_abundance = readfiles.read_interpolated_files_usingFilePath2(
        path = f"{path_without_file_name}/{file_names['abundance']}", 
        interpolation_type = "abundance",
    )
    gas_abundances = gas_abundances[["index"] + file_specific_columns_abundance]

    gas_temperature, file_specific_columns_temperature = readfiles.read_interpolated_files_usingFilePath2(
        path = f"{path_without_file_name}/{file_names['temperature']}", 
        interpolation_type = "temperature",
    )
    gas_temperature = gas_temperature[["index"] + file_specific_columns_temperature]

    ## Convert index to string
    gas['index'] = gas['index'].astype(str)  # or .astype(int) if applicable
    gas_abundances['index'] = gas_abundances['index'].astype(str)
    gas_temperature['index'] = gas_temperature['index'].astype(str)

    # Merge gas particles with abundances 
    gas = gas.merge(gas_abundances, how='inner', on=['index'], validate='one_to_one')

    # Merge gas particles with temperature
    gas = gas.merge(gas_temperature, how='inner', on=['index'], validate='one_to_one')

    # Read star particles
    star = readfiles.read_comprehensive_star_particles(galaxy_name, galaxy_type, redshift, directory_name)

    # Read semi_analytical_average_sobolev_smoothingLength.txt
    semi_analytical = readfiles.read_semianalytical_file_2(
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift, 
        directory_name = directory_name, 
        file_name = file_names['semi_analytical'],
    )

    # Read sed file 
    sed = readfiles.read_skirt_sed_file(galaxy_name, galaxy_type, redshift, directory_name = directory_name, inclination = "0")

    ############ Calculate galactic properties 
    ## sfr
    sfr_instantaneous = np.sum(gas["sfr"]) 
    sfr_5Myr = sfr_calculator(star_df = star, within_how_many_Myr = 5)
    sfr_10Myr = sfr_calculator(star_df = star, within_how_many_Myr = 10)
    sfr_100Myr = sfr_calculator(star_df = star, within_how_many_Myr = 100)

    ## gas mass
    total_gas_mass = np.sum(gas["mass"])
    total_star_mass = np.sum(star["mass"])

    ###### Metalicity
    
    ## Filtering 
    
    # Find the indices of gas and star particles within 8 kpc from the center
    R_inner = 8e3 #kpc 
    indices_gas_inner_galaxy = np.where(np.sqrt(gas['x']**2 + gas['y']**2 + gas['z']**2) < R_inner)[0] 
    indices_star_inner_galaxy = np.where(np.sqrt(star['x']**2 + star['y']**2 + star['z']**2) < R_inner)[0] 
    gas_inner_galaxy = gas.iloc[indices_gas_inner_galaxy]
    star_inner_galaxy = star.iloc[indices_star_inner_galaxy]

    # Use only the ism particles. Put a condition on density.   
    gas_hden_higher_than_1eminus2_condition = gas['hden'] > 1e-2
    gas_hden_higher_than_1 = gas['hden'] > 1e0

    ## calculation
    # all gas particles
    average_gas_metallicity = mass_average(
        property_array = gas['metallicity'].to_numpy(),
        mass = gas['mass'].to_numpy(),
    )

    average_star_metallicity = mass_average(
        property_array = star['metallicity'].to_numpy() / constants.solar_metallicity,
        mass = star['mass'].to_numpy(),
    )    

    # inner galaxy
    gas_metallicity_inner_galaxy = mass_average(
        property_array = gas_inner_galaxy['metallicity'].to_numpy(),
        mass = gas_inner_galaxy['mass'].to_numpy(),
    )

    star_metallicity_inner_galaxy = mass_average(
        property_array = star_inner_galaxy['metallicity'].to_numpy() / constants.solar_metallicity,
        mass = star_inner_galaxy['mass'].to_numpy()
    )

    # hden condition
    filtered_gas = gas[gas_hden_higher_than_1eminus2_condition].copy()
    gas_metallicity_hden_higher_than_1eminus2 = mass_average(
        property_array = filtered_gas['metallicity'].to_numpy(),
        mass = filtered_gas['mass'].to_numpy(),
    )

    filtered_gas = gas[gas_hden_higher_than_1].copy()
    gas_metallicity_hden_higher_than_1 = mass_average(
        property_array = filtered_gas['metallicity'].to_numpy(),
        mass = filtered_gas['mass'].to_numpy(),
    )

    ## line_luminosities
    # Krumholz - semi_analytical 
    semi_analytical_Lco = np.sum(semi_analytical['L_co'])
    # Cloudy
    total_line_luminosities = {}
    for line in lines: 
        total_line_luminosities[line] = sum(gas[line])
    

    ## Molecular gas mass
    h2_mass_semi_anaytical = np.sum(semi_analytical["Mh2"]) # if h2 mass is NaN probably it is zero already.
    h2_mass_cloudy = np.sum(gas['fh2'] * gas['mass'])

    # CO mass
    gas['mass_co'] = gas['fCO'] * gas['mass']
    co_mass_cloudy = np.sum(gas['mass_co'].to_numpy())

    ## halo mass
    try:
        Mhalo = halo_mass_calculator(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift
        )
    except Exception as e:
        Mhalo = get_halo_mass_from_previously_calculated_file(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift,
        )

    # NaN indices
    nan_indices = np.where(np.isnan(gas["L_co_10"]))[0]    

    #### Calculate infrared luminosity from sed data 
    # Define the bands 
    bands = {
        "all": {
            "min_wavelength": 1e-6,
            "max_wavelength": 1e6,    
        },
        "fir": {
            "min_wavelength": 3,
            "max_wavelength": 1000,
        }
    }    

    ## SED
    distance = 10 # Mpc 


    L_fir = functions.calculate_luminosity_from_sed(
        sed = sed, 
        min_wavelength = bands['fir']['min_wavelength'], 
        max_wavelenght = bands['fir']['max_wavelength'], 
        distance_in_Mpc = distance
    ) # erg / s

    L_skirt_full_wavelength_range = functions.calculate_luminosity_from_sed(
        sed = sed, 
        min_wavelength = bands['all']['min_wavelength'], 
        max_wavelenght = bands['all']['max_wavelength'], 
        distance_in_Mpc = distance
    ) # erg / s    


    # Get the half light radius
    half_light_wavelength_bands = {
        "visible": {
            "max": 700 * units.nm, # nm 
            "min": 400 * units.nm, # nm
        },
        "infrared": {
            "max": 1 * units.mm,    # mm
            "min": 780 * units.nm,  # nm
        },
    }


    R_half_light_visible = galactic_prop_funcs.half_light_radius(band = half_light_wavelength_bands['visible'], file_names = file_names, plot = False) # pc  
    R_half_light_stellar_mass = galactic_prop_funcs.half_mass_radius(particles_df = star) # pc

    ## Define the conditions for filtering the particles 
    conditions = {
        "20kpc": {
            "r_max" : 20e3, # pc
            "z_max" : 2e3, # pc 
        },
        "10kpc": {
            "r_max": 10e3, # pc
            "z_max": 2e3, # pc 
        },
        "r50visible_z500pc": {
            "r_max": R_half_light_visible, # pc
            "z_max": 500, # pc
        }
    }    


    ### Get the properties at the half light radius
    gas_half_light_visible_z500pc = galactic_prop_funcs.filter_particles(
        particles_df = gas.copy(),
        condition = conditions['r50visible_z500pc'],
    )  
    star_half_light_visible_z500pc = galactic_prop_funcs.filter_particles(
        particles_df = star.copy(),
        condition = conditions['r50visible_z500pc'],
    )
    

    # Calculate the pressures 
    Pgas_10kpc, Pstar_10kpc, Ptotal_10kpc, R_half_light_visible = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas, 
        star_particles = star, 
        filtering_condition = conditions['10kpc'], 
        R_half_light = R_half_light_visible,
        is_plot = False, 
    )    

    Pgas_20kpc, Pstar_20kpc, Ptotal_20kpc, R_half_light_visible  = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas, 
        star_particles = star, 
        filtering_condition = conditions['20kpc'], 
        R_half_light = R_half_light_visible,
        is_plot = False, 
    )        

    Pgas_half_light_visible_z500pc, Pstar_half_light_visible_z500pc, Ptotal_half_light_visible_z500pc, R_half_light_visible = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas,
        star_particles = star,
        filtering_condition = conditions['r50visible_z500pc'],
        R_half_light = R_half_light_visible,
        is_plot = False,
    )

    ### Calculate the mass averaged properties within a certain condition
    conditions2 = {
        "disk_half_light_visible": {
            "r_max": R_half_light_visible, # pc 
            "z_max": 500, # pc
        }
    }

    # Calculate the mass-weighted average for metallicities 
    metallicity_gas_half_light_visible_z500pc = calculate_mass_weighted_average_within_condition(
        particles = gas,
        condition = conditions2['disk_half_light_visible'],
        property_name = "metallicity"
    )

    metalicity_star_half_light_visible_z500pc = calculate_mass_weighted_average_within_condition(
        particles = star,
        condition = conditions2['disk_half_light_visible'],
        property_name = "metallicity" 
    ) / constants.solar_metallicity

    # Calculate mass averaged hden
    hden_mass_average_half_light_visible_z500pc = calculate_mass_weighted_average_within_condition(
        particles = gas,
        condition = conditions2['disk_half_light_visible'],
        property_name = "hden"
    )

    # Calculate mass averaged temperature
    fire_temperature_mass_average_half_light_visible_z500pc = calculate_mass_weighted_average_within_condition(
        particles = gas,
        condition = conditions2['disk_half_light_visible'],
        property_name = "temperature"
    )

    # Calculate the CO temperature
    co_weighted_temperature_mass_average_half_light_visible_z500pc = calculate_mass_weighted_average_within_condition(
        particles = gas,
        condition = conditions2['disk_half_light_visible'],
        property_name = "Tco"
    )

    # Calculate the H2 temperature
    h2_weighted_temperature_mass_average_half_light_visible_z500pc = calculate_mass_weighted_average_within_condition(
        particles = gas,
        condition = conditions2['disk_half_light_visible'],
        property_name = "Th2"
    )

    # Calculate the column density for certain gas particles determined by the condition
    column_density_half_light_visible_z500pc = calculate_column_density_within_condition(
        particles = gas,
        condition = conditions2['disk_half_light_visible']
    )

    # Calculate the line luminosities at the half light radius 
    total_line_luminosities_R50 = {}
    for line in lines:
        total_line_luminosities_R50[f"{line}_half_light_visible_z500pc"] = sum(gas_half_light_visible_z500pc[line])


    # Create a dictionary 
    other_properties = {
        "name": f"{galaxy_name}",
        "galaxy_type": f"{galaxy_type}",
        "redshift": f"{redshift}",
        "sfr": sfr_instantaneous, # Msolar / year
        "sfr_5Myr": sfr_5Myr, # Msolar / year
        "sfr_10Myr": sfr_10Myr, # Msolar / year
        "sfr_100Myr": sfr_100Myr, # Msolar / year
        "gas_mass": total_gas_mass, # Msolar
        "star_mass": total_star_mass, # Msolar
        "gas_average_metallicity": average_gas_metallicity, # Zsolar
        "star_average_metallicity": average_star_metallicity, # Zsolar
        "gas_metallicity_inner_galaxy": gas_metallicity_inner_galaxy, # Zsolar
        "star_metallicity_inner_galaxy": star_metallicity_inner_galaxy, # Zsolar
        "h2_mass_semi_analytical": h2_mass_semi_anaytical, # Msolar
        "h2_mass_cloudy": h2_mass_cloudy, # Msolar
        "co_mass_cloudy": co_mass_cloudy, # Msolar
        "halo_mass": Mhalo, # Msolar,
        "L_fir": L_fir, # erg / s
        "L_skirt_full_wavelength_range": L_skirt_full_wavelength_range, # erg / s
        "number_of_NaN_indices": len(nan_indices),
        "gas_metallicity_hden_higher_than_1eminus2": gas_metallicity_hden_higher_than_1eminus2, # Zsolar
        "gas_metallicity_hden_higher_than_1": gas_metallicity_hden_higher_than_1, # Zsolar
        "L_co_10_semi_analytical": semi_analytical_Lco, # K km s^-1 pc^2 
        "Pgas_10kpc": Pgas_10kpc, # K / cm3
        "Pstar_10kpc": Pstar_10kpc, # K / cm3
        "Ptotal_10kpc": Ptotal_10kpc, # K / cm3 
        "Pgas_20kpc": Pgas_20kpc, # K / cm3
        "Pstar_20kpc": Pstar_20kpc, # K / cm3
        "Ptotal_20kpc": Ptotal_20kpc, # K / cm3  
        "Pgas_half_light_visible_z500pc": Pgas_half_light_visible_z500pc, # K / cm3
        "Pstar_half_light_visible_z500pc": Pstar_half_light_visible_z500pc, # K / cm3
        "Ptotal_half_light_visible_z500pc": Ptotal_half_light_visible_z500pc, # K / cm3
        "R_half_light_visible": R_half_light_visible, # pc
        "R_half_light_stellar_mass": R_half_light_stellar_mass, # pc
        "metallicity_gas_half_light_visible_z500pc": metallicity_gas_half_light_visible_z500pc, # Zsolar
        "metalicity_star_half_light_visible_z500pc": metalicity_star_half_light_visible_z500pc, # Zsolar
        "hden_mass_average_half_light_visible_z500pc": hden_mass_average_half_light_visible_z500pc, # cm^-3
        "fire_temperature_mass_average_half_light_visible_z500pc": fire_temperature_mass_average_half_light_visible_z500pc, # K
        "co_weighted_temperature_mass_average_half_light_visible_z500pc": co_weighted_temperature_mass_average_half_light_visible_z500pc, # K
        "h2_weighted_temperature_mass_average_half_light_visible_z500pc": h2_weighted_temperature_mass_average_half_light_visible_z500pc, # K
        "column_density_half_light_visible_z500pc": column_density_half_light_visible_z500pc, # Msolar / pc2    
    }     

    # Merge two dictionaries
    galactic_properties = {**other_properties, **total_line_luminosities, **total_line_luminosities_R50}

    stop = time()
    
    print(f"For {galaxy_name}, it took {round((stop-start)/60, 3)} minutes")

    return galactic_properties


def run_for_single_galaxy(galaxy_name, galaxy_type, redshift, base_dir):
        try:
            return get_galactic_properties_for_a_galaxy(
                galaxy_name=galaxy_name, 
                galaxy_type=galaxy_type, 
                redshift=redshift, 
                base_dir=base_dir
            )
            
        except Exception as e:
            print(f"Exception occured for galaxy {galaxy_name}: \n{e}")    
            return None

def calculate_mass_weighted_average_within_condition(particles: pd.DataFrame, condition: dict, property_name: str) -> float:
    """
    Calculate the mass-weighted average for a property within a certain condition.
    The condition includes both radius and height.
    """
    r_max = condition.get('r_max', np.inf) # pc 
    z_max = condition.get('z_max', np.inf) # pc

    # Apply the condition
    within_condition = particles[
        (np.sqrt(particles['x']**2 + particles['y']**2) < r_max) &
        (np.abs(particles['z']) < z_max)
    ]

    if within_condition.empty:
        return np.nan

    # Calculate mass-weighted average
    property_array = within_condition[property_name].to_numpy()
    mass_array = within_condition["mass"].to_numpy()
    mass_weighted_average = np.sum(property_array * mass_array) / np.sum(mass_array)

    return mass_weighted_average

def calculate_column_density_within_condition(particles: pd.DataFrame, condition: dict) -> float:
    """
    Calculate the column density for certain gas particles determined by the condition.
    The condition includes both radius and height.
    """
    r_max = condition.get('r_max', np.inf) # pc 
    z_max = condition.get('z_max', np.inf) # pc

    # Apply the condition
    within_condition = particles[
        (np.sqrt(particles['x']**2 + particles['y']**2) < r_max) &
        (np.abs(particles['z']) < z_max)
    ]

    if within_condition.empty:
        return np.nan

    # Calculate column density
    area = np.pi * r_max**2
    column_density = np.sum(within_condition['mass']) / area # Msolar / pc2
  
    return column_density

if __name__ == "__main__":
    redshift = sys.argv[1]
    redshift = "{:.1f}".format(float(redshift))

    running_cluster = "niagara" # "cita" or "niagara"

    main(redshift, running_cluster=running_cluster)