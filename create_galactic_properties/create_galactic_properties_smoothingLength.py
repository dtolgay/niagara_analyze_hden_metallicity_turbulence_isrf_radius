import sys
import numpy as np
import pandas as pd 
from time import time

from tools_tolgay import functions, constants  # type: ignore 
from tools_tolgay import functions_readfiles as readfiles # type: ignore
# from tools_tolgay.functions_create_galactic_properties import sfr_calculator, halo_mass_calculator, mass_average # type: ignore
# from tools_tolgay.functions_create_galactic_properties import get_halo_mass_from_previously_calculated_file # type: ignore
import tools_tolgay.functions_create_galactic_properties as galactic_prop_funcs

import astropy.units as units
# import galactic_props_functions as galactic_prop_funcs

from concurrent.futures import ProcessPoolExecutor, as_completed

def main(redshift, max_workers):

    directory_name = "voronoi_1e5"
    all_args = []

    ###### particle_split ######
    all_args.extend([(name, "particle_split", redshift, directory_name) for name in [
        "m12i_r880_md"
    ]])

    ###### firebox ######
    # all_args.extend([(f"gal{i}", "firebox", redshift, directory_name) for i in range(0, 1)]) # TODO: Delete debugging purpose
    all_args.extend([(f"gal{i}", "firebox", redshift, directory_name) for i in range(int(1e3))])


    ###### zoom_in ######
    all_args.extend([(name, "zoom_in", redshift, directory_name) for name in [
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
        "m11q_r7100"
    ]])

    #########################

    galaxies_list = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(wrapper, args) for args in all_args]
        for future in as_completed(futures):
            result = future.result()
            if result is not None:
                galaxies_list.append(result)


    ####################################### Turn list to a pandas dataframe

    galaxies = pd.DataFrame(galaxies_list)
    print(f"galaxies: {galaxies}")
    # Order galaxies according to halo mass 
    try:
        galaxies = galaxies.sort_values(by='halo_mass', ascending=False).reset_index(drop=True)
    except Exception as e:
        print(f"Exception occurred while sorting the galaxies by halo mass:\n{e}")
        sys.exit(1)

    ####################################### Write to a file

    ## TODO: Save the file 
    save_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/data"
    file_path = f"{save_dir}/galactic_properties_smoothingLength_RBFInterpolator_z{int(float(redshift))}_using_luminosity_per_mass_values_{directory_name}.csv"

    # Append the DataFrame to the file
    galaxies.to_csv(file_path, mode='w', sep=',', index=False, float_format='%.5e')

    print(f"File saved to {file_path}")


    return 0

#### Functions defined below 

def wrapper(args):
    """Wrapper function to unpack arguments."""
    galaxy_name, galaxy_type, redshift, directory_name = args
    try:
        return get_galactic_properties_for_a_galaxy(
            galaxy_name=galaxy_name,
            galaxy_type=galaxy_type,
            redshift=redshift,
            directory_name=directory_name,
        )
    except Exception as e:
        print(f"Exception occurred for {galaxy_name} ({galaxy_type}):\n{e}")
        return None


def read_all_files(galaxy_name:str, galaxy_type:str, redshift:str, directory_name:str):

    base_fdir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
    read_dir = f"{base_fdir}/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"

    file_names = {
        # "interpolated_lines": "L_line_smoothingLength_hybridInterpolator_flux2Luminosity.txt",
        # "interpolated_lines": "line_emissions_RBFInterpolator_smoothingLength.txt",
        "interpolated_lines": "luminosity_from_luminosity_per_mass_RBFInterpolator_smoothingLength.txt",
        # "otherProperties": "otherProperties_smoothingLength_hybridInterpolator.txt",
        "abundance": "abundance_RBFInterpolator_smoothingLength.txt",
        "temperature": "temperature_RBFInterpolator_smoothingLength.txt",
        "semi_analytical": "semi_analytical_smoothingLength_cf_1.txt",
        "i0_total_fits": f"{read_dir}/{galaxy_name}_i0_total.fits",
        "skirt_wavelengths": f"{read_dir}/{galaxy_name}_grid_radiationField_wavelengths.dat",
        "other_properties": f"other_properties.csv", # Other properties extracted from the vanilla FIRE snapshots. 
    }

    
    # Read gas particles
    base_dir="/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
    path_without_file_name = f'{base_dir}/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}'

    gas, lines = readfiles.read_interpolated_files_usingFilePath2(
        path = f"{path_without_file_name}/{file_names['interpolated_lines']}",
        interpolation_type = "line_emissions",
    )

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

    # Read the mass fractions for the metals 
    gas_other_properties = readfiles.read_csv_files(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift, 
            directory_name = directory_name,
            file_name = file_names['other_properties'], 
            base_fdir = base_fdir,
        )   
    # Drop the metallicity and mass column
    gas_other_properties = gas_other_properties.drop(columns=['metallicity', 'mass'], errors='ignore')

    # Merge gas particles with the other properties 
    gas = readfiles.merge_on_rounded_coords(gas, gas_other_properties, ndigits=0)

    # Add the oxygen and hydrogen mass
    gas['oxygen_mass'] = gas['O_mass_fraction'] * gas['mass'] # Msolar
    gas['hydrogen_mass'] = ( 1 - (gas['metallicity'] * constants.solar_metallicity) ) * gas['mass']

    # Add the volume to gas particles
    gas['volume'] = gas['radius']**3 * (4/3) * np.pi # pc^3

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

    return gas, star, semi_analytical, sed, lines, file_names


def metallicity_dictionary(gas: pd.DataFrame, star: pd.DataFrame, conditions: dict) -> dict:
    
    ################# Masking the gas and star particles based on the conditions
    ### R50 and Disk
    gas_half_light_visible_z500pc = galactic_prop_funcs.filter_particles(
        particles_df = gas.copy(),
        condition = conditions['r50visible_z500pc'],
    )  
    star_half_light_visible_z500pc = galactic_prop_funcs.filter_particles(
        particles_df = star.copy(),
        condition = conditions['r50visible_z500pc'],
    )


    ############ All Galaxy
    average_gas_metallicity = galactic_prop_funcs.mass_average(
        property_array = gas['metallicity'].to_numpy(),
        mass = gas['mass'].to_numpy(),
    )

    average_star_metallicity = galactic_prop_funcs.mass_average(
        property_array = star['metallicity'].to_numpy() / constants.solar_metallicity,
        mass = star['mass'].to_numpy(),
    )    


    ############ Inner Galaxy
    ## Filtering -- Find the indices of gas and star particles within 8 kpc from the center

    R_inner = 8e3 #kpc 
    indices_gas_inner_galaxy = np.where(np.sqrt(gas['x']**2 + gas['y']**2 + gas['z']**2) < R_inner)[0] 
    indices_star_inner_galaxy = np.where(np.sqrt(star['x']**2 + star['y']**2 + star['z']**2) < R_inner)[0] 
    gas_inner_galaxy = gas.iloc[indices_gas_inner_galaxy]
    star_inner_galaxy = star.iloc[indices_star_inner_galaxy]

    gas_metallicity_inner_galaxy = galactic_prop_funcs.mass_average(
        property_array = gas_inner_galaxy['metallicity'].to_numpy(),
        mass = gas_inner_galaxy['mass'].to_numpy(),
    )

    star_metallicity_inner_galaxy = galactic_prop_funcs.mass_average(
        property_array = star_inner_galaxy['metallicity'].to_numpy() / constants.solar_metallicity,
        mass = star_inner_galaxy['mass'].to_numpy()
    )

    ############ hden conditions
    # Use only the ism particles. Put a condition on density.   
    gas_hden_higher_than_1eminus2_condition = gas['hden'] > 1e-2
    gas_hden_higher_than_1 = gas['hden'] > 1e0

    filtered_gas = gas[gas_hden_higher_than_1eminus2_condition].copy()
    gas_metallicity_hden_higher_than_1eminus2 = galactic_prop_funcs.mass_average(
        property_array = filtered_gas['metallicity'].to_numpy(),
        mass = filtered_gas['mass'].to_numpy(),
    )

    filtered_gas = gas[gas_hden_higher_than_1].copy()
    gas_metallicity_hden_higher_than_1 = galactic_prop_funcs.mass_average(
        property_array = filtered_gas['metallicity'].to_numpy(),
        mass = filtered_gas['mass'].to_numpy(),
    )

    ############ Metallicity at the R50 and within Disk
    metallicity_gas_half_light_visible_z500pc = calculate_weighted_average_within_condition(
        particles = gas,
        condition = conditions['r50visible_z500pc'],
        property_name = "metallicity",
        weight_name="mass",
    )

    metallicity_star_half_light_visible_z500pc = calculate_weighted_average_within_condition(
        particles = star,
        condition = conditions['r50visible_z500pc'],
        property_name = "metallicity" ,
        weight_name="mass",
    ) / constants.solar_metallicity


    ############ Metallicity in terms of oxygen abundance
    #### All Galaxy
    twelve_plus_logOH = 12 + np.log10( ( sum(gas['oxygen_mass']) / sum(gas['hydrogen_mass']) ) / ( constants.mo_molecular_mass / constants.mh_molecular_mass )) 

    #### R50 and Disk
    # Metallicity at the half light radius and within disk in terms of oxygen abundance
    twelve_plus_logOH_half_light_visible_z500pc = 12 + np.log10(
        ( np.sum(gas_half_light_visible_z500pc['oxygen_mass']) / np.sum(gas_half_light_visible_z500pc['hydrogen_mass']) ) / ( constants.mo_molecular_mass / constants.mh_molecular_mass )
    )


    metallicities = {
        "average_gas_metallicity": average_gas_metallicity, # Zsolar
        "average_star_metallicity": average_star_metallicity, # Zsolar
        "gas_metallicity_inner_galaxy": gas_metallicity_inner_galaxy, # Zsolar
        "star_metallicity_inner_galaxy": star_metallicity_inner_galaxy, # Zsolar
        "gas_metallicity_hden_higher_than_1eminus2": gas_metallicity_hden_higher_than_1eminus2, # Zsolar
        "gas_metallicity_hden_higher_than_1": gas_metallicity_hden_higher_than_1, # Zsolar
        "metallicity_gas_half_light_visible_z500pc": metallicity_gas_half_light_visible_z500pc, # Zsolar
        "metallicity_star_half_light_visible_z500pc": metallicity_star_half_light_visible_z500pc, # Zsolar
        "twelve_plus_logOH": twelve_plus_logOH, # [Unitless]
        "twelve_plus_logOH_half_light_visible_z500pc": twelve_plus_logOH_half_light_visible_z500pc, # Unitless
    }


    return metallicities


def pressure_dictionary(gas: pd.DataFrame, star: pd.DataFrame, conditions: dict, file_names: dict) -> dict:

    """
    Calculate the disk pressures for gas and star particles within certain conditions.
    """

    R_half_light_visible = get_R50_visible(file_names=file_names) # pc

    # Calculate the pressures 
    Pgas_10kpc, Pstar_10kpc, Ptotal_10kpc = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas, 
        star_particles = star, 
        filtering_condition = conditions['10kpc'], 
        R_half_light = R_half_light_visible,
        is_plot = False, 
    )    

    Pgas_20kpc, Pstar_20kpc, Ptotal_20kpc  = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas, 
        star_particles = star, 
        filtering_condition = conditions['20kpc'], 
        R_half_light = R_half_light_visible,
        is_plot = False, 
    )        

    Pgas_half_light_visible_z500pc, Pstar_half_light_visible_z500pc, Ptotal_half_light_visible_z500pc = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas,
        star_particles = star,
        filtering_condition = conditions['r50visible_z500pc'],
        R_half_light = R_half_light_visible,
        is_plot = False,
    )

    pressures_dict = {
        "Pgas_10kpc": Pgas_10kpc,       # K / cm3
        "Pstar_10kpc": Pstar_10kpc,     # K / cm3
        "Ptotal_10kpc": Ptotal_10kpc,   # K / cm3 
        "Pgas_20kpc": Pgas_20kpc,       # K / cm3
        "Pstar_20kpc": Pstar_20kpc,     # K / cm3
        "Ptotal_20kpc": Ptotal_20kpc,   # K / cm3  
        "Pgas_half_light_visible_z500pc": Pgas_half_light_visible_z500pc,       # K / cm3
        "Pstar_half_light_visible_z500pc": Pstar_half_light_visible_z500pc,     # K / cm3
        "Ptotal_half_light_visible_z500pc": Ptotal_half_light_visible_z500pc,   # K / cm3
    }

    return pressures_dict


def temperature_dictionary(gas: pd.DataFrame, star: pd.DataFrame, conditions: dict) -> dict:

    # Calculate mass averaged temperature
    fire_temperature_mass_average_half_light_visible_z500pc = calculate_weighted_average_within_condition(
        particles = gas,
        condition = conditions['r50visible_z500pc'],
        property_name = "temperature",
        weight_name="mass",
    )

    # Calculate the CO temperature
    co_weighted_temperature_mass_average_half_light_visible_z500pc = calculate_weighted_average_within_condition(
        particles = gas,
        condition = conditions['r50visible_z500pc'],
        property_name = "Tco",
        weight_name="mass",
    )

    # Calculate the H2 temperature
    h2_weighted_temperature_mass_average_half_light_visible_z500pc = calculate_weighted_average_within_condition(
        particles = gas,
        condition = conditions['r50visible_z500pc'],
        property_name = "Th2",
        weight_name="mass",
    )

    temperature_dict = {
        "fire_temperature_mass_average_half_light_visible_z500pc": fire_temperature_mass_average_half_light_visible_z500pc, # K
        "co_weighted_temperature_mass_average_half_light_visible_z500pc": co_weighted_temperature_mass_average_half_light_visible_z500pc, # K
        "h2_weighted_temperature_mass_average_half_light_visible_z500pc": h2_weighted_temperature_mass_average_half_light_visible_z500pc, # K
    }

    return temperature_dict


def sfr_dictionary(gas: pd.DataFrame, star: pd.DataFrame, file_names: dict) -> dict:
    """
    Calculate the star formation rates for the galaxy.
    """
    # SFRs
    sfr_instantaneous = np.sum(gas["sfr"]) 
    sfr_5Myr = galactic_prop_funcs.sfr_calculator(star_df = star, within_how_many_Myr = 5)
    sfr_10Myr = galactic_prop_funcs.sfr_calculator(star_df = star, within_how_many_Myr = 10)
    sfr_100Myr = galactic_prop_funcs.sfr_calculator(star_df = star, within_how_many_Myr = 100)

    sfr_dict = {
        "sfr_instantaneous": sfr_instantaneous, # Msolar / yr
        "sfr_5Myr": sfr_5Myr,                   # Msolar / yr
        "sfr_10Myr": sfr_10Myr,                 # Msolar / yr
        "sfr_100Myr": sfr_100Myr,               # Msolar / yr
    }

    # Star formation rate surface density 
    R50_visible = get_R50_visible(file_names=file_names) # pc

    sfr_dict_surface_density = {}
    for key in sfr_dict.keys():
        sfr_dict_surface_density[f"{key}_surface_density_twoR50"] = sfr_dict[key] / (np.pi * (2 * R50_visible)**2)


    # Specific star formation rate (sSFR)
    specific_sfr = {}
    total_star_mass = np.sum(star["mass"]) # Msolar
    for key in sfr_dict.keys():
        specific_sfr[f"{key}_sSFR"] = sfr_dict[key] / total_star_mass


    # Update the sfr dictionary
    sfr_dict.update(sfr_dict_surface_density)
    sfr_dict.update(specific_sfr)

    return sfr_dict


def SFE_dictionary(gas: pd.DataFrame, sfr_dictionary: dict, conditions: dict) -> dict:

    """
    Calculate the star formation efficiencies for the galaxy.
    """
    sfe_dict = {}

    total_gas_mass = np.sum(gas["mass"]) # Msolar
    for key in sfr_dictionary.keys():
        sfe_dict[f"sfe_{key}_all_gas"] = sfr_dictionary[key] / total_gas_mass

    total_gas_mass_within_R50_and_disk = galactic_prop_funcs.filter_particles(
        particles_df = gas.copy(),
        condition = conditions['r50visible_z500pc'],
    )['mass'].sum() # Msolar
    for key in sfr_dictionary.keys():
        sfe_dict[f"sfe_{key}_R50_and_disk"] = sfr_dictionary[key] / total_gas_mass_within_R50_and_disk

    return sfe_dict


def mass_dictionary(gas: pd.DataFrame, star: pd.DataFrame, semi_analytical: pd.DataFrame, galaxy_name: str, galaxy_type: str, redshift: str, conditions: dict) -> dict:

    ## gas and star mass
    total_gas_mass = np.sum(gas["mass"])
    total_star_mass = np.sum(star["mass"])

    ## halo mass
    try:
        Mhalo = galactic_prop_funcs.halo_mass_calculator(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift
        )
    except Exception as e:
        print(f"Exception occurred while calculating halo mass for {galaxy_name} ({galaxy_type}):\n{e}")
        Mhalo = galactic_prop_funcs.get_halo_mass_from_previously_calculated_file(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift,
        )    

    ## Molecular gas mass
    h2_mass_semi_anaytical = np.sum(semi_analytical["Mh2"]) # if h2 mass is NaN probably it is zero already.
    h2_mass_cloudy = np.sum(gas['fh2'] * gas['mass'])

    # CO mass
    gas['mass_co'] = gas['fCO'] * gas['mass']
    co_mass_cloudy = np.sum(gas['mass_co'].to_numpy())

    ##### Gas fractions
    # All galaxy 
    gas_fraction_all_galaxy = total_gas_mass / (total_gas_mass + total_star_mass)
    # Gas fraction within R50 and disk
    total_gas_mass_within_R50_and_disk = galactic_prop_funcs.filter_particles(
        particles_df = gas.copy(),
        condition = conditions['r50visible_z500pc'],
    )['mass'].sum() # Msolar
    total_star_mass_within_R50_and_disk = galactic_prop_funcs.filter_particles(
        particles_df = star.copy(),
        condition = conditions['r50visible_z500pc'],
    )['mass'].sum() # Msolar
    gas_fraction_R50_and_disk = total_gas_mass_within_R50_and_disk / (total_gas_mass_within_R50_and_disk + total_star_mass_within_R50_and_disk)


    ### Stellar mass fractions 
    stellar_mass_fraction_all_galaxy = total_star_mass / (total_star_mass + total_gas_mass)
    stellar_mass_fraction_within_R50_and_disk = total_star_mass_within_R50_and_disk / (total_star_mass_within_R50_and_disk + total_gas_mass_within_R50_and_disk)

    mass_dict = {
        "gas_mass": total_gas_mass, # Msolar
        "star_mass": total_star_mass, # Msolar
        "total_gas_mass_within_R50_and_disk": total_gas_mass_within_R50_and_disk, # Msolar
        "total_star_mass_within_R50_and_disk": total_star_mass_within_R50_and_disk, # Msolar
        "gas_fraction_all_galaxy": gas_fraction_all_galaxy, # Unitless
        "gas_fraction_R50_and_disk": gas_fraction_R50_and_disk, # Unitless
        "stellar_mass_fraction_all_galaxy": stellar_mass_fraction_all_galaxy, # Unitless
        "stellar_mass_fraction_within_R50_and_disk": stellar_mass_fraction_within_R50_and_disk, # Unitless
        "h2_mass_semi_analytical": h2_mass_semi_anaytical, # Msolar
        "h2_mass_cloudy": h2_mass_cloudy, # Msolar
        "co_mass_cloudy": co_mass_cloudy, # Msolar
        "halo_mass": Mhalo, # Msolar
    }

    return mass_dict


def skirt_related_bolometric_luminosity_dictinary(sed: pd.DataFrame, file_names: dict) -> dict:

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
        },
        "habing": {
            "min_wavelength": 0.0912, # 13.6 eV
            "max_wavelength": 0.1240, # 10 eV
        },
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

    L_skirt_habing = functions.calculate_luminosity_from_sed(
        sed = sed,
        min_wavelength = bands['habing']['min_wavelength'],
        max_wavelenght = bands['habing']['max_wavelength'],
        distance_in_Mpc = distance
    ) # erg / s

    skirt_dict = {
        "L_fir": L_fir, # erg / s
        "L_skirt_full_wavelength_range": L_skirt_full_wavelength_range, # erg / s,
        "L_skirt_habing": L_skirt_habing, # erg / s
    }

    return skirt_dict


def calculate_averaged_radiation_field(gas: pd.DataFrame, conditions: dict) -> dict:
    """
    Calculating the mass and volume averaged radiation field for all galaxy and for certain conditions.
    """

    ### For all galaxy 
    mass_average_radiation_field = calculate_weighted_average_within_condition(
        particles = gas.copy(),
        condition = conditions['all_galaxy'],
        property_name = "isrf",
        weight_name="mass",
    ) # G0

    volume_weighted_average_radiation_field = calculate_weighted_average_within_condition(
        particles = gas.copy(),
        condition = conditions['all_galaxy'],
        property_name = "isrf",
        weight_name="volume", 
    ) # G0

    ### For R50 and Disk
    mass_average_radiation_field_R50_and_disk = calculate_weighted_average_within_condition(
        particles = gas.copy(),
        condition = conditions['r50visible_z500pc'],
        property_name = "isrf",
        weight_name="mass",
    ) # G0 

    volume_weighted_average_radiation_field_R50_and_disk = calculate_weighted_average_within_condition(
        particles = gas.copy(),
        condition = conditions['r50visible_z500pc'],
        property_name = "isrf",
        weight_name="volume",
    ) # G0

    average_radiation_field_dict = {
        "mass_average_radiation_field": mass_average_radiation_field, # G0
        "volume_weighted_average_radiation_field": volume_weighted_average_radiation_field, # G0
        "mass_average_radiation_field_R50_and_disk": mass_average_radiation_field_R50_and_disk, # G0
        "volume_weighted_average_radiation_field_R50_and_disk": volume_weighted_average_radiation_field_R50_and_disk, # G0
    }

    return average_radiation_field_dict


def get_R50_visible(file_names: dict) -> float:

    ############ Half light radius calculations
    # Determine the bands
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

    R_half_light_visible = galactic_prop_funcs.half_light_radius(
        band = half_light_wavelength_bands['visible'], 
        i0_total_fits_file_path = file_names['i0_total_fits'],
        skirt_wavelengths_file_path = file_names['skirt_wavelengths'],
        plot = False
    ) # pc

    return R_half_light_visible


def get_galactic_properties_for_a_galaxy(galaxy_name:str, galaxy_type:str, redshift:str, directory_name:str, debugging:bool=False) -> dict:
    
    print(f"--------------  {galaxy_name}  --------------")
    
    start = time()

    gas, star, semi_analytical, sed, lines, file_names = read_all_files(
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift, 
        directory_name = directory_name
    )


    # Get the half light radius
    R_half_light_visible = get_R50_visible(file_names) # pc
    R_half_mass_stellar_mass = galactic_prop_funcs.half_mass_radius(particles_df = star) # pc


    ############ Define the conditions for filtering the particles 
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
        }, 
        "all_galaxy": {
            "r_max": np.inf, # pc
            "z_max": np.inf, # pc
        }
    }      

    ############ Calculate galactic properties 
    ### Line Luminosities
    # Krumholz - semi_analytical 
    semi_analytical_Lco = np.sum(semi_analytical['L_co'])
    # Cloudy
    total_line_luminosities = {}
    for line in lines: 
        total_line_luminosities[line] = sum(gas[line])
    
    # Calculate the line luminosities at the half light radius 
    total_line_luminosities_R50 = {}
    gas_half_light_visible_z500pc = galactic_prop_funcs.filter_particles(
        particles_df = gas.copy(),
        condition = conditions['r50visible_z500pc'],
    )
    for line in lines:
        total_line_luminosities_R50[f"{line}_half_light_visible_z500pc"] = sum(gas_half_light_visible_z500pc[line])

    if debugging:
        print("Line luminosities calculated successfully.")

    ### NaN indices
    nan_indices = np.where(np.isnan(gas["L_co_10"]))[0]    


    ### Calculate the mass averaged properties within a certain condition
    # Calculate mass averaged hden
    hden_mass_average_half_light_visible_z500pc = calculate_weighted_average_within_condition(
        particles = gas,
        condition = conditions['r50visible_z500pc'],
        property_name = "hden",
        weight_name="mass",
    )


    # Calculate the column density for certain gas particles determined by the condition
    column_density_half_light_visible_z500pc = calculate_column_density_within_condition(
        particles = gas,
        condition = conditions['r50visible_z500pc']
    )

    if debugging:
        print("Mass averaged hden and column density calculated successfully.")


    ### Calculate the mass properties
    mass_dict = mass_dictionary(
        gas = gas, 
        star = star, 
        semi_analytical = semi_analytical, 
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift,
        conditions = conditions,
    )

    if debugging:
        print("Mass properties calculated successfully.")

    ### Calculate sfr
    sfr_dict = sfr_dictionary(gas=gas.copy(), star=star.copy(), file_names=file_names)    

    if debugging:
        print("Star formation rates calculated successfully.")

    ### Calculate SFE
    sfe_dict = SFE_dictionary(gas=gas.copy(), sfr_dictionary=sfr_dict, conditions=conditions)

    if debugging:
        print("Star formation efficiencies calculated successfully.")

    ### Calculate metallicities 
    metallicities_dict = metallicity_dictionary(gas=gas.copy(), star=star.copy(), conditions=conditions)

    if debugging:
        print("Metallicities calculated successfully.")

    ### Calculate Pressures 
    pressures_dict = pressure_dictionary(
        gas=gas.copy(), 
        star=star.copy(), 
        conditions=conditions, 
        file_names=file_names
    )

    if debugging:
        print("Pressures calculated successfully.")

    ### Calculate temperatures
    temperature_dict = temperature_dictionary(gas=gas.copy(), star=star.copy(), conditions=conditions)

    if debugging:
        print("Temperatures calculated successfully.")

    ### Get the bolometric luminosity from the skirt sed file
    skirt_bolometric_luminosity_dict = skirt_related_bolometric_luminosity_dictinary(
        sed = sed.copy(),
        file_names = file_names
    )

    if debugging:
        print("Skirt bolometric luminosity calculated successfully.")

    ### Calculate the averaged radiation field
    average_radiation_field_dict = calculate_averaged_radiation_field(
        gas = gas.copy(),
        conditions = conditions,
    )

    if debugging:
        print("Averaged radiation field calculated successfully.")


    # Create a dictionary 
    galactic_props1 = {
        "name": f"{galaxy_name}",
        "galaxy_type": f"{galaxy_type}",
        "redshift": f"{redshift}",
        "number_of_NaN_indices": len(nan_indices),
        "R_half_light_visible": R_half_light_visible, # pc
        "R_half_mass_stellar_mass": R_half_mass_stellar_mass, # pc
        "hden_mass_average_half_light_visible_z500pc": hden_mass_average_half_light_visible_z500pc, # cm^-3
        "column_density_half_light_visible_z500pc": column_density_half_light_visible_z500pc, # Msolar / pc2,
        "L_co_10_semi_analytical": semi_analytical_Lco, # K km s^-1 pc^2 
    }     

    # Merge two dictionaries
    galactic_properties = {
        **galactic_props1, 
        **mass_dict,
        **sfr_dict, 
        **sfe_dict,
        **metallicities_dict, 
        **pressures_dict, 
        **temperature_dict,
        **skirt_bolometric_luminosity_dict,
        **average_radiation_field_dict,
        **total_line_luminosities, 
        **total_line_luminosities_R50, 
    }

    stop = time()
    
    print(f"For {galaxy_name}, it took {round((stop-start)/60, 3)} minutes")

    return galactic_properties


def calculate_weighted_average_within_condition(particles: pd.DataFrame, condition: dict, property_name: str, weight_name: str,) -> float:
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
    max_workers = int(sys.argv[2])

    main(redshift=redshift, max_workers=max_workers)