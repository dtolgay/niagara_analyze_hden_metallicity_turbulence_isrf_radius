import sys
sys.path.append("/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs")

import numpy as np
import pandas as pd 
from time import time

from tools import functions, constants  # type: ignore 
from tools import functions_readfiles as readfiles # type: ignore
from tools.functions_create_galactic_properties import sfr_calculator, halo_mass, mass_average # type: ignore

import astropy.units as units
import galactic_props_functions as galactic_prop_funcs

def main(redshift):

    # base file name to start searching files
    # base_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
    directory_name = "voronoi_1e6"
    galaxies_list = []

    # ####################################### For particle split 
    # galaxy_names = [
    #     "m12i_r880_md"
    # ]

    # for galaxy_name in galaxy_names:
    #     try: 
    #         galaxies_list.append(get_galactic_properties_for_a_galaxy(
    #             galaxy_name=galaxy_name, 
    #             galaxy_type="particle_split", 
    #             redshift=redshift, 
    #             directory_name=directory_name
    #         ))
    #     except Exception as e:
    #         print(f"Exception occured: \n{e}")    

    # ####################################### For firebox 

    # for i in range(int(1e3)):
    # for i in range(51):
    for i in [0]:
        try:
            galaxies_list.append(get_galactic_properties_for_a_galaxy(
                galaxy_name=f"gal{i}", 
                galaxy_type="firebox", 
                redshift=redshift, 
                directory_name=directory_name
            ))
            
        except Exception as e:
            print(f"Exception occured: \n{e}") 

    # ####################################### For zoom_in
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
    #     try:
    #         galaxies_list.append(get_galactic_properties_for_a_galaxy(
    #             galaxy_name=galaxy_name, 
    #             galaxy_type="zoom_in", 
    #             redshift=redshift, 
    #             directory_name=directory_name
    #         ))
            
    #     except Exception as e:
    #         print(f"Exception occured: \n{e}")   

    ####################################### Turn list to a pandas dataframe
    galaxies = pd.DataFrame(galaxies_list)


    ####################################### Write to a file

    save_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/data"
    file_path = f"{save_dir}/galactic_properties_smoothingLength_hybridInterpolator_z{int(float(redshift))}_usingIvalues.csv"
    # file_path = f"{save_dir}/delete.txt"

    # Append the DataFrame to the file
    galaxies.to_csv(file_path, mode='w', sep=',', index=False, float_format='%.5e')

    print(f"File saved to {file_path}")


    return 0

#### Functions defined below 


def get_galactic_properties_for_a_galaxy(galaxy_name:str, galaxy_type:str, redshift:str, directory_name:str):
    
    print(f"--------------  {galaxy_name}  --------------")
    
    start = time()


    save_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/data"

    base_fdir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
    read_dir = f"{base_fdir}/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"

    file_names = {
        "interpolated_lines": "L_line_smoothingLength_hybridInterpolator_flux2Luminosity.txt",
        "otherProperties": "otherProperties_smoothingLength_hybridInterpolator.txt",
        "semi_analytical": "semi_analytical_smoothingLength_cf_1.txt",
        "i0_total_fits": f"{read_dir}/{galaxy_name}_i0_total.fits",
        "skirt_wavelengths": f"{read_dir}/{galaxy_name}_grid_radiationField_wavelengths.dat",
        "write_file_name": f"{save_dir}/galactic_properties_smoothingLength_hybridInterpolator_z{int(float(redshift))}_usingIvalues.csv"
    }

    
    # Read gas particles
    gas, lines = readfiles.read_interpolated_Lline(
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift, 
        directory_name = directory_name, 
        file_name = file_names['interpolated_lines']
    )

    gas_other_properties = readfiles.read_otherProperties(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift, 
            directory_name = directory_name, 
            file_name = file_names['otherProperties'],
        )   

    # Merge gas particles 
    merged_columns = ['index', 'fh2', 'fCO']
    gas = gas.merge(gas_other_properties[merged_columns], how='inner', on=['index'], validate='one_to_one')


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
    sed = readfiles.read_skirt_sed_file(galaxy_name, galaxy_type, redshift, directory_name = "voronoi_1e6", inclination = "0")
    
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
    co_mass_cloudy = np.sum(gas['fCO'] * gas['mass'])

    ## halo mass
    Mhalo = halo_mass(
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift
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


    # Calculate disk pressure 
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

    conditions = {
        "20kpc": {
            "r_max" : 20e3, # pc
            "z_max" : 2e3, # pc 
        },
        "10kpc": {
            "r_max": 10e3, # pc
            "z_max": 2e3, # pc 
        },
    }



    Pgas_10kpc, Pstar_10kpc, Ptotal_10kpc = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas, 
        star_particles = star, 
        filtering_condition = conditions['10kpc'], 
        band = half_light_wavelength_bands['visible'],
        file_names = file_names,
        is_plot = False, 
    )    

    Pgas_20kpc, Pstar_20kpc, Ptotal_20kpc = galactic_prop_funcs.disk_pressure_calculator(
        gas_particles = gas, 
        star_particles = star, 
        filtering_condition = conditions['20kpc'], 
        band = half_light_wavelength_bands['visible'],
        file_names = file_names,
        is_plot = False, 
    )        


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
        "Pstar_10kpc": Pgas_10kpc, # K / cm3
        "Ptotal_10kpc": Ptotal_10kpc, # K / cm3 
        "Pgas_20kpc": Pgas_20kpc, # K / cm3
        "Pstar_20kpc": Pgas_20kpc, # K / cm3
        "Ptotal_20kpc": Ptotal_20kpc, # K / cm3         
    }     

    # Merge two dictionaries
    galactic_properties = {**other_properties, **total_line_luminosities}

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

if __name__ == "__main__":
    redshift = sys.argv[1]
    main(redshift)