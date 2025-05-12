import sys
sys.path.append("/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs")

import numpy as np
import pandas as pd 
from time import time

from tools import readsnap, functions, constants  # type: ignore 
from tools import functions_readfiles as readfiles # type: ignore

from concurrent.futures import ThreadPoolExecutor, as_completed # parallelize


def main():

    # base file name to start searching files
    # base_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius"
    directory_name = "voronoi_1e6"
    redshift = sys.argv[1]
    # redshift = "3.0"
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
    for i in range(50):
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
    header = f"""
    #Column 0: name 
    #Column 1: galaxy_type
    #Column 2: redshift 
    #Column 3: sfr_instantaneous [Msolar/year]
    #Column 4: sfr_5Myr [Msolar/year]
    #Column 5: sfr_10Myr [Msolar/year]
    #Column 6: sfr_100Myr [Msolar/year]
    #Column 7: gas_mass [Msolar]
    #Column 8: star_mass [Msolar]
    #Column 9: gas_average_metallicity [Zsolar]
    #Column 10: star_average_metallicity [Zsolar]
    #Column 11: alpha_co [Msolar / (K km s^-1 pc^2)]
    #Column 12: halo_mass [Msolar]
    #Column 13: number_of_NaN_indices [1]
    #Column 14: L_ly_alpha [erg s^-1]
    #Column 15: L_h_alpha [erg s^-1]
    #Column 16: L_h_beta [erg s^-1]
    #Column 17: L_co_10 [K km s^-1 pc^2] 
    #Column 18: L_co_21 [K km s^-1 pc^2] 
    #Column 19: L_co_32 [K km s^-1 pc^2] 
    #Column 20: L_co_43 [K km s^-1 pc^2] 
    #Column 21: L_co_54 [K km s^-1 pc^2] 
    #Column 22: L_co_65 [K km s^-1 pc^2] 
    #Column 23: L_co_76 [K km s^-1 pc^2] 
    #Column 24: L_co_87 [K km s^-1 pc^2] 
    #Column 25: L_13co [K km s^-1 pc^2] 
    #Column 26: L_c2 [erg s^-1]
    #Column 27: L_o3_88 [erg s^-1]
    #Column 28: L_o3_5006 [erg s^-1]
    #Column 29: L_o3_4958 [erg s^-1]     
    """

    save_dir = "/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/python_files/analyze_hden_metallicity_turbulence_isrf_radius/data"
    file_path = f"{save_dir}/galactic_properties_smoothingLength_nearestNDInterpolator_z{int(float(redshift))}_usingIvalues.csv"
    # file_path = f"{save_dir}/delete.txt"

    # # Write the header to the file
    # with open(file_path, 'w') as f:
    #     f.write(header.strip() + '\n')

    # Append the DataFrame to the file
    galaxies.to_csv(file_path, mode='w', sep=',', index=False, float_format='%.5e')

    print(f"File saved to {file_path}")


    return 0

#### Functions defined below 


def get_galactic_properties_for_a_galaxy(galaxy_name:str, galaxy_type:str, redshift:str, directory_name:str):
    
    print(f"--------------  {galaxy_name}  --------------")
    
    start = time()

    
    # Read gas particles
    
    file_name_lines = "L_line_smoothingLength_nearestNDInterpolator_flux2Luminosity.txt"
    gas, lines = readfiles.read_interpolated_Lline(
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift, 
        directory_name = directory_name, 
        file_name = file_name_lines
    )

    other_properties_file_name = "otherProperties_smoothingLength_nearestNDInterpolator.txt"
    gas_other_properties = readfiles.read_otherProperties(
            galaxy_name = galaxy_name, 
            galaxy_type = galaxy_type, 
            redshift = redshift, 
            directory_name = directory_name, 
            file_name = other_properties_file_name,
        )   

    # Merge gas particles 
    ## TODO: In the future change x,y,z with the index and merge according to index.
    merged_columns = ['index', 'fh2', 'fCO']
    gas = gas.merge(gas_other_properties[merged_columns], how='inner', on=['index'], validate='one_to_one')
        

    # Read star particles
    star = readfiles.read_comprehensive_star_particles(galaxy_name, galaxy_type, redshift, directory_name)

    
    # Read semi_analytical_average_sobolev_smoothingLength.txt
    semi_analytical_file_name = "semi_analytical_smoothingLength_cf_1.txt"
    semi_analytical = readfiles.read_semianalytical_file_2(
        galaxy_name = galaxy_name, 
        galaxy_type = galaxy_type, 
        redshift = redshift, 
        directory_name = directory_name, 
        file_name = semi_analytical_file_name,
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
    main()