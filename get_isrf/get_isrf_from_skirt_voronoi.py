# Importing modules
import sys 
sys.path.append("/scratch/m/murray/dtolgay/")


import numpy as np
from scipy import integrate
from scipy.spatial import KDTree
import time
import pandas as pd
import os

from tools import constants


def main(galaxy_name:str, which_FIRE:str, redshift:float):

    time_start_of_the_code = time.time()

    # directory_name = "voronoi_1e6_2"
    # directory_name = "zstar_doubled_voronoi_1e6"
    # directory_name = "mstar_doubled_voronoi_1e6"
    directory_name = "voronoi_1e5"
    # directory_name = "seperated_firebox_galaxies"
    # directory_name = "40kpc_voronoi_1e5"
    # directory_name = "voronoi_1e6"
    # directory_name = "voronoi_3e6"
    # directory_name = "voronoi_3e5"
    # directory_name = "voronoi_1e4"


    ############### Reading files
    # Define the name of the galaxy and gas particles files.

    # directory_path = f"/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius/{which_FIRE}/z{redshift}/{galaxy_name}/trial1"   
    if redshift in ["2.0", "3.0"]: 
        directory_path = f"/gpfs/fs0/scratch/r/rbond/dongwooc/scratch_rwa/doga/runs_hden_radius/{which_FIRE}/z{redshift}/{galaxy_name}/{directory_name}"    
    else: 
        sys.exit("Code only works for z=2.0. Implement other redshifts.")

    skirt_galaxy_path_wout_extension = f"{directory_path}/{galaxy_name}"

    print(f"directory_path: {directory_path}")


    ### Check if file exits. If it exists do not runs the code: 
    if os.path.isfile(f"{directory_path}/isrf_gas.txt"):
        print("File exits. Returning nothing!")
        return 0
    else:
        print(f"{directory_path}/isrf_gas.txt doesn't exitst. Continuing...")
    


    #### Get the wavelengths 
    wavelengths_filename = skirt_galaxy_path_wout_extension +"_grid_radiationField_wavelengths.dat"
    print(f"\nwavelengths_filename: {wavelengths_filename}")
    wavelength_indices, wavelengths = find_indices_of_wavelengths(wavelengths_filename, 
                                                                  min_wavelength = eV_2_micron(13.6),
                                                                  max_wavelength = eV_2_micron(6), 
                                                                  debug = True) 


    #### Read cell properties 
    start = time.time()

    cell_properties_file_name = skirt_galaxy_path_wout_extension + "_properties_spatialCell_cellprops.dat"
    print(f"\nStarting to read {cell_properties_file_name} file.")

    cell_content = np.loadtxt(cell_properties_file_name, delimiter = ' ')

    print(f"Read {cell_properties_file_name} file.")

    end = time.time()
    time_took_reading_cell_properties_file = end - start
    print(f"Took {round(time_took_reading_cell_properties_file, 2)} seconds")
    ####

    #### Import FIRE gas data
    gas_particles_file_name = "comprehensive_gas.txt"      
    gas_particles_path = f"{directory_path}/{gas_particles_file_name}"
    print(f"\nStarting to read {gas_particles_path} file.")

    gas_particles = np.loadtxt(gas_particles_path)
    print(f"Read {gas_particles_path} file.")  

    ####

    #### Read the radiationField_J.dat file 
    start = time.time()

    wavelength_indices_in_radiationField_J_file = wavelength_indices + 1 
    indices = np.append(0, wavelength_indices_in_radiationField_J_file)  # Add zero to get the cell index 

    radiationField_J_filename = skirt_galaxy_path_wout_extension + "_grid_radiationField_J.dat"
    print(f"\nI am reading the {radiationField_J_filename} file")
    radiationField_content = np.loadtxt(radiationField_J_filename, delimiter=' ', usecols = indices)  # Only the habing radiation field indices are imported

    end = time.time()
    time_took_reading_radiationField_J_dat = end - start
    print(f"{radiationField_J_filename} read! Took {round(time_took_reading_radiationField_J_dat / 60, 2) } minutes")
    #### 


    ############### Linking particles and cells
    # Search closest cell and gas particles
    print("\nStarted to link cells and gas particles")

    # Start the timer
    start = time.time()

    # Generate kd-tree for cells  
    cell_tree = KDTree(cell_content[:,1:4])

    # Initialize empty lists to store the indices of the gas and cells that are closest to each other
    gas_cell_coverage = []

    # To debugging store if the the particle is inside the cell. If it is bigger than 1 then it means gas particle is actually 
    # outside of the cell, and if it is smaller than one, it means gas particle is inside the closest cell. 
    distance_over_r = []

    counter = 0 # This is to track the index of the gas particles
    for gas_particle in gas_particles:
        # ask the query and record the index 
        k = 1  # Number of closest cells
        distances, cell_indices = cell_tree.query(gas_particle[0:3], k = k)

        # Matching the gas particle with the closest voronoi tesellation center.
        gas_cell_coverage.append([counter, cell_indices])

        if (counter % int(1e5) == 0):
            print(f"I am in the {counter}.th loop. Left {len(gas_particles) - counter}")

        # Incerement the counter by one 
        counter += 1

    # End the timer
    end = time.time()
    time_search_closest_gas_and_cell_particles = end - start
    print(f"It took {round(time_search_closest_gas_and_cell_particles, 2)} seconds to search closest cell and gas particles")
    
    gas_cell_coverage = np.array(gas_cell_coverage) 


    print("Voronoi tesellation method is used. The gas particles are matched with the closest voronoi tesellation center.")


    ############### Calculating the isrf
    # Filter the cell's such that I will calculate the ISRF for the cells that are linked to the gas particles to save time and not bother with the other cells. 
    # There are some repeating indices for cells. I need to get rid of the repetating values.
    cells_that_are_linked_with_gas_particles = np.unique(gas_cell_coverage[:,1]) # Automatically sorted.

    # Filter the cells 
    radiationField_content = radiationField_content[cells_that_are_linked_with_gas_particles]
    cell_content = cell_content[cells_that_are_linked_with_gas_particles]
    print("\nradiationField_content and cell_content arrays are filtered. Only the cells linked with gas particles are left")   


    # Calculate isrf
    print(f"\nISRF calculation started.")
    # Calculate the ISRF 
    start = time.time()

    isrf_for_cells = []
    for row in radiationField_content:
        G_over_G0 = return_isrf_in_habing_units(data = row, wavelengths = wavelengths)
    #     isrf_data.append({"cell_index": int(row[0]),
    #                       "G_in_habing_units": G_over_G0})
        isrf_for_cells.append([int(row[0]), G_over_G0])
        # 0.th column is the cell index 
        # 1.th column is the G_over_G0 value. ISRF in Habing Units


    isrf_for_cells = np.array(isrf_for_cells)  # Converting it to numpy array
    end = time.time()
    time_isrf_data_calculator = end - start
    print(f"Time it tooks to calculate is isrf is: {round(time_isrf_data_calculator, 2)} seconds")    


    ############### Assign isrf to gas particles

    print("\nISRF started to assing to gas particles")
    start = time.time()

    # TODO: Delete below
    # gas_index_isrf = []
    # for gas_index, cell_index in gas_cell_coverage:
        
    #     index_isrf_for_cells = np.where(isrf_for_cells[:,0] == cell_index)[0]
    #     isrf_corresponding_to_that_index = isrf_for_cells[index_isrf_for_cells][0][1] # There is no problem with [0][1] 
    #     gas_index_isrf.append([gas_index, isrf_corresponding_to_that_index])

    # Convert arrays to pandas dataframes
    df_gas_cell_coverage = pd.DataFrame(gas_cell_coverage, columns=['gas_index', 'cell_index'])
    df_isrf_for_cells = pd.DataFrame(isrf_for_cells, columns=['cell_index', 'isrf'])

    # Start the timer
    start = time.time()

    # Perform a merge operation
    merged_df = df_gas_cell_coverage.merge(df_isrf_for_cells, on='cell_index')
    merged_df = merged_df.sort_values(by='gas_index')

    # Extract the required columns
    gas_index_isrf = merged_df[['gas_index', 'isrf']].values.tolist()        


    end = time.time()
    print(f"Time it tooks to assign isrf to gas particles is: {round(end-start, 2)} seconds")    


    ############### Writing to a file 
    write_file_name = "isrf_gas.txt"
    write_file_path = f"{directory_path}/{write_file_name}"

    header = f"""isrf for {galaxy_name} galaxy
    Column 0: gas_index (1)
    Column 1: isrf (G0)
    """

    print(f"Started to write to {write_file_path}")
    np.savetxt(fname=write_file_path, X=gas_index_isrf, fmt='%.8e', header=header)
    print(f"File written to {write_file_path}")



    return 0 


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# Functions 

def eV_2_micron(energy):
    
    '''Conversion from electron volts to the wavelength
    
    #######
    
    Parameters: 
    
    energy: in eV
    
    Returns: 
    
    wavelength: in micron
    
    '''
    
    wavelength = constants.h * constants.c / energy * 1e6 # in microns
    
    return wavelength 

#############################################################################################################################

def return_isrf_in_habing_units(data:np.ndarray, wavelengths:np.ndarray, debug:bool=False) -> float:
    
    '''This function is created to integrate the filtered data
    
    #####
    
    Parameters:
    ------------
    data: list 
        This is the data corresponding to a single line in the radiationField_J.dat file.
        
    wavelengths: list
        This is the wavelengths that I am going to use to integrate the data.
        
    indices: np.ndarray 
        indices of the wavelengths corresponding to the Habing radiation field. The indices of the radiationField_J.dat and 
        radiationField_wavelenghts.dat file not the same. First column of the radiationField_J.dat file is the cell index number
        thus, mean intensity corresponding to habing radition field wavelengths are the columns corresponding to the indices+1.
    
    #####
    
    Return:
    ------------
    G_over_G0: float
        This is the interstellar radiation field corresponding in Habing units. Eqn 12.6, Draine. 
    
    
    '''
    
#     print("I am in the function return_isrf_in_habing_units")
    
    # Only get the flux values corresponding to that wavelengths 
    # Oth column belongs to the cell identification number  
    # Checked! Works correctly
    mean_intensity = data[1:]
    
    # Integration: https://en.wikipedia.org/wiki/Simpson%27s_rule
    integrated_value = integrate.simpson(y=mean_intensity, x=wavelengths)  # W/m^2

    if debug: 
        # Do other integration techniques as well to compare 
        # Trapezoidal rule
        integrated_value_trapz = integrate.trapz(y=mean_intensity, x=wavelengths)  # W/m^2
        print(f"integrated_value_trapz: {integrated_value_trapz}")
        print(f"integrated_value: {integrated_value}")
        print(f"integrated_value - integrated_value_trapz: {integrated_value - integrated_value_trapz}")
        print(f"integrated_value / integrated_value_trapz: {integrated_value / integrated_value_trapz}") 

    # Calculating the total radiation density (J/m^3). Eqn 1.9 in Rybicki and Lightman
    u = 4 * np.pi * integrated_value / constants.c 
    
    # Find G parameter
    G = u/constants.u_HAB
    
    # Normalizing the total radiation density with habing radiation field to find the isrf in Habing units. Eqn 12.6, Draine
    G_over_G0 = G/constants.G0
    
    return G_over_G0

#############################################################################################################################
def find_indices_of_wavelengths(radiationField_wavelengths_filename: str, 
                                min_wavelength: float, 
                                max_wavelength: float,
                                debug: bool = False) -> np.ndarray:
    
    """
    Returns indices of wavelengths.
    
    Parameters:
    -----------
    radiationField_wavelengths_filename : str
        Path to the file containing wavelength data.
        
    min_wavelength: float
        Smallest wavelength of interest.
        micron

    max_wavelength: float
        Largest wavelength of interest.
        micron

    Returns:
    --------
    indices:
        Indices corresponding to wavelengths.
        Note: The first column of the `radiationField_J.dat` file represents the cell index number. 
        Thus, when using these indices on that file, remember to increment each index by 1.
        
    wavelengths_filtered: 
        Center wavelengths in each measurement. 

    Example:
    --------
    >>> indices, wavelengths = find_indices_of_wavelengths_corresponding_to_the_habing_radiation_field("radiationField_wavelengths.dat")
    >>> print(indices)
    [index1, index2, ...]
    
    """
    
    
    print("I am in the function find_indices_of_wavelengths")
        
    print(f"{radiationField_wavelengths_filename} is started to read.")
    wavelengths_content = np.genfromtxt(radiationField_wavelengths_filename, delimiter = ' ', skip_header = 1)
    print(f"{radiationField_wavelengths_filename} was read.")
    
    wavelengths = wavelengths_content[:,0]
    
    indices = np.where((wavelengths < max_wavelength) & (wavelengths > min_wavelength))[0]
    wavelengths_filtered = wavelengths[indices] 

    if debug:
        print(f"Maximum wavelength: {max(wavelengths_filtered)} micron")
        print(f"Minimum wavelength: {min(wavelengths_filtered)} micron")
    
    return indices, wavelengths_filtered

##############################################################################################################################
##############################################################################################################################

if __name__ == "__main__":
    
    # ################### firebox
    # # galaxy_names = ["gal0", "gal5", "gal10", "gal15"]
    # galaxy_names = []
    # for i in range (29):
    # # for i in [11]: # TODO: Delete Testing purposses
    #     galaxy_name = f"gal{i}"
    #     galaxy_names.append(galaxy_name)

    # which_fire = ["firebox"] * len(galaxy_names)
    # redshift = ["0.0"] * len(galaxy_names)

    
    # with concurrent.futures.ThreadPoolExecutor(max_workers=15) as executor:
    #     results = list(executor.map(main, galaxy_names, which_fire, redshift))  
    # ################### 

    # galaxy_name = "m12i_res7100_md"
    # which_fire = "zoom_in"
    # redshift = "0.0"

    galaxy_name = sys.argv[1]
    which_fire = sys.argv[2]
    redshift = sys.argv[3]

    # Call main with the arguments
    main(galaxy_name, which_fire, redshift)    