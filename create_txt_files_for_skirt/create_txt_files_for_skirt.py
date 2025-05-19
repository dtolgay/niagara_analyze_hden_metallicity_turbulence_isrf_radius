#!/usr/bin/env python
# coding: utf-8

# Importing necessary librarires
import os
import sys
sys.path.append("/scratch/m/murray/dtolgay/")

import h5py
import numpy as np
from numpy import linalg as LA

import matplotlib
import matplotlib.pyplot as plt
# plt.style.use('seaborn-poster')

from math import pi

from tools import functions, readsnap, readsnap_FIREBox, readsnap_tripleLatte, constants  # tools directory is in the appended system directory
from tools.filter_rotate_galaxy import filter_rotate_galaxy

import random 
import sys # To run the batch script

from time import time

from tools.meshoid import Meshoid



def main(gas_particles, star_particles, header_info, galaxy_type, galaxy_name, snapshot_number=None):

    print(f"---------------------------------------------------- {galaxy_name} ----------------------------------------------------")

    write_base_fdir = "/gpfs/fs0/scratch/r/rbond/dongwooc/scratch_rwa/doga/runs_hden_radius"

    print(gas_particles.keys())

    hubble      = header_info['hubble']
    redshift    = header_info['redshift']   
    header_time = header_info['time']   

    gas_particles_df, star_particles_df = filter_rotate_galaxy(
        galaxy_name,    
        galaxy_type,    
        header_info, 
        gas_particles, 
        star_particles,   
        snapshot_number, # Not important if the galaxy is firebox    
    )    

    #################################### Plot
    R_max = 20 # kpc

    print("Plotting galaxy")
    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))  # Adjust figsize as needed

    # First subplot: x-y axis plot
    axs[0].hist2d(
        x=gas_particles_df["x"],
        y=gas_particles_df["y"],
        bins=500,
        norm=matplotlib.colors.LogNorm(),
        range=[[-R_max, R_max], [-R_max, R_max]]
    )
    axs[0].set_title('x-y axis plot')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('y')

    # Second subplot: z-y axis plot
    axs[1].hist2d(
        x=gas_particles_df["y"],  # Use the 'y' column for the x-axis
        y=gas_particles_df["z"],  # Keep the 'z' column for the y-axis
        bins=500,
        norm=matplotlib.colors.LogNorm(),
        range=[[-R_max, R_max], [-R_max, R_max]]
    )
    axs[1].set_title('z-y axis plot')
    axs[1].set_xlabel('y')
    axs[1].set_ylabel('z')

    # Show the figure
    plt.tight_layout()

    if galaxy_type == "zoom_in_tolgay":
        plt.savefig(f"{write_base_fdir}/zoom_in/z{round(redshift,1)}/images2/{galaxy_name}.png")
        print(f"Figure saved to: {write_base_fdir}/zoom_in/z{round(redshift,1)}/images2/{galaxy_name}.png")          
    else:
        plt.savefig(f"{write_base_fdir}/{galaxy_type}/z{round(redshift,1)}/images2/{galaxy_name}.png")
        print(f"Figure saved to: {write_base_fdir}/{galaxy_type}/z{round(redshift,1)}/images2/{galaxy_name}.png")    


    """
    Calculate the star age and smoothing length of star particles
    Append the mass resolution
    """

    #################### Calculating the stellar age
    age_star = functions.calculate_stellar_age(scale_factor_star=star_particles_df["scale_factor"], time=header_time)   # in Gyr

    if galaxy_type in ["zoom_in", "zoom_in_tolgay"]:
        star_particles_df["smoothing_length"] = np.ones(len(star_particles_df))*4*1e-3        # 4 pc = 0.004 kpc
        star_particles_df["mass_resolution"] = np.ones(len(star_particles_df))*7070*1e-10     # 1e10 Msolar
    elif galaxy_type == "particle_split":
        star_particles_df["smoothing_length"] = np.ones(len(star_particles_df))*2*1e-3        # 2 pc = 0.002 kpc
        star_particles_df["mass_resolution"] = np.ones(len(star_particles_df))*880*1e-10      # 1e10 Msolar
    elif galaxy_type == "firebox":
        star_particles_df["smoothing_length"] = np.ones(len(star_particles_df))*12*1e-3       # 12 pc = 0.012 kpc
        star_particles_df["mass_resolution"] = np.ones(len(star_particles_df))*6.26e4*1e-10   # 1e10 Msolar
    else:
        print("Wrong galaxy type! Exiting with 1.")
        exit(1)


    star_particles_df["age"] = age_star # in Gyr

    print("star age and smoothing length of the star particles are calculated.")


    #################### Gas Temperature
    # Calculating the gas temperature

    gas_particles_df["temperature"] = functions.gas_temperature_calculator(
        He_mass_fraction = gas_particles_df["He_mass_fraction"],
        electron_abundace_gas = gas_particles_df["electron_abundance"],
        internal_energy_gas = gas_particles_df["internal_energy"]
    )  
    # temperature_gas: K 

    print("Gas temperature is calculated.")


    #################### nH 
    # Calculating the hden 
    Y = gas_particles_df["He_mass_fraction"] / (4 * (1 - gas_particles_df["He_mass_fraction"]))
    mu_theoretical = (1 + 4 * Y) / (1 + Y + gas_particles_df["electron_abundance"]) # From the InternalEnergy of the Gizmo documentation
    ## Reference: http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html#snaps-reading
    gas_particles_df["mu_theoretical"] = mu_theoretical 

    print("varying mu is used!")
    mu = mu_theoretical
    print(f"average mu: {np.sum(mu)/len(mu)}")

    # # TODO: Check out this n_H calculation
    density_gas_unit_converted = gas_particles_df["density"] * 1e10 * constants.Msolar2kg / (constants.kpc2cm)**3 # kg/cm^3
    n_H = density_gas_unit_converted / (mu * constants.proton_mass)  # cm^-3
    gas_particles_df["hden"] = n_H

    print("hden calculated")


    #################### Sobolev length
    time_begin = time()
    print("Sobolev length calculation started.")
    # Initialize an instance of the Meshoid class
    gas_particles_meshoid = Meshoid(
        pos = gas_particles_df[['x', 'y', 'z']].values, 
        m = gas_particles_df['mass'].values, 
        kernel_radius = gas_particles_df['smoothing_length'].values, 
        verbose=True,
        n_jobs = 1,
    )

    # Calculate the density gradient using the D method
    # Assuming you want to calculate the first-order density gradient
    density_gradient = gas_particles_meshoid.D(
        f = gas_particles_df['density'].values, 
        order=1, 
        weighted=True
    )
    print("meshoid calculation for the density gradient finished.")

    # meshoid calculation
    density_gradient_norm = LA.norm(density_gradient, axis=1) 
    gas_particles_df["sobolev_length"] = gas_particles_df["density"] / density_gradient_norm
    gas_particles_df["average_sobolev_smoothingLength"] = (gas_particles_df["sobolev_length"] + gas_particles_df["smoothing_length"]) / 2
    print("Sobolev length calculation finished.")

    time_end = time()
    print(f"Time took to calculate sobolev length is: {time_end - time_begin} seconds.")

    #################### turbulence 
    print("Turbulence calculation started.")
    v_cm = np.array([
        gas_particles_meshoid.KernelAverageDoga(f = gas_particles_df['vx'].values),
        gas_particles_meshoid.KernelAverageDoga(f = gas_particles_df['vy'].values),
        gas_particles_meshoid.KernelAverageDoga(f = gas_particles_df['vz'].values)
    ]).T

    #############################################################
    
    ### Previous method is below
    # # Get the 32 closest neighbour
    # neighbouring_particles_indices = gas_particles_meshoid.NearestNeighbors()
    # neighbouring_particles_distances = gas_particles_meshoid.NeighborDistance()

    # # Calculate the center of mass
    # v_cm = []
    # for row, gas in gas_particles_df.iterrows():
    #     # Calculate the kernel 
    #     q = neighbouring_particles_distances[row] / gas["smoothing_length"]
    #     kernel = [Kernel(qi) for qi in q]
    #     normalized_kernel = np.array(kernel / sum(kernel))
        
    #     # Calculate the center of mass velocity
    #     neighbouring_gas_particles_velocities = gas_particles_df.loc[neighbouring_particles_indices[row], ["vx", "vy", "vz"]] 
        
    #     v_cm.append([
    #         sum(neighbouring_gas_particles_velocities["vx"] * normalized_kernel),
    #         sum(neighbouring_gas_particles_velocities["vy"] * normalized_kernel),
    #         sum(neighbouring_gas_particles_velocities["vz"] * normalized_kernel),
    #     ])

    #     if (row % 1e5 == 1):
    #         print(f"{row} finished. Left {len(gas_particles_df) - row}")
        
    # v_cm = np.array(v_cm)

    #############################################################

    velocity_difference = gas_particles_df[["vx", "vy", "vz"]].to_numpy() - v_cm

    # Compute the Euclidean norm (magnitude of velocity difference) for each particle
    # axis=1 computes the norm along the rows (i.e., for each particle)
    gas_particles_df["turbulence"] = np.linalg.norm(velocity_difference, axis=1)
    print("Turbulence calculation finished.")


    # TODO: Change the sobolev_length and turbulence calculation 
    # gas_particles_df["sobolev_length"] = np.ones(len(gas_particles_df))
    # gas_particles_df["average_sobolev_smoothingLength"] = np.ones(len(gas_particles_df))
    # gas_particles_df["turbulence"] = np.ones(len(gas_particles_df))

    #################### Radius
    radius = ((gas_particles_df["mass"]/gas_particles_df["density"])/(4/3 * np.pi))**(1/3) #kpc

    gas_particles_df["radius"] = radius #kpc

    print("radius calculated")


    ######### Unit conversions
    # Prevent multiple unit conversion 

    # Gas particles
    gas_particles_df[["x", "y", "z"]] *= constants.kpc2pc   # pc
    gas_particles_df["mass"] *= 1e10                        # Msolar
    gas_particles_df["density"] *= (1e10 * constants.M_sun2gr / (constants.kpc2cm)**3)  # gr/cm^3
    gas_particles_df[["smoothing_length", "radius", "sobolev_length", "average_sobolev_smoothingLength"]] *= constants.kpc2pc  # pc

    # Star particles 
    star_particles_df[["x", "y", "z"]] *= constants.kpc2pc  # pc
    star_particles_df[["mass", "mass_resolution"]] *= 1e10  # Msolar
    star_particles_df[["smoothing_length"]] *= constants.kpc2pc # pc
    star_particles_df["age"] *= 1e3 # Milion year



    ######### Saving to a file

    directory_name = "voronoi_1e5"


    # Save zoom_in and zoom_in_dtolgay to zoom_in directory

    if galaxy_type == "zoom_in_tolgay":
        galaxy_type = "zoom_in"

    path_to_save = f"""{write_base_fdir}/{galaxy_type}/z{round(redshift,1)}/{galaxy_name}/{directory_name}"""

    # Check if path exists 
    # Check if the directory exists
    if not os.path.exists(path_to_save):
        # Create the directory if it doesn't exist
        os.makedirs(path_to_save)
        print(f"Directory '{path_to_save}' is created.")
    else:
        print(f"Directory '{path_to_save}' already exists.")


    ### SKIRT Files
    ## Gas
    skirt_gas_header = f"""Gas particles for {galaxy_name} galaxy
    Column 1: x-coordinate (pc)
    Column 2: y-coordinate (pc)
    Column 3: z-coordinate (pc)
    Column 4: smoothing length (pc)
    Column 5: gas mass (Msun)
    Column 6: metallicity (1)
    Column 7: temperature (K)
    Column 8: vx (km/s)
    Column 9: vy (km/s)
    Column 10: vz (km/s)
    """


    skirt_gas_df = gas_particles_df[
        ["x",
        "y",
        "z",
        "smoothing_length",
        "mass",
        "metallicity",
        "temperature",
        "vx",
        "vy", 
        "vz"
        ]
    ]


    filename = f"{path_to_save}/skirt_gas.txt"
    np.savetxt(fname=filename, X = skirt_gas_df, fmt='%.18e', header=skirt_gas_header)
    print(f"File saved as {filename}/")


    ## Star
    skirt_star_header = f"""Star particles for {galaxy_name} galaxy
    Column 1: x-coordinate (pc)
    Column 2: y-coordinate (pc)
    Column 3: z-coordinate (pc)
    Column 4: smoothing length (pc)
    Column 5: vx (km/s)
    Column 6: vy (km/s)
    Column 7: vz (km/s)
    Column 8: initial mass (Msun)
    Column 9: metallicity (1)
    Column 10: age (Myr)
    """ 

    skirt_star_df = star_particles_df[
        ["x",
        "y",
        "z",
        "smoothing_length",
        "vx",
        "vy",
        "vz",
        "mass_resolution",
        "metallicity",
        "age"
        ]
    ]

    filename = f"{path_to_save}/skirt_star.txt"
    np.savetxt(fname=filename, X = skirt_star_df, fmt='%.18e', header=skirt_star_header)
    print(f"File saved as {filename}/")


    ### Comprehensive gas particles 
    ## Gas
    comprehensive_gas_header = f"""Gas particles for {galaxy_name} galaxy
    Column 0: x-coordinate (pc)
    Column 1: y-coordinate (pc)
    Column 2: z-coordinate (pc)
    Column 3: smoothing length (pc)
    Column 4: mass (Msolar)
    Column 5: metallicity (1)
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
    """

    comprehensive_gas_df = gas_particles_df[
        ["x",
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
        "star_formation_rate",
        "turbulence",
        "density",
        "mu_theoretical",
        "average_sobolev_smoothingLength"
        ]
    ]

    filename = f"{path_to_save}/comprehensive_gas.txt"
    np.savetxt(fname=filename, X = comprehensive_gas_df, fmt='%.18e', header=comprehensive_gas_header)
    print(f"File saved as {filename}/")


    ## Star
    comprehensive_star_header = f"""Star particles for {galaxy_name} galaxy
    Column 0: x-coordinate (pc)
    Column 1: y-coordinate (pc)
    Column 2: z-coordinate (pc)
    Column 3: vx (km/s)
    Column 4: vy (km/s)
    Column 5: vz (km/s)
    Column 6: metallicity (1)
    Column 7: mass (Msolar)
    Column 8: age (Myr)
    """

    comprehensive_star_df = star_particles_df[
        ["x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
        "metallicity",
        "mass",
        "age"
        ]
    ]

    filename = f"{path_to_save}/comprehensive_star.txt"
    np.savetxt(fname=filename, X=comprehensive_star_df, fmt='%.18e', header=comprehensive_star_header)
    print(f"File saved as {filename}/")



    return 0



######################## Functions  

def read_firebox(
    snap_dir_file_path, 
    snapshot_number,
    simulation_name,
    file_number
    ):

    # Reading the header
    print('Reading header')
    header_info = readsnap_FIREBox.readsnap(
        sdir=snap_dir_file_path, 
        snum=snapshot_number, 
        simulation_name=simulation_name, 
        file_number=file_number, 
        ptype=0, 
        header_only=1
    )

    hubble      = header_info['hubble']
    redshift    = header_info['redshift']   
    time        = header_info['time']  

    # Reading the gas_particles
    print('\nReading gas_particles')
    gas_particles    = readsnap_FIREBox.readsnap(
        sdir=snap_dir_file_path, 
        snum=snapshot_number, 
        simulation_name=simulation_name, 
        file_number=file_number, 
        ptype=0, 
        cosmological=1
    )
    # Reading the star particles
    print('\nReading star_particles')    
    star_particles   = readsnap_FIREBox.readsnap(
        sdir=snap_dir_file_path, 
        snum=snapshot_number, 
        simulation_name=simulation_name, 
        file_number=file_number, 
        ptype=4, 
        cosmological=1
    )   

    return gas_particles, star_particles, header_info


def read_firebox_seperated_galaxies(base_dir, galaxy_number):

    gas_particles = {}
    star_particles = {}
    header_info = {}

    with h5py.File(f'{base_dir}/gal_{galaxy_number}.hdf5', 'r') as f:
        for key in f['gas']:
            gas_particles[key] = f['gas'][key][:]
        for key in f['star']:
            star_particles[key] = f['star'][key][:]    
        for key in f['header']:
            dataset = f['header'][key]
            header_info[key] = dataset[()] if dataset.shape == () else dataset[:]

    return gas_particles, star_particles, header_info


def Kernel(q):
    """
    Un-normalized cubic-spline kernel function

    Arguments:
        q - array containing radii at which to evaluate the kernel,
        scaled to the kernel support radius (between 0 and 1)
    """
    if q <= 0.5:
        return 1 - 6 * q**2 + 6 * q**3
    elif q <= 1.0:
        return 2 * (1 - q) ** 3
    else:
        return 0.0



if __name__ == "__main__":

    galaxy_name = sys.argv[1] # Not used in Firebox. Give a dummy value
    galaxy_type = sys.argv[2] 
    file_number = int(sys.argv[3]) # Not used is zoom-in and particle_split. Give a dummy value It is from 0 to 999. integer. Which Firebox galaxy...
    redshift = str(sys.argv[4])

    # redshift = "2.0"  # Never compare float numbers...
    # Determine the snapshot number 
    if galaxy_type in ["zoom_in", "zoom_in_tolgay", "particle_split"]:
        # Determine snapshot_number
        if redshift == "0.0":
            snapshot_number = 600     # z = 0.0
        elif redshift == "0.5":
            snapshot_number = 381     # z = 0.5
        elif redshift == "1.0":
            snapshot_number = 277     # z = 1.0
        elif redshift == "2.0":
            snapshot_number = 172     # z = 2.0
        elif redshift == "3.0":
            snapshot_number = 120     # z = 3.0
        else:
            print(f"Exiting... Redshift is wrong. The given galaxy type is {redshift}")
            sys.exit(2)
        # Determine snap_dir_file_path
        if galaxy_type == "zoom_in":
            snap_dir_file_path = f'/mnt/raid-project/murray/FIRE/FIRE_2/Fei_analysis/md/{galaxy_name}/output'
        elif galaxy_type == "zoom_in_tolgay":
            snap_dir_file_path = f'/fs/lustre/project/murray/scratch/tolgay/metal_diffusion/{galaxy_name}/output'
        elif galaxy_type == "particle_split":
            snap_dir_file_path = f"/fs/lustre/project/murray/FIRE/FIRE_2/{galaxy_name}/output"
        else:
            print(f"Exiting... Galaxy type is wrong. The given galaxy type is {galaxy_type}")

    # TODO: Fix the path to file directories
    elif galaxy_type in ["firebox"]:
        snap_dir_file_path = f'/scratch/m/murray/dtolgay/post_processing_fire_outputs/firebox_halo_finder/z{redshift}'
        if redshift == "0.0":
            snapshot_number = 1200     # z = 0.0
        elif redshift == "0.5":
            print(f"Exiting... Currently there are no z=0.5 galaxies... {redshift}")
            sys.exit(2)                
        elif redshift == "1.0":
            snapshot_number = 554     # z = 1.0
        elif redshift == "2.0":
            snapshot_number = 344     # z = 2.0
        elif redshift == "3.0":
            snapshot_number = 240     # z = 3.0
        else:
            print(f"Exiting... Redshift is wrong. The given galaxy type is {redshift}")
            sys.exit(2)        
    else:
        print(f"Exiting... Galaxy type is wrong. The given galaxy type is {galaxy_type}")
        sys.exit(1)

    if galaxy_type in ["zoom_in", "zoom_in_tolgay", "particle_split"]:

        snapshot_number = '%03d' %snapshot_number    

        print(f"\n\n------------------ For {galaxy_name} ------------------\n")
        # Reading particles
        print('Reading gas particles')
        gas_particles  = readsnap.readsnap(snap_dir_file_path, snapshot_number, 0, cosmological=1)
        print('Reading star particles')
        star_particles = readsnap.readsnap(snap_dir_file_path, snapshot_number, 4, cosmological=1)

        # Reading Header to get scale_factor, hubble_run, and redshift information 
        print('Reading header')
        header_info = readsnap.readsnap(snap_dir_file_path, snapshot_number, 0, header_only=1)

        main(gas_particles, star_particles, header_info, galaxy_type, galaxy_name, snapshot_number)


    elif galaxy_type == "firebox":

        galaxy_name = f"gal{file_number}"        
    
        snapshot_number = '%03d' %snapshot_number
        # simulation_name = "FB15N1024"  
        simulation_name = "snapshot" # This is the new file name in Lichen's FIRE_CO directory

        # gas_particles, star_particles, header_info = read_firebox(
        #     snap_dir_file_path, 
        #     snapshot_number,
        #     simulation_name,
        #     file_number
        # )            
        

        gas_particles, star_particles, header_info = read_firebox_seperated_galaxies(
            base_dir = snap_dir_file_path,
            galaxy_number = file_number
        )

    
        main(gas_particles, star_particles, header_info, galaxy_type, galaxy_name, snapshot_number)    


    else: 
        print("There is a problem with the galaxy type!!")