#!/usr/bin/env python
# coding: utf-8

# Importing necessary librarires
import os
import sys

import h5py
import numpy as np
from numpy import linalg as LA


from math import pi

from tools_tolgay import functions, readsnap, readsnap_FIREBox, readsnap_tripleLatte, constants  
from tools_tolgay.filter_rotate_galaxy import filter_rotate_galaxy

import random 
import sys # To run the batch script

from time import time

from tools_tolgay.meshoid import Meshoid



def main(gas_particles, star_particles, header_info, galaxy_type, galaxy_name, snapshot_number, operating_cluster):

    print(f"---------------------------------------------------- {galaxy_name} ----------------------------------------------------")
    
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
        cluster_name=operating_cluster
    )


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
    gas_particles_df[["smoothing_length", "radius"]] *= constants.kpc2pc  # pc

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

    # Determine the path to save.
    if operating_cluster == "cita":
        path_to_save = f"/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{round(redshift,1)}/{galaxy_name}/{directory_name}"
    elif operating_cluster == "niagara":
        path_to_save = f"/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{round(redshift,1)}/{galaxy_name}/{directory_name}"
    elif operating_cluster == "trillium":
        path_to_save = f"/scratch/dtolgay/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{round(redshift,1)}/{galaxy_name}/{directory_name}"

    # Check if path exists 
    # Check if the directory exists
    if not os.path.exists(path_to_save):
        # Create the directory if it doesn't exist
        os.makedirs(path_to_save)
        print(f"Directory '{path_to_save}' is created.")
    else:
        print(f"Directory '{path_to_save}' already exists.")


    ### Saving the other properties of the particles to a csv
    header = """
    Other properties of the particles for {galaxy_name} galaxy
    Note: The id_fire is the particle ID in the FIRE simulation. It is used to identify particles across different snapshots.
    Column 1: x-coordinate (pc)
    Column 2: y-coordinate (pc)
    Column 3: z-coordinate (pc)
    Column 4: neutral hydrogen fraction (1)
    Column 5: electron abundance (1)
    Column 6: metallicity (1)
    Column 7: mass (Msolar)
    Column 8: He mass fraction (1)
    Column 9: C mass fraction (1)
    Column 10: N mass fraction (1)
    Column 11: O mass fraction (1)
    Column 12: Ne mass fraction (1)
    Column 13: Mg mass fraction (1)
    Column 14: Si mass fraction (1)
    Column 15: S mass fraction (1)
    Column 16: Ca mass fraction (1)
    Column 17: Fe mass fraction (1)
    """.format(galaxy_name=galaxy_name) # Not used in the code. Just to record the units for me. 
    
    other_properties_df = gas_particles_df[
        [
            'x', 'y', 'z', 
            'neutral_hydrogen_fraction', 'electron_abundance', 'metallicity', 'mass',
            'He_mass_fraction', 'C_mass_fraction', 'N_mass_fraction', 'O_mass_fraction',
            'Ne_mass_fraction', 'Mg_mass_fraction', 'Si_mass_fraction', 'S_mass_fraction',
            'Ca_mass_fraction', 'Fe_mass_fraction',
        ]].copy()

    # Save as csv file
    filename = f"{path_to_save}/other_properties.csv"
    other_properties_df.to_csv(
        filename,
        index=False,
        float_format="%.7e" 
    )

    # Print the file name and path
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
    operating_cluster = "trillium"  # cita, niagara or trillium. The cluster where the data is stored.

    # Get the snapshot directory and number
    snap_dir_file_path, snapshot_number = functions.return_snapshot_path_and_number(galaxy_type, galaxy_name, redshift, operating_cluster)

    print("snap_dir_file_path: ", snap_dir_file_path)
    print("snapshot_number: ", snapshot_number)

    if galaxy_type in ["zoom_in", "zoom_in_tolgay", "particle_split"]:

        snapshot_number = '%03d' %snapshot_number    

        print(f"\n\n------------------ For {galaxy_name} ------------------\n")
        # Reading particles
        print('Reading gas particles')
        gas_particles  = readsnap.readsnap(snap_dir_file_path, snapshot_number, 0, cosmological=1, loud=1)

        print('Reading star particles')
        star_particles = readsnap.readsnap(snap_dir_file_path, snapshot_number, 4, cosmological=1)

        # Reading Header to get scale_factor, hubble_run, and redshift information 
        print('Reading header')
        header_info = readsnap.readsnap(snap_dir_file_path, snapshot_number, 0, header_only=1)

        main(gas_particles, star_particles, header_info, galaxy_type, galaxy_name, snapshot_number, operating_cluster)


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
        

        if redshift in ["0.5", "1.0", "2.0", "3.0"]:
            gas_particles, star_particles, header_info = read_firebox_seperated_galaxies(
                base_dir = snap_dir_file_path,
                galaxy_number = file_number
            )
        elif redshift in ["0.0"]: 
            gas_particles, star_particles, header_info = read_firebox(
                snap_dir_file_path, 
                snapshot_number,
                simulation_name,
                file_number
            )
        else: 
            print(f"Exiting... Redshift is wrong. The given galaxy type is {redshift}")
            sys.exit(2)

        main(gas_particles, star_particles, header_info, galaxy_type, galaxy_name, snapshot_number, operating_cluster)    


    else: 
        print("There is a problem with the galaxy type!!")