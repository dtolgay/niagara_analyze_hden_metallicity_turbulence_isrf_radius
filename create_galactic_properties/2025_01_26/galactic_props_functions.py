import sys 
sys.path.append("/mnt/raid-cita/dtolgay/FIRE/post_processing_fire_outputs")
from tools import constants # type: ignore 

import astropy.units as units
from astropy.io import fits
import numpy as np
import pandas as pd 
from pprint import pprint

import matplotlib
import matplotlib.pyplot as plt

def calculate_coordinates_from_fits_file(hdul):
    
    # Read the header and the data 
    fits_data = hdul[0].data    
    header = hdul[0].header
    
#     pprint(header)
    
    ## Find the coordinates 
    # Extract spatial axis parameters from the header
    crpix1, crval1, cdelt1 = header['CRPIX1'], header['CRVAL1'], header["CDELT1"]  # X-axis
    crpix2, crval2, cdelt2 = header['CRPIX2'], header['CRVAL2'], header["CDELT2"]  # Y-axis
    
    # Generate X and Y coordinates
    x_coords = (np.arange(fits_data.shape[1]) + 1 - crpix1) * cdelt1 + crval1
    y_coords = (np.arange(fits_data.shape[2]) + 1 - crpix2) * cdelt2 + crval2
    #### Convert coordinates to pc
    # arcseconds -> radians 
    arcseconds2radians = 4.84814e-6
    x_coords *= arcseconds2radians
    y_coords *= arcseconds2radians
    # radians -> pc 
    distance = 10 # Mpc
    Mpc2parsec = 1e6 # Mpc -> pc
    x_coords *= distance * Mpc2parsec  # pc
    y_coords *= distance * Mpc2parsec  # pc
    
    # Meshgrid to calculate the coordinate of each pixel.
    X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')
    X = X.ravel() # pc
    Y = Y.ravel() # pc
    
    return X, Y

def seperate_into_annuli(data):
    
    Rgal = np.sqrt(data['x']**2 + data['y']**2)
    bins = np.linspace(0, 20e3, num=101) # The boundaries for 100 annuli's are created
    
    # Calculate digitized bins 
    digitized_indices_for_gas_bins = np.digitize(Rgal, bins)
    
    # Calculate center Rgal for every annuli
    annuli_center_Rgas = (bins[:-1] + bins[1:])/2  
    
    return digitized_indices_for_gas_bins, annuli_center_Rgas 

def half_light_radius(band, file_names, plot=False):
    
    # Read the fits emission file 
    hdul = fits.open(file_names['i0_total_fits'])
    fits_data = hdul[0].data
    
    # Read the wavelengths file 
    wavelengths_unfiltered = np.loadtxt(file_names['skirt_wavelengths'])[:,0] * units.micron
    
    if plot:
        plt.hist(np.log10(wavelengths_unfiltered.to(units.nm).value), bins=100)
        plt.axvline(np.log10(band['max'].value), color='red', linestyle='--', label='Max')
        plt.axvline(np.log10(band['min'].value), color='orange', linestyle='--', label='Min')
        plt.legend()
        plt.show()

    # Find the indices of the wavelengths inside the band 
    condition1 = wavelengths_unfiltered <= band['max']
    condition2 = wavelengths_unfiltered >= band['min']
    indices = np.where(condition1 & condition2)[0]
    # Filter the wavelengths 
    wavelengths = wavelengths_unfiltered[indices].value
    
    # Filter the fits file 
    fits_data = fits_data[indices,:,:]
    
    # Integrate the emission data 
    integrated_data = np.trapz(fits_data, wavelengths, axis=0) # W/m²/arcsec²
    log_integrated_data = np.log10(integrated_data)

    # Create df
    X, Y = calculate_coordinates_from_fits_file(hdul = hdul)
    size_of_a_rectangle = (Y[1] - Y[0])**2 * constants.pc2m**2 # m2 -- The size is the same 
    data = pd.DataFrame(np.array([X,Y,integrated_data.ravel()]).T, columns=["x", "y", "flux"]) 
    data["luminosity"] = data["flux"] * size_of_a_rectangle # Watts
    # Luminosity bulunmadan da sadece flux ile half light radius bulunulabilirdi. Her pixel icin alan sabit oldugunda, flux 
    # sabit bir sayi ile carpiliyor. Dolayisiyla, flux veya luminosity kullanildiginda degisen bir sey olmuyor.
    
    if plot:
        plt.hist2d(data['x'], data['y'], bins=[1024, 1024], weights=np.log10(data['flux']), cmap='inferno')
        plt.colorbar(label='(W/m²/arcsec²)')
        plt.xlabel('X (pc)')
        plt.ylabel('Y (pc)')
        plt.show()    
    
    
    # Data is created. Now seperate data into annuli
    digitized_indices_for_gas_bins, annuli_center_Rgas = seperate_into_annuli(data = data)
    
    # Calculate flux in each annuli
    luminosity_in_annuli = []
    for annuli_number in range(min(digitized_indices_for_gas_bins), max(digitized_indices_for_gas_bins)):
        indices = np.where(digitized_indices_for_gas_bins == annuli_number)
        luminosity_in_annuli.append(sum(data.iloc[indices]['luminosity']))
    
    if plot: 
        plt.scatter(annuli_center_Rgas, luminosity_in_annuli)
        plt.xlabel("Rgal [pc]")
        plt.ylabel("Luminosity [Watts]")
        plt.show()
        
    # Find half light radius
    total_luminosity = sum(data["luminosity"])
    integrated_luminosity_upto_certain_r = 0
    for i in range(len(luminosity_in_annuli)):
        integrated_luminosity_upto_certain_r += luminosity_in_annuli[i]
        if integrated_luminosity_upto_certain_r > total_luminosity/2:
            R_half_life = annuli_center_Rgas[i] # [pc] half light radius is where the half of the total luminosity is reached.
            break
    
    return R_half_life 

def filter_particles(particles_df, condition):
    
    r = np.sqrt(particles_df['x']**2 + particles_df['y']**2)
    conditions = (r < condition['r_max']) & (abs(particles_df['z']) < condition['z_max'])
    
    return particles_df[conditions].copy()

def disk_pressure_calculator(gas_particles, star_particles, filtering_condition, band, file_names, is_plot=False):
        
    # Filter the particles
    gas_particles = filter_particles(particles_df = gas_particles, condition = filtering_condition)
    star_particles = filter_particles(particles_df = star_particles, condition = filtering_condition)
    
    
    if is_plot:
        particles = star_particles.copy()
        Rmax = 20e3
        plt.hist2d(
            x = particles['x'],
            y = particles['y'],
            bins = 500,
            weights = particles['mass'],
            norm=matplotlib.colors.LogNorm(),        
            range = [[-Rmax, Rmax], [-Rmax, Rmax]]
        )
        plt.xlabel('x [pc]')
        plt.ylabel('y [pc]')
        plt.show()

        plt.hist2d(
            x = particles['x'],
            y = particles['z'],
            bins = 500,
            weights = particles['mass'],
            norm=matplotlib.colors.LogNorm(),        
            range = [[-Rmax, Rmax], [-Rmax, Rmax]]
        )    
        plt.xlabel('x [pc]')
        plt.ylabel('z [pc]')
        plt.show()
    
    
    #### Calculate the pressure 
    area = np.pi * filtering_condition['r_max']**2 # pc2
#     volume = area * filtering_condition['z_max'] #pc3 # Not used now. Ellison's method is used
    ## Column Density calculation
    sigma_gas = sum(gas_particles['mass']) / area # Msolar / pc2
    ## Star Density calculation
#     rho_star_dtolgay = sum(star_particles['mass']) / volume # Msolar / pc3 # Not used now. Ellison's method is used
    sigma_star = sum(star_particles['mass']/ area) # Msolar / pc2
    # Calculate the half light radius 
    R_half_light = half_light_radius(band, file_names, plot=is_plot) # pc   
    # Rstar eqn. 4. in Ellison
    R_star = R_half_light / 1.68 # [pc]
    rho_star = sigma_star / (0.54 * R_star) # [Msolar / pc2] eqn. 3 in Ellison
    
    ### Disk pressure calculation
    ## Convering to SI units 
    sigma_gas_unit_converted = sigma_gas * constants.Msolar2kg / constants.pc2m**2
    rho_star_unit_converted = rho_star * constants.Msolar2kg / constants.pc2m**3
    ## Gas pressure 
    Pgas = np.pi * constants.gravitational_constant * sigma_gas_unit_converted**2 / 2 # [N/m2 = J/m3]
    Pgas *= 1 / constants.kb / constants.m2cm**3 # [K/cm3]
    
    ## Star pressure 
    # Assuming constant σ_z 
    velocity_dispersion_in_z_direction = 11*1e3 # 11 km/s -- Ellison 2024
    Pstar = sigma_gas_unit_converted * np.sqrt(2 * constants.gravitational_constant * rho_star_unit_converted) * velocity_dispersion_in_z_direction # [N/m2]
    Pstar *= 1 / constants.kb / constants.m2cm**3 # [K/cm3] 
    
    ## Total pressure 
    Ptotal = Pgas + Pstar # [K/cm3]
    
    return Pgas, Pstar, Ptotal
