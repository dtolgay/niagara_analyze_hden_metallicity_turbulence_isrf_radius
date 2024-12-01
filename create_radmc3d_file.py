import sys
sys.path.append("/scratch/m/murray/dtolgay")
from tools import functions_readfiles as readFiles

import numpy as np
import pandas as pd 

def main(galaxy_name, galaxy_type, redshift, directory_name):

    base_dir = "/home/m/murray/dtolgay/scratch/post_processing_fire_outputs/skirt/runs_hden_radius"

    ## Read gas particles 
    gas = readFiles.read_otherProperties(
        base_dir, 
        galaxy_name, 
        galaxy_type, 
        redshift, 
        directory_name, 
        file_name="otherProperties_hybridInterpolator.txt"
    )

    ## Read the gas temperature 

    ## Merge two dataframes together.

    # Calculate Mh2 
    gas['h2_number_density'] = gas['hden'] * gas['fh2'] # 1/cm3

    ## Change the units 
    # velocity 
    gas['vx'] *= 1e5 # cm/sec
    gas['vy'] *= 1e5 # cm/sec
    gas['vz'] *= 1e5 # cm/sec


    ## Write to a file 

    return None 


if __name__ == "__main__":
    galaxy_name = sys.argv[1]
    galaxy_type = sys.argv[2]
    redshift = sys.argv[3]
    directory_name = "voronoi_1e6"

    main(galaxy_name, galaxy_type, redshift, directory_name)