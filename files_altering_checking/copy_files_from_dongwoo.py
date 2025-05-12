

import os


def copy_from_dongwoo_to_my_scratch(redshift, directory_name, filename):


    # Firebox
    galaxy_type = "firebox" 
    for gal_number in range(1000):
        galaxy_name = f"gal{gal_number}"
        source_dir = f"/gpfs/fs0/scratch/r/rbond/dongwooc/scratch_rwa/doga/runs_hden_radius/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}/{filename}"

        # Destination directory
        destination_dir = f"/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/runs_hden_radius/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}/"        
        # Check if folder exists
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)
            print(f"Directory {destination_dir} created")

        # Copy file to destination directory
        destination_file = os.path.join(destination_dir, filename)
        os.system(f"scp {source_dir} {destination_file}")
        print(f"{source_dir} -> {destination_file}")

    return None 


def copy_skirt_files_from_dongwoo_to_migrate_files_from_cita(redshift, directory_name):


    # Firebox
    galaxy_type = "firebox" 
    for gal_number in range(1000):
        galaxy_name = f"gal{gal_number}"

        skirt_file_names = {
            "total_fits": {
                "destination_directory": "_i0_total",
                "file_name_without_galaxy_name": "_i0_total.fits",
                "extension": "fits",
            },
            "wavelengths": {
                "destination_directory": "_grid_radiationField_wavelengths",
                "file_name_without_galaxy_name": "_grid_radiationField_wavelengths.dat",
                "extension": "dat",
            },
            "sed":{
                "destination_directory": "_i0_sed",
                "file_name_without_galaxy_name": "_i0_sed.dat",
                "extension": "dat",
            }
        }


        for skirt_file in skirt_file_names.keys():
            source_dir = f"/gpfs/fs0/scratch/r/rbond/dongwooc/scratch_rwa/doga/runs_hden_radius/{galaxy_type}/z{redshift}/{galaxy_name}/{directory_name}"
            source_file_name = f"{galaxy_name}{skirt_file_names[skirt_file]['file_name_without_galaxy_name']}"
            
            
            destination_directory = skirt_file_names[skirt_file]["destination_directory"]   
            destination_filename = f"{galaxy_name}.{skirt_file_names[skirt_file]['extension']}"
            destination_dir = f"/scratch/m/murray/dtolgay/post_processing_fire_outputs/skirt/migrated_files_from_cita/{destination_directory}"


            # Copy file to destination directory
            soruce_file = os.path.join(source_dir, source_file_name)
            destination_file = os.path.join(destination_dir, destination_filename)
            os.system(f"scp {soruce_file} {destination_file}")
            print(f"{soruce_file} -> {destination_file}")    


    return None


if __name__ == "__main__":

    redshift = "2.0"
    directory_name = "voronoi_1e5"
    # filename = "cloudy_gas_particles.txt"
    # filename = "comprehensive_star.txt"
    
    # copy_from_dongwoo_to_my_scratch(
    #     redshift=redshift, 
    #     directory_name=directory_name, 
    #     filename=filename
    # )

    copy_skirt_files_from_dongwoo_to_migrate_files_from_cita(
        redshift=redshift, 
        directory_name=directory_name, 
    )    