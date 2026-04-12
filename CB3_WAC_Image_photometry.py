import numpy as np
import matplotlib.pyplot as plt
import csv
import os

import functs
import Get_image_info

image_files = Get_image_info.CB3_WAC_image_files
planet_cassini_distance = Get_image_info.CB3_WAC_planet_cassini_distance
exposure_time = Get_image_info.CB3_WAC_exposure
image_gain = Get_image_info.CB3_WAC_gain
phase_angle = Get_image_info.CB3_WAC_phase_angle
OBStime = Get_image_info.CB3_WAC_OBStime

image_files_length = len(image_files)

camera_filter = "CB3"
camera = "WAC" # !!! Note: always capitalize
pixel_angular_size = 59.749e-6 #[rad per pix]

#------------------------------------------------------------

path = f"Titan_images/{camera_filter}_{camera}/"

if not os.path.exists(path + "apertures"):
    print("No apertures dir. found")
    os.makedirs(path + "apertures")
    print(path + "apertures CREATED")
    
if not os.path.exists(path + "aperture_rejects"):
    print("No aperture_rejects dir. found")
    os.makedirs(path + "aperture_rejects")
    print(path + "aperture_rejects CREATED")
    
if not os.path.exists(path + "images"):
    print("No images dir. found")
    os.makedirs(path + "images")
    print(path + "images CREATED")
    
if not os.path.exists(path + "Canny_edge"):
    print("No Canny_edge dir. found")
    os.makedirs(path + "Canny_edge")
    print(path + "Canny_edge CREATED")

#------------------------------------------------------------

n_images_with_different_gain = 0
images_with_different_gain = []

n_failed_called_images = 0
failed_called_images = []

total_rejected_images = 0
total_processed_images = 0

no_centers_found = 0
no_center_images_found = []

fields = ['Image file', 'Observation time [MJD]', 'Phase angle [deg]', 'Planet-satellite distance [km]', 'Aperture count', 'Aperture area', 'Sky median', 'Sky std', 'Source count', 'Exposure time [s]', 'Gain', 'Target Radius [pix]']

with open(path + f"{camera_filter}_{camera}_photometry_data.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)
    
for n in range(image_files_length):
    print(f"IMG {n+1}/{image_files_length}: {image_files[n]}")
    
    if (image_gain[n] == '95'):
        gain = 95
        satellite_distance = planet_cassini_distance[n] * 2
    elif (image_gain[n] == '29'):
        gain = 29
        satellite_distance = planet_cassini_distance[n]
    elif (image_gain[n] == '12'):
        gain = 12
        satellite_distance = planet_cassini_distance[n]
    elif (image_gain[n] == '215'):
        gain = 215
        satellite_distance = planet_cassini_distance[n] * 4
    else: 
        images_with_different_gain.append(str(image_files[n]))
        n_images_with_different_gain += 1
        
        gain = int(image_gain[n])
        satellite_distance = planet_cassini_distance[n]

    try:
        image_data = functs.get_image_data(image_files[n], camera_filter, camera)
    except Exception as e:
        print(f"{e}\n")
        failed_called_images.append(str(image_files[n]))
        n_failed_called_images += 1
        
    functs.save_image_as_png(image_files[n], camera_filter, camera)
    
    try:
        aperture_center, Canny_edge_radius = functs.find_aperture_center(image_files[n], camera_filter, camera, 100, 10)
        
    except Exception as e:
        print(f"{e}")
        # "No center could be found"
        no_center_images_found.append(str(image_files[n]))
        no_centers_found += 1
        total_rejected_images += 1
        
        plt.imsave(f"Titan_images/{camera_filter}_{camera}/Canny_edge/Failed_CannyEdge_{image_files[n]}.png",
                   image_data, cmap='gray')
        
    else:
        print(f'Canny edge center: {aperture_center}')
        print(f'Canny edge radius: {Canny_edge_radius}')
        extra_aperture_radius = 10
        annulus_width = 10
        try:
            aperture, annulus = functs.get_titan_aperture(image_data, 
                                                          satellite_distance, 
                                                          pixel_angular_size,
                                                          aperture_center, 
                                                          extra_aperture_radius, 
                                                          annulus_width)
            
        except Exception as e:
            print(f"{e}")
            plt.savefig(path + f"aperture_rejects/{image_files[n]}.png", dpi = 120, bbox_inches='tight')
            plt.show()
            
            total_rejected_images += 1
                
        else:
            plt.savefig(path + f"apertures/{image_files[n]}.png", dpi = 120, bbox_inches='tight')
            plt.show()
            
            total_processed_images += 1
            
            planet_radius = 2675 #km
            target_radius = functs.CCD_Target_radius(planet_radius, planet_cassini_distance[n], pixel_angular_size)
            
            aperture_count, aperture_area, sky_median, sky_std, source_count = functs.image_photometry(image_data, aperture, annulus)
            photometric_data = np.array([image_files[n], OBStime[n], phase_angle[n], planet_cassini_distance[n], aperture_count, aperture_area, sky_median, sky_std, source_count, exposure_time[n], gain, target_radius])
        
            with open(path + f"{camera_filter}_{camera}_photometry_data.csv", "a") as csvfile:
                csvwriter = csv.writer(csvfile)
                csvwriter.writerow(photometric_data)
        
print(f'\nTotal no. of images: {image_files_length}')
print(f'No. of rejected images: {total_rejected_images}')
print(f'No. of processed images: {total_processed_images}')

print(f'No. of failed called images: {n_failed_called_images}')
print(f'No. of images that connot find a center: {no_centers_found}')
