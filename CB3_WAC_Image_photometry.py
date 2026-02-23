import numpy as np
import matplotlib.pyplot as plt
import csv

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

xaxis_aperture = "improved"
yaxis_aperture = "improved"
annulus_width = 30
#------------------------------------------------------------

total_rejected_images = 0
image_index_rejected = []
total_processed_images = 0
#image_index_processed = []

num_of_images_with_different_gain = 0
images_with_different_gain = []

failed_image_call = 0
failed_image_call_file = []

fields = ['Image file', 'Observation time [MJD]', 'Phase angle [deg]', 'Planet-satellite distance [km]', 'Aperture count', 'Aperture area', 'Sky median', 'Sky std', 'Source count', 'Exposure time [s]', 'Gain']

path = f"Titan_images/{camera_filter}_{camera}/"

xaxis_aperture = "improved"
yaxis_aperture = "improved"
annulus_width = 15

with open(path + f"{camera_filter}_{camera}_photometry_data.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)
    
# Inital iteration
for n in range(image_files_length):
    print(f"IMG {n+1}: {image_files[n]}")
    
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
        print("Image has gain different than stated")   
        images_with_different_gain.append(str(image_files[n]))
        num_of_images_with_different_gain += 1
        
        gain = int(image_gain[n])
        satellite_distance = planet_cassini_distance[n]

    try:
        image_data = functs.get_image_data(image_files[n], camera_filter, camera)
    except Exception as e:
        print(f"{e}\n")
        failed_image_call_file.append(str(image_files[n]))
        failed_image_call += 1
        
    image_size = image_data.shape[0]
    
    try:
        xaxis_shift_factor = 1
        yaxis_shift_factor = 0.6
        
        if (image_size >= 1000):
            extra_aperture_radius = 30
            annulus_width = 40
        else:
            extra_aperture_radius = 20
            annulus_width = 35
            
        aperture, annulus = functs.get_titan_aperture(image_data, satellite_distance, pixel_angular_size, 
                                                              xaxis_shift_factor, yaxis_shift_factor,
                                                              extra_aperture_radius, annulus_width)
    except Exception as e:
        print(f"{e}")
        plt.savefig(path + f"aperture_rejects/{image_files[n]}.png", dpi = 120, bbox_inches='tight')
        plt.show()
        
        total_rejected_images += 1
        image_index_rejected.append(n)
    else:
        plt.savefig(path + f"apertures/{image_files[n]}.png", dpi = 120, bbox_inches='tight')
        plt.show()
        
        total_processed_images += 1
            
        aperture_count, aperture_area, sky_median, sky_std, source_count = functs.image_photometry(image_data, aperture, annulus)
        photometric_data = np.array([[image_files[n], OBStime[n], phase_angle[n], planet_cassini_distance[n], aperture_count, aperture_area, sky_median, sky_std, source_count, exposure_time[n], gain]])

        with open(path + f"{camera_filter}_{camera}_photometry_data.csv", "a") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(photometric_data)
        
print(f'Total no. of images: {image_files_length}')
print(f'No. of rejected images: {total_rejected_images}')
print(f'No. of processed images: {total_processed_images}')

