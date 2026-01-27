import numpy as np
import matplotlib.pyplot as plt
import csv
import traceback

import functs
import Get_image_info

image_files = Get_image_info.CB2_WAC_image_files
image_numbers = Get_image_info.CB2_WAC_image_numbers
planet_cassini_distance = Get_image_info.CB2_WAC_planet_cassini_distance
exposure_time = Get_image_info.CB2_WAC_exposure
image_gain = Get_image_info.CB2_WAC_gain
phase_angle = Get_image_info.CB2_WAC_phase_angle
image_files_length = len(image_files)

# !!! Check folder names
camera_filter = "CB2"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix

total_rejected_images = 0
total_processed_images = 0

num_of_images_with_different_gain = 0
images_with_different_gain = []

fields = ['Image file', 'Satellite distance [km]', 'Exposure time [s]', 'Gain [e per DN]',
          'Aperture count [DN]', 'Aperture area', 'Annulus radius', 'Sky median [DN]',
          'Aperture area sky count [DN]', 'Source count [DN]', 'Source intenisty [e per s per m2]', 'Phase angle [deg]']

csvpath = f"Titan_images/{camera_filter}_{camera}/"

with open(csvpath + 'photometry_data.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)

for n in range(image_files_length):
    print(f"IMG {n+1}: {image_files[n]}")
    
    if (image_gain[n] == '95'):
        gain = 95
        satellite_distance = planet_cassini_distance[n] * 2
    elif (image_gain[n] == '29'):
        gain = 29
        satellite_distance = planet_cassini_distance[n]
    else:
        print("!!! Found image with different gain number other than 29 or 95 !!!")
        images_with_different_gain.append(str(image_files[n]))
        satellite_distance = planet_cassini_distance[n]
        num_of_images_with_different_gain += 1

    image_data = functs.get_image_data(image_files[n], camera_filter, camera)
    
    try:
        aperture, annulus = functs.get_titan_aperture(image_data, satellite_distance, pixel_angular_size)
        plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_aperture.png", dpi = 120, bbox_inches='tight')
        plt.show()
    except Exception:
        traceback.print_exc()
        plt.savefig(f"Titan_images/{camera_filter}_{camera}/aperture_rejects/{image_files[n]}_aperture.png", dpi = 120, bbox_inches='tight')
        plt.show()
        total_rejected_images += 1
    else:
        print("Titan found inside the CCD image.\n")
        photometric_data = functs.image_photometry(image_data, satellite_distance, exposure_time[n], gain, aperture, annulus)
        photometric_data = [image_files[n]] + photometric_data + [phase_angle[n]]
        photometric_data = np.array([photometric_data])
        
        with open(csvpath + 'photometry_data.csv', 'a') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerows(photometric_data)
        
        total_processed_images += 1
        
print(f"Total number of images: {image_files_length}")
print(f"Number of rejected images: {total_rejected_images}")
print(f"Number of processed images: {total_processed_images}")
