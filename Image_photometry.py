
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry

import functs
import Get_image_info

image_files = Get_image_info.CB2_wide_image_files
image_numbers = Get_image_info.CB2_wide_image_numbers
titan_distance = Get_image_info.CB2_wide_planet_cassini_distance
exposure_time = Get_image_info.CB2_wide_exposure
image_files_length = len(image_files)
print(f"Total number of image files: {image_files_length}\n")

# !!! Check folder names
camera_filter = "CB2"
camera = "Wide" # !!! Note: always capitalize

pixel_angular_size_wide_camera = 59.749e-6 #rad per pix

total_rejected_images = 0
total_processed_images = 0

processed_image_ids = []
processed_image_skyreduced_apature_count = []
processed_image_count_intenisty = []

for n in range(image_files_length):
    print(f"IMG {n+1}: {image_files[n]}")
    
    try:
        aperture, annulus_aperture = functs.get_titan_aperture(image_files[n], image_numbers[n], titan_distance[n], pixel_angular_size_wide_camera)
    except:
        print("Titan is outside of the CCD image. Cannot produce an aperture for this image.\n")
        total_rejected_images += 1
    else:
        print("Titan found inside the CCD image.\n")
        image_number, titan_aperture_sky_reduced, titan_intensity = functs.image_photometry(image_files[n], image_numbers[n], titan_distance[n], exposure_time[n], aperture, annulus_aperture)
        processed_image_ids.append(image_number)
        processed_image_skyreduced_apature_count.append(titan_aperture_sky_reduced)
        processed_image_count_intenisty.append(titan_intensity)
        total_processed_images += 1
        
print(f"Number of rejected images: {total_rejected_images}")
print(f"Number of processed images: {total_processed_images}")

