#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pdr
from astropy.io import fits

import functs
import Get_image_info


image_files = Get_image_info.CB2_narrow_image_files
image_files_length = len(image_files)
print(f"Total number of image files: {image_files_length}\n")

# !!! Check folder names
camera_filter = "CB2"
camera = "NAC" # !!! Note: always capitalize

for n in range(image_files_length):
    print(f"IMG {n+1}: {image_files[n]}")
    functs.pdr_to_fits_images(image_files[n], camera_filter, camera)
    functs.plot_fits_images(image_files[n], camera_filter, camera)
    functs.plot_number_of_rejected_pixels(image_files[n], camera_filter, camera)
    
print(f"{camera_filter} {camera} Camera images converted successfully!")    
    
    