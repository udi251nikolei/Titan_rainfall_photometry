import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry
from astropy.visualization import simple_norm
from astropy.visualization import ZScaleInterval

import functs
import Get_image_info

image_files = Get_image_info.CB2_NAC_image_files
image_numbers = Get_image_info.CB2_NAC_image_numbers
planet_cassini_distance = Get_image_info.CB2_NAC_planet_cassini_distance
exposure_time = Get_image_info.CB2_NAC_exposure
image_gain = Get_image_info.CB2_NAC_gain
image_files_length = len(image_files)
print(f"Total number of image files: {image_files_length}\n")

# !!! Check folder names
camera_filter = "CB2"
camera = "NAC" # !!! Note: always capitalize

pixel_angular_size = 5.9907e-6 #rad per pix

total_rejected_images = 0
total_processed_images = 0

processed_image_ids = []
processed_image_skyreduced_apature_count = []
processed_image_count_intenisty = []

for n in range(1, 11):
    print(f"IMG {n+1}: {image_files[n]}")
    
    if (image_gain[n] == '95'):
        satellite_distance = planet_cassini_distance[n] * 2
    elif (image_gain[n] == '29'):
        satellite_distance = planet_cassini_distance[n]
    else:
        print("!!! Found image with different gain number other than 29 or 95 !!!")
        break
    
    image_data = functs.get_image_data(image_files[n], camera_filter, camera)
    
    functs.get_titan_aperture(image_data, satellite_distance, pixel_angular_size)
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_method1_cumulativebound_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    functs.get_titan_aperture_2(image_data, satellite_distance, pixel_angular_size)
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_method2_maxbound_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    zscale_interval = ZScaleInterval(contrast=0.05)
    zscale_image_data = zscale_interval(image_data)
    
    snorm = simple_norm(zscale_image_data, 'log', log_a=10)
    norm_image = snorm(zscale_image_data)
    functs.draw_titan_aperture(norm_image, satellite_distance, pixel_angular_size)
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_method3_log10normalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    snorm = simple_norm(zscale_image_data, 'sqrt')
    norm_image = snorm(zscale_image_data)
    functs.draw_titan_aperture(norm_image, satellite_distance, pixel_angular_size)
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_method4_sqrtnormalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    snorm = simple_norm(zscale_image_data, 'power', power=2)
    norm_image = snorm(zscale_image_data)
    functs.draw_titan_aperture(norm_image, satellite_distance, pixel_angular_size)
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_method5_squarednormalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    snorm = simple_norm(zscale_image_data, 'log', log_a=100)
    norm_image = snorm(zscale_image_data)
    functs.draw_titan_aperture(norm_image, satellite_distance, pixel_angular_size)
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_method6_log100normalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    """
    try:
        aperture, annulus = functs.get_titan_aperture(image_files[n], camera_filter, camera, satellite_distance, pixel_angular_size)
    except:
        print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
        total_rejected_images += 1
    else:
        print("Titan found inside the CCD image.\n")
        #image_number, titan_aperture_sky_reduced, titan_intensity = functs.image_photometry(image_files[n], image_numbers[n], camera_filter, camera, satellite_distance, exposure_time[n], aperture, annulus)
        total_processed_images += 1
    """
        
print(f"Number of rejected images: {total_rejected_images}")
print(f"Number of processed images: {total_processed_images}")

