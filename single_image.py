import numpy as np
import matplotlib.pyplot as plt

import functs


image_files = "W1575527942_1_CALIB"
planet_cassini_distance = 127284.345
image_gain = "95"
exposure_time = 15


camera_filter = "CB3"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix
#pixel_angular_size = 5.9907e-6 #rad per pix

xaxis_aperture = "improved"
yaxis_aperture = "improved"
annulus_width = 30

#--------------------------------------------------------

if (image_gain == '95'):
    gain = 95
    satellite_distance = planet_cassini_distance * 2
elif (image_gain == '29'):
    gain = 29
    satellite_distance = planet_cassini_distance
else:
    print("Found image with different gain number other than 29 or 95")   
        
    gain = int(image_gain)
    satellite_distance = planet_cassini_distance

try:
    image_data = functs.get_image_data(image_files, camera_filter, camera)
except Exception as e:
    print(f"{e}\n")

image_size = image_data.shape[0]

planet_radius = 2875 #km
moon_pixel_radius = functs.CCD_body_radius(planet_radius, planet_cassini_distance, pixel_angular_size)

try:
    xaxis_shift_factor = 0.8
    yaxis_shift_factor = 0.8
    
    if (image_size > 1000):
        extra_aperture_radius = 60
    else:
        extra_aperture_radius = 40
    
    aperture, annulus = functs.get_titan_aperture(image_data, 'None', satellite_distance, pixel_angular_size, 
                                                          xaxis_shift_factor, yaxis_shift_factor,
                                                          xaxis_aperture, yaxis_aperture, 
                                                          extra_aperture_radius, annulus_width)
except Exception as e:
    print(f"{e}")
    plt.show()
            
else:
    plt.show()
    aperture_count, aperture_count_error, aperture_area, sky_mean, sky_median, annulus_area, source_count, source_count_error, SNR, sky_std = functs.image_photometry(image_data, aperture, annulus)

        
print(f'Aperture count: {aperture_count}')
print(f'Aperture count error: {aperture_count_error}')
print(f'Aperture area: {aperture_area}')
print(f'Annulus area: {annulus_area}')
print(f'Area ratio: {aperture_area/annulus_area}\n')
print(f'Sky mean: {sky_mean}')
print(f'Sky median: {sky_median}')
print(f'Sky std: {sky_std}')

tmp = (sky_mean - sky_median) / sky_std
if (tmp <= 0.3):
    print(tmp)
    mode = 2.5*sky_median - 1.5*sky_mean
else:
    print(tmp)
    mode = sky_median
    
print(f'Mode: {mode}')

print(f'SNR: {SNR}')
print(f'Source count: {source_count}')
print(f'Source error: {source_count_error}\n')


total_e_count = (source_count * gain)

km_to_m = 1e3
intensity = (source_count * gain) / (exposure_time * 4*np.pi * (planet_cassini_distance * km_to_m)**2)

print(f'Intenisty [e- per s per m^2]: {intensity}')
