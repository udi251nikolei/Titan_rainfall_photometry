import numpy as np
import matplotlib.pyplot as plt

import functs


image_files = "W1683616251_1_CALIB"
planet_cassini_distance = 150076.613
image_gain = "29"
exposure_time = 0.82

'''
image_files = "W1646541673_1_CALIB"
planet_cassini_distance = 1960168.965
image_gain = "29"
exposure_time = 15
'''

camera_filter = "CB2"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix
#pixel_angular_size = 5.9907e-6 #rad per pix


xaxis_aperture = "improved"
yaxis_aperture = "improved"
annulus_width = 20

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

range_of_aperture_sizes = np.arange(-20, 65, 5)
aperture_sizes = []
aperture_counts = []
sky_medians = []
source_counts = []
SNR_list = []

for i in range(len(range_of_aperture_sizes)):
    try:
        xaxis_shift_factor = 0.6
        yaxis_shift_factor = 1.2
        aperture, annulus = functs.get_titan_aperture(image_data, 'None', satellite_distance, pixel_angular_size, 
                                                          xaxis_shift_factor, yaxis_shift_factor,
                                                          xaxis_aperture, yaxis_aperture, 
                                                          range_of_aperture_sizes[i], annulus_width)
    except Exception as e:
        print(f"{e}")
        plt.show()
            
    else:
        plt.show()
        aperture_count, aperture_area, sky_median, annulus_area, source_count, SNR, sky_std = functs.image_photometry(image_data, aperture, annulus)
        
        aperture_counts.append(aperture_count)
        sky_medians.append(sky_median)
        source_counts.append(source_count)
        aperture_sizes.append(range_of_aperture_sizes[i])
        
        km_to_m = 1e3
        #intensity = (source_count * gain) / (exposure_time * 4*np.pi * (planet_cassini_distance * km_to_m)**2)
        
        print(f'Aperture area: {aperture_area}')
        print(f'Annulus area: {annulus_area}')
        print(f'Sky: {sky_median}')
        print(f'Sky std: {sky_std}')
        print(f'SNR: {SNR}')
        print(f'Source count: {source_count}\n')

plt.plot(aperture_sizes, aperture_counts, marker='o')
plt.title('Aperture counts')
plt.show()

plt.plot(aperture_sizes, sky_medians, marker='o')
plt.title('Sky median')
plt.show()

plt.plot(aperture_sizes, source_counts, marker='o')
plt.title('Source counts')
plt.show()
        
        
xaxis_FWHM_left, xaxis_FWHM_right, yaxis_FWHM_left, yaxis_FWHM_right = functs.FWHM(image_data)
print(f'xaxis FWHM left seperation: {xaxis_FWHM_left}')
print(f'xaxis FWHM right seperation: {xaxis_FWHM_right}')
print(f'avg xaxis FWHM radius: {(xaxis_FWHM_left + xaxis_FWHM_right)/2}')
print(f'Optimal xaxis aperture radius: {(xaxis_FWHM_left + xaxis_FWHM_right)}')
print(f'yaxis FWHM left seperation: {yaxis_FWHM_left}')
print(f'yaxis FWHM right seperation: {yaxis_FWHM_right}')
print(f'avg yaxis FWHM radius: {(yaxis_FWHM_left + yaxis_FWHM_right)/2}')
print(f'Optimal yaxis aperture radius: {(yaxis_FWHM_left + yaxis_FWHM_right)}')
