import numpy as np
import matplotlib.pyplot as plt

import functs


image_files = "W1683616251_1_CALIB"
planet_cassini_distance = 150076.613
image_gain = "29"
exposure_time = 0.82

camera_filter = "CB2"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix
#pixel_angular_size = 5.9907e-6 #rad per pix


xaxis_aperture = "improved"
yaxis_aperture = "improved"
annulus_width = 15

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
if (image_size == 254):
    satellite_distance = planet_cassini_distance * 4

range_of_aperture_sizes = [0] #np.arange(-50, 50, 5)
aperture_sizes = []
aperture_counts = []
sky_medians = []
source_counts = []

for i in range(len(range_of_aperture_sizes)):
    try:
        xaxis_shift_factor = 0.6
        yaxis_shift_factor = 1
        aperture, annulus = functs.get_titan_aperture(image_data, 'None', satellite_distance, pixel_angular_size, 
                                                      xaxis_shift_factor, yaxis_shift_factor,
                                                      xaxis_aperture, yaxis_aperture, 
                                                      0, annulus_width)
        
        print(f'Aperture size: {range_of_aperture_sizes[i]}')
        
    except Exception as e:
        print(f"{e}")
        plt.show()
        break
        
    else:
        plt.show()
        
        aperture_count, sky_median, source_count = functs.image_photometry(image_data, aperture, annulus)
        aperture_counts.append(aperture_count)
        sky_medians.append(sky_median)
        source_counts.append(source_count)
        aperture_sizes.append(0)
        
        print(f'Source count w/out sky: {source_count}\n')
        
        
plt.plot(aperture_sizes, aperture_counts, marker='o')
plt.title('Aperture counts')
plt.show()

plt.plot(aperture_sizes, sky_medians, marker='o')
plt.title('Sky median')
plt.show()

plt.plot(aperture_sizes, source_counts, marker='o')
plt.title('Source counts')
plt.show()
