import numpy as np
import matplotlib.pyplot as plt
import csv

import functs

#image_file = "N1642241415_1_CALIB"
#satellite_distance = 1212494.476 * 2
#exposure_time = 5.6

#image_file = "N1577095185_1_CALIB"
#satellite_distance = 1283698.306
#exposure_time = 8.2

#image_file = "W1624421052_1_CALIB"
#satellite_distance = 177548.592
#exposure_time = 0.82

image_file = "W1683616251_1_CALIB"
satellite_distance = 150076.613
exposure_time = 0.82
gain = 29

# !!! Check folder names
camera_filter = "CB2"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix for Wide angle
#pixel_angular_size = 5.9907e-6 #rad per pix for Narrow angle

fields = ['Image file', 'Satellite distance [km]', 'Exposure time [s]', 'Gain [e per DN]',
          'Aperture count [DN]', 'Aperture area', 'Annulus radius', 'Sky median [DN]',
          'Aperture area sky count [DN]', 'Source count [DN]', 'Source intenisty [e per s per m2]']

with open("test_output.csv", 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)

image_data = functs.get_image_data(image_file, camera_filter, camera)

try:
    aperture, annulus = functs.get_titan_aperture(image_data, satellite_distance, pixel_angular_size)
    #plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_cumulative_aperture.png", dpi = 120, bbox_inches='tight')
    plt.show()
except:
    print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
    #plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_files[n]}_cumulative_aperture.png", dpi = 120, bbox_inches='tight')
    plt.show()
else:
    print("Titan found inside the CCD image.\n")
    photometric_data = functs.image_photometry(image_data, image_file, satellite_distance, exposure_time, gain, aperture, annulus)
    photometric_data = np.array([photometric_data])
    
    with open("test_output.csv", 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(photometric_data)
    
