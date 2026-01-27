import numpy as np
import matplotlib.pyplot as plt

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

# !!! Check folder names
camera_filter = "CB2"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix for Wide angle
#pixel_angular_size = 5.9907e-6 #rad per pix for Narrow angle

total_rejected_images = 0
total_processed_images = 0

method = []
processed_image_total_apature_count = []

image_data = functs.get_image_data(image_file, camera_filter, camera)


print("Method 1: Finding 50% cumulative middle bound")
try:
    aperture, annulus = functs.get_titan_aperture(image_data, satellite_distance, pixel_angular_size)
except:
    print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
    plt.show()
    total_rejected_images += 1
else:
    print("Titan found inside the CCD image.\n")
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_file}_method1_cumulativebound_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    method.append('Finding 50%\ncumulative')
    total_aperture_count_sky_reduced = functs.image_photometry(image_data, satellite_distance, exposure_time, aperture, annulus)
    processed_image_total_apature_count.append(total_aperture_count_sky_reduced)
    total_processed_images += 1
        
print("Method 2: Finding max pixel bound")
try:
    aperture, annulus = functs.get_titan_aperture_2(image_data, satellite_distance, pixel_angular_size)
except:
    print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
    plt.show()
    total_rejected_images += 1
else:
    print("Titan found inside the CCD image.\n")
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_file}_method2_maxbound_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    method.append('Finding max\npixel margin')
    total_aperture_count_sky_reduced = functs.image_photometry(image_data, satellite_distance, exposure_time, aperture, annulus)
    processed_image_total_apature_count.append(total_aperture_count_sky_reduced)
    total_processed_images += 1
        
print("Method 3: zscaled-log normalized image")
try:
    aperture, annulus = functs.draw_titan_aperture(image_data, satellite_distance, pixel_angular_size, 'log10')
except:
    print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
    plt.show()
    total_rejected_images += 1
else:
    print("Titan found inside the CCD image.\n")
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_file}_method3_log10normalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    method.append('log10 norm')
    total_aperture_count_sky_reduced = functs.image_photometry(image_data, satellite_distance, exposure_time, aperture, annulus)
    processed_image_total_apature_count.append(total_aperture_count_sky_reduced)
    total_processed_images += 1
        
print("Method 4: zscaled-sqrt normalized image")
try:
    aperture, annulus = functs.draw_titan_aperture(image_data, satellite_distance, pixel_angular_size, 'sqrt')
except:
    print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
    plt.show()
    total_rejected_images += 1
else:
    print("Titan found inside the CCD image.\n")
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_file}_method4_sqrtnormalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    method.append('Sqrt norm')
    total_aperture_count_sky_reduced = functs.image_photometry(image_data, satellite_distance, exposure_time, aperture, annulus)
    processed_image_total_apature_count.append(total_aperture_count_sky_reduced)
    total_processed_images += 1
        
print("Method 5: zscaled-squared normalized image")
try:
    aperture, annulus = functs.draw_titan_aperture(image_data, satellite_distance, pixel_angular_size, 'squared')
except:
    print("Titan is outside of the CCD image. Cannot produce an aperture for this image and proceed to image photometry.\n")
    plt.show()
    total_rejected_images += 1
else:
    print("Titan found inside the CCD image.\n")
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/apertures/{image_file}_method5_squarednormalizedIMG_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    method.append('Square norm')
    total_aperture_count_sky_reduced = functs.image_photometry(image_data, satellite_distance, exposure_time, aperture, annulus)
    processed_image_total_apature_count.append(total_aperture_count_sky_reduced)
    total_processed_images += 1

print(f"Number of rejected images: {total_rejected_images}")
print(f"Number of processed images: {total_processed_images}")

plt.plot(method, processed_image_total_apature_count, marker='o')
plt.xticks(rotation=75)
plt.xlabel("Methods")
plt.ylabel("Aperture counts")
plt.show()
    
