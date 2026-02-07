import numpy as np
import matplotlib.pyplot as plt
import csv

import functs
import Get_image_info

image_files = Get_image_info.CB2_WAC_image_files
#image_numbers = Get_image_info.CB2_WAC_image_numbers
planet_cassini_distance = Get_image_info.CB2_WAC_planet_cassini_distance
exposure_time = Get_image_info.CB2_WAC_exposure
image_gain = Get_image_info.CB2_WAC_gain
phase_angle = Get_image_info.CB2_WAC_phase_angle
image_files_length = len(image_files)

camera_filter = "CB2"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix

#------------------------------------------------------------

total_rejected_images = 0
image_index_rejected = []
total_processed_images = 0
image_index_processed = []

num_of_images_with_different_gain = 0
images_with_different_gain = []

failed_image_call = 0
failed_image_call_file = []

fields = ['Image file', 'Phase angle [deg]', 'Aperture size', 'Source count [DN]', 'Planet-satellite distance [km]', 'Exposure time [s]', 'Gain', 'Intensity']

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
    else:
        print("Image has gain different than 29 or 95")   
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
        
    extra_aperture_radius = np.arange(-50, 55, 5)
    aperture_sizes = []
    aperture_counts = []
    sky_medians = []
    source_counts = []
    
    for i in range(len(extra_aperture_radius)):
        try:
            xaxis_shift_factor = 0.8
            yaxis_shift_factor = 0.8
            aperture, annulus = functs.get_titan_aperture(image_data, 'log100', satellite_distance, pixel_angular_size,  
                                                          xaxis_shift_factor, yaxis_shift_factor,
                                                          xaxis_aperture, yaxis_aperture, 
                                                          extra_aperture_radius[i], annulus_width)
        except Exception as e:
            print(f"{e}")
            plt.savefig(path + f"aperture_rejects/{image_files[n]}_aperture{extra_aperture_radius[i]}.png", dpi = 120, bbox_inches='tight')
            plt.show()
            break
        else:
            plt.savefig(path + f"apertures/{image_files[n]}_aperture{extra_aperture_radius[i]}.png", dpi = 120, bbox_inches='tight')
            plt.show()
            
            aperture_count, sky_median, source_count = functs.image_photometry(image_data, aperture, annulus)
            aperture_counts.append(aperture_count)
            sky_medians.append(sky_median)
            source_counts.append(source_count)
            aperture_sizes.append(extra_aperture_radius[i])
            
    if (len(aperture_counts) == 0):
        print('Reject image')
        total_rejected_images += 1
    else:
        plt.plot(aperture_sizes, aperture_counts, marker='o')
        plt.title('Aperture count')
        plt.savefig(path + f"apertures/{image_files[n]}_aperture_count.png", dpi = 120, bbox_inches='tight')
        plt.show()
        
        plt.plot(aperture_sizes, sky_medians, marker='o')
        plt.title('Sky median')
        plt.savefig(path + f"apertures/{image_files[n]}_sky_median.png", dpi = 120, bbox_inches='tight')
        plt.show()
    
        plt.plot(aperture_sizes, source_counts, marker='o')
        plt.title('Source counts')
        plt.savefig(path + f"apertures/{image_files[n]}_source_count.png", dpi = 120, bbox_inches='tight')
        plt.show()
         
        km_to_m = 1e3
        intensity = (source_counts[-1] * gain) / (exposure_time[n] * 4*np.pi * (planet_cassini_distance[n] * km_to_m)**2)
                
        photometric_data = [image_files[n]] + [phase_angle[n]] + [extra_aperture_radius[-1]] + [source_counts[-1]] + [planet_cassini_distance[n]] + [exposure_time[n]] + [gain] + [intensity]
        photometric_data = np.array([photometric_data])
        '''
        with open(path + f"{camera_filter}_{camera}_photometry_data.csv", "a") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerows(photometric_data)
        '''
        
        total_processed_images += 1
   
print(f"Total images: {image_files_length}")
print(f"Total rejected images: {total_rejected_images}")
print(f"Total processed images: {total_processed_images}")

'''
elif (len(aperture_counts) < 5) and (len(aperture_counts) > 0):
    plt.plot(aperture_sizes, sky_medians, marker='o')
    plt.title('Sky median')
    plt.savefig(path + f"apertures/{image_files[n]}_aperture_plot.png", dpi = 120, bbox_inches='tight')
    plt.show()

    plt.plot(aperture_sizes, source_counts, marker='o')
    plt.title('Source counts')
    plt.savefig(path + f"apertures/{image_files[n]}_aperture_plot.png", dpi = 120, bbox_inches='tight')
    plt.show()
         
    km_to_m = 1e3
    intensity = (source_counts[-1] * gain) / (exposure_time[n] * 4*np.pi * (planet_cassini_distance[n] * km_to_m)**2)
            
    photometric_data = [image_files[n]] + [phase_angle[n]] + [extra_aperture_radius[-1]] + [source_counts[-1]] + [planet_cassini_distance[n]] + [exposure_time[n]] + [gain] + [intensity]
    photometric_data = np.array([photometric_data])
            
    with open(path + f"{camera_filter}_{camera}_photometry_data.csv", "a") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerows(photometric_data)
        
    total_processed_images += 1
'''
