import numpy as np
import matplotlib.pyplot as plt
import cv2
from astropy.io import fits

import functs

image_file = "W1557735069_1"
planet_cassini_distance = 241851.914
image_gain = "29"

camera_filter = "CB3"
camera = "WAC" # !!! Note: always capitalize

pixel_angular_size = 59.749e-6 #rad per pix
#pixel_angular_size = 5.9907e-6 #rad per pix

#------------------------------------------------------

def enhance_image_contrast(image_file, camera_filter, camera):

    image_data = fits.getdata(f"Titan_images/{camera_filter}_{camera}/{image_file}.fits") 
    image_data = image_data[2:-2, 2:-2]
    rejected_values = (image_data <= 0)
    image_data[rejected_values] = 0
    
    plt.imshow(image_data, cmap='gray')
    plt.show()
    
    image_data = image_data*1e3
    show_image = (image_data >= 150)
    image_data[show_image] = 150

    plt.imshow(image_data, cmap='gray')
    plt.show()
    plt.imsave(f"Titan_images/img_dump/{image_file}.png", image_data, cmap='gray')

    return None

def check_Canny_edge_detection(image_file, camera_filter, camera, p1, p2):

    image = cv2.imread(f"Titan_images/img_dump/{image_file}.png")
    
    image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    image_blur = cv2.medianBlur(image_gray, 5)
        
    n_rows = image_blur.shape[0]
    
    circles = cv2.HoughCircles(
        image_blur,
        cv2.HOUGH_GRADIENT,
        dp=1,
        minDist=n_rows,
        param1=p1,
        param2=p2,
        minRadius=5,
        maxRadius=0)
            
    if circles is not None:
        circles = np.uint16(np.around(circles))
        for i in circles[0, :]:
            center = (i[0], i[1])
            print(f'Center: {center}')
            # circle center; gray
            cv2.circle(image, center, 1, (0, 100, 100), 3)
            # circle outline; green
            radius = i[2]
            cv2.circle(image, center, radius, (0, 255, 0), 1)
        
        plt.imshow(image, cmap='gray')
        plt.show()
                
        #plt.imsave(f"Titan_images/img_dump/{image_file}.png", image)
            
    else:
        return None

    return circles

def check_Canny_edge_detection_parameters(image_file, camera_filter, camera):

    image = cv2.imread(f"Titan_images/img_dump/{image_file}.png")
    
    image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    image_blur = cv2.medianBlur(image_gray, 5)
        
    n_rows = image_blur.shape[0]
    
    p1 = np.arange(100, 1100, 100)
    p2 = np.arange(10, 110, 10)
    
    result = []
    
    for a in range(len(p1)):
        for b in range(len(p2)):
            #print(f'({p1[a]}, {p2[b]})')
            circles = cv2.HoughCircles(
                image_blur,
                cv2.HOUGH_GRADIENT,
                dp=1,
                minDist=n_rows,
                param1=p1[a],
                param2=p2[b],
                minRadius=5,
                maxRadius=0)
            
            if circles is not None:
                circles = np.uint16(np.around(circles))
                for i in circles[0, :]:
                    center = (i[0], i[1])
                    #print(f'Center: {center}')
                    # circle center; gray
                    cv2.circle(image, center, 1, (0, 100, 100), 3)
                    # circle outline; green
                    radius = i[2]
                    cv2.circle(image, center, radius, (0, 255, 0), 1)
                    
                #plt.imshow(image, cmap='gray')
                #plt.show()
                
                result.append(f'({p1[a]}, {p2[b]})')
                
                #plt.imsave(f"Titan_images/img_dump/{image_file}.png", image)
            
            else:
                #raise Exception("No centers can be found in the image.")
                continue
            
    return result

#--------------------------------------------------------

if (image_gain == '95'):
    gain = 95
    satellite_distance = planet_cassini_distance * 2
elif (image_gain == '29'):
    gain = 29
    satellite_distance = planet_cassini_distance

try:
    image_data = functs.get_image_data(image_file, camera_filter, camera)
except Exception as e:
    print(f"1. {e}\n")
    
enhance_image_contrast(image_file, camera_filter, camera)
CannyEdge = check_Canny_edge_detection(image_file, camera_filter, camera, 100, 10)

if CannyEdge is None:
    print('No centers could be found.')
    
else:
    aperture_center = CannyEdge[0, 0, :2]
    CannyEdge_radius = CannyEdge[0, 0, -1]
    
    extra_aperture_radius = np.arange(-50, 60, 5)
    annulus_width = 10
    
    aperture_sizes = []
    source_counts = []
    for i in range(len(extra_aperture_radius)):
        print(f'\nExtra aperture: {extra_aperture_radius[i]}')
        print(f'Canny Radius: {CannyEdge_radius}')
        try:
            aperture, annulus = functs.get_titan_aperture(image_data, 
                                                          satellite_distance, 
                                                          pixel_angular_size,
                                                          aperture_center, 
                                                          extra_aperture_radius[i], 
                                                          annulus_width)
            
        except Exception as e:
           print(f"3. {e}")
           plt.show()
           
           break
               
        else:
           plt.show()
           
           aperture_count, aperture_area, sky_median, sky_std, source_count = functs.image_photometry(image_data, aperture, annulus)
           aperture_sizes.append(extra_aperture_radius[i])
           source_counts.append(source_count)
   
    plt.plot(aperture_sizes, source_counts, marker='o')
    plt.show()

