import numpy as np
import matplotlib.pyplot as plt
import pdr
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry
from astropy.visualization import simple_norm
from astropy.visualization import ZScaleInterval


def proper_image_name(image_numbers, camera):
    # Input: an 1-line array of image numbers
    # Returns: an array of proper image file names to use 
    
    input_length = len(image_numbers)
    
    image_files_list = []
    
    if (camera == "WAC"):
        for n in range(input_length):
            file_name_format = "W{}_1_CALIB"
            image_name = file_name_format.format(image_numbers[n])
            image_files_list.append(image_name)
    elif (camera == "NAC"):
        for n in range(input_length):
            file_name_format = "N{}_1_CALIB"
            image_name = file_name_format.format(image_numbers[n])
            image_files_list.append(image_name)
    else:
        raise Exception("Either 'Wide' or 'Narrow' has to be specified")
        
    return np.array(image_files_list)


# For Save_IMG_as_FITS.py
#-------------------------------------------------------------

def pdr_to_fits_images(image_file, camera_filter, camera):
    
    IMG_name = image_file
    IMG_path = f"Titan_images/{camera_filter}_{camera}/{IMG_name}.IMG" 
    data = pdr.read(IMG_path)
    image_data = np.array(data['IMAGE'], dtype=float)
    
    convert_to_fits = fits.PrimaryHDU(data=image_data)
    convert_to_fits.writeto(f'Titan_images/{camera_filter}_{camera}/{IMG_name}.fits', overwrite=True)
        
    return print(f"{IMG_name}.fits created\n")


def plot_fits_images(image_file, camera_filter, camera):

    image_data = fits.getdata(f"Titan_images/{camera_filter}_{camera}/{image_file}.fits") 
        # If wrong {camera_filter} and {camera}, that folder is not present and an error will occure.
    value_adjust = (image_data >= 1.0) | (image_data <= 0) #reject pixel values under x and over x
    image_data[value_adjust] = 0
        
    plt.imshow(image_data, cmap='gray')
    plt.gca().invert_yaxis()
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/fits_images/{image_file}.png", dpi = 120, bbox_inches='tight')
    plt.show()
        
    return None

def plot_number_of_rejected_pixels(image_file, camera_filter, camera):
    
    image_data = fits.getdata(f"Titan_images/{camera_filter}_{camera}/{image_file}.fits") 
    
    upper_lim = np.arange(0.1, 3, 0.01)
    total_rejection = np.zeros_like(upper_lim)

    for i in range(len(upper_lim)):
        value_adjust = (image_data >= upper_lim[i]) | (image_data <= 0)
        rej_no = value_adjust.sum()
        total_rejection[i] = rej_no
        
    plt.plot(upper_lim, total_rejection)
    plt.title("No. of rejection pixel values in image above diff. upper limits")
    plt.xlabel("Upper limit")
    plt.ylabel("No. of rej. pixels")
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/fits_images/{image_file}_rejected_pix_plot.png", dpi = 120, bbox_inches='tight')
    plt.show() 
    
    return None

# Image_photometry
#-------------------------------------------------------------

def get_image_data(image_file, camera_filter, camera):
    image_data = fits.getdata(f"Titan_images/{camera_filter}_{camera}/" + image_file + ".fits")
    image_data = image_data[1:-1, 1:-1]
    rejected_values = (image_data >= 1.0) | (image_data <= 0) #reject pixel values below and above
    image_data[rejected_values] = 0
    
    return image_data

def get_titan_aperture(image_data, planet_cassini_distance, pixel_angular_size): #!!!
    
    total_image_value = image_data.sum()
    print(f"Total image value: {total_image_value}")
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    print(f"Square length of array: {image_pixel_length}\n")
    
    CCD_range = np.arange(0, image_pixel_length)
    
    cumulate_middle_xaxis = (np.cumsum(xaxis_marginal) <= 0.5)
    middle_bound_xaxis = cumulate_middle_xaxis.sum()
    #print(f"Middle cumulative bound from left xaxis: {middle_bound_xaxis}")
    find_max_xaxis_margin = (max(xaxis_marginal) == xaxis_marginal)
    max_xaxis_margin = CCD_range[find_max_xaxis_margin]
    max_xaxis_margin = max_xaxis_margin[0]
    xaxis_horizontal_limit = max(xaxis_marginal) / 6
    xaxis_left_bound_booleen = (xaxis_marginal[:max_xaxis_margin] <= xaxis_horizontal_limit)
    xaxis_left_bound = xaxis_left_bound_booleen.sum()
    xaxis_left_bound_to_maxvalue_seperation = abs(max_xaxis_margin - xaxis_left_bound)
    xaxis_right_bound_booleen = (xaxis_marginal[max_xaxis_margin:][::-1] <= xaxis_horizontal_limit)
    xaxis_right_bound = len(xaxis_marginal) - xaxis_right_bound_booleen.sum()
    xaxis_right_bound_to_maxvalue_seperation = abs(max_xaxis_margin - xaxis_right_bound)
    xaxis_center_shift = abs(xaxis_left_bound_to_maxvalue_seperation - xaxis_right_bound_to_maxvalue_seperation) / 2
    if (xaxis_right_bound_to_maxvalue_seperation >= xaxis_left_bound_to_maxvalue_seperation):
        improved_center_xaxis = middle_bound_xaxis + xaxis_center_shift
    else:
        improved_center_xaxis = middle_bound_xaxis - xaxis_center_shift
    
    cumulate_middle_yaxis = (np.cumsum(yaxis_marginal) <= 0.5)
    middle_bound_yaxis = cumulate_middle_yaxis.sum()
    #print(f"Middle cumulative bound from left yaxis: {middle_bound_yaxis}")
    find_max_yaxis_margin = (max(yaxis_marginal) == yaxis_marginal)
    max_yaxis_margin = CCD_range[find_max_yaxis_margin]
    max_yaxis_margin = max_yaxis_margin[0]
    yaxis_horizontal_limit = max(yaxis_marginal) / 6
    yaxis_left_bound_booleen = (yaxis_marginal[:max_yaxis_margin] <= yaxis_horizontal_limit)
    yaxis_left_bound = yaxis_left_bound_booleen.sum()
    yaxis_left_bound_to_maxvalue_seperation = abs(max_yaxis_margin - yaxis_left_bound)
    yaxis_right_bound_booleen = (yaxis_marginal[max_yaxis_margin:][::-1] <= yaxis_horizontal_limit)
    yaxis_right_bound = len(yaxis_marginal) - yaxis_right_bound_booleen.sum()
    yaxis_right_bound_to_maxvalue_seperation = abs(max_yaxis_margin - yaxis_right_bound)
    yaxis_center_shift = abs(yaxis_left_bound_to_maxvalue_seperation - yaxis_right_bound_to_maxvalue_seperation) / 2
    if (yaxis_right_bound_to_maxvalue_seperation >= yaxis_left_bound_to_maxvalue_seperation):
        improved_center_yaxis = middle_bound_yaxis + yaxis_center_shift
    else:
        improved_center_yaxis = middle_bound_yaxis - yaxis_center_shift
    
    planet_radius = 2975 #km
    planet_pixel_radius = np.arctan(planet_radius/planet_cassini_distance) / pixel_angular_size
    print(f"Titan raidus in pixels: {planet_pixel_radius}")
    
    if (planet_pixel_radius >= 100):
        x_apature = improved_center_xaxis
        y_apature = improved_center_yaxis
        
        x_apature_imporved = middle_bound_xaxis
        y_apature_imporved = middle_bound_yaxis
        
        extra_apature_radius = 20
        aperture_radius = planet_pixel_radius + extra_apature_radius
        print(f"Aperture raidus in pixels: {aperture_radius}\n")
    else:
        x_apature = middle_bound_xaxis
        y_apature = middle_bound_yaxis
        
        x_apature_imporved = middle_bound_xaxis
        y_apature_imporved = middle_bound_yaxis
        
        extra_apature_radius = 20
        aperture_radius = planet_pixel_radius + extra_apature_radius
        print(f"Aperture radius in pixels: {aperture_radius}\n")

    r_apature = aperture_radius
    
    aperture = CircularAperture((x_apature, y_apature), r=r_apature)
    
    inner_annulus = r_apature + 2
    outer_annulus = inner_annulus + 10
    annulus = CircularAnnulus((x_apature, y_apature), r_in=inner_annulus, r_out=outer_annulus)
    
    plot_image_aperture(image_data, x_apature, y_apature, r_apature, aperture, annulus, x_apature_imporved, y_apature_imporved, xaxis_horizontal_limit, yaxis_horizontal_limit)
    
    if ((x_apature + r_apature) >= image_pixel_length):
        raise Exception("Titan is outside image.")
    elif ((x_apature - r_apature) <= 0):
        raise Exception("Titan is outside image.")
        
    if ((y_apature + r_apature) >= image_pixel_length):
        raise Exception("Titan is outside image.")
    elif ((y_apature - r_apature) <= 0):
        raise Exception("Titan is outside image.")
    
    return (aperture, annulus)

def get_titan_aperture_CB3_WAC_second_try(image_data, planet_cassini_distance, pixel_angular_size): #!!!
    
    total_image_value = image_data.sum()
    print(f"Total image value: {total_image_value}")
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    print(f"Square length of array: {image_pixel_length}\n")
    
    CCD_range = np.arange(0, image_pixel_length)
    
    cumulate_middle_xaxis = (np.cumsum(xaxis_marginal) <= 0.5)
    middle_bound_xaxis = cumulate_middle_xaxis.sum()
    #print(f"Middle cumulative bound from left xaxis: {middle_bound_xaxis}")
    find_max_xaxis_margin = (max(xaxis_marginal) == xaxis_marginal)
    max_xaxis_margin = CCD_range[find_max_xaxis_margin]
    max_xaxis_margin = max_xaxis_margin[0]
    xaxis_horizontal_limit = max(xaxis_marginal) / 6
    xaxis_left_bound_booleen = (xaxis_marginal[:max_xaxis_margin] <= xaxis_horizontal_limit)
    xaxis_left_bound = xaxis_left_bound_booleen.sum()
    xaxis_left_bound_to_maxvalue_seperation = abs(max_xaxis_margin - xaxis_left_bound)
    xaxis_right_bound_booleen = (xaxis_marginal[max_xaxis_margin:][::-1] <= xaxis_horizontal_limit)
    xaxis_right_bound = len(xaxis_marginal) - xaxis_right_bound_booleen.sum()
    xaxis_right_bound_to_maxvalue_seperation = abs(max_xaxis_margin - xaxis_right_bound)
    xaxis_center_shift = abs(xaxis_left_bound_to_maxvalue_seperation - xaxis_right_bound_to_maxvalue_seperation) / 2
    if (xaxis_right_bound_to_maxvalue_seperation >= xaxis_left_bound_to_maxvalue_seperation):
        improved_center_xaxis = middle_bound_xaxis + xaxis_center_shift
    else:
        improved_center_xaxis = middle_bound_xaxis - xaxis_center_shift
    
    cumulate_middle_yaxis = (np.cumsum(yaxis_marginal) <= 0.5)
    middle_bound_yaxis = cumulate_middle_yaxis.sum()
    #print(f"Middle cumulative bound from left yaxis: {middle_bound_yaxis}")
    find_max_yaxis_margin = (max(yaxis_marginal) == yaxis_marginal)
    max_yaxis_margin = CCD_range[find_max_yaxis_margin]
    max_yaxis_margin = max_yaxis_margin[0]
    yaxis_horizontal_limit = max(yaxis_marginal) / 6
    yaxis_left_bound_booleen = (yaxis_marginal[:max_yaxis_margin] <= yaxis_horizontal_limit)
    yaxis_left_bound = yaxis_left_bound_booleen.sum()
    yaxis_left_bound_to_maxvalue_seperation = abs(max_yaxis_margin - yaxis_left_bound)
    yaxis_right_bound_booleen = (yaxis_marginal[max_yaxis_margin:][::-1] <= yaxis_horizontal_limit)
    yaxis_right_bound = len(yaxis_marginal) - yaxis_right_bound_booleen.sum()
    yaxis_right_bound_to_maxvalue_seperation = abs(max_yaxis_margin - yaxis_right_bound)
    yaxis_center_shift = abs(yaxis_left_bound_to_maxvalue_seperation - yaxis_right_bound_to_maxvalue_seperation) / 2
    if (yaxis_right_bound_to_maxvalue_seperation >= yaxis_left_bound_to_maxvalue_seperation):
        improved_center_yaxis = middle_bound_yaxis + yaxis_center_shift
    else:
        improved_center_yaxis = middle_bound_yaxis - yaxis_center_shift
    
    planet_radius = 2975 #km
    planet_pixel_radius = np.arctan(planet_radius/planet_cassini_distance) / pixel_angular_size
    print(f"Titan raidus in pixels: {planet_pixel_radius}")
    
    if (planet_pixel_radius >= 100):
        x_apature = middle_bound_xaxis
        y_apature = improved_center_yaxis
        
        x_apature_imporved = middle_bound_xaxis
        y_apature_imporved = middle_bound_yaxis
        
        extra_apature_radius = 10
        aperture_radius = planet_pixel_radius + extra_apature_radius
        print(f"Aperture raidus in pixels: {aperture_radius}\n")
    else:
        x_apature = middle_bound_xaxis
        y_apature = middle_bound_yaxis
        
        x_apature_imporved = middle_bound_xaxis
        y_apature_imporved = middle_bound_yaxis
        
        extra_apature_radius = 20
        aperture_radius = planet_pixel_radius + extra_apature_radius
        print(f"Aperture radius in pixels: {aperture_radius}\n")

    r_apature = aperture_radius
    
    aperture = CircularAperture((x_apature, y_apature), r=r_apature)
    
    inner_annulus = r_apature + 2
    outer_annulus = inner_annulus + 10
    annulus = CircularAnnulus((x_apature, y_apature), r_in=inner_annulus, r_out=outer_annulus)
    
    plot_image_aperture(image_data, x_apature, y_apature, r_apature, aperture, annulus, x_apature_imporved, y_apature_imporved, xaxis_horizontal_limit, yaxis_horizontal_limit)
    
    if ((x_apature + r_apature) >= image_pixel_length):
        raise Exception("Titan is outside image.")
    elif ((x_apature - r_apature) <= 0):
        raise Exception("Titan is outside image.")
        
    if ((y_apature + r_apature) >= image_pixel_length):
        raise Exception("Titan is outside image.")
    elif ((y_apature - r_apature) <= 0):
        raise Exception("Titan is outside image.")
    
    return (aperture, annulus)

def plot_image_aperture(image_data, x_apature, y_apature, r_apature, aperture, annulus, x_apature2, y_apature2, xaxis_horizontal, yaxis_horizontal):
    
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.imshow(image_data, cmap='gray')
    aperture.plot(color = "white", lw = 1.5, linestyle = "dashed")
    annulus.plot(color = "green", lw = 1.5)

    ax.set_xlabel("CCD's x-axis")
    ax.set_ylabel("CCD's y-axis")
    
    total_image_value = image_data.sum()
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
        
    divider = make_axes_locatable(ax)
    ax_histx = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)

    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    
    CCD_size = np.arange(0, image_pixel_length)

    ax_histx.plot(CCD_size, xaxis_marginal, color="tab:gray", alpha=0.2)
    ax_histx.set_xlim(0, len(xaxis_marginal))
    ax_histx.set_ylim(0)
    ax_histx.fill_between(CCD_size, 0, xaxis_marginal, color="tab:gray", alpha=0.5)
    ax_histx.axvline(x_apature, color="black", linestyle="dashed")
    ax_histx.axvline(x_apature - r_apature, color="tab:gray", linestyle="dashed")
    ax_histx.axvline(x_apature + r_apature, color="tab:gray", linestyle="dashed")
    ax_histx.axvline(x_apature2, color="blue")
    ax_histx.axvline(x_apature2 - r_apature, color="blue", linestyle="dashed")
    ax_histx.axvline(x_apature2 + r_apature, color="blue", linestyle="dashed")
    ax_histx.axhline(xaxis_horizontal, color="black", linestyle="dashed")

    ax_histy.plot(yaxis_marginal, CCD_size, color="tab:gray", alpha=0.2)
    ax_histy.set_xlim(0)
    ax_histy.set_ylim(0, len(yaxis_marginal))
    ax_histy.fill_between(yaxis_marginal, 0, CCD_size, color="tab:gray", alpha=0.5)
    ax_histy.axhline(y_apature, color="black", linestyle="dashed")
    ax_histy.axhline(y_apature - r_apature, color="tab:gray", linestyle="dashed")
    ax_histy.axhline(y_apature + r_apature, color="tab:gray", linestyle="dashed")
    ax_histy.axhline(y_apature2, color="blue")
    ax_histy.axhline(y_apature2 - r_apature, color="blue", linestyle="dashed")
    ax_histy.axhline(y_apature2 + r_apature, color="blue", linestyle="dashed")
    ax_histy.axvline(yaxis_horizontal, color="black", linestyle="dashed")
    
    return None

def draw_aperture_for_norm_image(image_data, method, planet_cassini_distance, pixel_angular_size): #!!!
    
    zscale_interval = ZScaleInterval(contrast=0.05)
    zscale_image_data = zscale_interval(image_data)
    
    if (method == 'log10'):
        snorm = simple_norm(zscale_image_data, 'log', log_a=10)
        zscaled_normalized_image = snorm(zscale_image_data) 
    elif (method == 'log100'):
        snorm = simple_norm(zscale_image_data, 'log', log_a=100)
        zscaled_normalized_image = snorm(zscale_image_data)
    elif (method == 'sqrt'):
        snorm = simple_norm(zscale_image_data, 'sqrt')
        zscaled_normalized_image = snorm(zscale_image_data)
    elif (method == 'squared'):
        snorm = simple_norm(zscale_image_data, 'power', power=2)
        zscaled_normalized_image = snorm(zscale_image_data)
        
    planet_radius = 2925 #km
    planet_pixel_radius = np.arctan(planet_radius/planet_cassini_distance) / pixel_angular_size
    print(f"Titan raidus in pixels: {planet_pixel_radius}")
       
    if (planet_pixel_radius >= 200):
        x_aperture2, y_aperture2, xaxis_horizontal_limit, yaxis_horizontal_limit, x_aperture, y_aperture = aperture_centers_from_max_pix(zscaled_normalized_image, 2)
        extra_apature_radius = 0
        r_aperture = planet_pixel_radius + extra_apature_radius
        print(f"Aperture raidus in pixels: {r_aperture}\n")
    else:
        x_aperture, y_aperture, xaxis_horizontal_limit, yaxis_horizontal_limit, x_aperture2, y_aperture2 = aperture_centers_from_max_pix(zscaled_normalized_image, 1.2)
        extra_apature_radius = 0
        r_aperture = planet_pixel_radius + extra_apature_radius
        print(f"Aperture radius in pixels: {r_aperture}\n")
    
    aperture = CircularAperture((x_aperture, y_aperture), r=r_aperture)
    
    inner_annulus = r_aperture + 2
    outer_annulus = inner_annulus + 10
    annulus = CircularAnnulus((x_aperture, y_aperture), r_in=inner_annulus, r_out=outer_annulus)
    
    plot_aperture_for_norm_image(image_data, method, x_aperture, y_aperture, r_aperture, aperture, annulus, x_aperture2, y_aperture2, xaxis_horizontal_limit, yaxis_horizontal_limit)
    
    """
    if ((x_aperture + r_aperture) >= image_pixel_length):
        raise Exception("Titan is outside image.")
    elif ((x_aperture - r_aperture) <= 0):
        raise Exception("Titan is outside image.")
        
    if ((y_aperture + r_aperture) >= image_pixel_length):
        raise Exception("Titan is outside image.")
    elif ((y_aperture - r_aperture) <= 0):
        raise Exception("Titan is outside image.")
    
    if ((x_aperture + outer_annulus) >= image_pixel_length):
        raise Exception("Titan outside image along the xaxis.")
    elif ((x_aperture - outer_annulus) <= 0):
        raise Exception("Titan outside image along the xaxis.")
        
    if ((y_aperture + outer_annulus) >= image_pixel_length):
        raise Exception("Titan outside image along the yaxis.")
    elif ((y_aperture - outer_annulus) <= 0):
        raise Exception("Titan outside image along the yaxis.")
    """
    
    return (aperture, annulus)

def aperture_centers_from_max_pix(image_data, horizontal_limit):
    
    total_image_value = image_data.sum()
    print(f"Total image value: {total_image_value}")
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    print(f"Square length of array: {image_pixel_length}\n")
    
    CCD_range = np.arange(0, image_pixel_length)
    
    find_max_xaxis_margin = (max(xaxis_marginal) == xaxis_marginal)
    max_xaxis_margin = CCD_range[find_max_xaxis_margin]
    max_xaxis_margin = max_xaxis_margin[0]
    print(f"xaxis center: {max_xaxis_margin}")
    xaxis_horizontal_limit = max(xaxis_marginal) / horizontal_limit
    xaxis_left_bound_booleen = (xaxis_marginal[:max_xaxis_margin] <= xaxis_horizontal_limit)
    xaxis_left_bound = xaxis_left_bound_booleen.sum()
    xaxis_left_bound_to_maxvalue_bound_distance = abs(max_xaxis_margin - xaxis_left_bound)
    xaxis_right_bound_booleen = (xaxis_marginal[max_xaxis_margin:][::-1] <= xaxis_horizontal_limit)
    xaxis_right_bound = len(xaxis_marginal) - xaxis_right_bound_booleen.sum()
    xaxis_right_bound_to_maxvalue_bound_distance = abs(max_xaxis_margin - xaxis_right_bound)
    xaxis_center_shift = abs(xaxis_left_bound_to_maxvalue_bound_distance - xaxis_right_bound_to_maxvalue_bound_distance) / 2
    if (xaxis_right_bound_to_maxvalue_bound_distance >= xaxis_left_bound_to_maxvalue_bound_distance):
        improved_center_xaxis = max_xaxis_margin + xaxis_center_shift
    else:
        improved_center_xaxis = max_xaxis_margin - xaxis_center_shift
    
    find_max_yaxis_margin = (max(yaxis_marginal) == yaxis_marginal)
    max_yaxis_margin = CCD_range[find_max_yaxis_margin]
    max_yaxis_margin = max_yaxis_margin[0]
    print(f"yaxis center: {max_yaxis_margin}")
    yaxis_horizontal_limit = max(yaxis_marginal) / horizontal_limit
    yaxis_left_bound_booleen = (yaxis_marginal[:max_yaxis_margin] <= yaxis_horizontal_limit)
    yaxis_left_bound = yaxis_left_bound_booleen.sum()
    yaxis_left_bound_to_maxvalue_bound_distance = abs(max_yaxis_margin - yaxis_left_bound)
    yaxis_right_bound_booleen = (yaxis_marginal[max_yaxis_margin:][::-1] <= yaxis_horizontal_limit)
    yaxis_right_bound = len(yaxis_marginal) - yaxis_right_bound_booleen.sum()
    yaxis_right_bound_to_maxvalue_bound_distance = abs(max_yaxis_margin - yaxis_right_bound)
    yaxis_center_shift = abs(yaxis_left_bound_to_maxvalue_bound_distance - yaxis_right_bound_to_maxvalue_bound_distance) / 2
    if (yaxis_right_bound_to_maxvalue_bound_distance >= yaxis_left_bound_to_maxvalue_bound_distance):
        improved_center_yaxis = max_yaxis_margin + yaxis_center_shift
    else:
        improved_center_yaxis = max_yaxis_margin - yaxis_center_shift
    
    return (max_xaxis_margin, max_yaxis_margin, xaxis_horizontal_limit, yaxis_horizontal_limit, improved_center_xaxis, improved_center_yaxis)

def plot_aperture_for_norm_image(image_data, method, x_aperture, y_aperture, r_aperture, aperture, annulus, x_aperture2, y_aperture2, xaxis_horizontal_limit, yaxis_horizontal_limit):
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    zscale_interval = ZScaleInterval(contrast=0.05)
    zscale_image_data = zscale_interval(image_data)
    
    if (method == 'log10'):
        snorm = simple_norm(zscale_image_data, 'log', log_a=10)
        zscaled_normalized_image = snorm(zscale_image_data) 
    elif (method == 'log100'):
        snorm = simple_norm(zscale_image_data, 'log', log_a=100)
        zscaled_normalized_image = snorm(zscale_image_data)
    elif (method == 'sqrt'):
        snorm = simple_norm(zscale_image_data, 'sqrt')
        zscaled_normalized_image = snorm(zscale_image_data)
    elif (method == 'squared'):
        snorm = simple_norm(zscale_image_data, 'power', power=2)
        zscaled_normalized_image = snorm(zscale_image_data)

    ax.imshow(zscaled_normalized_image, origin="lower", cmap="gray")
    aperture.plot(color = "white", lw = 1.5, linestyle = "dashed")
    #annulus.plot(color = "green", lw = 1.5)

    ax.set_xlabel("CCD's x-axis")
    ax.set_ylabel("CCD's y-axis")
    
    total_image_value = zscaled_normalized_image.sum()
    yaxis_marginal = zscaled_normalized_image.sum(axis=1)/total_image_value
    xaxis_marginal = zscaled_normalized_image.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
        
    divider = make_axes_locatable(ax)
    ax_histx = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)

    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    
    CCD_size = np.arange(0, image_pixel_length)

    ax_histx.plot(CCD_size, xaxis_marginal, color="tab:gray", alpha=0.2)
    ax_histx.set_xlim(0, len(xaxis_marginal))
    ax_histx.set_ylim(0)
    ax_histx.fill_between(CCD_size, 0, xaxis_marginal, color="tab:gray", alpha=0.5)
    ax_histx.axvline(x_aperture, color="black", linestyle="dashed")
    ax_histx.axvline(x_aperture + r_aperture, color="tab:gray", linestyle="dashed")
    ax_histx.axvline(x_aperture - r_aperture, color="tab:gray", linestyle="dashed")
    ax_histx.axhline(xaxis_horizontal_limit, color="black", linestyle="dashed")
    ax_histx.axvline(x_aperture2, color="blue")
    ax_histx.axvline(x_aperture2 + r_aperture, color="blue", linestyle="dashed")
    ax_histx.axvline(x_aperture2 - r_aperture, color="blue", linestyle="dashed")
    
    ax_histy.plot(yaxis_marginal, CCD_size, color="tab:gray", alpha=0.2)
    ax_histy.set_xlim(0)
    ax_histy.set_ylim(0, len(yaxis_marginal))
    ax_histy.fill_between(yaxis_marginal, 0, CCD_size, color="tab:gray", alpha=0.5)
    ax_histy.axhline(y_aperture, color="black", linestyle="dashed")
    ax_histy.axhline(y_aperture + r_aperture, color="tab:gray", linestyle="dashed")
    ax_histy.axhline(y_aperture - r_aperture, color="tab:gray", linestyle="dashed")
    ax_histy.axvline(yaxis_horizontal_limit, color="black", linestyle="dashed")
    ax_histy.axhline(y_aperture2, color="blue")
    ax_histy.axhline(y_aperture2 + r_aperture, color="blue", linestyle="dashed")
    ax_histy.axhline(y_aperture2 - r_aperture, color="blue", linestyle="dashed")
    
    return None

def image_photometry(image_data, planet_distance_in_km, exposure_time, gain, aperture, annulus):
    
    km_to_m = 1e3
    planet_distance = planet_distance_in_km * km_to_m
    
    aperstats = ApertureStats(image_data, annulus)
    annulus_area = aperstats.sum_aper_area
    sky_median = aperstats.median
    aperture_area = aperture.area_overlap(image_data)
    total_sky = sky_median * aperture_area
    
    image_photometry = aperture_photometry(image_data, aperture)
    image_photometry["aperture area"] = aperture_area
    image_photometry["annulus median sky value"] = sky_median
    image_photometry["total sky in aperture"] = total_sky
    aperture_count_skysub = image_photometry["aperture_sum"] - total_sky
    image_photometry["aperture_sum_skysub"] = aperture_count_skysub
    source_intenisty = image_photometry["aperture_sum_skysub"] * gain / (exposure_time * 4*np.pi * planet_distance**2)
    image_photometry["source intenisty"] = source_intenisty
    
    total_aperture_count = float(image_photometry["aperture_sum"][0])
    aperture_area = float(image_photometry["aperture area"][0])
    annulus_area = annulus.r_out - annulus.r_in
    sky_median = float(image_photometry["annulus median sky value"][0])
    aperture_sky_count = float(image_photometry["total sky in aperture"][0])
    source_count_skysub = float(image_photometry["aperture_sum_skysub"][0])
    intenisty = float(image_photometry["source intenisty"][0])
    
    photometry_data = [planet_distance_in_km, exposure_time, gain, 
                       total_aperture_count, aperture_area, annulus_area, sky_median,
                       aperture_sky_count, source_count_skysub, intenisty]
    
    return photometry_data
