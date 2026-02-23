import numpy as np
import matplotlib.pyplot as plt
import pdr
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry
from astropy.visualization import simple_norm
from astropy.visualization import ZScaleInterval
from astropy.stats import SigmaClip


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
    IMG_path = f"Titan_images/{camera_filter}_{camera}/{IMG_name}_CALIB.IMG"
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

def get_titan_aperture(image_data, planet_cassini_distance, pixel_angular_size, xaxis_shift_factor, yaxis_shift_factor, extra_apature_radius, extra_annulus_radius, xaxis_aperture_center='improved_xaxis_center', yaxis_aperture_center='improved_yaxis_center'): #!!!
    
    image_size = image_data.shape[0]
    
    cumulated_center_xaxis, cumulated_center_yaxis, improved_center_xaxis, improved_center_yaxis = find_aperture_center(image_data, xaxis_shift_factor, yaxis_shift_factor)
    
    planet_radius = 2875 #km
    planet_pixel_radius = CCD_Target_radius(planet_radius, planet_cassini_distance, pixel_angular_size)
    print(f"Titan raidus in pixels: {planet_pixel_radius}")
    
    if (xaxis_aperture_center == 'improved_xaxis_center'):
        x_apature = improved_center_xaxis
    else:
        x_apature = cumulated_center_xaxis
        
    if (yaxis_aperture_center == 'improved_yaxis_center'):
        y_apature = improved_center_yaxis
    else:
        y_apature = cumulated_center_yaxis
    
    if (planet_pixel_radius <= 100):
       x_apature = cumulated_center_xaxis
       y_apature = cumulated_center_yaxis
        
       extra_apature_radius = 10
       aperture_radius = planet_pixel_radius + extra_apature_radius
       
    aperture_radius = planet_pixel_radius + extra_apature_radius
    r_apature = aperture_radius
    
    aperture = CircularAperture((x_apature, y_apature), r=r_apature)
    
    inner_annulus = r_apature + 15
    outer_annulus = inner_annulus + extra_annulus_radius
    annulus = CircularAnnulus((x_apature, y_apature), r_in=inner_annulus, r_out=outer_annulus)
    
    plot_image_aperture(image_data, x_apature, y_apature, r_apature, aperture, annulus, cumulated_center_xaxis, cumulated_center_yaxis)
    
    if ((x_apature + r_apature) >= image_size):
        raise Exception("Titan is outside image.")
    elif ((x_apature - r_apature) <= 0):
        raise Exception("Titan is outside image.")
        
    if ((y_apature + r_apature) >= image_size):
        raise Exception("Titan is outside image.")
    elif ((y_apature - r_apature) <= 0):
        raise Exception("Titan is outside image.")
    
    return (aperture, annulus)

def CCD_Target_radius(planet_radius, planet_cassini_distance, pixel_angular_size):
    return (np.arctan(planet_radius/planet_cassini_distance) / pixel_angular_size)


def find_aperture_center(image_data, xaxis_shift_factor, yaxis_shift_factor, horizontal_limit=6):

    total_image_value = image_data.sum()
    #print(f"Total image value: {total_image_value}")
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    #print(f"Square length of array: {image_pixel_length}\n")
    
    CCD_range = np.arange(0, image_pixel_length)
    
    half_cumulated_center_xaxis = (np.cumsum(xaxis_marginal) <= 0.5)
    cumulated_center_xaxis = half_cumulated_center_xaxis.sum()
    #print(f"Middle cumulative bound from left xaxis: {middle_bound_xaxis}")
    find_max_xaxis_margin = (max(xaxis_marginal) == xaxis_marginal)
    max_xaxis_margin = CCD_range[find_max_xaxis_margin]
    max_xaxis_margin = max_xaxis_margin[0]
    xaxis_horizontal_limit = max(xaxis_marginal) / horizontal_limit
    xaxis_left_bound_booleen = (xaxis_marginal[:max_xaxis_margin] <= xaxis_horizontal_limit)
    xaxis_left_bound = xaxis_left_bound_booleen.sum()
    xaxis_left_bound_to_maxvalue_seperation = abs(cumulated_center_xaxis - xaxis_left_bound)
    xaxis_right_bound_booleen = (xaxis_marginal[max_xaxis_margin:][::-1] <= xaxis_horizontal_limit)
    xaxis_right_bound = len(xaxis_marginal) - xaxis_right_bound_booleen.sum()
    xaxis_right_bound_to_maxvalue_seperation = abs(cumulated_center_xaxis - xaxis_right_bound)
    xaxis_center_shift = abs(xaxis_left_bound_to_maxvalue_seperation - xaxis_right_bound_to_maxvalue_seperation) / xaxis_shift_factor
    if (xaxis_right_bound_to_maxvalue_seperation >= xaxis_left_bound_to_maxvalue_seperation):
        improved_center_xaxis = cumulated_center_xaxis + xaxis_center_shift
    else:
        improved_center_xaxis = cumulated_center_xaxis - xaxis_center_shift
    
    half_cumulated_center_yaxis = (np.cumsum(yaxis_marginal) <= 0.5)
    cumulated_center_yaxis = half_cumulated_center_yaxis.sum()
    #print(f"Middle cumulative bound from left yaxis: {middle_bound_yaxis}")
    find_max_yaxis_margin = (max(yaxis_marginal) == yaxis_marginal)
    max_yaxis_margin = CCD_range[find_max_yaxis_margin]
    max_yaxis_margin = max_yaxis_margin[0]
    yaxis_horizontal_limit = max(yaxis_marginal) / horizontal_limit
    yaxis_left_bound_booleen = (yaxis_marginal[:max_yaxis_margin] <= yaxis_horizontal_limit)
    yaxis_left_bound = yaxis_left_bound_booleen.sum()
    yaxis_left_bound_to_maxvalue_seperation = abs(cumulated_center_yaxis - yaxis_left_bound)
    yaxis_right_bound_booleen = (yaxis_marginal[max_yaxis_margin:][::-1] <= yaxis_horizontal_limit)
    yaxis_right_bound = len(yaxis_marginal) - yaxis_right_bound_booleen.sum()
    yaxis_right_bound_to_maxvalue_seperation = abs(cumulated_center_yaxis - yaxis_right_bound)
    yaxis_center_shift = abs(yaxis_left_bound_to_maxvalue_seperation - yaxis_right_bound_to_maxvalue_seperation) / yaxis_shift_factor
    if (yaxis_right_bound_to_maxvalue_seperation >= yaxis_left_bound_to_maxvalue_seperation):
        improved_center_yaxis = cumulated_center_yaxis + yaxis_center_shift
    else:
        improved_center_yaxis = cumulated_center_yaxis - yaxis_center_shift

    return (cumulated_center_xaxis, cumulated_center_yaxis, improved_center_xaxis, improved_center_yaxis)

def plot_image_aperture(image_data, x_apature, y_apature, r_apature, aperture, annulus, x_apature2, y_apature2, horizontal_limit=6):
    
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
    
    xaxis_horizontal = max(xaxis_marginal) / horizontal_limit
    yaxis_horizontal = max(yaxis_marginal) / horizontal_limit
        
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

def image_photometry(image_data, aperture, annulus):
    
    aperture_stats = ApertureStats(image_data, aperture, sigma_clip=None)
    aperture_area = aperture_stats.sum_aper_area.value
    aperture_count = aperture_stats.sum
    
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    annulus_stats = ApertureStats(image_data, annulus, sigma_clip=sigclip)
    sky_median = annulus_stats.median
    sky_std = annulus_stats.std
    total_sky = sky_median * aperture_area
    
    source_count = aperture_count - total_sky
    
    return (aperture_count, aperture_area, sky_median, sky_std, source_count)
