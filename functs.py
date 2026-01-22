import numpy as np
import matplotlib.pyplot as plt
import pdr
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats, aperture_photometry

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
    plt.savefig(f"Titan_images/{camera_filter}_{camera}/{image_file}.png", dpi = 120, bbox_inches='tight')
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

def draw_titan_aperture(image_data, planet_cassini_distance, pixel_angular_size):
    
    total_image_value = image_data.sum()
    print(f"Total image value: {total_image_value}")
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    print(f"Square length of array: {image_pixel_length}\n")
    
    cumulate_middle_xaxis = (np.cumsum(xaxis_marginal) <= 0.5)
    middle_bound_xaxis = cumulate_middle_xaxis.sum()
    print(f"0.5 cumulative bound along xaxis columns: {middle_bound_xaxis}")
    
    cumulate_middle_yaxis = (np.cumsum(yaxis_marginal) <= 0.5)
    middle_bound_yaxis = cumulate_middle_yaxis.sum()
    print(f"0.5 cumulative bound along yaxis rows: {middle_bound_yaxis}")
    
    planet_radius = 2975 #2575km
    planet_pixel_radius = np.arctan(planet_radius/planet_cassini_distance) / pixel_angular_size
    print(f"Planet raidus in pixels: {planet_pixel_radius}\n")
    
    #------------------------------------------------------------------
    
    #aperture = CircularAperture((x_apature, y_apature), r=r_apature)
    
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.imshow(image_data, cmap='gray')
    #aperture.plot(color = "white", lw = 1.5,  linestyle = "dashed")

    ax.set_xlabel("CCD's x-axis")
    ax.set_ylabel("CCD's y-axis")
        
    divider = make_axes_locatable(ax)
    ax_histx = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)

    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    
    square_image_range = np.arange(0, image_pixel_length)

    ax_histx.plot(square_image_range, xaxis_marginal, color="tab:gray", alpha=0.2)
    #ax_histx.set_xlim(0, len(xaxis_marginal))
    #ax_histx.set_ylim(0)
    ax_histx.fill_between(square_image_range, 0, xaxis_marginal, color="tab:gray", alpha=0.5)
    ax_histx.axvline(middle_bound_xaxis, color="black", linestyle="dashed")
    ax_histx.axvline(middle_bound_xaxis + planet_pixel_radius, color="tab:blue", linestyle="dashed")
    ax_histx.axvline(middle_bound_xaxis - planet_pixel_radius, color="tab:blue", linestyle="dashed")
    
    ax_histy.plot(yaxis_marginal, square_image_range, color="tab:gray", alpha=0.2)
    #ax_histy.set_xlim(0)
    #ax_histy.set_ylim(0, len(yaxis_marginal))
    ax_histy.fill_between(yaxis_marginal, 0, square_image_range, color="tab:gray", alpha=0.5)
    ax_histy.axhline(middle_bound_yaxis, color="black", linestyle="dashed")
    ax_histy.axhline(middle_bound_yaxis + planet_pixel_radius, color="tab:blue", linestyle="dashed")
    ax_histy.axhline(middle_bound_yaxis - planet_pixel_radius, color="tab:blue", linestyle="dashed")
    
    plt.show()
    
    return None
    
    
    

def get_titan_aperture(image_file, image_number, camera_filter, camera, planet_cassini_distance, pixel_angular_size):
    
    image_data = get_image_data(image_file, camera_filter, camera)
    
    total_image_value = image_data.sum()
    print(f"Total image value: {total_image_value}")
    
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    print(f"Square length of array: {image_pixel_length}\n")
    
    cumulate_middle_yaxis = (np.cumsum(yaxis_marginal) <= 0.5)
    middle_bound_yaxis = cumulate_middle_yaxis.sum()
    print(f"Middle cumulative bound from left yaxis: {middle_bound_yaxis}")
    cumulate_quarter_yaxis_left = (np.cumsum(yaxis_marginal) <= 0.05)
    quarter_bound_yaxis_left = cumulate_quarter_yaxis_left.sum()
    quarter_middle_bound_yaxis_left_diff = abs(middle_bound_yaxis - quarter_bound_yaxis_left)
    print(f"Quarter cumulative bound from left yaxis: {quarter_bound_yaxis_left}")
    cumulate_quarter_yaxis_right = (np.cumsum(yaxis_marginal[::-1]) <= 0.05)
    quarter_bound_yaxis_right = len(yaxis_marginal) - cumulate_quarter_yaxis_right.sum()
    quarter_middle_bound_yaxis_right_diff = abs(middle_bound_yaxis - quarter_bound_yaxis_right)
    print(f"Quarter cumulative from right yaxis: {quarter_bound_yaxis_right}\n")
    middle_bound_yaxis_shift = abs(quarter_middle_bound_yaxis_right_diff - quarter_middle_bound_yaxis_left_diff)
    
    if (quarter_middle_bound_yaxis_right_diff >= quarter_middle_bound_yaxis_left_diff):
        middle_bound_yaxis_improved  = middle_bound_yaxis + middle_bound_yaxis_shift
    else:
        middle_bound_yaxis_improved = middle_bound_yaxis - middle_bound_yaxis_shift
    
    cumulate_middle_xaxis = (np.cumsum(xaxis_marginal) <= 0.5)
    middle_bound_xaxis = cumulate_middle_xaxis.sum()
    print(f"Middle cumulative bound from left xaxis: {middle_bound_xaxis}")
    cumulate_quarter_xaxis_left = (np.cumsum(xaxis_marginal) <= 0.05)
    quarter_bound_xaxis_left = cumulate_quarter_xaxis_left.sum()
    quarter_middle_bound_xaxis_left_diff = abs(middle_bound_xaxis - quarter_bound_xaxis_left)
    print(f"Quarter cumulative bound from left xaxis: {quarter_bound_xaxis_left}")
    cumulate_quarter_xaxis_right = (np.cumsum(xaxis_marginal[::-1]) <= 0.05)
    quarter_bound_xaxis_right = len(xaxis_marginal) - cumulate_quarter_xaxis_right.sum()
    quarter_middle_bound_xaxis_right_diff = abs(middle_bound_xaxis - quarter_bound_xaxis_right)
    print(f"Quarter cumulative from right xaxis: {quarter_bound_xaxis_right}\n")
    middle_bound_xaxis_shift = abs(quarter_middle_bound_xaxis_right_diff - quarter_middle_bound_xaxis_left_diff)
    
    if (quarter_middle_bound_xaxis_right_diff >= quarter_middle_bound_xaxis_left_diff):
        middle_bound_xaxis_improved  = middle_bound_xaxis + middle_bound_xaxis_shift
    else:
        middle_bound_xaxis_improved = middle_bound_xaxis - middle_bound_xaxis_shift
    
    planet_radius = 2575 #km
    planet_pixel_radius = np.arctan(planet_radius/planet_cassini_distance) / pixel_angular_size
    print(f"Planet raidus in pixels: {planet_pixel_radius}\n")
    
    if (planet_pixel_radius >= 50):
        titan_center_xcoord = middle_bound_xaxis_improved
        titan_center_ycoord = middle_bound_yaxis_improved
        extra_apature_radius = 60 #pix
        titan_apature_radius = planet_pixel_radius + extra_apature_radius
    else:
        titan_center_xcoord = middle_bound_xaxis
        titan_center_ycoord = middle_bound_yaxis
        extra_apature_radius = 40 #pix
        titan_apature_radius = planet_pixel_radius + extra_apature_radius
    
    x_apature = titan_center_xcoord
    y_apature = titan_center_ycoord
    r_apature = titan_apature_radius
    
    aperture = CircularAperture((x_apature, y_apature), r=r_apature)
    annulus = CircularAnnulus((x_apature, y_apature), r_in = (r_apature + 50), r_out = (r_apature + 100))
    
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.imshow(image_data, cmap='gray')
    aperture.plot(color = "white", lw = 1.5,  linestyle = "dashed")
    annulus.plot(color = "tab:blue", lw = 1.5)

    ax.set_xlabel("CCD's x-axis")
    ax.set_ylabel("CCD's y-axis")
    
    total_image_value = image_data.sum()
    yaxis_marginal = image_data.sum(axis=1)/total_image_value
    xaxis_marginal = image_data.sum(axis=0)/total_image_value
    image_pixel_length = len(yaxis_marginal)
    
    if (len(yaxis_marginal) >= 1000):
        ax.text(750, 980, f"W{image_number}", color="white")
    else:
        ax.text(350, 480, f"W{image_number}", color="white")
        
    divider = make_axes_locatable(ax)
    ax_histx = divider.append_axes("top", 1.2, pad=0.1, sharex=ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)

    ax_histx.xaxis.set_tick_params(labelbottom=False)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    
    CCD_range = np.arange(0, image_pixel_length)

    ax_histx.plot(CCD_range, xaxis_marginal, color="tab:gray", alpha=0.2)
    ax_histx.set_xlim(0, len(xaxis_marginal))
    ax_histx.set_ylim(0)
    ax_histx.fill_between(CCD_range, 0, xaxis_marginal, color="tab:gray", alpha=0.5)
    ax_histx.axvline(middle_bound_xaxis, color="tab:blue", linestyle="dashed")
    ax_histx.axvline(quarter_bound_xaxis_left, color="tab:blue")
    ax_histx.axvline(quarter_bound_xaxis_right, color="tab:blue")
    ax_histx.axvline(middle_bound_xaxis - planet_pixel_radius, color="tab:green", linestyle="dashed")
    ax_histx.axvline(middle_bound_xaxis + planet_pixel_radius, color="tab:green", linestyle="dashed")
    ax_histx.axvline(x_apature, color="black", linestyle="dashed")
    ax_histx.axvline(x_apature - r_apature, color="tab:gray", linestyle="dashed")
    ax_histx.axvline(x_apature + r_apature, color="tab:gray", linestyle="dashed")

    ax_histy.plot(yaxis_marginal, CCD_range, color="tab:gray", alpha=0.2)
    ax_histy.set_xlim(0)
    ax_histy.set_ylim(0, len(yaxis_marginal))
    ax_histy.fill_between(yaxis_marginal, 0, CCD_range, color="tab:gray", alpha=0.5)
    ax_histy.axhline(middle_bound_yaxis, color="tab:blue", linestyle="dashed")
    ax_histy.axhline(quarter_bound_yaxis_left, color="tab:blue")
    ax_histy.axhline(quarter_bound_yaxis_right, color="tab:blue")
    ax_histy.axhline(middle_bound_yaxis - planet_pixel_radius, color="tab:green", linestyle="dashed")
    ax_histy.axhline(middle_bound_yaxis + planet_pixel_radius, color="tab:green", linestyle="dashed")
    ax_histy.axhline(y_apature, color="black", linestyle="dashed")
    ax_histy.axhline(y_apature - r_apature, color="tab:gray", linestyle="dashed")
    ax_histy.axhline(y_apature + r_apature, color="tab:gray", linestyle="dashed")
    
    #plt.savefig(f"Titan_images/CB2_Wide/apatures/{image_file}_annulus.png", dpi = 120, bbox_inches='tight')
    plt.show()
    
    return (aperture, annulus)

def image_photometry(image_file, image_number, planet_distance_in_km, exposure_time, aperture, annulus_aperture):
    
    image_data = get_image_data(image_file)
    
    km_to_m = 1e3
    planet_distance = planet_distance_in_km * km_to_m
    
    aperstats = ApertureStats(image_data, annulus_aperture)
    sky_median = aperstats.median
    aperture_area = aperture.area_overlap(image_data)
    total_sky = sky_median * aperture_area
    
    print(f"The aperture area: {aperture_area} pix^2")
    
    image_photometry = aperture_photometry(image_data, aperture)
    image_photometry["id"] = image_number
    image_photometry["total_sky"] = total_sky
    titan_aperture_sky_reduced = image_photometry["aperture_sum"] - total_sky
    image_photometry["aperture_sum_skysub"] = titan_aperture_sky_reduced
    titan_intensity = image_photometry["aperture_sum_skysub"] / (exposure_time * 4*np.pi * planet_distance**2)
    image_photometry["intenisty"] = titan_intensity
    #for col in image_photometry.colnames:
    #    image_photometry[col].info.format = '%.8g'
        
    print(image_photometry)
    
    with open(f"Titan_images/CB2_Wide/apatures/W{image_number}_photometry.txt", "w") as out_file:
        line = str(image_photometry)
        out_file.write(line)
        
    return (image_number, titan_aperture_sky_reduced, titan_intensity)    
    
