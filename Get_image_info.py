#!/usr/bin/python

import numpy as np
from astropy.table import Table

import functs

# Get image names from CB2 Wide Camera
#---------------------------------------------
"""
# !!! Check folder names
camera_filter = "CB2"
camera = "Wide" # !!! Note: always capitalize

CB2_wide_table = Table.read('Titan_images/CB2_Wide/data.csv', format='ascii')

CB2_wide_image_numbers = np.array(CB2_wide_table['Image Number [Cassini ISS]'])
CB2_wide_image_files = functs.proper_image_name(CB2_wide_image_numbers)
CB2_wide_image_size = np.array(CB2_wide_table['Greater Size in Pixels'])
CB2_wide_planet_cassini_distance = np.array(CB2_wide_table['Body Center Distance (Min) [Titan] (km)'])
CB2_wide_phase_angle = np.array(CB2_wide_table['Phase Angle at Body Center (Min) [Titan] (degrees)'])
CB2_wide_exposure = np.array(CB2_wide_table['Exposure Duration (secs)'])
CB2_wide_gain = np.array(CB2_wide_table['Gain Mode [Cassini ISS]'])
"""

# Get image names from CB2 Narrow Camera
#---------------------------------------------

CB2_narrow_table = Table.read('Titan_images/CB2_NAC/data.csv', format='ascii')

CB2_narrow_image_numbers = np.array(CB2_narrow_table['Image Number [Cassini ISS]'])
CB2_narrow_image_files = functs.proper_image_name(CB2_narrow_image_numbers, "NAC")
CB2_narrow_image_size = np.array(CB2_narrow_table['Greater Size in Pixels'])
CB2_narrow_planet_cassini_distance = np.array(CB2_narrow_table['Body Center Distance (Max) [Titan] (km)'])
CB2_narrow_phase_angle = np.array(CB2_narrow_table['Phase Angle at Body Center (Min) [Titan] (degrees)'])
CB2_narrow_exposure = np.array(CB2_narrow_table['Exposure Duration (secs)'])
CB2_narrow_gain = np.array(CB2_narrow_table['Gain Mode [Cassini ISS]'])

# Get image names from CB3 Wide Camera
#---------------------------------------------


# Get image names from CB3 Narrow Camera
#---------------------------------------------


