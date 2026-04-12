import numpy as np
from astropy.table import Table
from astropy.time import Time

# Get image names from CB2 Wide Camera
#---------------------------------------------
CB2_WAC_table = Table.read('Titan_images/CB2_WAC/CB2WACdata.csv', format='ascii')

CB2_WAC_image_files = np.array(CB2_WAC_table['Primary File Spec'])
CB2_WAC_gain = np.array(CB2_WAC_table['Gain Mode [Cassini ISS]'])
for n in range(len(CB2_WAC_image_files)):
    CB2_WAC_image_files[n] = CB2_WAC_image_files[n][-17:-4]
    CB2_WAC_gain[n] = CB2_WAC_gain[n][:2]
CB2_WAC_planet_cassini_distance = np.array(CB2_WAC_table['Body Center Distance (Min) [Titan] (km)'])
CB2_WAC_phase_angle = np.array(CB2_WAC_table['Phase Angle at Body Center (Min) [Titan] (degrees)'])
CB2_WAC_exposure = np.array(CB2_WAC_table['Exposure Duration (secs)'])
CB2_WAC_OBStime = np.array(CB2_WAC_table['Observation Start Time (YMDhms)'])
time = Time(CB2_WAC_OBStime, format='isot', scale='utc')
CB2_WAC_OBStime = time.mjd

# Get image names from CB2 Narrow Camera
#---------------------------------------------
CB2_NAC_table = Table.read('Titan_images/CB2_NAC/CB2NACdata.csv', format='ascii')

CB2_NAC_image_files = np.array(CB2_NAC_table['Primary File Spec'])
CB2_NAC_gain = np.array(CB2_NAC_table['Gain Mode [Cassini ISS]'])
for n in range(len(CB2_NAC_image_files)):
    CB2_NAC_image_files[n] = CB2_NAC_image_files[n][-17:-4]
    CB2_NAC_gain[n] = CB2_NAC_gain[n][:2]
CB2_NAC_planet_cassini_distance = np.array(CB2_NAC_table['Body Center Distance (Min) [Titan] (km)'])
CB2_NAC_phase_angle = np.array(CB2_NAC_table['Phase Angle at Body Center (Min) [Titan] (degrees)'])
CB2_NAC_exposure = np.array(CB2_NAC_table['Exposure Duration (secs)'])
CB2_NAC_OBStime = np.array(CB2_NAC_table['Observation Start Time (YMDhms)'])
time = Time(CB2_NAC_OBStime, format='isot', scale='utc')
CB2_NAC_OBStime = time.mjd

# Get image names from CB3 Wide Camera
#---------------------------------------------
CB3_WAC_table = Table.read('Titan_images/CB3_WAC/CB3WACdata.csv', format='ascii')

CB3_WAC_image_files = np.array(CB3_WAC_table['Primary File Spec'])
CB3_WAC_gain = np.array(CB3_WAC_table['Gain Mode [Cassini ISS]'])
for n in range(len(CB3_WAC_image_files)):
    CB3_WAC_image_files[n] = CB3_WAC_image_files[n][-17:-4]
    CB3_WAC_gain[n] = CB3_WAC_gain[n][:2]
CB3_WAC_planet_cassini_distance = np.array(CB3_WAC_table['Body Center Distance (Min) [Titan] (km)'])
CB3_WAC_phase_angle = np.array(CB3_WAC_table['Phase Angle at Body Center (Min) [Titan] (degrees)'])
CB3_WAC_exposure = np.array(CB3_WAC_table['Exposure Duration (secs)'])
CB3_WAC_OBStime = np.array(CB3_WAC_table['Observation Start Time (YMDhms)'])
time = Time(CB3_WAC_OBStime, format='isot', scale='utc')
CB3_WAC_OBStime = time.mjd

# Get image names from CB3 Narrow Camera
#---------------------------------------------

CB3_NAC_table = Table.read('Titan_images/CB3_NAC/CB3NACdata.csv', format='ascii')

CB3_NAC_image_files = np.array(CB3_NAC_table['Primary File Spec'])
CB3_NAC_gain = np.array(CB3_NAC_table['Gain Mode [Cassini ISS]'])
for n in range(len(CB3_NAC_image_files)):
    CB3_NAC_image_files[n] = CB3_NAC_image_files[n][-17:-4]
    CB3_NAC_gain[n] = CB3_NAC_gain[n][:2]
CB3_NAC_planet_cassini_distance = np.array(CB3_NAC_table['Body Center Distance (Min) [Titan] (km)'])
CB3_NAC_phase_angle = np.array(CB3_NAC_table['Phase Angle at Body Center (Min) [Titan] (degrees)'])
CB3_NAC_exposure = np.array(CB3_NAC_table['Exposure Duration (secs)'])
CB3_NAC_OBStime = np.array(CB3_NAC_table['Observation Start Time (YMDhms)'])
time = Time(CB3_NAC_OBStime, format='isot', scale='utc')
CB3_NAC_OBStime = time.mjd


