import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import ListedColormap
import csv

import functs

CB2_NAC_photometry = Table.read('Titan_images/CB2_NAC/CB2_NAC_photometry_data.csv', format='ascii')

CB2_NAC_images = np.array(CB2_NAC_photometry['Image file'])
CB2_NAC_time = np.array(CB2_NAC_photometry['Observation time [MJD]'])
CB2_NAC_phase_angle = np.array(CB2_NAC_photometry['Phase angle [deg]'])
CB2_NAC_distance = np.array(CB2_NAC_photometry['Planet-satellite distance [km]'])
CB2_NAC_source_count = np.array(CB2_NAC_photometry['Source count'])
CB2_NAC_exposure = np.array(CB2_NAC_photometry['Exposure time [s]'])
CB2_NAC_gain = np.array(CB2_NAC_photometry['Gain'])
#CB2_NAC_target_radius = np.array(CB2_NAC_photometry['Target Radius [pix]'])
print(f'Number of images: {len(CB2_NAC_images)}')

'''
a = 1122
b = 1125
image_1 = CB2_NAC_images[a]
image_2 = CB2_NAC_images[b]
print(f'Image1: {image_1}')
print(f'Image2: {image_2}')
d1 = CB2_NAC_distance[a]
d2 = CB2_NAC_distance[b]
print(f'D1: {d1}')
print(f'D2: {d2}')
print(f'Distance ratio: {d1/d2}')
g_1 = CB2_NAC_gain[a]
g_2 = CB2_NAC_gain[b]
print(f'Gain1: {g_1}')
print(f'Gain2: {g_2}')
source_1 = CB2_NAC_source_count[a]
source_2 = CB2_NAC_source_count[b]
print(f'source 1: {source_1}')
print(f'source 2: {source_2}')
print(f'Source count factor: {source_2/source_1}')
i1 = (source_1 * 4*np.pi * CB2_NAC_distance[a]**2)
i2 = (source_2 * 4*np.pi * CB2_NAC_distance[b]**2)
print(f'Intensity1: {i1}')
print(f'Intensity2: {i2}')
print(f'Result ratio: {i1/i2}')
size1 = CB2_NAC_target_radius[a] * 2
size2 = CB2_NAC_target_radius[b]
print(f'Target_size1: {size1}')
print(f'Target_size2: {size2}')
print(f'Size ratio: {size1/size2}')
i1 = source_1 / (np.pi * size1**2)
i2 = source_2 / (np.pi * size2**2)
print(f'alt Intensity1: {i1}')
print(f'alt Intensity2: {i2}')
print(f'alt Result ratio: {i1/i2}\n')
'''

sort_OBStime = np.argsort(CB2_NAC_time)
CB2_NAC_images = CB2_NAC_images[sort_OBStime]
CB2_NAC_time = CB2_NAC_time[sort_OBStime]
CB2_NAC_phase_angle = CB2_NAC_phase_angle[sort_OBStime]
CB2_NAC_distance = CB2_NAC_distance[sort_OBStime]
CB2_NAC_source_count = CB2_NAC_source_count[sort_OBStime]
CB2_NAC_exposure = CB2_NAC_exposure[sort_OBStime]
CB2_NAC_gain = CB2_NAC_gain[sort_OBStime]
#CB2_NAC_target_radius = CB2_NAC_target_radius[sort_OBStime]

gain29 = 0
gain95 = 0
gain12 = 0
gain215 = 0
gain_other = 0

for i in range(len(CB2_NAC_gain)):
    if (CB2_NAC_gain[i] == 95):
        gain95 += 1
    elif (CB2_NAC_gain[i] == 29):
        gain29 += 1
    elif (CB2_NAC_gain[i] == 12):
        gain12 += 1
    elif (CB2_NAC_gain[i] == 215):
        gain215 += 1
    else:
        gain_other += 1
        
print(f'gain12, gain29, gain95, gain215, gain_other: {gain12}, {gain29}, {gain95}, {gain215}, {gain_other}')

CB2_NAC_gain95 = (CB2_NAC_gain == 95)
print(f'No. of gain95: {CB2_NAC_gain95.sum()}')
CB2_NAC_source_count[CB2_NAC_gain95] = CB2_NAC_source_count[CB2_NAC_gain95] * 4

planet_radius = 2675 #km
pixel_angular_size = 5.9907e-6 #rad per pix
CB2_NAC_target_radius = functs.CCD_Target_radius(planet_radius, CB2_NAC_distance, pixel_angular_size)


CB2_NAC_intensity = CB2_NAC_source_count / (np.pi * CB2_NAC_target_radius**2)
#CB2_NAC_intensity = CB2_NAC_source_count * 4*np.pi * CB2_NAC_distance**2

CB2_NAC_phase_angle = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_phase_angle, 3.4722e-3)
CB2_NAC_intensity = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_intensity, 3.4722e-3)
CB2_NAC_exposure = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_exposure, 3.4722e-3)
CB2_NAC_gain = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_gain, 3.4722e-3)
CB2_NAC_images = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_images, 3.4722e-3)
CB2_NAC_distance = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_distance, 3.4722e-3)
CB2_NAC_source_count = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_source_count, 3.4722e-3)
CB2_NAC_target_radius = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_target_radius, 3.4722e-3)
CB2_NAC_time = functs.adjust_Satellite_Ramping(CB2_NAC_time, CB2_NAC_time, 3.4722e-3)

#limit_intensity = (CB2_NAC_intensity >= 0.82e16)
limit_intensity = (CB2_NAC_intensity >= 0)
NAC_phase_angle = CB2_NAC_phase_angle[limit_intensity]
NAC_intensity = CB2_NAC_intensity[limit_intensity]
NAC_exposure = CB2_NAC_exposure[limit_intensity]
NAC_gain = CB2_NAC_gain[limit_intensity]
NAC_time = CB2_NAC_time[limit_intensity]
NAC_images = CB2_NAC_images[limit_intensity]
NAC_distance = CB2_NAC_distance[limit_intensity]
NAC_target_radius = CB2_NAC_target_radius[limit_intensity]
NAC_source_count = CB2_NAC_source_count[limit_intensity]
print(limit_intensity.sum())

#----------------------------------------------------
# Phase angle vs. Intensity vs. OBStime
'''
plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.xlim(-2, 82)
#plt.xticks(np.linspace(30, 50, 9))
#plt.ylim(1.61e18, 1.77e18)
#plt.axvline(38.3625, linestyle='dashed', color='black')
#plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/CB2_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# OBStime vs. Phase angle vs. Intensity
'''
plt.scatter(NAC_time, NAC_phase_angle, c=NAC_intensity, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.axhline(38.3625, linestyle='dashed', color='black')
plt.axhline(49.0875, linestyle='dashed', color='black')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Intensity', rotation=90)
plt.gca().invert_yaxis()
plt.ylabel('Phase angle [deg]')
plt.xlabel('OBS time [MJD]')
plt.savefig("Titan_images/CB2_NAC_OBStime_vs_Phase_angle_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# OBStime vs. Intensity vs. Phase angle
'''
plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
#plt.xlim(54000, 58000)
#plt.ylim(1.61e18, 1.77e18)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB2_NAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# Distance vs. Intensity vs. OBStime
'''
plt.scatter(NAC_distance, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime', rotation=90)
plt.gca().invert_xaxis()
plt.ylabel('Intensity')
plt.xlabel('Distance [km]')
plt.title('Distance vs. Intensity vs. Phase angle')
plt.savefig("Titan_images/CB2_NAC_Distance_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# Distance vs. OBStime vs. Intensity
'''
plt.scatter(NAC_distance, NAC_time, c=NAC_intensity, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime', rotation=90)
plt.gca().invert_xaxis()
plt.xlabel('Distance')
plt.ylabel('OBStime')
plt.title('Distance vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB2_NAC_Distance_vs_OBStime_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''

