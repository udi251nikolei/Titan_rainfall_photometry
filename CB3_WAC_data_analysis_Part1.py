import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import ListedColormap
import csv

import functs

CB3_WAC_photometry = Table.read('Titan_images/CB3_WAC/CB3_WAC_photometry_data.csv', format='ascii')

CB3_WAC_images = np.array(CB3_WAC_photometry['Image file'])
CB3_WAC_time = np.array(CB3_WAC_photometry['Observation time [MJD]'])
CB3_WAC_phase_angle = np.array(CB3_WAC_photometry['Phase angle [deg]'])
CB3_WAC_distance = np.array(CB3_WAC_photometry['Planet-satellite distance [km]'])
CB3_WAC_source_count = np.array(CB3_WAC_photometry['Source count'])
CB3_WAC_exposure = np.array(CB3_WAC_photometry['Exposure time [s]'])
CB3_WAC_gain = np.array(CB3_WAC_photometry['Gain'])
#CB3_WAC_target_radius = np.array(CB3_WAC_photometry['Target Radius [pix]'])
print(f'Number of images: {len(CB3_WAC_images)}')
'''
a = 1122
b = 1125
image_1 = CB3_WAC_images[a]
image_2 = CB3_WAC_images[b]
print(f'Image1: {image_1}')
print(f'Image2: {image_2}')
d1 = CB3_WAC_distance[a]
d2 = CB3_WAC_distance[b]
print(f'D1: {d1}')
print(f'D2: {d2}')
print(f'Distance ratio: {d1/d2}')
g_1 = CB3_WAC_gain[a]
g_2 = CB3_WAC_gain[b]
print(f'Gain1: {g_1}')
print(f'Gain2: {g_2}')
source_1 = CB3_WAC_source_count[a]
source_2 = CB3_WAC_source_count[b]
print(f'source 1: {source_1}')
print(f'source 2: {source_2}')
print(f'Source count factor: {source_2/source_1}')
i1 = (source_1 * 4*np.pi * CB3_WAC_distance[a]**2)
i2 = (source_2 * 4*np.pi * CB3_WAC_distance[b]**2)
print(f'Intensity1: {i1}')
print(f'Intensity2: {i2}')
print(f'Result ratio: {i1/i2}')
size1 = CB3_WAC_target_radius[a] * 2
size2 = CB3_WAC_target_radius[b]
print(f'Target_size1: {size1}')
print(f'Target_size2: {size2}')
print(f'Size ratio: {size1/size2}')
i1 = source_1 / (np.pi * size1**2)
i2 = source_2 / (np.pi * size2**2)
print(f'alt Intensity1: {i1}')
print(f'alt Intensity2: {i2}')
print(f'alt Result ratio: {i1/i2}\n')
'''
sort_OBStime = np.argsort(CB3_WAC_time)
CB3_WAC_images = CB3_WAC_images[sort_OBStime]
CB3_WAC_time = CB3_WAC_time[sort_OBStime]
CB3_WAC_phase_angle = CB3_WAC_phase_angle[sort_OBStime]
CB3_WAC_distance = CB3_WAC_distance[sort_OBStime]
CB3_WAC_source_count = CB3_WAC_source_count[sort_OBStime]
CB3_WAC_exposure = CB3_WAC_exposure[sort_OBStime]
CB3_WAC_gain = CB3_WAC_gain[sort_OBStime]
#CB3_WAC_target_radius = CB3_WAC_target_radius[sort_OBStime]

gain29 = 0
gain95 = 0
gain12 = 0
gain215 = 0
gain_other = 0

for i in range(len(CB3_WAC_gain)):
    if (CB3_WAC_gain[i] == 95):
        gain95 += 1
    elif (CB3_WAC_gain[i] == 29):
        gain29 += 1
    elif (CB3_WAC_gain[i] == 12):
        gain12 += 1
    elif (CB3_WAC_gain[i] == 215):
        gain215 += 1
    else:
        gain_other += 1
        
print(f'gain12, gain29, gain95, gain215, gain_other: {gain12}, {gain29}, {gain95}, {gain215}, {gain_other}')

CB3_WAC_gain95 = (CB3_WAC_gain == 95)
print(f'No. of gain95: {CB3_WAC_gain95.sum()}')
CB3_WAC_source_count[CB3_WAC_gain95] = CB3_WAC_source_count[CB3_WAC_gain95] * 3.68

planet_radius = 2675 #km
pixel_angular_size = 59.749e-6 #rad per pix
CB3_WAC_target_radius = functs.CCD_Target_radius(planet_radius, CB3_WAC_distance, pixel_angular_size)

CB3_WAC_intensity = CB3_WAC_source_count / (np.pi * CB3_WAC_target_radius**2)
#B3_WAC_intensity = CB3_WAC_source_count * 4*np.pi * CB3_WAC_distance**2


CB3_WAC_phase_angle = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_phase_angle, 3.4722e-3)
CB3_WAC_intensity = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_intensity, 3.4722e-3)
CB3_WAC_exposure = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_exposure, 3.4722e-3)
CB3_WAC_gain = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_gain, 3.4722e-3)
CB3_WAC_images = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_images, 3.4722e-3)
CB3_WAC_distance = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_distance, 3.4722e-3)
CB3_WAC_source_count = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_source_count, 3.4722e-3)
CB3_WAC_target_radius = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_target_radius, 3.4722e-3)
CB3_WAC_time = functs.adjust_Satellite_Ramping(CB3_WAC_time, CB3_WAC_time, 3.4722e-3)

#limit_intensity = (CB3_WAC_intensity >= 0.82e16)
limit_intensity = (CB3_WAC_intensity >= 0.105)
WAC_phase_angle = CB3_WAC_phase_angle[limit_intensity]
WAC_intensity = CB3_WAC_intensity[limit_intensity]
WAC_exposure = CB3_WAC_exposure[limit_intensity]
WAC_gain = CB3_WAC_gain[limit_intensity]
WAC_time = CB3_WAC_time[limit_intensity]
WAC_images = CB3_WAC_images[limit_intensity]
WAC_distance = CB3_WAC_distance[limit_intensity]
WAC_target_radius = CB3_WAC_target_radius[limit_intensity]
WAC_source_count = CB3_WAC_source_count[limit_intensity]
print(limit_intensity.sum())

#----------------------------------------------------
# Phase angle vs. Intensity vs. OBStime
'''
plt.scatter(WAC_phase_angle, WAC_intensity, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.xlim(-2, 82)
#plt.xticks(np.linspace(30, 50, 9))
#plt.ylim(1.61e18, 1.77e18)
#plt.axvline(38.3625, linestyle='dashed', color='black')
#plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/CB3_WAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# OBStime vs. Phase angle vs. Intensity
'''
plt.scatter(WAC_time, WAC_phase_angle, c=WAC_intensity, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.axhline(38.3625, linestyle='dashed', color='black')
plt.axhline(49.0875, linestyle='dashed', color='black')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Intensity', rotation=90)
plt.gca().invert_yaxis()
plt.ylabel('Phase angle [deg]')
plt.xlabel('OBS time [MJD]')
plt.savefig("Titan_images/CB3_WAC_OBStime_vs_Phase_angle_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# OBStime vs. Intensity vs. Phase angle
'''
plt.scatter(WAC_time, WAC_intensity, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
#plt.xlim(54000, 58000)
#plt.ylim(1.61e18, 1.77e18)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB3_WAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# Distance vs. Intensity vs. OBStime
'''
plt.scatter(WAC_distance, WAC_intensity, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime', rotation=90)
plt.gca().invert_xaxis()
plt.ylabel('Intensity')
plt.xlabel('Distance [km]')
plt.title('Distance vs. Intensity vs. Phase angle')
plt.savefig("Titan_images/CB3_WAC_Distance_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# Distance vs. OBStime vs. Intensity
'''
plt.scatter(WAC_distance, WAC_time, c=WAC_intensity, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime', rotation=90)
plt.gca().invert_xaxis()
plt.xlabel('Distance')
plt.ylabel('OBStime')
plt.title('Distance vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_WAC_Distance_vs_OBStime_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
