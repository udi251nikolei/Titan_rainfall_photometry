import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import ListedColormap
import csv

import functs
import methane_rainbow_model

CB3_NAC_photometry = Table.read('Titan_images/CB3_NAC/CB3_NAC_photometry_data.csv', format='ascii')

CB3_NAC_images = np.array(CB3_NAC_photometry['Image file'])
CB3_NAC_time = np.array(CB3_NAC_photometry['Observation time [MJD]'])
CB3_NAC_phase_angle = np.array(CB3_NAC_photometry['Phase angle [deg]'])
CB3_NAC_distance = np.array(CB3_NAC_photometry['Planet-satellite distance [km]'])
CB3_NAC_source_count = np.array(CB3_NAC_photometry['Source count'])
CB3_NAC_exposure = np.array(CB3_NAC_photometry['Exposure time [s]'])
CB3_NAC_gain = np.array(CB3_NAC_photometry['Gain'])
print(f'Number of images: {len(CB3_NAC_images)}')

sort_OBStime = np.argsort(CB3_NAC_time)
CB3_NAC_images = CB3_NAC_images[sort_OBStime]
CB3_NAC_time = CB3_NAC_time[sort_OBStime]
CB3_NAC_phase_angle = CB3_NAC_phase_angle[sort_OBStime]
CB3_NAC_distance = CB3_NAC_distance[sort_OBStime]
CB3_NAC_source_count = CB3_NAC_source_count[sort_OBStime]
CB3_NAC_exposure = CB3_NAC_exposure[sort_OBStime]
CB3_NAC_gain = CB3_NAC_gain[sort_OBStime]

# Saturn in the background of images
images_with_Saturn = (CB3_NAC_images == 'N1684686307_1') | (CB3_NAC_images == 'N1691555298_1') | (CB3_NAC_images == 'N1691555408_1') | (CB3_NAC_images == 'N1691555518_1') | (CB3_NAC_images == 'N1710198856_1') | (CB3_NAC_images == 'N1710198966_1') | (CB3_NAC_images == 'N1710199076_1')
CB3_NAC_images = CB3_NAC_images[~images_with_Saturn]
CB3_NAC_time = CB3_NAC_time[~images_with_Saturn]
CB3_NAC_phase_angle = CB3_NAC_phase_angle[~images_with_Saturn]
CB3_NAC_distance = CB3_NAC_distance[~images_with_Saturn]
CB3_NAC_source_count = CB3_NAC_source_count[~images_with_Saturn]
CB3_NAC_exposure = CB3_NAC_exposure[~images_with_Saturn]
CB3_NAC_gain = CB3_NAC_gain[~images_with_Saturn]
print(f'Number of images: {len(CB3_NAC_images)}')

'''
a = 383
b = 384
image_1 = CB3_NAC_images[a]
image_2 = CB3_NAC_images[b]
print(f'Image1: {image_1}')
print(f'Image2: {image_2}')
d1 = CB3_NAC_distance[a]
d2 = CB3_NAC_distance[b]
print(f'D1: {d1}')
print(f'D2: {d2}')
print(f'Distance ratio: {d1/d2}')
g_1 = CB3_NAC_gain[a]
g_2 = CB3_NAC_gain[b]
print(f'Gain1: {g_1}')
print(f'Gain2: {g_2}')
source_1 = CB3_NAC_source_count[a] * 4
source_2 = CB3_NAC_source_count[b]
print(f'source 1: {source_1}')
print(f'source 2: {source_2}')
print(f'Source count ratio: {source_1/source_2}')
i1 = (source_1 * 4*np.pi * CB3_NAC_distance[a]**2)
i2 = (source_2 * 4*np.pi * CB3_NAC_distance[b]**2)
print(f'Intensity1: {i1}')
print(f'Intensity2: {i2}')
print(f'Result ratio: {i1/i2}')
size1 = CB3_NAC_target_size[a] * 2
size2 = CB3_NAC_target_size[b]
print(f'Target_size1: {size1}')
print(f'Target_size2: {size2}')
print(f'Size ratio: {size1/size2}')
i1 = source_1 / (np.pi * size1**2)
i2 = source_2 / (np.pi * size2**2)
print(f'alt Intensity1: {i1}')
print(f'alt Intensity2: {i2}')
print(f'alt Result ratio: {i1/i2}\n')
'''

gain29 = 0
gain95 = 0
gain12 = 0
gain215 = 0
gain_other = 0

for i in range(len(CB3_NAC_gain)):
    if (CB3_NAC_gain[i] == 95):
        gain95 += 1
    elif (CB3_NAC_gain[i] == 29):
        gain29 += 1
    elif (CB3_NAC_gain[i] == 12):
        gain12 += 1
    elif (CB3_NAC_gain[i] == 215):
        gain215 += 1
    else:
        gain_other += 1
        
print(f'gain12, gain29, gain95, gain215, gain_other: {gain12}, {gain29}, {gain95}, {gain215}, {gain_other}')

phaseAngle_0_to_30 = (CB3_NAC_phase_angle < 30)
phaseAngle_30_to_50 = (CB3_NAC_phase_angle >= 30) & (CB3_NAC_phase_angle <= 50)
phaseAngle_50_to_80 = (CB3_NAC_phase_angle > 50)
print(f'Phase angles: below 30 = {phaseAngle_0_to_30.sum()}, 30-50 = {phaseAngle_30_to_50.sum()}, above 50 = {phaseAngle_50_to_80.sum()}')

CB3_NAC_gain95 = (CB3_NAC_gain == 95)
print(f'No. of gain95: {CB3_NAC_gain95.sum()}')
CB3_NAC_source_count[CB3_NAC_gain95] = CB3_NAC_source_count[CB3_NAC_gain95] * 3.965

planet_radius = 2675 #km
pixel_angular_size = 5.9907e-6 #rad per pix
CB3_NAC_target_radius = functs.CCD_Target_radius(planet_radius, CB3_NAC_distance, pixel_angular_size)

CB3_NAC_intensity = CB3_NAC_source_count / (np.pi * CB3_NAC_target_radius**2)
#CB3_NAC_intensity = CB3_NAC_source_count * 4*np.pi * CB3_NAC_distance**2


CB3_NAC_phase_angle = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_phase_angle, 3.4722e-3)
CB3_NAC_intensity = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_intensity, 3.4722e-3)
CB3_NAC_exposure = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_exposure, 3.4722e-3)
CB3_NAC_gain = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_gain, 3.4722e-3)
CB3_NAC_images = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_images, 3.4722e-3)
CB3_NAC_distance = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_distance, 3.4722e-3)
CB3_NAC_source_count = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_source_count, 3.4722e-3)
CB3_NAC_target_radius = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_target_radius, 3.4722e-3)
CB3_NAC_time = functs.adjust_Satellite_Ramping(CB3_NAC_time, CB3_NAC_time, 3.4722e-3)

'''
#limit_intensity = (CB3_NAC_intensity >= 0.8e18) & (CB3_NAC_intensity <= 1.6e18)
limit_intensity = (CB3_NAC_intensity >= 0) #(CB3_NAC_intensity >= 0.12) & (CB3_NAC_intensity <= 0.21)
NAC_phase_angle = CB3_NAC_phase_angle[limit_intensity]
NAC_intensity = CB3_NAC_intensity[limit_intensity]
NAC_exposure = CB3_NAC_exposure[limit_intensity]
NAC_gain = CB3_NAC_gain[limit_intensity]
NAC_time = CB3_NAC_time[limit_intensity]
NAC_images = CB3_NAC_images[limit_intensity]
NAC_distance = CB3_NAC_distance[limit_intensity]
NAC_source_count = CB3_NAC_source_count[limit_intensity]
NAC_target_radius = CB3_NAC_target_radius[limit_intensity]
print(limit_intensity.sum())
'''
NAC_phase_angle = CB3_NAC_phase_angle
NAC_intensity = CB3_NAC_intensity
NAC_exposure = CB3_NAC_exposure
NAC_gain = CB3_NAC_gain
NAC_time = CB3_NAC_time
NAC_images = CB3_NAC_images
NAC_distance = CB3_NAC_distance
NAC_source_count = CB3_NAC_source_count
NAC_target_radius = CB3_NAC_target_radius

'''
outliers = (CB3_NAC_intensity >= 0) #& (CB3_NAC_intensity <= 0.12)
outlier_images = CB3_NAC_images[outliers]
outlier_intensity = CB3_NAC_intensity[outliers]
outlier_phase_angle = CB3_NAC_phase_angle[outliers]
outlier_time = CB3_NAC_time[outliers]
#print(list(zip(outlier_images, outlier_phase_angle, outlier_intensity)))
'''

#limit_angle = (CB3_NAC_phase_angle >= 30) & (CB3_NAC_phase_angle <= 50)
#phase_angle_30_to_50 = CB3_NAC_phase_angle[limit_angle]
#intensity_30_to_50 = CB3_NAC_intensity[limit_angle]
#time_30_to_50 = CB3_NAC_time[limit_angle]

#-----------------------------------------------------------------

NAC_phaseAngle_continuousOBS, NAC_phaseAngle_singleOBS = functs.find_Satellites_continuous_OBS(NAC_time, NAC_phase_angle, 1)
NAC_intensity_continuousOBS, NAC_intensity_singleOBS = functs.find_Satellites_continuous_OBS(NAC_time, NAC_intensity, 1)
NAC_time_continuousOBS, NAC_time_singleOBS = functs.find_Satellites_continuous_OBS(NAC_time, NAC_time, 1)

#-----------------------------------------------------------------
# Phase angle vs. Intensity vs. OBStime
'''
plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.plot(x_angles, lambertian_fit, linestyle='dashed', color='tab:blue')
plt.xlim(-2, 82)
#plt.xticks(np.linspace(30, 50, 9))
plt.ylim(0.09, 0.25)
#plt.axvline(38.3625, linestyle='dashed', color='black')
#plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity.png", dpi = 120, bbox_inches='tight')
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
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Phase_angle_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------
# OBStime vs. Intensity vs. Phase angle
'''
plt.scatter(NAC_time, NAC_intensity, c=NAC_phase_angle, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
#plt.xlim(54000, 58000)
#plt.ylim(0.12, 0.22)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
'''
plt.scatter(time_30_to_50, intensity_30_to_50, c=phase_angle_30_to_50, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
#plt.xlim(54000, 58000)
plt.ylim(0.12, 0.22)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------

angles_model = methane_rainbow_model.angles
intensity_939nm = methane_rainbow_model.intensity_939nm

#plt.plot(angles_model, intensity_939nm)
#plt.show()
'''
def linear_fit(x, m, c):
    return (x*m + c)

angle_limit = (angles_model <= 31.6125)
x_known = angles_model[angle_limit]
y_known = intensity_939nm[angle_limit]
param_linear = np.polyfit(x_known, y_known, 1)
x_interpolate = np.arange(20, 30, 0.1)
y_interpolate = linear_fit(x_interpolate, *param_linear)
x_zerolvl = (y_interpolate < 0.00045)
print(list(zip(x_interpolate[x_zerolvl], y_interpolate[x_zerolvl])))


angle_limit = (angles_model >= 49.0875)
x_known = angles_model[angle_limit]
y_known = intensity_939nm[angle_limit] #0.00014385
#print(list(zip(x_known, y_known)))
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
plt.savefig("Titan_images/CB3_NAC_Distance_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
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
plt.savefig("Titan_images/CB3_NAC_Distance_vs_OBStime_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''

#----------------------------------------------------
'''
NAC_intensityTmp = NAC_intensity / NAC_PhaseAngle_factors

plt.scatter(NAC_phase_angle, NAC_intensityTmp, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.ylim(1.61e18, 1.77e18)
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

plt.scatter(NAC_time, NAC_intensityTmp, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(54000, 58000)
#plt.ylim(1.61e18, 1.77e18)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.title('OBStime vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''


'''
plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(54250, MissionPhase1_avgIntensity, s=40, marker='+', color='black')
plt.scatter(54250, MissionPhase1_midIntensity, s=40, marker='+', color='red')
plt.scatter(54800, MissionPhase2_avgIntensity, s=40, marker='+', color='black')
plt.scatter(54800, MissionPhase2_midIntensity, s=40, marker='+', color='red')
plt.scatter(55400, MissionPhase3_avgIntensity, s=40, marker='+', color='black')
plt.scatter(55400, MissionPhase3_midIntensity, s=40, marker='+', color='red')
plt.scatter(55900, MissionPhase4_avgIntensity, s=40, marker='+', color='black')
plt.scatter(55900, MissionPhase4_midIntensity, s=40, marker='+', color='red')
plt.scatter(56250, MissionPhase5_avgIntensity, s=40, marker='+', color='black')
plt.scatter(56250, MissionPhase5_midIntensity, s=40, marker='+', color='red')
plt.scatter(56800, MissionPhase6_avgIntensity, s=40, marker='+', color='black')
plt.scatter(56800, MissionPhase6_midIntensity, s=40, marker='+', color='red')
plt.scatter(57500, MissionPhase7_avgIntensity, s=40, marker='+', color='black')
plt.scatter(57500, MissionPhase7_midIntensity, s=40, marker='+', color='red')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

NAC_intensity[MissionPhase1] = NAC_intensity[MissionPhase1] - (MissionPhase1_midIntensity - MissionPhase4_avgIntensity)
NAC_intensity[MissionPhase2] = NAC_intensity[MissionPhase2] - (MissionPhase2_midIntensity - MissionPhase4_midIntensity)
NAC_intensity[MissionPhase3] = NAC_intensity[MissionPhase3] - (MissionPhase3_midIntensity - MissionPhase4_midIntensity) + 0.03e18
NAC_intensity[MissionPhase4] = NAC_intensity[MissionPhase4] - (MissionPhase4_avgIntensity - MissionPhase4_avgIntensity)
NAC_intensity[MissionPhase5] = NAC_intensity[MissionPhase5] - (MissionPhase5_avgIntensity - MissionPhase4_avgIntensity)
NAC_intensity[MissionPhase6] = NAC_intensity[MissionPhase6] - (MissionPhase6_avgIntensity - MissionPhase4_avgIntensity)
NAC_intensity[MissionPhase7] = NAC_intensity[MissionPhase7] - (MissionPhase7_avgIntensity - MissionPhase4_avgIntensity)

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.show()
'''

'''
yaxis_align = (NAC_time > 55000) & (NAC_time <= 56300)
NAC_intensity_align = NAC_intensity[yaxis_align]
NAC_time_align = NAC_time[yaxis_align]
horizontal_alignment = np.mean(NAC_intensity_align)
print(horizontal_alignment)


MissionPhase1_to_2 = (NAC_time <= 55000)
NAC_intensity_MissionPhase1_to_2 = NAC_intensity[MissionPhase1_to_2]
NAC_time_MissionPhase1_to_2 = NAC_time[MissionPhase1_to_2]
params = np.polyfit(NAC_time_MissionPhase1_to_2, NAC_intensity_MissionPhase1_to_2, 1)
yfit_MissionPhase1_to_2 = linear_fit(NAC_time_MissionPhase1_to_2, *params)
align_MissionPhase1 = yfit_MissionPhase1_to_2 - horizontal_alignment
#MissionPhase1 = (NAC_time <= 54600)
#NAC_intensity[MissionPhase1_to_2] = NAC_intensity[MissionPhase1_to_2] - align_MissionPhase1

MissionPhase2_to_3 = (NAC_time > 55000) & (NAC_time <= 55500)
NAC_intensity_MissionPhase2_to_3 = NAC_intensity[MissionPhase2_to_3]
NAC_time_MissionPhase2_to_3 = NAC_time[MissionPhase2_to_3]
params = np.polyfit(NAC_time_MissionPhase2_to_3, NAC_intensity_MissionPhase2_to_3, 1)
yfit_MissionPhase2_to_3 = linear_fit(NAC_time_MissionPhase2_to_3, *params)
align_MissionPhase2_to_3 = yfit_MissionPhase2_to_3 - horizontal_alignment
#NAC_intensity[MissionPhase2_to_3] = NAC_intensity[MissionPhase2_to_3] - align_MissionPhase2_to_3

MissionPhase3_to_4 = (NAC_time > 55000) & (NAC_time <= 56300)
NAC_intensity_MissionPhase3_to_4 = NAC_intensity[MissionPhase3_to_4]
NAC_time_MissionPhase3_to_4 = NAC_time[MissionPhase3_to_4]
params = np.polyfit(NAC_time_MissionPhase3_to_4, NAC_intensity_MissionPhase3_to_4, 1)
yfit_MissionPhase3_to_4 = linear_fit(NAC_time_MissionPhase3_to_4, *params)
align_MissionPhase3_to_4 = yfit_MissionPhase3_to_4 - horizontal_alignment
#NAC_intensity[MissionPhase3_to_4] = NAC_intensity[MissionPhase3_to_4] - align_MissionPhase3_to_4

MissionPhase4_to_5 = (NAC_time > 55500) & (NAC_time <= 56600)
NAC_intensity_MissionPhase4_to_5 = NAC_intensity[MissionPhase4_to_5]
NAC_time_MissionPhase4_to_5 = NAC_time[MissionPhase4_to_5]
params = np.polyfit(NAC_time_MissionPhase4_to_5, NAC_intensity_MissionPhase4_to_5, 1)
yfit_MissionPhase4_to_5 = linear_fit(NAC_time_MissionPhase4_to_5, *params)
align_MissionPhase4_to_5 = yfit_MissionPhase4_to_5 - horizontal_alignment
#NAC_intensity[MissionPhase4_to_5] = NAC_intensity[MissionPhase4_to_5] - align_MissionPhase4_to_5

MissionPhase6_to_7 = (NAC_time > 56600)
NAC_intensity_MissionPhase6_to_7 = NAC_intensity[MissionPhase6_to_7]
NAC_time_MissionPhase6_to_7 = NAC_time[MissionPhase6_to_7]
params = np.polyfit(NAC_time_MissionPhase6_to_7, NAC_intensity_MissionPhase6_to_7, 1)
yfit_MissionPhase6_to_7 = linear_fit(NAC_time_MissionPhase6_to_7, *params)
align_MissionPhase6_to_7 = yfit_MissionPhase6_to_7 - horizontal_alignment
#NAC_intensity[MissionPhase6_to_7] = NAC_intensity[MissionPhase6_to_7] - align_MissionPhase6_to_7

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(NAC_time_MissionPhase1_to_2, yfit_MissionPhase1_to_2, s=40, marker='+', color='black')
plt.scatter(NAC_time_MissionPhase2_to_3, yfit_MissionPhase2_to_3, s=40, marker='+', color='black')
plt.scatter(NAC_time_MissionPhase3_to_4, yfit_MissionPhase3_to_4, s=40, marker='+', color='black')
plt.scatter(NAC_time_MissionPhase4_to_5, yfit_MissionPhase4_to_5, s=40, marker='+', color='black')
plt.scatter(NAC_time_MissionPhase6_to_7, yfit_MissionPhase6_to_7, s=40, marker='+', color='black')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()



plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(NAC_time_MissionPhase1_to_5, yfit_MissionPhase1_to_5, s=40, marker='+', color='black')
plt.scatter(NAC_time_MissionPhase6_to_7, yfit_MissionPhase6_to_7, s=40, marker='+', color='black')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

NAC_intensity[MissionPhase1_to_5] = NAC_intensity[MissionPhase1_to_5] - align_MissionPhase1_to_5
NAC_intensity[MissionPhase6_to_7] = NAC_intensity[MissionPhase6_to_7] - align_MissionPhase6_to_7

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()






yfit_laterMission_shift_amounts = yfit_laterMission - yfit_align
NAC_intensity[laterMission_OBS] = NAC_intensity[laterMission_OBS] - yfit_laterMission_shift_amounts
yfit_laterMission = (yfit_laterMission - yfit_laterMission_shift_amounts)


plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(NAC_time_earlyMission, yfit_earlyMission, s=40, marker='+', color='black')
plt.scatter(NAC_time_laterMission, yfit_laterMission, s=40, marker='+', color='black')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

#----------------------------------------------------

plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')

plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(54000, 58000)
#plt.ylim(0.175, 0.205)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.title('OBStime vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------

plt.scatter(NAC_distance, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg', rotation=90)
plt.gca().invert_xaxis()
plt.ylabel('Intensity')
plt.xlabel('Distance [km]')
plt.title('Distance vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_NAC_Distance_vs_Intensity_vs_phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------

plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

params = np.polyfit(NAC_phase_angle, NAC_intensity, 1)
yfit_normalization = linear_fit(NAC_phase_angle, *params)

#NAC_intensity = (NAC_intensity / yfit_normalization)

plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()


time_limit = (NAC_time >= 55700) & (NAC_time <= 56200)
NAC_phase_angle = NAC_phase_angle[time_limit]
NAC_intensity = NAC_intensity[time_limit]
NAC_time = NAC_time[time_limit]
NAC_images = NAC_images[time_limit]
NAC_distance = NAC_distance[time_limit]
NAC_gain = NAC_gain[time_limit]
NAC_exposure = NAC_exposure[time_limit]
NAC_source_count = NAC_source_count[time_limit]
NAC_target_size = NAC_target_size[time_limit]

print(min(NAC_intensity), max(NAC_intensity))

plt.scatter(NAC_time, NAC_intensity, c=NAC_gain, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(54000, 58000)
plt.ylim(0.175, 0.205)
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

fields = ['Image file', 'Observation time [MJD]', 'Phase angle [deg]', 'Planet-satellite distance [km]', 'Source count', 'Exposure time [s]', 'Gain', 'Intensity']
path = "Titan_images/"
with open(path + "CB3_NAC_data.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)
    
for i in range(time_limit.sum()):
    photometric_data = np.array([NAC_images[i], NAC_time[i], NAC_phase_angle[i], NAC_distance[i], NAC_source_count[i], NAC_exposure[i], NAC_gain[i], NAC_intensity[i]])

    with open(path + "CB3_NAC_data.csv", "a") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(photometric_data)


#-----------------------------------------------------------


fields = ['Image file', 'Observation time [MJD]', 'Phase angle [deg]', 'Planet-satellite distance [km]', 'Source count', 'Exposure time [s]', 'Gain', 'Intensity']
path = "Titan_images/"
with open(path + "CB3_NAC_data.csv", "w") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(fields)
    
for i in range(len(NAC_time)):
    photometric_data = np.array([NAC_images[i], NAC_time[i], NAC_phase_angle[i], NAC_distance[i], NAC_source_count[i], NAC_exposure[i], NAC_gain[i], NAC_intensity[i]])

    with open(path + "CB3_NAC_data.csv", "a") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(photometric_data)


'''



