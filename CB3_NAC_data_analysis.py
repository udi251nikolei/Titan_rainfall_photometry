import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import ListedColormap
import csv

import functs

CB3_NAC_photometry = Table.read('Titan_images/CB3_NAC/CB3_NAC_photometry_data.csv', format='ascii')

CB3_NAC_images = np.array(CB3_NAC_photometry['Image file'])
CB3_NAC_time = np.array(CB3_NAC_photometry['Observation time [MJD]'])
CB3_NAC_phase_angle = np.array(CB3_NAC_photometry['Phase angle [deg]'])
CB3_NAC_distance = np.array(CB3_NAC_photometry['Planet-satellite distance [km]'])
CB3_NAC_source_count = np.array(CB3_NAC_photometry['Source count'])
CB3_NAC_exposure = np.array(CB3_NAC_photometry['Exposure time [s]'])
CB3_NAC_gain = np.array(CB3_NAC_photometry['Gain'])
CB3_NAC_target_size = np.array(CB3_NAC_photometry['Target size [pix]'])

print(f'Number of images: {len(CB3_NAC_images)}')

sort_OBStime = np.argsort(CB3_NAC_time)
CB3_NAC_images = CB3_NAC_images[sort_OBStime]
CB3_NAC_time = CB3_NAC_time[sort_OBStime]
CB3_NAC_phase_angle = CB3_NAC_phase_angle[sort_OBStime]
CB3_NAC_distance = CB3_NAC_distance[sort_OBStime]
CB3_NAC_source_count = CB3_NAC_source_count[sort_OBStime]
CB3_NAC_exposure = CB3_NAC_exposure[sort_OBStime]
CB3_NAC_gain = CB3_NAC_gain[sort_OBStime]
CB3_NAC_target_size = CB3_NAC_target_size[sort_OBStime]

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

CB3_NAC_gain95 = (CB3_NAC_gain == 95)
print(f'No. of gain95: {CB3_NAC_gain95.sum()}')
CB3_NAC_source_count[CB3_NAC_gain95] = CB3_NAC_source_count[CB3_NAC_gain95] * 4

CB3_NAC_target_size[CB3_NAC_gain95] = CB3_NAC_target_size[CB3_NAC_gain95] * 2

CB3_NAC_intensity = CB3_NAC_source_count / (np.pi * CB3_NAC_target_size**2)
#CB3_NAC_intensity = CB3_NAC_source_count * 4*np.pi * CB3_NAC_distance**2

CB3_NAC_phase_angle_compensate = functs.Lambertian_reflector_phase_function(CB3_NAC_phase_angle)
angle_range = np.linspace(0, 90, 9001)
lambertian_curve = 0.25*functs.Lambertian_reflector_phase_function(angle_range)
simple_illumination_curve = functs.illumination_factor(angle_range)

#CB3_NAC_intensity = CB3_NAC_intensity / CB3_NAC_phase_angle_compensate

#limit_intensity = (CB3_NAC_intensity >= 1e18) & (CB3_NAC_intensity <= 2.2e18)
limit_intensity = (CB3_NAC_intensity >= 0) & (CB3_NAC_intensity <= 0.2)
NAC_phase_angle = CB3_NAC_phase_angle[limit_intensity]
NAC_intensity = CB3_NAC_intensity[limit_intensity]
NAC_exposure = CB3_NAC_exposure[limit_intensity]
NAC_gain = CB3_NAC_gain[limit_intensity]
NAC_time = CB3_NAC_time[limit_intensity]
NAC_images = CB3_NAC_images[limit_intensity]
NAC_distance = CB3_NAC_distance[limit_intensity]
NAC_source_count = CB3_NAC_source_count[limit_intensity]
NAC_target_size = CB3_NAC_target_size[limit_intensity]
print(limit_intensity.sum())

NAC_phase_angle = functs.adjust_Satellite_Ramping(NAC_time, NAC_phase_angle, 6.944e-3)
NAC_intensity = functs.adjust_Satellite_Ramping(NAC_time, NAC_intensity, 6.944e-3)
NAC_exposure = functs.adjust_Satellite_Ramping(NAC_time, NAC_exposure, 6.944e-3)
NAC_gain = functs.adjust_Satellite_Ramping(NAC_time, NAC_gain, 6.944e-3)
NAC_images = functs.adjust_Satellite_Ramping(NAC_time, NAC_images, 6.944e-3)
NAC_distance = functs.adjust_Satellite_Ramping(NAC_time, NAC_distance, 6.944e-3)
NAC_source_count = functs.adjust_Satellite_Ramping(NAC_time, NAC_source_count, 6.944e-3)
NAC_target_size = functs.adjust_Satellite_Ramping(NAC_time, NAC_target_size, 6.944e-3)
NAC_time = functs.adjust_Satellite_Ramping(NAC_time, NAC_time, 6.944e-3)
'''
NAC_phase_angle_continuous_OBS, NAC_phase_angle = functs.find_Satellites_continuous_OBS(NAC_time, NAC_phase_angle, 1)
NAC_intensity_continuous_OBS, NAC_intensity = functs.find_Satellites_continuous_OBS(NAC_time, NAC_intensity, 1)
NAC_exposure_continuous_OBS, NAC_exposure = functs.find_Satellites_continuous_OBS(NAC_time, NAC_exposure, 1)
NAC_gain_continuous_OBS, NAC_gain = functs.find_Satellites_continuous_OBS(NAC_time, NAC_gain, 1)
NAC_images_continuous_OBS, NAC_images = functs.find_Satellites_continuous_OBS(NAC_time, NAC_images, 1)
NAC_distance_continuous_OBS, NAC_distance = functs.find_Satellites_continuous_OBS(NAC_time, NAC_distance, 1)
NAC_source_count_continuous_OBS, NAC_source_count = functs.find_Satellites_continuous_OBS(NAC_time, NAC_source_count, 1)
NAC_target_size_continuous_OBS, NAC_target_size = functs.find_Satellites_continuous_OBS(NAC_time, NAC_target_size, 1)
NAC_time_continuous_OBS, NAC_time = functs.find_Satellites_continuous_OBS(NAC_time, NAC_time, 1)
'''
#print(sorted(zip(NAC_time, NAC_phase_angle, NAC_intensity)))

#----------------------------------------------------
'''
plt.scatter(NAC_phase_angle_continuous_OBS, NAC_intensity_continuous_OBS, c=NAC_time_continuous_OBS, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
#plt.ylim(0.175, 0.205)
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

functs.plot_Satellites_continuous_OBS(NAC_time_continuous_OBS, NAC_phase_angle_continuous_OBS, NAC_intensity_continuous_OBS)

#----------------------------------------------------

plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.ylim(0.175, 0.205)
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''

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
'''
plt.scatter(NAC_time, NAC_phase_angle, c=NAC_intensity, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.axhline(38.3625, linestyle='dashed', color='black')
plt.axhline(49.0875, linestyle='dashed', color='black')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Intensity', rotation=90)
plt.gca().invert_yaxis()
plt.ylabel('Phase angle [deg]')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Phase_angle_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#----------------------------------------------------

plt.scatter(NAC_time, NAC_intensity, c=NAC_phase_angle, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
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
cbar.ax.set_ylabel('OBStime', rotation=90)
plt.ylabel('Intensity')
plt.xlabel('Distance [km]')
plt.title('Distance vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_NAC_Distance_vs_Intensity_vs_phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------

def linear_fit(x, m, c):
    return (x*m + c)

def two_orders_poly_fitting(x, a1, a2, c):
    return (a1*x**2 + a2*x + c)

def three_orders_poly_fitting(x, a1, a2, a3, c):
    return (a1*x**3 + a2*x**2 + a3*x + c)

def fourth_order_poly_fitting(x, a1, a2, a3, a4, c):
    return (a1*x**4 + a2*x**3 + a3*x**2 + a4*x + c)

data_time_limit = (NAC_time < 56600)
NAC_phase_angle_earlyMission = NAC_phase_angle[data_time_limit]
NAC_intensity_earlyMission = NAC_intensity[data_time_limit]
NAC_exposure_earlyMission = NAC_exposure[data_time_limit]
NAC_gain_earlyMission = NAC_gain[data_time_limit]
NAC_time_earlyMission = NAC_time[data_time_limit]
NAC_images_earlyMission = NAC_images[data_time_limit]
NAC_distance_earlyMission = NAC_distance[data_time_limit]
NAC_target_size_earlyMission = NAC_target_size[data_time_limit]

params = np.polyfit(NAC_time_earlyMission, NAC_intensity_earlyMission, 1)
yfit_earlyMission = linear_fit(NAC_time_earlyMission, *params)

data_time_limit = (NAC_time >= 56600)
NAC_phase_angle_laterMission = NAC_phase_angle[data_time_limit]
NAC_intensity_laterMission = NAC_intensity[data_time_limit]
NAC_exposure_laterMission = NAC_exposure[data_time_limit]
NAC_gain_laterMission = NAC_gain[data_time_limit]
NAC_time_laterMission = NAC_time[data_time_limit]
NAC_images_laterMission = NAC_images[data_time_limit]
NAC_distance_laterMission = NAC_distance[data_time_limit]
NAC_target_size_laterMission = NAC_target_size[data_time_limit]

yfit_align = linear_fit(NAC_time_laterMission, *params) #This are the data points I want to shift the later mission part to the early mission linear fit; using the previous params
params = np.polyfit(NAC_time_laterMission, NAC_intensity_laterMission, 1)
yfit_laterMission = linear_fit(NAC_time_laterMission, *params)

yfit_laterMission_shift_amounts = yfit_laterMission - yfit_align

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(NAC_time_earlyMission, yfit_earlyMission, s=40, marker='+', color='black')
plt.scatter(NAC_time_laterMission, yfit_laterMission, s=40, marker='+', color='black')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB3_NAC_OBStime_diff_adjust.png", dpi = 120, bbox_inches='tight')
plt.show()

align_laterMission_to_earlyMission_intensity = (NAC_intensity_laterMission - yfit_laterMission_shift_amounts)
yfit_laterMission = (yfit_laterMission - yfit_laterMission_shift_amounts)

NAC_intensity = np.hstack((norm_intensity_earlyMission, norm_intensity_laterMission))
NAC_time = np.hstack((NAC_time_earlyMission, NAC_time_laterMission))
NAC_phase_angle = np.hstack((NAC_phase_angle_earlyMission, NAC_phase_angle_laterMission))
NAC_exposure = np.hstack((NAC_exposure_earlyMission, NAC_exposure_laterMission))
NAC_gain = np.hstack((NAC_gain_earlyMission, NAC_gain_laterMission))
NAC_images = np.hstack((NAC_images_earlyMission, NAC_images_laterMission))
NAC_distance = np.hstack((NAC_distance_earlyMission, NAC_distance_laterMission))
NAC_target_size = np.hstack((NAC_target_size_earlyMission, NAC_target_size_laterMission))

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(NAC_time_earlyMission, yfit_earlyMission, s=40, marker='+', color='black')
plt.scatter(NAC_time_laterMission, yfit_laterMission, s=40, marker='+', color='black')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/CB3_NAC_OBStime_diff_adjust.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------
'''
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
'''
#----------------------------------------------------
'''
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
'''
'''
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
'''

#-----------------------------------------------------------

'''
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
