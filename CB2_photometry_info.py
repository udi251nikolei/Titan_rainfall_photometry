import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import csv
from matplotlib.colors import ListedColormap

import functs

CB2_WAC_photometry = Table.read('Titan_images/CB2_WAC/CB2_WAC_photometry_data.csv', format='ascii')

CB2_WAC_images = np.array(CB2_WAC_photometry['Image file'])
CB2_WAC_time = np.array(CB2_WAC_photometry['Observation time [MJD]'])
CB2_WAC_phase_angle = np.array(CB2_WAC_photometry['Phase angle [deg]'])
CB2_WAC_distance = np.array(CB2_WAC_photometry['Planet-satellite distance [km]'])
CB2_WAC_source_count = np.array(CB2_WAC_photometry['Source count'])
CB2_WAC_exposure = np.array(CB2_WAC_photometry['Exposure time [s]'])
CB2_WAC_gain = np.array(CB2_WAC_photometry['Gain'])

CB2_WAC_gain95 = (CB2_WAC_gain == 95)
CB2_WAC_source_count[CB2_WAC_gain95] = CB2_WAC_source_count[CB2_WAC_gain95] * 4

CB2_WAC_phase_angle_compensate = functs.Lambertian_reflector_phase_function(CB2_WAC_phase_angle)
CB2_WAC_source_count = CB2_WAC_source_count / CB2_WAC_phase_angle_compensate

#planet_radius = 2875 #km
#CB2_WAC_target_pixel_radius = functs.CCD_Target_radius(planet_radius, CB2_WAC_distance, 59.749e-6)

#CB2_WAC_intensity = CB2_WAC_source_count / (np.pi * CB2_WAC_target_pixel_radius**2)
CB2_WAC_intensity = CB2_WAC_source_count * 4*np.pi * CB2_WAC_distance**2

phase_angle35_55 = (CB2_WAC_phase_angle >= 30.0125) & (CB2_WAC_phase_angle <= 49.9875)
WAC_phase_angle = CB2_WAC_phase_angle[phase_angle35_55]
WAC_intensity = CB2_WAC_intensity[phase_angle35_55]
WAC_exposure = CB2_WAC_exposure[phase_angle35_55]
WAC_gain = CB2_WAC_gain[phase_angle35_55]
WAC_time = CB2_WAC_time[phase_angle35_55]
WAC_images = CB2_WAC_images[phase_angle35_55]
WAC_phase_angle_compensate = CB2_WAC_phase_angle_compensate[phase_angle35_55]

WAC_phase_angle = functs.Satellite_Ramping_adjust(WAC_time, WAC_phase_angle, 6.944e-3)
WAC_intensity = functs.Satellite_Ramping_adjust(WAC_time, WAC_intensity, 6.944e-3)
WAC_phase_angle_compensate = functs.Satellite_Ramping_adjust(WAC_time, WAC_phase_angle_compensate, 6.944e-3)
WAC_time = functs.Satellite_Ramping_adjust(WAC_time, WAC_time, 6.944e-3)

colourmap = ListedColormap(['red', 'orange', 'yellow', 'green', 'royalblue', 'darkviolet', 'deeppink'])
plt.scatter(WAC_phase_angle, WAC_intensity, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.colorbar()
plt.axvline(48.7625, linestyle='dashed', color='black')
plt.axvline(38.9625, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('CB2 WAC, ~752 nm')
plt.savefig("Titan_images/CB2_WAC_angle_vs_intensity_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

colourmap = ListedColormap(['red', 'orange', 'yellow', 'green', 'royalblue', 'darkviolet', 'deeppink'])
plt.scatter(WAC_time, WAC_phase_angle, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.axhline(48.7625, linestyle='dashed', color='black')
plt.axhline(38.9625, linestyle='dashed', color='black')
plt.colorbar()
plt.gca().invert_yaxis()
plt.ylabel('Phase angle [deg]')
plt.xlabel('OBS time [MJD]')
plt.title('CB2 WAC, ~752 nm')
plt.savefig("Titan_images/CB2_WAC_OBStime_vs_angle_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

plt.scatter(WAC_time, WAC_intensity, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.title('CB2 WAC, ~752 nm')
plt.savefig("Titan_images/CB2_WAC_OBStime_vs_Intensity_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#--------------------------------

CB2_NAC_photometry = Table.read('Titan_images/CB2_NAC/CB2_NAC_photometry_data.csv', format='ascii')

CB2_NAC_images = np.array(CB2_NAC_photometry['Image file'])
CB2_NAC_time = np.array(CB2_NAC_photometry['Observation time [MJD]'])
CB2_NAC_phase_angle = np.array(CB2_NAC_photometry['Phase angle [deg]'])
CB2_NAC_distance = np.array(CB2_NAC_photometry['Planet-satellite distance [km]'])
CB2_NAC_source_count = np.array(CB2_NAC_photometry['Source count'])
CB2_NAC_exposure = np.array(CB2_NAC_photometry['Exposure time [s]'])
CB2_NAC_gain = np.array(CB2_NAC_photometry['Gain'])

"""
a = 884
b = 885
phase_1 = CB2_NAC_phase_angle[a]
phase_2 = CB2_NAC_phase_angle[b]
print(f'Phase1: {phase_1}')
print(f'Phase2: {phase_2}')
d1 = CB2_NAC_distance[a]
d2 = CB2_NAC_distance[b]
print(f'D1: {d1}')
print(f'D2: {d2}')
print(f'Distance ratio: {d1/d2}')
source_1 = CB2_NAC_source_count[a]
source_2 = CB2_NAC_source_count[b] * 4
print(f'source 1: {source_1}')
print(f'source 2: {source_2}')
print(f'Source count ratio: {source_1/source_2}')
i1 = (source_1 * 4*np.pi * CB2_NAC_distance[a]**2)
i2 = (source_2 * 4*np.pi * CB2_NAC_distance[b]**2)
print(f'Intensity1: {i1}')
print(f'Intensity2: {i2}')
print(f'Result ratio: {i1/i2}\n')
"""

CB2_NAC_gain95 = (CB2_NAC_gain == 95)
CB2_NAC_source_count[CB2_NAC_gain95] = CB2_NAC_source_count[CB2_NAC_gain95] * 4

CB2_NAC_phase_angle_compensate = functs.Lambertian_reflector_phase_function(CB2_NAC_phase_angle)
CB2_NAC_source_count = CB2_NAC_source_count / CB2_NAC_phase_angle_compensate

#planet_radius = 2875 #km
#CB2_NAC_target_pixel_radius = functs.CCD_Target_radius(planet_radius, CB2_NAC_distance, 5.9907e-6)

#CB2_NAC_intensity = CB2_NAC_source_count / (np.pi * CB2_NAC_target_pixel_radius**2)
CB2_NAC_intensity = CB2_NAC_source_count * 4*np.pi * CB2_NAC_distance**2

phase_angle35_55 = (CB2_NAC_phase_angle >= 30.0125) & (CB2_NAC_phase_angle <= 49.9875)
NAC_phase_angle = CB2_NAC_phase_angle[phase_angle35_55]
NAC_intensity = CB2_NAC_intensity[phase_angle35_55]
NAC_exposure = CB2_NAC_exposure[phase_angle35_55]
NAC_gain = CB2_NAC_gain[phase_angle35_55]
NAC_time = CB2_NAC_time[phase_angle35_55]
NAC_images = CB2_NAC_images[phase_angle35_55]
NAC_distance = CB2_NAC_distance[phase_angle35_55]
NAC_source_count = CB2_NAC_source_count[phase_angle35_55]
NAC_phase_angle_compensate = CB2_NAC_phase_angle_compensate[phase_angle35_55]

NAC_phase_angle = functs.Satellite_Ramping_adjust(NAC_time, NAC_phase_angle, 6.944e-3)
NAC_intensity = functs.Satellite_Ramping_adjust(NAC_time, NAC_intensity, 6.944e-3)
NAC_phase_angle_compensate = functs.Satellite_Ramping_adjust(NAC_time, NAC_phase_angle_compensate, 6.944e-3)
NAC_time = functs.Satellite_Ramping_adjust(NAC_time, NAC_time, 6.944e-3)

colourmap = ListedColormap(['red', 'orange', 'yellow', 'green', 'royalblue', 'darkviolet', 'deeppink'])
plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.colorbar()
plt.axvline(48.7625, linestyle='dashed', color='black')
plt.axvline(38.9625, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('CB2 NAC, ~752 nm')
plt.savefig("Titan_images/CB2_NAC_angle_vs_intensity_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

colourmap = ListedColormap(['red', 'orange', 'yellow', 'green', 'royalblue', 'darkviolet', 'deeppink'])
plt.scatter(NAC_time, NAC_phase_angle, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.axhline(48.7625, linestyle='dashed', color='black')
plt.axhline(38.9625, linestyle='dashed', color='black')
plt.colorbar()
plt.gca().invert_yaxis()
plt.ylabel('Phase angle [deg]')
plt.xlabel('OBS time [MJD]')
plt.title('CB2 NAC, ~752 nm')
plt.savefig("Titan_images/CB2_NAC_OBStime_vs_angle_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.colorbar()
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.title('CB2 NAC, ~752 nm')
plt.savefig("Titan_images/CB2_NAC_OBStime_vs_Intensity_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#--------------------------------
'''
plt.scatter(WAC_phase_angle, WAC_intensity, marker='o', s=20, linewidths=0.5, edgecolors='black', color='tab:blue', label='WAC')
plt.ylim(0.15, 0.25)
plt.title('CB2 (933.6 to 943.5 nm)')
plt.legend()
plt.show()

plt.scatter(NAC_phase_angle, NAC_intensity, marker='o', s=20, linewidths=0.5, edgecolors='black', color='tab:orange', label='NAC')
plt.title('CB2 (933.6 to 943.5 nm)')
plt.legend()
plt.show()
'''
'''
exp = WAC_exposure
plt.scatter(WAC_phase_angle, WAC_intensity, c=exp, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.colorbar()
plt.title('Exposure comparison')
plt.show()

gain = WAC_gain
plt.scatter(WAC_phase_angle, WAC_intensity, c=gain, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.colorbar()
plt.title('Gain comparison')
plt.show()
'''
'''
time = WAC_time
plt.scatter(WAC_phase_angle, WAC_intensity, c=time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.colorbar()
plt.title('OBS time comparison')
plt.show()
'''
'''
exp = NAC_exposure
plt.scatter(NAC_phase_angle, NAC_intensity, c=exp, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.title('Exposure comparison')
plt.colorbar()
plt.show()

gain = NAC_gain
plt.scatter(NAC_phase_angle, NAC_intensity, c=gain, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.colorbar()
plt.title('Gain comparison')
plt.show()

from matplotlib.colors import ListedColormap
colourmap = ListedColormap(['red', 'orange', 'yellow', 'green', 'royalblue', 'darkviolet'])
time = NAC_time
plt.scatter(NAC_phase_angle, NAC_intensity, c=time, cmap=colourmap, marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.colorbar()
plt.xlabel('Phase angle [deg]')
plt.xlabel('Intensity')
plt.axvline(34)
#plt.axhline(i1)
plt.axvline(37.5)
#plt.axhline(i2)
plt.title('OBS time comparison')
plt.savefig("Titan_images/CB2_NAC_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''

'''
time = NAC_time
plt.scatter(NAC_time, NAC_phase_angle, c=time, cmap='Set1', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
#plt.scatter(WAC_time, WAC_phase_angle, marker='o', s=20, linewidths=0.5, edgecolors='black', color='tab:blue', label='WAC')
plt.colorbar()
plt.ylim(29, 51)
plt.xlim(54000, 58000)
plt.xlabel('Obs. time')
plt.ylabel('Phase angle')
plt.title('OBS time')
plt.legend()
#plt.savefig("Titan_images/pictures/NAC_OBS_time.png", dpi = 120, bbox_inches='tight')
plt.show()
'''

#---------------------------------------------------------

def linear_fit(x, m, c):
    return (x*m + c)
'''
params = np.polyfit(NAC_phase_angle, NAC_intensity, 1)
print(params)
y_norm = linear_fit(NAC_phase_angle, *params)
NAC_intensity_norm = (NAC_intensity / y_norm) - 1

time = NAC_time
plt.scatter(NAC_phase_angle, NAC_intensity, c=time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.plot(NAC_phase_angle, y_norm)
plt.colorbar()
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.title('OBS time comparison')
plt.legend()
#plt.savefig("Titan_images/pictures/NAC_linear_fit2.png", dpi = 120, bbox_inches='tight')
plt.show()


params = np.polyfit(WAC_phase_angle, WAC_intensity, 1)
print(params)
y_norm = linear_fit(WAC_phase_angle, *params)
WAC_intensity_norm = (WAC_intensity / y_norm) - 1


time = WAC_time
plt.scatter(WAC_phase_angle, WAC_intensity, c=time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
plt.plot(WAC_phase_angle, y_norm)
plt.colorbar()
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
plt.title('OBS time comparison')
plt.legend()
#plt.savefig("Titan_images/pictures/WAC_linear_fit2.png", dpi = 120, bbox_inches='tight')
plt.show()

plt.scatter(WAC_phase_angle, WAC_intensity_norm, marker='o', s=20, linewidths=0.5, edgecolors='black', color='tab:blue', label='WAC')
plt.scatter(NAC_phase_angle, NAC_intensity_norm, marker='o', s=20, linewidths=0.5, edgecolors='black', color='tab:orange', label='NAC')
plt.xlabel('Phase angle')
plt.ylabel('Intensity')
plt.legend()
#plt.savefig("Titan_images/pictures/Normalized_result2.png", dpi = 120, bbox_inches='tight')
plt.show()
'''

