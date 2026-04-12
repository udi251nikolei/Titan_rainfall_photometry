import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import ListedColormap
import csv

import functs
import CB3_NAC_data_analysis_Part2
import CB3_WAC_data_analysis_Part2
import methane_rainbow_model

NAC_time = CB3_NAC_data_analysis_Part2.NAC_time
NAC_phase_angle = CB3_NAC_data_analysis_Part2.NAC_phase_angle
NAC_intensity = CB3_NAC_data_analysis_Part2.NAC_intensity
NAC_distance = CB3_NAC_data_analysis_Part2.NAC_distance

WAC_time = CB3_WAC_data_analysis_Part2.WAC_time
WAC_phase_angle = CB3_WAC_data_analysis_Part2.WAC_phase_angle
WAC_intensity = CB3_WAC_data_analysis_Part2.WAC_intensity
WAC_distance = CB3_WAC_data_analysis_Part2.WAC_distance

intensity = np.hstack((NAC_intensity, WAC_intensity))
phase_angle = np.hstack((NAC_phase_angle, WAC_phase_angle))
time = np.hstack((NAC_time, WAC_time))

#--------------------------------------------------------
'''
plt.scatter(phase_angle, intensity, c=time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
#plt.ylim(1.61e18, 1.77e18)
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashed', color='black')
plt.axvline(49.0875, linestyle='dashed', color='black')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/CB3_NAC_Phase_angle_vs_intensity_vs_OBStime.png", dpi = 120, bbox_inches='tight')
plt.show()
'''
#--------------------------------------------------------
theta_intervals = np.linspace(30, 50, 11)
avg_intensity, PhaseAngle_xaxis = functs.avg_intensity_btw_PhaseAngle_interval(phase_angle, intensity, theta_intervals)

angles_model = methane_rainbow_model.angles
intensity_939nm = methane_rainbow_model.intensity_939nm + 1

plt.plot(angles_model, intensity_939nm, linestyle='solid', alpha=0.7, color='tab:blue', label='Rainbow model @ 939 nm')
plt.scatter(NAC_phase_angle, NAC_intensity, marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.7, label='NAC')
plt.scatter(WAC_phase_angle, WAC_intensity, marker='x', s=20, linewidths=0.5, color='black', alpha=1, label='WAC')
plt.plot(PhaseAngle_xaxis, avg_intensity, '-o', color='black', alpha=0.9)
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
#plt.ylim(1.61e18, 1.77e18)
#for i in theta_intervals:
#    plt.axvline(i, linestyle='dashed', color='black', alpha=0.3)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.legend(loc='upper left')
plt.savefig("Titan_images/CB3_avg_Intensity.png", dpi = 120, bbox_inches='tight')
plt.show()

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

#----------------------------------------------------
# OBStime vs. Intensity vs. Phase angle

plt.scatter(NAC_time, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(54000, 58000)
#plt.ylim(1.61e18, 1.77e18)
cbar = plt.colorbar()
cbar.ax.set_ylabel('Phase angle [deg]', rotation=90)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.title('OBStime vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_NAC_OBStime_vs_Intensity_vs_Phase_angle_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------
# Distance vs. Intensity vs. OBStime

plt.scatter(NAC_distance, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime', rotation=90)
plt.gca().invert_xaxis()
plt.ylabel('Intensity')
plt.xlabel('Distance [km]')
plt.title('Distance vs. Intensity vs. Phase angle')
#plt.savefig("Titan_images/CB3_NAC_Distance_vs_Intensity_vs_OBStime_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()

#----------------------------------------------------
# Distance vs. OBStime vs. Intensity

plt.scatter(NAC_distance, NAC_time, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Intensity', rotation=90)
plt.gca().invert_xaxis()
plt.xlabel('Distance')
plt.ylabel('OBStime')
plt.title('Distance vs. OBStime vs. Intensity')
#plt.savefig("Titan_images/CB3_NAC_Distance_vs_OBStime_vs_Intensity_colourmap.png", dpi = 120, bbox_inches='tight')
plt.show()
'''


