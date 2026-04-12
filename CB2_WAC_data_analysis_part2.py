import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import functs
import CB2_WAC_data_analysis_part1
import methane_rainbow_model

WAC_images = CB2_WAC_data_analysis_part1.WAC_images
WAC_time = CB2_WAC_data_analysis_part1.WAC_time
WAC_phase_angle = CB2_WAC_data_analysis_part1.WAC_phase_angle
WAC_intensity = CB2_WAC_data_analysis_part1.WAC_intensity
WAC_distance = CB2_WAC_data_analysis_part1.WAC_distance
WAC_source_count = CB2_WAC_data_analysis_part1.WAC_source_count
WAC_exposure = CB2_WAC_data_analysis_part1.WAC_exposure
WAC_gain = CB2_WAC_data_analysis_part1.WAC_gain
WAC_target_radius = CB2_WAC_data_analysis_part1.WAC_target_radius

camera_filter = 'CB2'
camera = 'WAC'

#----------------------------------------------------
# Identified large intensity shift in trend over the course of Cassini's mission. 
# Here I adjust the later mission Intensity vs OBS time to fit the early mission trend

def Lamberitan_curve_fit(x, A):
    lambertian_curve = A * functs.Lambertian_reflector_phase_function(x)
    return lambertian_curve

def linear_fit(x, m, c):
    return (x*m + c)

def two_orders_poly_fit(x, a1, a2, c):
    return (a1*x**2 + a2*x + c)

def three_orders_poly_fitting(x, a1, a2, a3, c):
    return (a1*x**3 + a2*x**2 + a3*x + c)

def fourth_order_poly_fitting(x, a1, a2, a3, a4, c):
    return (a1*x**4 + a2*x**3 + a3*x**2 + a4*x + c)

MissionPhase1 = (WAC_time >= 53000) & (WAC_time < 54100)
WAC_images_MissionPhase1 = WAC_images[MissionPhase1]
WAC_intensity_MissionPhase1 = WAC_intensity[MissionPhase1]
WAC_time_MissionPhase1 = WAC_time[MissionPhase1]
WAC_PhaseAngle_MissionPhase1 = WAC_phase_angle[MissionPhase1]
#print(list(zip(WAC_images_MissionPhase1, WAC_time_MissionPhase1, WAC_PhaseAngle_MissionPhase1, WAC_intensity_MissionPhase1)))

MissionPhase2 = (WAC_time >= 54100) & (WAC_time < 54700)
WAC_images_MissionPhase2 = WAC_images[MissionPhase2]
WAC_intensity_MissionPhase2 = WAC_intensity[MissionPhase2]
WAC_time_MissionPhase2 = WAC_time[MissionPhase2]
WAC_PhaseAngle_MissionPhase2 = WAC_phase_angle[MissionPhase2]

MissionPhase3 = (WAC_time >= 54700) & (WAC_time < 56000)
WAC_images_MissionPhase3 = WAC_images[MissionPhase3]
WAC_intensity_MissionPhase3 = WAC_intensity[MissionPhase3]
WAC_time_MissionPhase3 = WAC_time[MissionPhase3]
WAC_PhaseAngle_MissionPhase3 = WAC_phase_angle[MissionPhase3]

MissionPhase4 = (WAC_time >= 56000) #& (WAC_time < 57500)
WAC_images_MissionPhase4 = WAC_images[MissionPhase4]
WAC_intensity_MissionPhase4 = WAC_intensity[MissionPhase4]
WAC_time_MissionPhase4 = WAC_time[MissionPhase4]
WAC_PhaseAngle_MissionPhase4 = WAC_phase_angle[MissionPhase4]

#-------------------------------------------------------

WAC_intensity[MissionPhase1] = WAC_intensity[MissionPhase1] - 0.0
WAC_intensity[MissionPhase2] = WAC_intensity[MissionPhase2] + 0.004
WAC_intensity[MissionPhase3] = WAC_intensity[MissionPhase3] + 0.0 
WAC_intensity[MissionPhase4] = WAC_intensity[MissionPhase4] - 0.008

x_angles = np.arange(0, 90, 0.1)
lambertian_fit = 0.268*functs.Lambertian_reflector_phase_function(x_angles)

DISR_aerosol_model_xaxis = np.arange(0, 140, 10)
DISR_aerosol_model_yaxis = np.array([0.2425, 0.2325, 0.22, 0.2, 0.1825, 0.1625, 0.145, 0.1275, 0.115, 0.1025, 0.095, 0.0875, 0.0825, 0.0825]) + 0.037
DISR_param = np.polyfit(DISR_aerosol_model_xaxis, DISR_aerosol_model_yaxis, 4)
DISR_aerosol_model = fourth_order_poly_fitting(np.arange(0, 90, 0.1), *DISR_param)\

zero_lvl = (WAC_phase_angle <= 25) | (WAC_phase_angle >= 49.1625)
param = np.polyfit(WAC_phase_angle[zero_lvl], WAC_intensity[zero_lvl], 3)
y_curvefit = three_orders_poly_fitting(x_angles, *param)

#----------------------------------------------------------------------

plt.plot(x_angles, lambertian_fit, linestyle='dashed', color='tab:blue')
#plt.plot(x_angles, DISR_aerosol_model, linestyle='dashed', color='tab:green')
plt.plot(x_angles, y_curvefit, linestyle='dashed', color='tab:green')
plt.scatter(WAC_phase_angle, WAC_intensity, marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.5)
#plt.scatter(WAC_phase_angle[MissionPhase1], WAC_intensity[MissionPhase1], marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.5)
#plt.scatter(WAC_phase_angle[MissionPhase3], WAC_intensity[MissionPhase3], marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.5)
#plt.scatter(WAC_phase_angle[MissionPhase4], WAC_intensity[MissionPhase4], marker='x', s=20, color='black')
#plt.axvline(38.3625, linestyle='dashed', color='black')
#plt.axvline(49.0875, linestyle='dashed', color='black')
#plt.ylim(0.1, 0.26)
plt.xlim(-2, 82)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.savefig("Titan_images/plots/24-CB2_WAC_angle_vs_intesnity.png", dpi = 120, bbox_inches='tight')
plt.show()

plt.scatter(WAC_time, WAC_intensity, c=WAC_phase_angle, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Time', rotation=90)
plt.scatter(WAC_time[MissionPhase4], WAC_intensity[MissionPhase4], marker='x', s=20, color='black')
#plt.ylim(0.08, 0.25)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

#----------------------------------------------------------------------

WAC_intensity[MissionPhase1] = WAC_intensity[MissionPhase1] / three_orders_poly_fitting(WAC_phase_angle[MissionPhase1], *param)
WAC_intensity[MissionPhase2] = WAC_intensity[MissionPhase2] / three_orders_poly_fitting(WAC_phase_angle[MissionPhase2], *param)
WAC_intensity[MissionPhase3] = WAC_intensity[MissionPhase3] / three_orders_poly_fitting(WAC_phase_angle[MissionPhase3], *param)
WAC_intensity[MissionPhase4] = WAC_intensity[MissionPhase4] / three_orders_poly_fitting(WAC_phase_angle[MissionPhase4], *param)

angles_model = methane_rainbow_model.angles
intensity_939nm = methane_rainbow_model.intensity_939nm + 1

plt.plot(angles_model, intensity_939nm, linestyle='dashed', color='black')
plt.scatter(WAC_phase_angle, WAC_intensity, c=WAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='WAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
#plt.scatter(WAC_phase_angle[MissionPhase3], WAC_intensity[MissionPhase3], marker='x', s=20, color='black')
#plt.plot(x_angles, DISR_aerosol_model, linestyle='dashed', color='tab:blue')
#plt.ylim(0.08, 0.25)
plt.xlim(30, 50)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.savefig("Titan_images/plots/25-CB2_WAC_normalized.png", dpi = 120, bbox_inches='tight')
plt.show()


