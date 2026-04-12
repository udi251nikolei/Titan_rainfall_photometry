import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import functs
import CB2_NAC_data_analysis_part1
import methane_rainbow_model

NAC_images = CB2_NAC_data_analysis_part1.NAC_images
NAC_time = CB2_NAC_data_analysis_part1.NAC_time
NAC_phase_angle = CB2_NAC_data_analysis_part1.NAC_phase_angle
NAC_intensity = CB2_NAC_data_analysis_part1.NAC_intensity
NAC_distance = CB2_NAC_data_analysis_part1.NAC_distance
NAC_source_count = CB2_NAC_data_analysis_part1.NAC_source_count
NAC_exposure = CB2_NAC_data_analysis_part1.NAC_exposure
NAC_gain = CB2_NAC_data_analysis_part1.NAC_gain
NAC_target_radius = CB2_NAC_data_analysis_part1.NAC_target_radius

camera_filter = 'CB2'
camera = 'NAC'

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

MissionPhase1 = (NAC_time >= 53000) & (NAC_time < 53900)
NAC_images_MissionPhase1 = NAC_images[MissionPhase1]
NAC_intensity_MissionPhase1 = NAC_intensity[MissionPhase1]
NAC_time_MissionPhase1 = NAC_time[MissionPhase1]
NAC_PhaseAngle_MissionPhase1 = NAC_phase_angle[MissionPhase1]
#print(list(zip(NAC_images_MissionPhase1, NAC_time_MissionPhase1, NAC_PhaseAngle_MissionPhase1, NAC_intensity_MissionPhase1)))

MissionPhase2 = (NAC_time >= 53900) & (NAC_time < 54900)
NAC_images_MissionPhase2 = NAC_images[MissionPhase2]
NAC_intensity_MissionPhase2 = NAC_intensity[MissionPhase2]
NAC_time_MissionPhase2 = NAC_time[MissionPhase2]
NAC_PhaseAngle_MissionPhase2 = NAC_phase_angle[MissionPhase2]

MissionPhase3 = (NAC_time >= 54900) & (NAC_time < 55800)
NAC_images_MissionPhase3 = NAC_images[MissionPhase3]
NAC_intensity_MissionPhase3 = NAC_intensity[MissionPhase3]
NAC_time_MissionPhase3 = NAC_time[MissionPhase3]
NAC_PhaseAngle_MissionPhase3 = NAC_phase_angle[MissionPhase3]

MissionPhase4 = (NAC_time >= 55800) & (NAC_time < 57000)
NAC_images_MissionPhase4 = NAC_images[MissionPhase4]
NAC_intensity_MissionPhase4 = NAC_intensity[MissionPhase4]
NAC_time_MissionPhase4 = NAC_time[MissionPhase4]
NAC_PhaseAngle_MissionPhase4 = NAC_phase_angle[MissionPhase4]

MissionPhase5 = (NAC_time >= 57000) & (NAC_time < 57500)
NAC_images_MissionPhase5 = NAC_images[MissionPhase5]
NAC_intensity_MissionPhase5 = NAC_intensity[MissionPhase5]
NAC_time_MissionPhase5 = NAC_time[MissionPhase5]
NAC_PhaseAngle_MissionPhase5 = NAC_phase_angle[MissionPhase5]

MissionPhase6 = (NAC_time >= 57500)
NAC_images_MissionPhase6 = NAC_images[MissionPhase6]
NAC_intensity_MissionPhase6 = NAC_intensity[MissionPhase6]
NAC_time_MissionPhase6 = NAC_time[MissionPhase6]
NAC_PhaseAngle_MissionPhase6 = NAC_phase_angle[MissionPhase6]

#-------------------------------------------------------

NAC_intensity[MissionPhase1] = NAC_intensity[MissionPhase1] - 0.0
NAC_intensity[MissionPhase2] = NAC_intensity[MissionPhase2] + 0.004 #within angle-of-interest
NAC_intensity[MissionPhase3] = NAC_intensity[MissionPhase3] + 0.0 #within angle-of-interest
NAC_intensity[MissionPhase4] = NAC_intensity[MissionPhase4] + 0.001 #within angle-of-interest
NAC_intensity[MissionPhase5] = NAC_intensity[MissionPhase5] - 0.004 #within angle-of-interest
NAC_intensity[MissionPhase6] = NAC_intensity[MissionPhase6] - 0.006

x_angles = np.arange(0, 90, 0.1)
lambertian_fit = 0.25*functs.Lambertian_reflector_phase_function(x_angles)

DISR_aerosol_model_xaxis = np.arange(0, 140, 10)
DISR_aerosol_model_yaxis = np.array([0.2425, 0.2325, 0.22, 0.2, 0.1825, 0.1625, 0.145, 0.1275, 0.115, 0.1025, 0.095, 0.0875, 0.0825, 0.0825]) + 0.019
DISR_param = np.polyfit(DISR_aerosol_model_xaxis, DISR_aerosol_model_yaxis, 4)
DISR_aerosol_model = fourth_order_poly_fitting(x_angles, *DISR_param)

MissionPhase3_zero = (NAC_PhaseAngle_MissionPhase3 <= 25) | (NAC_PhaseAngle_MissionPhase3 >= 49.1625)
param_MissionPhase3 = np.polyfit(NAC_PhaseAngle_MissionPhase3[MissionPhase3_zero], NAC_intensity_MissionPhase3[MissionPhase3_zero], 3)
param_MissionPhase3[-1] = param_MissionPhase3[-1] - 0.002
MissionPhase3_curvefit = three_orders_poly_fitting(np.arange(0, 90, 0.1), *param_MissionPhase3)

albedo_curvefit = (NAC_phase_angle <= 25) | (NAC_phase_angle >= 49.1625)
param = np.polyfit(NAC_phase_angle[albedo_curvefit], NAC_intensity[albedo_curvefit], 3)
y_curvefit = three_orders_poly_fitting(np.arange(0, 90, 0.1), *param)

angles_model = methane_rainbow_model.angles
intensity_939nm = methane_rainbow_model.intensity_939nm + 1
#print(list(zip(angles_model, intensity_939nm)))

#----------------------------------------------------------------------

plt.plot(x_angles, lambertian_fit, linestyle='dashed', color='tab:blue')
#plt.plot(x_angles, DISR_aerosol_model, linestyle='dashed', color='tab:orange')
plt.plot(x_angles, MissionPhase3_curvefit, linestyle='dashed', color='tab:green')
plt.scatter(NAC_phase_angle, NAC_intensity, marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.5)
#plt.scatter(NAC_phase_angle[MissionPhase6], NAC_intensity[MissionPhase6], marker='x', s=20, color='black')
#plt.axvline(38.3625, linestyle='dashed', color='black')
#plt.axvline(49.0875, linestyle='dashed', color='black')
plt.ylim(0.1, 0.28)
plt.xlim(-2, 82)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/plots/22-CB2_NAC_angle_vs_intesnity.png", dpi = 120, bbox_inches='tight')
plt.show()

plt.scatter(NAC_time, NAC_intensity, c=NAC_phase_angle, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Time', rotation=90)
plt.scatter(NAC_time[MissionPhase6], NAC_intensity[MissionPhase6], marker='x', s=20, color='black')
#plt.ylim(0.08, 0.25)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
plt.show()

#----------------------------------------------------------------------

lambertian_fit = 0.25*functs.Lambertian_reflector_phase_function(NAC_phase_angle)
DISR_aerosol_model = fourth_order_poly_fitting(NAC_phase_angle, *DISR_param)
MissionPhase3_curvefit = three_orders_poly_fitting(NAC_phase_angle, *param_MissionPhase3)

NAC_intensity = NAC_intensity / MissionPhase3_curvefit

'''
NAC_intensity[MissionPhase1] = NAC_intensity[MissionPhase1] / three_orders_poly_fitting(NAC_phase_angle[MissionPhase1], *param_MissionPhase3)
NAC_intensity[MissionPhase2] = NAC_intensity[MissionPhase2] / three_orders_poly_fitting(NAC_phase_angle[MissionPhase2], *param_MissionPhase3)
NAC_intensity[MissionPhase3] = NAC_intensity[MissionPhase3] / three_orders_poly_fitting(NAC_phase_angle[MissionPhase3], *param_MissionPhase3)
NAC_intensity[MissionPhase4] = NAC_intensity[MissionPhase4] / three_orders_poly_fitting(NAC_phase_angle[MissionPhase4], *param_MissionPhase3)
NAC_intensity[MissionPhase5] = NAC_intensity[MissionPhase5] / three_orders_poly_fitting(NAC_phase_angle[MissionPhase5], *param_MissionPhase3)
NAC_intensity[MissionPhase6] = NAC_intensity[MissionPhase6] / three_orders_poly_fitting(NAC_phase_angle[MissionPhase6], *param_MissionPhase3)
'''

angles_model = methane_rainbow_model.angles
intensity_939nm = methane_rainbow_model.intensity_939nm + 1

plt.plot(angles_model, intensity_939nm, linestyle='dashed', color='black')
plt.scatter(NAC_phase_angle, NAC_intensity, c=NAC_time, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
#plt.scatter(NAC_phase_angle[MissionPhase3], NAC_intensity[MissionPhase3], marker='x', s=20, color='black')
#plt.plot(x_angles, DISR_aerosol_model, linestyle='dashed', color='tab:blue')
#plt.ylim(0.08, 0.25)
plt.xlim(30, 50)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/plots/23-CB2_NAC_normalized.png", dpi = 120, bbox_inches='tight')
plt.show()
