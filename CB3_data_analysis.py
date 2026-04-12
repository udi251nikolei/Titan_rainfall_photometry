import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import functs
import CB3_NAC_data_analysis_Part1
import CB3_WAC_data_analysis_Part1
import methane_rainbow_model

NAC_images = CB3_NAC_data_analysis_Part1.NAC_images
NAC_time = CB3_NAC_data_analysis_Part1.NAC_time
NAC_phase_angle = CB3_NAC_data_analysis_Part1.NAC_phase_angle
NAC_intensity = CB3_NAC_data_analysis_Part1.NAC_intensity
NAC_distance = CB3_NAC_data_analysis_Part1.NAC_distance
NAC_source_count = CB3_NAC_data_analysis_Part1.NAC_source_count
NAC_exposure = CB3_NAC_data_analysis_Part1.NAC_exposure
NAC_gain = CB3_NAC_data_analysis_Part1.NAC_gain
NAC_target_radius = CB3_NAC_data_analysis_Part1.NAC_target_radius

NAC_intensity = NAC_intensity + 0.014

WAC_images = CB3_WAC_data_analysis_Part1.WAC_images
WAC_time = CB3_WAC_data_analysis_Part1.WAC_time
WAC_phase_angle = CB3_WAC_data_analysis_Part1.WAC_phase_angle
WAC_intensity = CB3_WAC_data_analysis_Part1.WAC_intensity
WAC_distance = CB3_WAC_data_analysis_Part1.WAC_distance
WAC_source_count = CB3_WAC_data_analysis_Part1.WAC_source_count
WAC_exposure = CB3_WAC_data_analysis_Part1.WAC_exposure
WAC_gain = CB3_WAC_data_analysis_Part1.WAC_gain
WAC_target_radius = CB3_WAC_data_analysis_Part1.WAC_target_radius

WAC_intensity = WAC_intensity - 0.009

#----------------------------------------------------
# Identified large intensity shift in trend over the course of Cassini's mission. 
# Here I adjust the later mission Intensity vs OBS time to fit the early mission trend

def linear_fit(x, m, c):
    return (x*m + c)

def two_orders_poly_fit(x, a1, a2, c):
    return (a1*x**2 + a2*x + c)

def three_orders_poly_fit(x, a1, a2, a3, c):
    return (a1*x**3 + a2*x**2 + a3*x + c)

def four_orders_poly_fit(x, a1, a2, a3, a4, c):
    return (a1*x**4 + a2*x**3 + a3*x**2 + a4*x + c)

def five_orders_poly_fit(x, a1, a2, a3, a4, a5, c):
    return (a1*x**5 + a2*x**4 + a3*x**3 + a4*x**2 + a5*x + c)

#--------------------------------------------------------------------
# For NAC
NAC_Phase1 = (NAC_time < 53550)
NAC_Phase2 = (NAC_time >= 53550) & (NAC_time < 54300)
NAC_Phase3 = (NAC_time >= 54300) & (NAC_time < 54600)
NAC_Phase4 = (NAC_time >= 54600) & (NAC_time < 54800)
NAC_Phase5 = (NAC_time >= 54800) & (NAC_time < 55100)
NAC_Phase6 = (NAC_time >= 55100) & (NAC_time < 55600)
NAC_Phase7 = (NAC_time >= 55600) & (NAC_time < 55800)
NAC_Phase8 = (NAC_time >= 55800) & (NAC_time < 56000)
NAC_Phase9 = (NAC_time >= 56000) & (NAC_time < 56500)
NAC_Phase10 = (NAC_time >= 56500) & (NAC_time < 57200)
NAC_Phase11 = (NAC_time >= 57200) & (NAC_time < 58000)

NAC_intensity[NAC_Phase1] = NAC_intensity[NAC_Phase1] - 0.006
NAC_intensity[NAC_Phase2] = NAC_intensity[NAC_Phase2] + 0.005
NAC_intensity[NAC_Phase3] = NAC_intensity[NAC_Phase3] + 0.001
NAC_intensity[NAC_Phase4] = NAC_intensity[NAC_Phase4] + 0.001
NAC_intensity[NAC_Phase5] = NAC_intensity[NAC_Phase5] - 0.003
NAC_intensity[NAC_Phase6] = NAC_intensity[NAC_Phase6] + 0.00
NAC_intensity[NAC_Phase7] = NAC_intensity[NAC_Phase7] + 0.002
NAC_intensity[NAC_Phase8] = NAC_intensity[NAC_Phase8] + 0.001
NAC_intensity[NAC_Phase9] = NAC_intensity[NAC_Phase9] - 0.00
NAC_intensity[NAC_Phase10] = NAC_intensity[NAC_Phase10] - 0.005
NAC_intensity[NAC_Phase11] = NAC_intensity[NAC_Phase11] - 0.004

#RoI = (NAC_phase_angle[NAC_Phase6] >= 24.1) & (NAC_phase_angle[NAC_Phase6] <= 49.1625)
#print(list(zip(np.around(NAC_time[NAC_Phase6][RoI], 3), NAC_phase_angle[NAC_Phase6][RoI], np.around(NAC_target_radius[NAC_Phase6][RoI], 3), NAC_exposure[NAC_Phase6][RoI], NAC_gain[NAC_Phase6][RoI], np.around(NAC_intensity[NAC_Phase6][RoI], 5))))

# For WAC
WAC_Phase1 = (WAC_time >= 53300) & (WAC_time < 53900)
WAC_Phase2 = (WAC_time >= 53900) & (WAC_time < 54300)
WAC_Phase3 = (WAC_time >= 54300) & (WAC_time < 54600)
WAC_prePhase4 = (WAC_time >= 54600) & (WAC_time < 54700)
WAC_Phase4 = (WAC_time >= 54700) & (WAC_time < 54800)
WAC_Phase5 = (WAC_time >= 54800) & (WAC_time < 55100)
WAC_Phase6 = (WAC_time >= 55100) & (WAC_time < 55500)
WAC_Phase7 = (WAC_time >= 55500) & (WAC_time < 55800)
WAC_Phase8 = (WAC_time >= 55800) & (WAC_time < 56400)
WAC_Phase9 = (WAC_time >= 56400) & (WAC_time < 57000)                                      
WAC_Phase10 = (WAC_time >= 57000) & (WAC_time < 57500)
WAC_Phase11 = (WAC_time >= 57700)

WAC_intensity[WAC_Phase1] = WAC_intensity[WAC_Phase1] + 0.00
WAC_intensity[WAC_Phase2] = WAC_intensity[WAC_Phase2] + 0.01
WAC_intensity[WAC_Phase3] = WAC_intensity[WAC_Phase3] + 0.009
WAC_intensity[WAC_prePhase4] = WAC_intensity[WAC_prePhase4] + 0.01
WAC_intensity[WAC_Phase4] = WAC_intensity[WAC_Phase4] + 0.015
WAC_intensity[WAC_Phase5] = WAC_intensity[WAC_Phase5] + 0.011
WAC_intensity[WAC_Phase6] = WAC_intensity[WAC_Phase6] + 0.01
WAC_intensity[WAC_Phase7] = WAC_intensity[WAC_Phase7] + 0.009
WAC_intensity[WAC_Phase8] = WAC_intensity[WAC_Phase8] + 0.012
WAC_intensity[WAC_Phase9] = WAC_intensity[WAC_Phase9] + 0.0
WAC_intensity[WAC_Phase10] = WAC_intensity[WAC_Phase10] + 0.00
WAC_intensity[WAC_Phase11] = WAC_intensity[WAC_Phase11] - 0.004

#RoI = (WAC_phase_angle[WAC_Phase2] >= 24.1) & (WAC_phase_angle[WAC_Phase2] <= 49.1625)
#print(list(zip(np.around(WAC_time[WAC_Phase2][RoI], 3), WAC_phase_angle[WAC_Phase2][RoI], np.around(WAC_target_radius[WAC_Phase2][RoI], 3), WAC_exposure[WAC_Phase2][RoI], WAC_gain[WAC_Phase2][RoI], np.around(WAC_intensity[WAC_Phase2][RoI], 5))))

#--------------------------------------------------------------------

x_angles = np.linspace(0, 90, 9001)
lambertian_fit = functs.Lambertian_reflector_phase_function(x_angles, 0.2275)

DISR_aerosol_model_xaxis = np.arange(0, 140, 10)
DISR_aerosol_model_yaxis = np.array([0.2425, 0.2325, 0.22, 0.2, 0.1825, 0.1625, 0.145, 0.1275, 0.115, 0.1025, 0.095, 0.0875, 0.0825, 0.0825])
DISR_param = np.polyfit(DISR_aerosol_model_xaxis, DISR_aerosol_model_yaxis, 4)
DISR_aerosol_model = four_orders_poly_fit(x_angles, *DISR_param)

# NAC Curvefits
Phase3_startfit = functs.Lambertian_reflector_phase_function(np.arange(0, 10, 1), 0.2275) + 0.002
intensity_NAC_Phase2_3 = np.hstack((Phase3_startfit, NAC_intensity[NAC_Phase3], four_orders_poly_fit(np.arange(55, 85, 1), *DISR_param)))
PhaseAngle_NAC_Phase2_3 = np.hstack((np.arange(0, 10, 1), NAC_phase_angle[NAC_Phase3], np.arange(55, 85, 1)))
surface_albedo_NAC_Phase2_3 = (PhaseAngle_NAC_Phase2_3 <= 24.1) | (PhaseAngle_NAC_Phase2_3 >= 49.1625)
param_NAC_Phase2_3 = np.polyfit(PhaseAngle_NAC_Phase2_3[surface_albedo_NAC_Phase2_3], intensity_NAC_Phase2_3[surface_albedo_NAC_Phase2_3], 3)
curvefit_NAC_Phase2_3 = three_orders_poly_fit(x_angles, *param_NAC_Phase2_3)

intensity_NAC_Phase4_5 = np.hstack((NAC_intensity[NAC_Phase4], NAC_intensity[NAC_Phase5]))
PhaseAngle_NAC_Phase4_5 = np.hstack((NAC_phase_angle[NAC_Phase4], NAC_phase_angle[NAC_Phase5]))
surface_albedo_NAC_Phase4_5 = (PhaseAngle_NAC_Phase4_5 <= 24.1) | (PhaseAngle_NAC_Phase4_5 >= 49.1625)
param_NAC_Phase4_5 = np.polyfit(PhaseAngle_NAC_Phase4_5[surface_albedo_NAC_Phase4_5], intensity_NAC_Phase4_5[surface_albedo_NAC_Phase4_5], 3)
curvefit_NAC_Phase4_5 = three_orders_poly_fit(x_angles, *param_NAC_Phase4_5)

Phase6_startfit = three_orders_poly_fit(np.arange(0, 20, 2), *param_NAC_Phase4_5)
intensity_NAC_Phase6 = np.hstack((Phase6_startfit, NAC_intensity[NAC_Phase6]))
PhaseAngle_NAC_Phase6 = np.hstack((np.arange(0, 20, 2), NAC_phase_angle[NAC_Phase6]))
surface_albedo_NAC_Phase6 = (PhaseAngle_NAC_Phase6 <= 24.1) | (PhaseAngle_NAC_Phase6 >= 49.1625)
param_NAC_Phase6 = np.polyfit(PhaseAngle_NAC_Phase6[surface_albedo_NAC_Phase6], intensity_NAC_Phase6[surface_albedo_NAC_Phase6], 3)
curvefit_NAC_Phase6 = three_orders_poly_fit(x_angles, *param_NAC_Phase6)

intensity_NAC_Phase7 = NAC_intensity[NAC_Phase7]
PhaseAngle_NAC_Phase7 = NAC_phase_angle[NAC_Phase7]
surface_albedo_NAC_Phase7 = (PhaseAngle_NAC_Phase7 <= 24.1) | (PhaseAngle_NAC_Phase7 >= 49.1625)
param_NAC_Phase7 = np.polyfit(PhaseAngle_NAC_Phase7[surface_albedo_NAC_Phase7], intensity_NAC_Phase7[surface_albedo_NAC_Phase7], 3)
curvefit_NAC_Phase7 = three_orders_poly_fit(x_angles, *param_NAC_Phase7)

intensity_NAC_Phase8 = NAC_intensity[NAC_Phase8]
PhaseAngle_NAC_Phase8 = NAC_phase_angle[NAC_Phase8]
surface_albedo_NAC_Phase8 = (PhaseAngle_NAC_Phase8 <= 24.1) | (PhaseAngle_NAC_Phase8 >= 49.1625)
param_NAC_Phase8 = np.polyfit(PhaseAngle_NAC_Phase8[surface_albedo_NAC_Phase8], intensity_NAC_Phase8[surface_albedo_NAC_Phase8], 3)
curvefit_NAC_Phase8 = three_orders_poly_fit(x_angles, *param_NAC_Phase8)

intensity_NAC_Phase9 = NAC_intensity[NAC_Phase9]
PhaseAngle_NAC_Phase9 = NAC_phase_angle[NAC_Phase9]
surface_albedo_NAC_Phase9 = (PhaseAngle_NAC_Phase9 <= 24.1) | (PhaseAngle_NAC_Phase9 >= 49.1625)
param_NAC_Phase9 = np.polyfit(PhaseAngle_NAC_Phase9[surface_albedo_NAC_Phase9], intensity_NAC_Phase9[surface_albedo_NAC_Phase9], 3)
curvefit_NAC_Phase9 = three_orders_poly_fit(x_angles, *param_NAC_Phase9)

intensity_NAC_Phase10 = NAC_intensity[NAC_Phase10]
PhaseAngle_NAC_Phase10 = NAC_phase_angle[NAC_Phase10]
surface_albedo_NAC_Phase10 = (PhaseAngle_NAC_Phase10 <= 24.1) | (PhaseAngle_NAC_Phase10 >= 49.1625)
param_NAC_Phase10 = np.polyfit(PhaseAngle_NAC_Phase10[surface_albedo_NAC_Phase10], intensity_NAC_Phase10[surface_albedo_NAC_Phase10], 3)
curvefit_NAC_Phase10 = three_orders_poly_fit(x_angles, *param_NAC_Phase10)

intensity_NAC_Phase11 = NAC_intensity[NAC_Phase11]
PhaseAngle_NAC_Phase11 = NAC_phase_angle[NAC_Phase11]
surface_albedo_NAC_Phase11 = (PhaseAngle_NAC_Phase11 <= 24.1) | (PhaseAngle_NAC_Phase11 >= 49.1625)
param_NAC_Phase11 = np.polyfit(PhaseAngle_NAC_Phase11[surface_albedo_NAC_Phase11], intensity_NAC_Phase11[surface_albedo_NAC_Phase11], 3)
curvefit_NAC_Phase11 = three_orders_poly_fit(x_angles, *param_NAC_Phase11)

intensity_NAC_Phase10_11 = np.hstack((NAC_intensity[NAC_Phase10], NAC_intensity[NAC_Phase11]))
PhaseAngle_NAC_Phase10_11 = np.hstack((NAC_phase_angle[NAC_Phase10], NAC_phase_angle[NAC_Phase11]))
surface_albedo_NAC_Phase10_11 = (PhaseAngle_NAC_Phase10_11 <= 24.1) | (PhaseAngle_NAC_Phase10_11 >= 49.1625)
param_NAC_Phase10_11 = np.polyfit(PhaseAngle_NAC_Phase10_11[surface_albedo_NAC_Phase10_11], intensity_NAC_Phase10_11[surface_albedo_NAC_Phase10_11], 3)
curvefit_NAC_Phase10_11 = three_orders_poly_fit(x_angles, *param_NAC_Phase10_11)


# WAC Curvefits
intensity_NAC_WAC_Phase1 = np.hstack((NAC_intensity[NAC_Phase1], WAC_intensity[WAC_Phase1]))
PhaseAngle_NAC_WAC_Phase1 = np.hstack((NAC_phase_angle[NAC_Phase1], WAC_phase_angle[WAC_Phase1]))
surface_albedo_NAC_WAC_Phase1 = (PhaseAngle_NAC_WAC_Phase1 <= 24.1) | (PhaseAngle_NAC_WAC_Phase1 >= 49.1625)
param_NAC_WAC_Phase1 = np.polyfit(PhaseAngle_NAC_WAC_Phase1[surface_albedo_NAC_WAC_Phase1], intensity_NAC_WAC_Phase1[surface_albedo_NAC_WAC_Phase1], 3)
curvefit_NAC_WAC_Phase1 = three_orders_poly_fit(x_angles, *param_NAC_WAC_Phase1)

phase4_5_6_endfit = four_orders_poly_fit(np.arange(85, 90, 1), *DISR_param)
intensity_WAC_Phase4_5_6 = np.hstack((WAC_intensity[WAC_Phase4], WAC_intensity[WAC_Phase5], WAC_intensity[WAC_Phase6], phase4_5_6_endfit))
PhaseAngle_WAC_Phase4_5_6 = np.hstack((WAC_phase_angle[WAC_Phase4], WAC_phase_angle[WAC_Phase5], WAC_phase_angle[WAC_Phase6], np.arange(85, 90, 1)))
surface_albedo_WAC_Phase4_5_6 = (PhaseAngle_WAC_Phase4_5_6 <= 24.1) | (PhaseAngle_WAC_Phase4_5_6 >= 49.1625)
param_WAC_Phase4_5_6 = np.polyfit(PhaseAngle_WAC_Phase4_5_6[surface_albedo_WAC_Phase4_5_6], intensity_WAC_Phase4_5_6[surface_albedo_WAC_Phase4_5_6], 3)
curvefit_WAC_Phase4_5_6 = three_orders_poly_fit(x_angles, *param_WAC_Phase4_5_6)

phase7_9_startfit = four_orders_poly_fit(np.arange(0, 16, 3), *DISR_param) - 0.005
intensity_WAC_Phase7_9 = np.hstack((phase7_9_startfit, WAC_intensity[WAC_Phase7][0:18], WAC_intensity[WAC_Phase7][27], WAC_intensity[WAC_Phase7][29:-11], WAC_intensity[WAC_Phase7][-7:-1], WAC_intensity[WAC_Phase9]))
PhaseAngle_WAC_Phase7_9 = np.hstack((np.arange(0, 16, 3), WAC_phase_angle[WAC_Phase7][0:18], WAC_phase_angle[WAC_Phase7][27], WAC_phase_angle[WAC_Phase7][29:-11], WAC_phase_angle[WAC_Phase7][-7:-1], WAC_phase_angle[WAC_Phase9]))
surface_albedo_WAC_Phase7_9 = (PhaseAngle_WAC_Phase7_9 <= 24.1) | (PhaseAngle_WAC_Phase7_9 >= 49.1625)
param_WAC_Phase7_9 = np.polyfit(PhaseAngle_WAC_Phase7_9[surface_albedo_WAC_Phase7_9], intensity_WAC_Phase7_9[surface_albedo_WAC_Phase7_9], 3)
curvefit_WAC_Phase7_9 = three_orders_poly_fit(x_angles, *param_WAC_Phase7_9)

intensity_WAC_Phase2_3_8_10 = np.hstack((WAC_intensity[WAC_Phase2][0:-8], WAC_intensity[WAC_Phase2][-6:-1], WAC_intensity[WAC_Phase3], WAC_intensity[WAC_Phase8][0:10], WAC_intensity[WAC_Phase8][21:27], WAC_intensity[WAC_Phase8][29:-12], WAC_intensity[WAC_Phase8][-8:-1], WAC_intensity[WAC_Phase10]))
PhaseAngle_WAC_Phase2_3_8_10 = np.hstack((WAC_phase_angle[WAC_Phase2][0:-8], WAC_phase_angle[WAC_Phase2][-6:-1], WAC_phase_angle[WAC_Phase3], WAC_phase_angle[WAC_Phase8][0:10], WAC_phase_angle[WAC_Phase8][21:27], WAC_phase_angle[WAC_Phase8][29:-12], WAC_phase_angle[WAC_Phase8][-8:-1], WAC_phase_angle[WAC_Phase10]))
surface_albedo_WAC_Phase2_3_8_10 = (PhaseAngle_WAC_Phase2_3_8_10 <= 24.1) | (PhaseAngle_WAC_Phase2_3_8_10 >= 49.1625)
param_WAC_Phase2_3_8_10 = np.polyfit(PhaseAngle_WAC_Phase2_3_8_10[surface_albedo_WAC_Phase2_3_8_10], intensity_WAC_Phase2_3_8_10[surface_albedo_WAC_Phase2_3_8_10], 3)
curvefit_WAC_Phase2_3_8_10 = three_orders_poly_fit(x_angles, *param_WAC_Phase2_3_8_10)

phase11_startfit = four_orders_poly_fit(np.arange(0, 18, 3), *DISR_param) + 0.006
intensity_WAC_Phase11 = np.hstack((phase11_startfit, WAC_intensity[WAC_Phase11][:15], WAC_intensity[WAC_Phase11][15:-26], WAC_intensity[WAC_Phase11][-24:-5]))
PhaseAngle_WAC_Phase11 = np.hstack((np.arange(0, 18, 3), WAC_phase_angle[WAC_Phase11][:15], WAC_phase_angle[WAC_Phase11][15:-26], WAC_phase_angle[WAC_Phase11][-24:-5]))
surface_albedo_WAC_Phase11 = (PhaseAngle_WAC_Phase11 <= 24.1) | (PhaseAngle_WAC_Phase11 >= 49.1625)
param_WAC_Phase11 = np.polyfit(PhaseAngle_WAC_Phase11[surface_albedo_WAC_Phase11], intensity_WAC_Phase11[surface_albedo_WAC_Phase11], 3)
curvefit_WAC_Phase11 = three_orders_poly_fit(x_angles, *param_WAC_Phase11)

#-----------------------------------------------------------------------------------------
# Identifying data points whose day of detection aligns with the day that clouds were detected on Titan by ISS and/or VIMS

time_cloud_detection = [54113, 54129, 54153, 54156, 54169, 
                        54185, 54200, 54216, 54232, 54233, 
                        54248, 54264, 54280, 54281, 54300, 
                        54343, 54345, 54375, 54423, 54439, 
                        54454, 54470, 54485, 54518, 54551, 
                        54553, 54582, 54598, 54614, 54615, 
                        54678, 54693, 54739, 54753, 54773, 
                        54783, 54789, 54805, 54806, 54821, 
                        54869, 54877, 54916, 54917, 54925, 
                        54941, 54949, 54956, 54972, 54979, 
                        54988, 55004, 55005, 55020, 55021, 
                        55068, 55097, 55116, 55150, 55177,
                        55178, 55193, 55195, 55196, 55208,
                        55210, 55224, 55243, 55260, 55277, 
                        55291, 55352, 55360, 55368, 55384, 
                        55425, 55449, 55463, 55466, 55483, 
                        55487, 55498, 55511, 55623, 55689, 
                        55732, 55756, 55908, 55956, 55976, 
                        55988, 55995, 56014, 56170, 56183, 
                        56196, 56202, 56244, 56435, 56499, 
                        56547, 56794, 56826, 56858, 56859, 
                        56890, 56903, 56954, 57097, 57149, 
                        57210, 57293, 57339, 57386, 57403, 
                        57419, 57434, 57453, 57468, 57482, 
                        57497, 57499, 57501, 57505, 57508, 
                        57514, 57529, 57546, 57547, 57578, 
                        57580, 57585, 57594, 57610, 57658, 
                        57688, 57690, 57691, 57706, 57721, 
                        57740, 57752, 57753, 57767, 57801, 
                        57816, 57817, 57829, 57832, 57833, 
                        57846, 57861, 57865, 57880, 57912, 
                        57945, 57958, 57960, 57963, 57976, 
                        57995, 57997, 58007, 58009]

NAC_cloud_detection_Phase1 = np.isin(NAC_time[NAC_Phase1].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase2 = np.isin(NAC_time[NAC_Phase2].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase3 = np.isin(NAC_time[NAC_Phase3].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase4 = np.isin(NAC_time[NAC_Phase4].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase5 = np.isin(NAC_time[NAC_Phase5].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase6 = np.isin(NAC_time[NAC_Phase6].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase7 = np.isin(NAC_time[NAC_Phase7].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase8 = np.isin(NAC_time[NAC_Phase8].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase9 = np.isin(NAC_time[NAC_Phase9].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase10 = np.isin(NAC_time[NAC_Phase10].astype(int), time_cloud_detection)
NAC_cloud_detection_Phase11 = np.isin(NAC_time[NAC_Phase11].astype(int), time_cloud_detection)

WAC_cloud_detection_Phase1 = np.isin(WAC_time[WAC_Phase1].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase2 = np.isin(WAC_time[WAC_Phase2].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase3 = np.isin(WAC_time[WAC_Phase3].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase4 = np.isin(WAC_time[WAC_Phase4].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase5 = np.isin(WAC_time[WAC_Phase5].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase6 = np.isin(WAC_time[WAC_Phase6].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase7 = np.isin(WAC_time[WAC_Phase7].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase8 = np.isin(WAC_time[WAC_Phase8].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase9 = np.isin(WAC_time[WAC_Phase9].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase10 = np.isin(WAC_time[WAC_Phase10].astype(int), time_cloud_detection)
WAC_cloud_detection_Phase11 = np.isin(WAC_time[WAC_Phase11].astype(int), time_cloud_detection)

#-----------------------------------------------------------------------------------------

plt.plot(x_angles, lambertian_fit, linestyle='dashed', color='tab:blue')
plt.plot(x_angles, DISR_aerosol_model, linestyle='dashed', color='tab:orange')
plt.plot(x_angles, curvefit_NAC_Phase6, linestyle='dashed', color='tab:red')
plt.scatter(NAC_phase_angle[NAC_Phase6], NAC_intensity[NAC_Phase6], marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.5)
#plt.scatter(NAC_phase_angle[NAC_Phase3], NAC_intensity[NAC_Phase3], marker='o', s=20, linewidths=0.5, edgecolors='black', color='gray', alpha=0.5)
plt.axvline(24.1, linestyle='dashed', color='black')
plt.axvline(49.1625, linestyle='dashed', color='black')
plt.ylim(0.08, 0.27)
plt.xlim(-2, 92)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/plots/NAC/PhaseAngle-vs-Intenisty_NAC_Phase8.png", dpi = 120, bbox_inches='tight')
plt.show()

#plt.scatter(NAC_time, NAC_intensity, c=NAC_phase_angle, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.scatter(NAC_time, NAC_intensity, c=NAC_phase_angle, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
cbar = plt.colorbar()
cbar.ax.set_ylabel('Time', rotation=90)
plt.scatter(NAC_time[NAC_Phase6], NAC_intensity[NAC_Phase6], marker='x', s=20, color='black')
#plt.scatter(NAC_time[NAC_Phase3], NAC_intensity[NAC_Phase3], marker='x', s=20, color='black')
plt.ylim(0.08, 0.26)
plt.ylabel('Intensity')
plt.xlabel('OBS time [MJD]')
#plt.savefig("Titan_images/plots/NAC/OBStime-vs-Intenisty_NAC_Phase8.png", dpi = 120, bbox_inches='tight')
plt.show()

#-----------------------------------------------------------------------------------------

NAC_intensity[NAC_Phase1] = NAC_intensity[NAC_Phase1] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase1], *param_NAC_WAC_Phase1)
NAC_intensity[NAC_Phase2] = NAC_intensity[NAC_Phase2] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase2], *param_NAC_Phase2_3)
NAC_intensity[NAC_Phase3] = NAC_intensity[NAC_Phase3] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase3], *param_NAC_Phase2_3)
NAC_intensity[NAC_Phase4] = NAC_intensity[NAC_Phase4] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase4], *param_NAC_Phase4_5)
NAC_intensity[NAC_Phase5] = NAC_intensity[NAC_Phase5] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase5], *param_NAC_Phase4_5)
NAC_intensity[NAC_Phase6] = NAC_intensity[NAC_Phase6] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase6], *param_NAC_Phase6)
NAC_intensity[NAC_Phase7] = NAC_intensity[NAC_Phase7] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase7], *param_NAC_Phase7)
NAC_intensity[NAC_Phase8] = NAC_intensity[NAC_Phase8] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase8], *param_NAC_Phase8)
NAC_intensity[NAC_Phase9] = NAC_intensity[NAC_Phase9] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase9], *param_NAC_Phase9)
NAC_intensity[NAC_Phase10] = NAC_intensity[NAC_Phase10] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase10], *param_NAC_Phase10_11)
NAC_intensity[NAC_Phase11] = NAC_intensity[NAC_Phase11] / three_orders_poly_fit(NAC_phase_angle[NAC_Phase11], *param_NAC_Phase10_11)

WAC_intensity[WAC_Phase1] = WAC_intensity[WAC_Phase1] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase1], *param_NAC_WAC_Phase1)
WAC_intensity[WAC_Phase2] = WAC_intensity[WAC_Phase2] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase2], *param_WAC_Phase2_3_8_10)
WAC_intensity[WAC_Phase3] = WAC_intensity[WAC_Phase3] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase3], *param_WAC_Phase2_3_8_10)
WAC_intensity[WAC_prePhase4] = WAC_intensity[WAC_prePhase4] / three_orders_poly_fit(WAC_phase_angle[WAC_prePhase4], *param_WAC_Phase4_5_6)
WAC_intensity[WAC_Phase4] = WAC_intensity[WAC_Phase4] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase4], *param_WAC_Phase4_5_6)
WAC_intensity[WAC_Phase5] = WAC_intensity[WAC_Phase5] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase5], *param_WAC_Phase4_5_6)
WAC_intensity[WAC_Phase6] = WAC_intensity[WAC_Phase6] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase6], *param_WAC_Phase4_5_6)
WAC_intensity[WAC_Phase7] = WAC_intensity[WAC_Phase7] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase7], *param_WAC_Phase7_9)
WAC_intensity[WAC_Phase8] = WAC_intensity[WAC_Phase8] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase8], *param_WAC_Phase2_3_8_10)
WAC_intensity[WAC_Phase9] = WAC_intensity[WAC_Phase9] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase9], *param_WAC_Phase7_9)
WAC_intensity[WAC_Phase10] = WAC_intensity[WAC_Phase10] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase10], *param_WAC_Phase2_3_8_10)
WAC_intensity[WAC_Phase11] = WAC_intensity[WAC_Phase11] / three_orders_poly_fit(WAC_phase_angle[WAC_Phase11], *param_WAC_Phase11)

angles_model = methane_rainbow_model.angles
intensity_939nm = methane_rainbow_model.intensity_939nm + 1

plt.plot(angles_model, intensity_939nm, linestyle='solid', color='black')
plt.axhline(1, linestyle='dashed', linewidth=0.5, color='black')
plt.scatter(NAC_phase_angle[NAC_Phase6], NAC_intensity[NAC_Phase6], c=NAC_time[NAC_Phase6], cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
#plt.scatter(NAC_phase_angle[NAC_Phase3], NAC_intensity[NAC_Phase3], c=NAC_time[NAC_Phase3], cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
#plt.scatter(NAC_phase_angle[NAC_Phase6][RoI][4:21], NAC_intensity[NAC_Phase6][RoI][4:21], marker='x', s=20, edgecolors='black', color='black')
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.ylim(0.96, 1.07)
plt.xlim(24.1, 51)
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
#plt.savefig("Titan_images/plots/NAC/Normalized_NAC_Phase8.png", dpi = 120, bbox_inches='tight')
plt.show()

