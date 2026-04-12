import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.colors import ListedColormap
import csv

import functs
import CB3_NAC_data_analysis_Part1

NAC_time = CB3_NAC_data_analysis_Part1.NAC_time
NAC_phase_angle = CB3_NAC_data_analysis_Part1.NAC_phase_angle
NAC_intensity = CB3_NAC_data_analysis_Part1.NAC_intensity

NAC_phase_angle_continuous_OBS, NAC_phase_angle = functs.find_Satellites_continuous_OBS(NAC_time, NAC_phase_angle, 1)
NAC_intensity_continuous_OBS, NAC_intensity = functs.find_Satellites_continuous_OBS(NAC_time, NAC_intensity, 1)
NAC_time_continuous_OBS, NAC_time = functs.find_Satellites_continuous_OBS(NAC_time, NAC_time, 1)

NAC_PhaseAngle_factors = 1.675e18*functs.Lambertian_reflector_phase_function(NAC_phase_angle_continuous_OBS)
NAC_intensity_continuous_OBS = NAC_intensity_continuous_OBS / NAC_PhaseAngle_factors

plt.scatter(NAC_phase_angle_continuous_OBS, NAC_intensity_continuous_OBS, c=NAC_time_continuous_OBS, cmap='gist_rainbow', marker='o', s=20, linewidths=0.5, edgecolors='black', label='NAC')
plt.xlim(29, 51)
plt.xticks(np.linspace(30, 50, 9))
#plt.ylim(1.61e18, 1.77e18)
cbar = plt.colorbar()
cbar.ax.set_ylabel('OBStime [MJD]', rotation=90)
plt.axvline(38.3625, linestyle='dashdot', color='grey')
plt.axvline(49.0875, linestyle='dashdot', color='grey')
plt.xlabel('Phase angle [deg]')
plt.ylabel('Intensity')
plt.title('Phase angle vs. Intensity vs. OBStime')
plt.savefig("Titan_images/CB3_NAC_continuous_OBS.png", dpi = 120, bbox_inches='tight')
plt.show()

#--------------------------------------------------------

functs.plot_Satellites_continuous_OBS(NAC_time_continuous_OBS, NAC_phase_angle_continuous_OBS, NAC_intensity_continuous_OBS)
plt.savefig("Titan_images/CB3_NAC_continuous_plot.png", dpi = 120, bbox_inches='tight')
plt.show()
