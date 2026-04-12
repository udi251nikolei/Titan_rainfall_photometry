import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

primary_methane_rainbow_table = Table.read('Gustav_rainbow_data/rainbow_methan_Titan_primary_deviation_angle.csv', format='ascii')
primary_angles = np.array(primary_methane_rainbow_table['angle (deg)'], dtype=float)
primary_752nm_intensity = np.array(primary_methane_rainbow_table['intensity, 752 nm, 0.0526 deg'], dtype=float)
primary_939nm_intensity = np.array(primary_methane_rainbow_table['intensity, 939 nm, 0.0526 deg'], dtype=float)

secondary_methane_rainbow_table = Table.read('Gustav_rainbow_data/rainbow_methan_Titan_secondary_deviation_angle.csv', format='ascii')
secondary_angles = (np.array(secondary_methane_rainbow_table['angle (deg)'], dtype=float) * -1)[::-1]
secondary_752nm_intensity = (np.array(secondary_methane_rainbow_table['intensity, 752 nm, 0.0526 deg'], dtype=float))[::-1]
secondary_939nm_intensity = (np.array(secondary_methane_rainbow_table['intensity, 939 nm, 0.0526 deg'], dtype=float))[::-1]
'''
plt.plot(primary_angles, primary_752nm_intensity, linestyle='solid', color='tab:orange', label='752 nm')
plt.plot(secondary_angles, secondary_752nm_intensity, linestyle='solid', color='tab:orange')
plt.plot(primary_angles, primary_939nm_intensity, linestyle='solid', color='tab:blue', label='939 nm')
plt.plot(secondary_angles, secondary_939nm_intensity, linestyle='solid', color='tab:blue')
plt.xlim(35, 55)
plt.ylim(-0.003, 0.07)
plt.ylabel('Intensity')
plt.xlabel('Scattering angle (deg)')
plt.grid(True)
plt.legend(loc='upper left')
#plt.savefig("Titan_images/plots/rainbow_model.png", dpi = 250, bbox_inches='tight')
plt.show()
'''
#-----------------------------------------------------------

first_secondary_angle = min(secondary_angles)
tail_primary_angle = max(primary_angles)

find_primary_start = (first_secondary_angle <= primary_angles)
find_secondary_tail = (tail_primary_angle >= secondary_angles)

angles = primary_angles[find_primary_start]
primary_752nm_intensity = primary_752nm_intensity[find_primary_start]
primary_939nm_intensity = primary_939nm_intensity[find_primary_start]
secondary_752nm_intensity = secondary_752nm_intensity[find_secondary_tail]
secondary_939nm_intensity = secondary_939nm_intensity[find_secondary_tail]

intensity_752nm = primary_752nm_intensity + secondary_752nm_intensity
intensity_939nm = primary_939nm_intensity + secondary_939nm_intensity

#-----------------------------------------------------------

find_primary_939nm_peak = (intensity_939nm == max(intensity_939nm))
primary_939nm_angle_at_peak = angles[find_primary_939nm_peak][0]

find_seconday_939nm = (angles <= 40)
find_seconday_939nm_peak = (intensity_939nm[find_seconday_939nm] == max(intensity_939nm[find_seconday_939nm]))
secondary_peak_angles = angles[find_seconday_939nm]
secondary_939nm_angle_at_peak = secondary_peak_angles[find_seconday_939nm_peak][0]

#print(f'Primary 939nm peak at {primary_939nm_angle_at_peak} deg')
#print(f'Secondary 939nm peak at {secondary_939nm_angle_at_peak} deg\n')

find_primary_752nm_peak = (intensity_752nm == max(intensity_752nm))
primary_752nm_angle_at_peak = angles[find_primary_752nm_peak][0]

find_seconday_752nm = (angles <= 40)
find_seconday_752nm_peak = (intensity_752nm[find_seconday_752nm] == max(intensity_752nm[find_seconday_752nm]))
secondary_peak_angles = angles[find_seconday_752nm]
secondary_752nm_angle_at_peak = secondary_peak_angles[find_seconday_752nm_peak][0]

#print(f'Primary 752nm peak at {primary_752nm_angle_at_peak} deg')
#print(f'Secondary 752nm peak at {secondary_752nm_angle_at_peak} deg')

#-----------------------------------------------------------
'''
plt.plot(angles, intensity_752nm, label='752 nm')
plt.plot(angles, intensity_939nm, label='939 nm')
plt.xlabel('Angle [deg]')
plt.ylabel('Rel. intensity')
plt.legend()
#plt.savefig("Titan_images/plots/1-Rainbow-feature.png", dpi = 120, bbox_inches='tight')
plt.show()
'''