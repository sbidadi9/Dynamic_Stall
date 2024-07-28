import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Akima1DInterpolator

# Inputs
frequency_edge = 2.84 #1.84
frequency_flap = 2.84

edge_amp_m = 2.63
flap_amp_m = 2.63

rot_avg_deg = 90
rot_amp_deg = 10.0
frequency_rot = 0.75

time_step_size = 0.0002667

# --- For steady state simulation --
# Static
num_iter_steady = 10 #50 #4000
steady_state_time = num_iter_steady*time_step_size

time_array = np.linspace(0.0, steady_state_time, num_iter_steady+1)
edge_offset_array = np.linspace(0.0, 0.0, num_iter_steady+1)
flap_offset_array = np.linspace(0.0, 0.0, num_iter_steady+1)
theta_deg_array = np.linspace(0.0, 0.0, num_iter_steady+1)

# --- For edgewise and flapwise sinusoidal motion --
# Dynamic
sin_num_iter = 1000 #13200 #19000
sin_duration = sin_num_iter*time_step_size

# time
time_sin = np.linspace(1.0e-6, sin_duration, sin_num_iter+1)

# pitching motion
sin_rotation_deg = rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_sin)

# edgwise motion
sin_edge_m = edge_amp_m * np.sin(2 * np.pi * frequency_edge * time_sin)

# flapwise motion
sin_flap_m = flap_amp_m * np.sin(2 * np.pi * frequency_flap * time_sin)

# --- Combined array ---
time_array = np.concatenate((time_array[:-1], time_array[-1]+time_sin))

edge_offset_array = np.concatenate((edge_offset_array, sin_edge_m[:-1]))
delta_edge_offset_array = [j - i for i, j in zip(edge_offset_array[: -1], edge_offset_array[1 :])]

flap_offset_array = np.concatenate((flap_offset_array, sin_flap_m[:-1]))
delta_flap_offset_array = [j - i for i, j in zip(flap_offset_array[: -1], flap_offset_array[1 :])]

theta_deg_array = np.concatenate((theta_deg_array, sin_rotation_deg[:-1]))
delta_theta_deg_array = [j - i for i, j in zip(theta_deg_array[: -1], theta_deg_array[1 :])]

# Plot to verify
plt.figure()
#plt.plot(time_array, edge_offset_array, "g", label="Edge offset (m)")
plt.plot(time_array, theta_deg_array, "g")
plt.xlabel("time (s)")
plt.ylabel("Mixed Displacement (m)")
plt.legend()
plt.show()

# Now Create the output definition
print("  mesh_motion:\n    - name: arbitrary_motion_airfoil\n      mesh_parts:\n      - fluid-HEX\n      frame: non_inertial\n      motion:")

# Rotations
for i in range(len(time_array) - 1):
    print(f"      - type: rotation\n        angle: {delta_theta_deg_array[i]}\n        start_time: {time_array[i]+1e-6}\n        end_time: {time_array[i+1]}\n        axis: [0.0, 0.0, 1.0]\n        origin: [{delta_flap_offset_array[i]}, {delta_edge_offset_array[i]}, 0.0]")

## Displacements
for i in range(len(time_array) - 1):
    print(f"      - type: translation\n        start_time: {time_array[i]+1e-6}\n        end_time: {time_array[i+1]}\n        displacement: [{delta_flap_offset_array[i]}, {delta_edge_offset_array[i]}, 0.0]")
