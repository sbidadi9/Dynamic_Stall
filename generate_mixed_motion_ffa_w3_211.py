import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Akima1DInterpolator
import sys

# Inputs
frequency_edge = 1.99 #1.84
frequency_flap = 1.99 #1.59

edge_amp_m = -1.009 #1.6
flap_amp_m = -0.14 #5.35

#rot_avg_deg = 90
rot_amp_deg = 3.69 #10.0
frequency_rot = 1.99 #2.38

time_step_size = 0.0002667

T=0.5

# --- For steady state simulation --
# Static
num_iter_steady = 10
steady_state_time = num_iter_steady*time_step_size

time_array = np.linspace(0.0, steady_state_time, num_iter_steady+1)

# Dynamic
sin_num_iter = 9390
sin_duration = sin_num_iter*time_step_size

# time
time_sin = np.linspace(1.0e-6, sin_duration, sin_num_iter+1)
time_array = np.concatenate((time_array[:-1], time_array[-1]+time_sin))

def generate_mixed_motion(rot_avg_deg, theta_deg_array, x_offset_array, y_offset_array):

    # edgwise motion
    sin_edge_m = edge_amp_m * np.sin(2 * np.pi * frequency_edge * time_sin)
    x_sin_edge_m = sin_edge_m * np.cos(np.radians(rot_avg_deg))
    y_sin_edge_m = sin_edge_m * np.sin(np.radians(rot_avg_deg)) 

    # flapwise motion
    sin_flap_m = flap_amp_m * np.sin(2 * np.pi * frequency_flap * time_sin)
    x_sin_flap_m = sin_flap_m * np.cos(np.radians(rot_avg_deg))
    y_sin_flap_m = sin_flap_m * np.sin(np.radians(rot_avg_deg)) 

    x_offset_array = np.concatenate((x_offset_array, x_sin_edge_m[:-1] + x_sin_flap_m[:-1]))
    delta_x_offset_array = [j - i for i, j in zip(x_offset_array[: -1], x_offset_array[1 :])]

    y_offset_array = np.concatenate((y_offset_array, y_sin_edge_m[:-1] + y_sin_flap_m[:-1]))
    delta_y_offset_array = [j - i for i, j in zip(y_offset_array[: -1], y_offset_array[1 :])]

    # pitching motion
    sin_rotation_deg = rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_sin)

    theta_deg_array = np.concatenate((theta_deg_array, sin_rotation_deg[:-1]))
    delta_theta_deg_array = [j - i for i, j in zip(theta_deg_array[: -1], theta_deg_array[1 :])]

    # Plot to verify
    plt.figure()
    plt.plot(time_array/T, x_offset_array, "r", label="x offset (m)")
    plt.plot(time_array/T, y_offset_array, "b", label="y offset (m)")
    plt.plot(time_array/T, theta_deg_array, "g", label="theta offset (deg)")
    plt.xlabel("t/T")
    plt.ylabel("Mixed Displacement [211]")
    plt.ylim([-4.1, 4.1])
    plt.grid()
    plt.legend()
    plt.show()

    # Now Create the output definition
    print("  mesh_motion:\n    - name: arbitrary_motion_airfoil\n      mesh_parts:\n      - fluid-HEX\n      frame: non_inertial\n      motion:")

    # Rotations
    for i in range(len(time_array) - 1):
        print(f"      - type: rotation\n        angle: {delta_theta_deg_array[i]}\n        start_time: {time_array[i]+1e-6}\n        end_time: {time_array[i+1]}\n        axis: [0.0, 0.0, 1.0]\n        origin: [{delta_x_offset_array[i]}, {delta_y_offset_array[i]}, 0.0]")

    ## Displacements
    for i in range(len(time_array) - 1):
        print(f"      - type: translation\n        start_time: {time_array[i]+1e-6}\n        end_time: {time_array[i+1]}\n        displacement: [{delta_x_offset_array[i]}, {delta_y_offset_array[i]}, 0.0]")


if __name__=="__main__":
    rot_avg_deg = sys.argv[1]

    theta_deg_array = np.linspace(0.0, 0.0, num_iter_steady+1)
    x_offset_array = np.linspace(0.0, 0.0, num_iter_steady+1)
    y_offset_array = np.linspace(0.0, 0.0, num_iter_steady+1)

    generate_mixed_motion(int(rot_avg_deg), theta_deg_array, x_offset_array, y_offset_array)
