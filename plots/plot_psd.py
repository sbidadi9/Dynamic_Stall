# coding: utf-8
import numpy as np
from scipy.fftpack import fft, fftfreq
import scipy.signal
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

rho = 1.2
u_infty = 75.0
dyn_pres = 0.5 * rho * (u_infty ** 2)
#aoa = [-180, -33, -30, -27, -24, -21, -18, -15]
#aoa = [15, 18, 21, 24, 27, 30, 33, 175, 180]    

aoa = [-180, -175, -170, -165,
                -160, -152, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -70, -60, -50, -40,
                -33, -30, -27, -24, -21, -18,
                -15, -12, -9, -6, -3,
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 110, 120,
                128, 136, 144, 152, 160,
                165, 170, 175, 180]

aoa_rad = np.multiply(aoa, np.pi/180.0)
chord_lengthn = np.sin(aoa_rad)
  
N = 2500
T = 0.000266667 # sample spacing
fs = 1/T
sf = 1

nalu_runs = '/lustre/orion/cfd116/scratch/sbidadi/ffa_w3/nalu_runs'

def get_case_data_u_75():
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/nalu_runs/ffa_w3_301/aoa_90/ffa_w3_301_90_combined.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

def get_case_data_u_15():
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/nalu_runs/ffa_w3_301/aoa_90_u_15/ffa_w3_301_90_combined.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

with PdfPages('psd_ffa_w3_301.pdf') as pfpgs:
    fig = plt.figure()
    cdata_75 = get_case_data_u_75()
    cdata_15 = get_case_data_u_15()
    cl_15 = (cdata_15.iloc[-N:]["Fpy"] + cdata_15.iloc[-N:]["Fvy"])/(0.5 * rho * (15 ** 2))/4.0
    cl_75 = (cdata_75.iloc[-N:]["Fpy"] + cdata_75.iloc[-N:]["Fvy"])/(0.5 * rho * (75 ** 2))/4.0
    time_15 = (cdata_15.iloc[-N:]["Time"])
    time_75 = (cdata_75.iloc[-N:]["Time"])

    sampling_freq = fftfreq(N, T)[1:N//2] 


#    cl_fft_15 = fft(np.array(cl_15) - np.average(cl_15))
#    cl_fft_75 = fft(np.array(cl_75) - np.average(cl_75))

    cl_fft_15 = fft(np.array(cl_15))
    cl_fft_75 = fft(np.array(cl_75))



    plt.plot(sampling_freq, 2.0/N * np.abs(cl_fft_15[1:N//2]), label='u=15m/s')
    plt.plot(sampling_freq, 2.0/N * np.abs(cl_fft_75[1:N//2]), label='u=75m/s')

    plt.xlim([0.0, 40])       

#       plt.ylim([1.e-10,5.0])
#       plt.title('IDDES: ' + g + ' Second order upwind')

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('PSD of $C_l$')
    plt.legend(loc=0)
    plt.show()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
