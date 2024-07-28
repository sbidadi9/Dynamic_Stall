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
  
N = 8000
T = 0.00013335 # sample spacing
fs = 1/T
sf = 5

nalu_runs = '/lustre/orion/cfd116/scratch/sbidadi/ffa_w3/nalu_runs'

def get_avg_data(af_name, aoa):
    """Get average lift and drag data for a given grid"""
    case_data = [pd.read_csv(af_name+'/'+'data_files/'+af_name+'_'+str(aoa[i])+'.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float).iloc[-N:] for i, c in enumerate(aoa)]
    case_cl = [ np.average((c["Fpy"] + c["Fvy"])/dyn_pres/4.0) for c in case_data]
    case_cd = [ np.average((c["Fpx"] + c["Fvx"])/dyn_pres/4.0) for c in case_data]
    return case_cl, case_cd

def get_case_data(af_name, aoa_one_airfoil):
    print(af_name+'/'+'data_files/'+af_name+'_'+str(aoa_one_airfoil)+'.dat')
    return pd.read_csv(af_name+'/'+'data_files/'+af_name+'_'+str(aoa_one_airfoil)+'.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

af_name='ffa_w3_211'




with PdfPages('flat_plate_fft_cl.pdf') as pfpgs:
    fig = plt.figure()
    for i,c in enumerate(aoa):
        
        if c == 33 or c == 18 or c == 21 or c==24 or c == 27 or c == 30 or c == -12 or c == -15 or c==-18 or c == -21 or c == -24:
           cdata = get_case_data(af_name, c)
           cl = (cdata.iloc[-N:]["Fpy"] + cdata.iloc[-N:]["Fvy"])/dyn_pres/4.0
           time = (cdata.iloc[-N:]["Time"])
           sampling_freq = fftfreq(N, T)[1:N//2]
           st = np.abs(sampling_freq*chord_lengthn[i]/u_infty) 
           cl_fft = fft(np.array(cl) - np.average(cl))
#        plt.loglog(sampling_freq, 2.0/N * np.abs(cl_fft[1:N//2]), label=c)
#           plt.semilogy(sampling_freq, 2.0/N * np.abs(cl_fft[1:N//2]), label=c)
           plt.plot(st, 2.0/N * np.abs(cl_fft[1:N//2]), label=c)



        plt.xlim([0, 1])
        
#       plt.ylim([1.e-10,5.0])
#       plt.title('IDDES: ' + g + ' Second order upwind')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('PSD of $C_l$')
    plt.legend(loc=0)
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
