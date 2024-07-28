# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl 
from cycler import cycler

prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
print(colors)


mpl.rcParams['lines.linewidth'] = 2 
mpl.rcParams['axes.titlesize'] = 30
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['legend.fontsize'] = 8
#mpl.rcParams['figure.figsize'] = (6.328, 5.328)
mpl.rcParams['figure.figsize'] = (7.328, 5.328)
mpl.rcParams["figure.autolayout"] = True

rho = 1.2
u_infty = 75.0
dyn_pres = 0.5 * rho * (u_infty ** 2)
N = 100
T = 0.35

#aoa_array = [-180, -175, -170, -165,
#                -160, -152, -144, -136, -128,
#                -120, -110, -100, 
#                -90, -80, -70, -60, -50, -40,
#                -33, -30, -27, -24, -21, -18,
#                -15, -12, -9, -6, -3,
#                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
#                40, 50, 60, 70, 80, 90,
#                100, 110, 120,
#                128, 136, 144, 152, 160,
#                165, 170, 175, 180]

# FFA_W3_211:
#af_name='ffa_w3_211'
aoa_211 = [-180, -175, -170, -165,
                -160, -144, -128,
                -120, -110, -100, 
                -90, -80, -50, -40,
                -33, -30, -27, -24, -21, -18,
                -15, -12, -9, -6, -3,
                0, 3, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 110, 120,
                136, 144, 152, 160,
                165, 170, 175, 180]

#af_name='ffa_w3_241'
aoa_241 = [-180, -175, -170, -165,
                -160, -152, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -70, -60, -50, -40,
                -33, -30, -27, -24, -21, -18,
                -15, -12, -6, -3,
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 110, 120,
                128, 136, 144, 152, 160,
                165, 170, 175, 180]

#af_name='ffa_w3_270'
aoa_270 = [-180, -175, -170, -165,
                -160, -152, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -70, -60, -40,
                -33, -30, -27, -24, -21, -18,
                -15, -12, -9, -6, -3,
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 110, 120,
                128, 136, 144, 152, 160,
                165, 170, 175, 180]

#af_name='ffa_w3_301'
aoa_301 = [-180, -175, -170, -165,
                -160, -152, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -70, -60, -50, -40,
                -30, -27, -24, -21, -18,
                -15, -12, -9, -6,
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80,
                100, 110, 120,
                128, 136, 144, 152, 160,
                165, 170, 175, 180]

#af_name='ffa_w3_330'
aoa_330 = [-175, -170, -165,
                -160, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -70, -60, -50, -40,
                -33, -30, -27, -24, -21, -18,
                -15, -12, -9, -6, -3,
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 110, 120,
                136, 144, 152, 160,
                165, 170, 175, 180]

# FFA_W3_360:
#af_name='ffa_w3_360'
aoa_360 = [-180, -175, -170, -165,
                -160, -152, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -60, -50, -40,
                -33, -27, -24, -21, -18,
                -15, -12, -9, -6, -3,
                0, 3, 6, 12, 15, 18, 21, 24, 27, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 120,
                128, 136, 144, 152, 160,
                165, 170, 175, 180]

# FFA_W3_500:
#af_name='ffa_w3_500'
aoa_500 = [-180, -175, -170, -165,
                -160, -152, -144, -136, -128,
                -110, -100, 
                -80, -70, -60, -50, -40,
                -33, -27, -24, -21, -18,
                -15, -12, -9, -6, -3,
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33,
                50, 60, 70, 80, 90,
                100, 110, 120,
                136, 144, 152, 160,
                165, 170, 180]

# Flat Plate:
aoa_fp = [-180, -175, -170, -165,
                -160, -144, -136, -128,
                -120, -110, -100, 
                -90, -80, -70, -60, -50, -40,
                -33, -30, -27, -24, -21, -18,
                -12, -9, -6, -3,
                0, 3, 6, 9, 12, 15, 18, 24, 30, 33,
                40, 50, 60, 70, 80, 90,
                100, 110, 120,
                128, 136, 144, 152, 160,
                165, 170, 175, 180]

aoa_airfoils = [aoa_211, aoa_241, aoa_270, aoa_301, aoa_330, aoa_360, aoa_500, aoa_fp]
aoa_rad_airfoils = [np.multiply(aoa_211, np.pi/180.0), np.multiply(aoa_241, np.pi/180.0), 
                    np.multiply(aoa_270, np.pi/180.0), np.multiply(aoa_301, np.pi/180.0),
                    np.multiply(aoa_330, np.pi/180.0), np.multiply(aoa_360, np.pi/180.0),
                    np.multiply(aoa_500, np.pi/180.0), np.multiply(aoa_fp, np.pi/180.0)]

af_type = ['ffa_w3_211', 'ffa_w3_241', 'ffa_w3_270', 'ffa_w3_301', 'ffa_w3_330', 'ffa_w3_360', 'ffa_w3_500', 'flat_plate']
af_thick = ['211', '241', '270', '301', '330', '360', '500', 'fp']

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_211_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_211_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_241_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_241_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_270_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_270_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_301_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_301_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_330_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_330_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_360_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_360_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

#ref_data = {
#    'cl' : pd.read_csv('ref_data/cl_ffa_w3_500_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cl']),
#    'cd' : pd.read_csv('ref_data/cd_ffa_w3_500_ref_data.txt',header=None,sep='\t', lineterminator='\n',names=['aoa','cd']),
#}

def get_case_data_u_75():
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/nalu_runs/ffa_w3_301/aoa_90/ffa_w3_301_90_combined.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

def get_case_data_u_15():
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/nalu_runs/ffa_w3_301/aoa_90_u_15/ffa_w3_301_90_combined.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

with PdfPages('cl_aoa_ffa_w3_301.pdf') as pfpgs:
     fig = plt.figure()
     
     cdata_75 = get_case_data_u_75()
     cl_75 = (cdata_75["Fpy"] + cdata_75["Fvy"])/(0.5 * rho * (75 ** 2))/4.0
     cd_75 = (cdata_75["Fpx"] + cdata_75["Fvx"])/(0.5 * rho * (75 ** 2))/4.0
     time_75 = (cdata_75["Time"])/T
     print(np.average(cl_75), np.average(cd_75))

     cdata_15 = get_case_data_u_15()
     cl_15 = (cdata_15["Fpy"] + cdata_15["Fvy"])/(0.5 * rho * (15 ** 2))/4.0
     cd_15 = (cdata_15["Fpx"] + cdata_15["Fvx"])/(0.5 * rho * (15 ** 2))/4.0
     time_15 = (cdata_15["Time"])/T
     print(np.average(cl_15), np.average(cd_15))


     plt.subplot(1, 2, 1)
     plt.plot(time_15.iloc[1:], cl_15.iloc[1:], label='u=15m/s')
     plt.plot(time_75.iloc[1:], cl_75.iloc[1:], label='u=75m/s')
     plt.xlabel('t/T')
     plt.ylabel('$C_l$')
     plt.legend(loc='best', ncol=8)
#     plt.xlim(left=0.1)

     plt.subplot(1, 2, 2)
     plt.plot(time_15.iloc[1:], cd_15.iloc[1:], label='u=15m/s')
     plt.plot(time_75.iloc[1:], cd_75.iloc[1:], label='u=75m/s')
     plt.xlabel('t/T')
     plt.ylabel('$C_d$')
     plt.legend(loc='best', ncol=8)
#     plt.xlim(left=0.1)
     plt.ylim([1.5,5.0])

     plt.show()

     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)
