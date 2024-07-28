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

aoa = [-50, -40, -30, 30, 40, 50]
af_type = ['ffa_w3_211', 'ffa_w3_241', 'ffa_w3_270', 'ffa_w3_301', 'ffa_w3_330', 'ffa_w3_360', 'ffa_w3_500']
af_thick = ['211', '241', '270', '301', '330', '360', '500']

def get_case_data(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

def get_case_data_restart(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '_restart.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

with PdfPages('cl_vs_time_ffa_w3_211.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.5
     j = 0      
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
#         cdata_res = get_case_data_restart(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
#         cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
#         cl_tot = np.concatenate([cl.iloc[1:], cl_res[78:]])
#         time_tot = np.concatenate([cdata["Time"].iloc[1:],cdata_res["Time"].iloc[78:]])
         cl_tot = cl.iloc[1:]
         time_tot = cdata["Time"].iloc[1:]
         plt.plot(time_tot/T, cl_tot, label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (211)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

with PdfPages('cl_vs_time_ffa_w3_241.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.40
     j = 1      
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cdata_res = get_case_data_restart(i,j)

         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
         cl_tot = np.concatenate([cl.iloc[1:], cl_res[234:]])
         time_tot = np.concatenate([cdata["Time"].iloc[1:],cdata_res["Time"].iloc[234:]])
#         cl_tot = cl.iloc[1:]
#         time_tot = cdata["Time"].iloc[1:]
 
         plt.plot(time_tot/T, cl_tot, label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (241)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

with PdfPages('cl_vs_time_ffa_w3_270.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.36
     j = 2 
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         plt.plot(cdata["Time"].iloc[1:]/T, cl.iloc[1:], label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (270)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

with PdfPages('cl_vs_time_ffa_w3_301.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.32
     j = 3
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         plt.plot(cdata["Time"].iloc[1:]/T, cl.iloc[1:], label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (301)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

with PdfPages('cl_vs_time_ffa_w3_330.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.28
     j = 4
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         plt.plot(cdata["Time"].iloc[1:]/T, cl.iloc[1:], label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (330)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

with PdfPages('cl_vs_time_ffa_w3_360.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.25
     j = 5
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         plt.plot(cdata["Time"].iloc[1:]/T, cl.iloc[1:], label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (360)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

with PdfPages('cl_vs_time_ffa_w3_500.pdf') as pfpgs:
     fig = plt.figure()

     T = 0.25
     j = 6      
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         plt.plot(cdata["Time"].iloc[1:]/T, cl.iloc[1:], label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$ (500)')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)

