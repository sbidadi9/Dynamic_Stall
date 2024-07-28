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
T = 0.5

aoa = [-50, -40, -30, 30, 40, 50]
af_type = ['ffa_w3_211', 'ffa_w3_241', 'ffa_w3_270', 'ffa_w3_301', 'ffa_w3_330', 'ffa_w3_360', 'ffa_w3_500']
af_thick = ['211', '241', '270', '301', '330', '360', '500']

def get_case_data(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

with PdfPages('cl_vs_time_ffa_w3_211.pdf') as pfpgs:
     fig = plt.figure()

     j = 0      
     for i, c in enumerate(aoa):
         cdata = get_case_data(i,j)
         cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
         plt.plot(cdata["Time"].iloc[1:]/T, cl.iloc[1:], label=aoa[i])

     plt.xlabel('t/T')
     plt.ylabel(r'$C_l$')
     plt.legend(loc='lower right', ncol=8)
     plt.tight_layout()
     pfpgs.savefig()
     plt.close(fig)


