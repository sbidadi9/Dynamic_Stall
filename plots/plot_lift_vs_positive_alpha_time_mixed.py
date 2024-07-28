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
mpl.rcParams['figure.figsize'] = (7.328, 5.328)
mpl.rcParams["figure.autolayout"] = True

rho = 1.2
u_infty = 75.0
dyn_pres = 0.5 * rho * (u_infty ** 2)
N = 1880
T = 0.35

Ns=8000

rot_amp_deg = 3.69 
frequency_rot = 1.99

# Static
time_step_size = 0.0002667
num_iter_steady = 10
steady_state_time = num_iter_steady*time_step_size
time_array = np.linspace(0.0, steady_state_time, num_iter_steady+1)

aoa = [-50, -40, -30, 30, 40, 50]
aoa_p = [30, 40, 50]


af_type = ['ffa_w3_211', 'ffa_w3_241', 'ffa_w3_270', 'ffa_w3_301', 'ffa_w3_330', 'ffa_w3_360', 'ffa_w3_500']
af_thick = ['211', '241', '270', '301', '330', '360', '500']
N_iters = [1880, 1530, 1350, 1220, 1050, 948, 962]

def get_case_data(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

def get_case_data_restart(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '_restart.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

def get_avg_data(af_name, aoa):
    """Get average lift and drag data for a given grid"""
    case_data = [pd.read_csv(af_name+'/'+'data_files/'+af_name+'_'+str(aoa[i])+'.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float).iloc[-Ns:] for i, c in enumerate(aoa)]
    case_cl = [ np.average((c["Fpy"] + c["Fvy"])/dyn_pres/4.0) for c in case_data]
    case_cd = [ np.average((c["Fpx"] + c["Fvx"])/dyn_pres/4.0) for c in case_data]
    if af_name == "ffa_w3_500":
       case_cm = [ np.average((c["Mtz"]/dyn_pres/4.0)*0.316) for c in case_data]
    else:
       case_cm = [ np.average((c["Mtz"]/dyn_pres/4.0)/4.0) for c in case_data]

    return case_cl, case_cd, case_cm


with PdfPages('cl_ffa_w3_211_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 0
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = '^', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color='black')
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (211)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_ffa_w3_241_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 1
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = '^', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color='black')
 
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (241)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_ffa_w3_270_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 2
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = '^', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color='black')
 
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (270)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_ffa_w3_301_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 3
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = '^', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color=colors[k])
 
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (301)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_ffa_w3_330_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 4
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color='black')
 
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (330)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_ffa_w3_360_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 5
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = '^', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color='black')
 
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (360)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_ffa_w3_500_positive_aoa_vs_time.pdf') as pfpgs:
    fig = plt.figure()

    j = 6
    k = 0

    case_cl = np.zeros(np.size(aoa_p)) 
    case_aoa = np.zeros(np.size(aoa_p))

    case_cl_static, case_cd_static, case_cm_static = get_avg_data(af_type[j], aoa_p)

    for i, c in enumerate(aoa):

        if c > 28.0:

           case_aoa[k] = c
           cdata = get_case_data(i,j)
           time = cdata["Time"]
           cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0

           time_tot = time[-N_iters[j]:]
           cl_tot = cl.iloc[-N_iters[j]:]
           alpha_tot = c + rot_amp_deg * np.sin(2 * np.pi * frequency_rot * time_tot)

           case_cl[k] = np.average(cl_tot)
           plt.plot(case_aoa[k], case_cl_static[k], marker = '^', markersize=9, color='black', label=af_thick[j])
           plt.plot(case_aoa[k], case_cl[k], marker = 'o', markersize=9, color='black', label=af_thick[j])
#           plt.scatter(alpha_tot[::10], cl_tot[::10], marker = '+', color=colors[k])
           plt.plot(alpha_tot, cl_tot, linestyle='solid', linewidth=3, color='black')
 
           k = k + 1

    plt.plot(case_aoa, case_cl_static, linestyle='dashed', linewidth=3, color='black')

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$ (500)')
    plt.xlim([25.0,55.0])
    plt.ylim([0.7,2.8])
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()


