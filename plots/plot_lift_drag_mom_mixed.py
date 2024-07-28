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
N = 1000
T = 0.35

af_num = 5

#aoa = [-50, -40, -30, 30, 40, 50]
aoa = [-50, -40, -30, 30, 40, 50]

af_type = ['ffa_w3_211', 'ffa_w3_241', 'ffa_w3_270', 'ffa_w3_301', 'ffa_w3_330', 'ffa_w3_360', 'ffa_w3_500']
af_thick = ['211', '241', '270', '301', '330', '360', '500']
N_iters = [1880, 1530, 1350, 1220, 1050, 948, 962]

def get_case_data(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)

def get_case_data_restart(i,j):
    return pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '_restart.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float)


def get_avg_data(j):
    """Get average lift and drag data for a given grid"""
    case_data = [pd.read_csv('/scratch/sbidadi/Dynamic_Stall/final/nalu_runs/' + af_type[j] + '/aoa_' + str(aoa[i]) + '/' + af_type[j] + '_' + str(aoa[i]) + '.dat', sep="\s+", skiprows=1,header=None, names=[ "Time","Fpx","Fpy","Fpz","Fvx","Fvy","Fvz","Mtx","Mty","Mtz","Y+min","Y+max"],dtype=float).iloc[-N_iters[j]:] for i, c in enumerate(aoa)]

    case_cl = [ np.average((c["Fpy"] + c["Fvy"])/dyn_pres/4.0) for c in case_data]
    case_cd = [ np.average((c["Fpx"] + c["Fvx"])/dyn_pres/4.0) for c in case_data]
    case_cm = [ np.average((c["Mtz"]/dyn_pres/4.0)/4.0) for c in case_data] 

    return case_cl, case_cd, case_cm

with PdfPages('cm_all_airfoils_positive_aoa.pdf') as pfpgs:
    fig = plt.figure()

    for j, af_name in enumerate(af_type):

        case_cd = np.zeros(3) 
        case_aoa = np.zeros(3)
        k = 0

        for i, c in enumerate(aoa):

            if c > 28.0:
               
               case_aoa[k] = c
               cdata = get_case_data(i,j)

               cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
               cd = (cdata["Fpx"] + cdata["Fvx"])/dyn_pres/4.0
               cm = cdata["Mtz"]/dyn_pres/4.0

               if j == 1:
                  cdata_res = get_case_data_restart(i,j)
                  cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
                  cd_res = (cdata_res["Fpx"] + cdata_res["Fvx"])/dyn_pres/4.0
                  cm_res = cdata_res["Mtz"]/dyn_pres/4.0

                  cl_tot = np.concatenate([cl.iloc[N:], cl_res[234:]])
                  cd_tot = np.concatenate([cd.iloc[N:], cd_res[234:]])
                  cm_tot = np.concatenate([cm.iloc[N:], cm_res[234:]])

               else:
                  cl_tot = cl.iloc[N:]
                  cd_tot = cd.iloc[N:]
                  cm_tot = cm.iloc[N:]

               case_cd[k] = np.average(cd_tot)
               k = k + 1

#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

        plt.plot(case_aoa, case_cd, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_m$')
    plt.xlim([28.0,52.0])
#    plt.ylim([0.0,2.0])
    plt.legend() 
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cm_all_airfoils_negative_aoa.pdf') as pfpgs:
    fig = plt.figure()

    for j, af_name in enumerate(af_type):

        case_cd = np.zeros(3) 
        case_aoa = np.zeros(3)
        k = 0

        for i, c in enumerate(aoa):

            if c < -28.0:
               
               case_aoa[k] = c
               cdata = get_case_data(i,j)

               cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
               cd = (cdata["Fpx"] + cdata["Fvx"])/dyn_pres/4.0
               cm = cdata["Mtz"]/dyn_pres/4.0

               if j < 2:
                  cdata_res = get_case_data_restart(i,j)
                  cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
                  cd_res = (cdata_res["Fpx"] + cdata_res["Fvx"])/dyn_pres/4.0
                  cm_res = cdata_res["Mtz"]/dyn_pres/4.0

                  cl_tot = np.concatenate([cl.iloc[N:], cl_res[234:]])
                  cd_tot = np.concatenate([cd.iloc[N:], cd_res[234:]])
                  cm_tot = np.concatenate([cm.iloc[N:], cm_res[234:]])

               else:
                  cl_tot = cl.iloc[N:]
                  cd_tot = cd.iloc[N:]
                  cm_tot = cm.iloc[N:]

               case_cd[k] = np.average(cd_tot)
               k = k + 1

#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])
        plt.plot(case_aoa, case_cd, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_m$')
    plt.xlim([-52.0,-28.0])
#    plt.ylim([0.0,2.0])
    plt.legend() 
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cd_all_airfoils_positive_aoa.pdf') as pfpgs:
    fig = plt.figure()

    for j, af_name in enumerate(af_type):

        case_cd = np.zeros(3) 
        case_aoa = np.zeros(3)
        k = 0

        for i, c in enumerate(aoa):

            if c > 28.0:
               
               case_aoa[k] = c
               cdata = get_case_data(i,j)

               cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
               cd = (cdata["Fpx"] + cdata["Fvx"])/dyn_pres/4.0
               cm = cdata["Mtz"]/dyn_pres/4.0

               if j < 2:
                  cdata_res = get_case_data_restart(i,j)
                  cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
                  cd_res = (cdata_res["Fpx"] + cdata_res["Fvx"])/dyn_pres/4.0
                  cm_res = cdata_res["Mtz"]/dyn_pres/4.0

                  cl_tot = np.concatenate([cl.iloc[N:], cl_res[234:]])
                  cd_tot = np.concatenate([cd.iloc[N:], cd_res[234:]])
                  cm_tot = np.concatenate([cm.iloc[N:], cm_res[234:]])

               else:
                  cl_tot = cl.iloc[N:]
                  cd_tot = cd.iloc[N:]
                  cm_tot = cm.iloc[N:]

               case_cd[k] = np.average(cd_tot)
               k = k + 1

#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

        plt.plot(case_aoa, case_cd, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_d$')
    plt.xlim([28.0,52.0])
    plt.ylim([0.0,2.0])
    plt.legend() 
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cd_all_airfoils_negative_aoa.pdf') as pfpgs:
    fig = plt.figure()

    for j, af_name in enumerate(af_type):

        case_cd = np.zeros(3) 
        case_aoa = np.zeros(3)
        k = 0

        for i, c in enumerate(aoa):

            if c < -28.0:
               
               case_aoa[k] = c
               cdata = get_case_data(i,j)

               cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
               cd = (cdata["Fpx"] + cdata["Fvx"])/dyn_pres/4.0
               cm = cdata["Mtz"]/dyn_pres/4.0

               if j < 2:
                  cdata_res = get_case_data_restart(i,j)
                  cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
                  cd_res = (cdata_res["Fpx"] + cdata_res["Fvx"])/dyn_pres/4.0
                  cm_res = cdata_res["Mtz"]/dyn_pres/4.0

                  cl_tot = np.concatenate([cl.iloc[N:], cl_res[234:]])
                  cd_tot = np.concatenate([cd.iloc[N:], cd_res[234:]])
                  cm_tot = np.concatenate([cm.iloc[N:], cm_res[234:]])

               else:
                  cl_tot = cl.iloc[N:]
                  cd_tot = cd.iloc[N:]
                  cm_tot = cm.iloc[N:]

               case_cd[k] = np.average(cd_tot)
               k = k + 1

#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])
        plt.plot(case_aoa, case_cd, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_d$')
    plt.xlim([-52.0,-28.0])
    plt.ylim([0.0,2.0])
    plt.legend() 
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()


with PdfPages('cl_all_airfoils_positive_aoa.pdf') as pfpgs:
    fig = plt.figure()

    for j, af_name in enumerate(af_type):

        case_cl = np.zeros(3) 
        case_aoa = np.zeros(3)
        k = 0

        for i, c in enumerate(aoa):

            if c > 28.0:
               
               case_aoa[k] = c
               cdata = get_case_data(i,j)

               cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
               cd = (cdata["Fpx"] + cdata["Fvx"])/dyn_pres/4.0
               cm = cdata["Mtz"]/dyn_pres/4.0

               if j < 2:
                  cdata_res = get_case_data_restart(i,j)
                  cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
                  cd_res = (cdata_res["Fpx"] + cdata_res["Fvx"])/dyn_pres/4.0
                  cm_res = cdata_res["Mtz"]/dyn_pres/4.0

                  cl_tot = np.concatenate([cl.iloc[N:], cl_res[234:]])
                  cd_tot = np.concatenate([cd.iloc[N:], cd_res[234:]])
                  cm_tot = np.concatenate([cm.iloc[N:], cm_res[234:]])

               else:
                  cl_tot = cl.iloc[N:]
                  cd_tot = cd.iloc[N:]
                  cm_tot = cm.iloc[N:]

               case_cl[k] = np.average(cl_tot)
               k = k + 1

#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])
        plt.plot(case_aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$')
    plt.xlim([28.0,52.0])
    plt.ylim([0.75,2.0])
    plt.legend() 
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

with PdfPages('cl_all_airfoils_negative_aoa.pdf') as pfpgs:
    fig = plt.figure()

    for j, af_name in enumerate(af_type):

        case_cl = np.zeros(3) 
        case_aoa = np.zeros(3)
        k = 0

        for i, c in enumerate(aoa):

            if c < -28.0:
               
               case_aoa[k] = c
               cdata = get_case_data(i,j)

               cl = (cdata["Fpy"] + cdata["Fvy"])/dyn_pres/4.0
               cd = (cdata["Fpx"] + cdata["Fvx"])/dyn_pres/4.0
               cm = cdata["Mtz"]/dyn_pres/4.0

               if j < 2:
                  cdata_res = get_case_data_restart(i,j)
                  cl_res = (cdata_res["Fpy"] + cdata_res["Fvy"])/dyn_pres/4.0
                  cd_res = (cdata_res["Fpx"] + cdata_res["Fvx"])/dyn_pres/4.0
                  cm_res = cdata_res["Mtz"]/dyn_pres/4.0

                  cl_tot = np.concatenate([cl.iloc[N:], cl_res[234:]])
                  cd_tot = np.concatenate([cd.iloc[N:], cd_res[234:]])
                  cm_tot = np.concatenate([cm.iloc[N:], cm_res[234:]])

               else:
                  cl_tot = cl.iloc[N:]
                  cd_tot = cd.iloc[N:]
                  cm_tot = cm.iloc[N:]

               case_cl[k] = np.average(cl_tot)
               k = k + 1

#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

        plt.plot(case_aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])

    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_l$')
    plt.xlim([-52.0,-28.0])
    plt.ylim([-2.0,-0.75])
    plt.legend() 
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)
    plt.show()

#with PdfPages('cl_all_airfoils.pdf') as pfpgs:
#    fig = plt.figure()
#    for j, af_name in enumerate(af_type):
#        case_cl, case_cd, case_cm = get_avg_data(j)   
#        plt.plot(aoa, case_cl, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])
#
##    plt.plot(ref_data['cl']['aoa'], ref_data['cl']['cl'],'+-', label='Ref. Data', color='black')
#
#    plt.xlabel(r'$\alpha$')
#    plt.ylabel('$C_l$')
#    plt.xlim([-52.0,52.0])
#    plt.ylim([-1.6,2.5]) 
#    plt.minorticks_on()
#    plt.grid()
#    plt.tight_layout()
#    pfpgs.savefig()
#    plt.close(fig)
#    plt.show()

with PdfPages('cd_all_airfoils.pdf') as pfpgs:
    fig = plt.figure()
    for j, af_name in enumerate(af_type):
        case_cl, case_cd, case_cm = get_avg_data(j)
        plt.plot(aoa, case_cd, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])
#    plt.plot(ref_data['cd']['aoa'], ref_data['cd']['cd'],'+-', label='Ref. Data', color='black')
    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_d$')
    plt.xlim([-52.0,52.0]) 
    plt.ylim([-0.1,2.0])  
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)

with PdfPages('cm_all_airfoils.pdf') as pfpgs:
    fig = plt.figure()
    for j, af_name in enumerate(af_type):
        case_cl, case_cd, case_cm = get_avg_data(j)
        plt.plot(aoa, case_cm, marker = 'o', markersize=8, color=colors[j], linestyle = 'solid', label=af_thick[j])
    plt.xlabel(r'$\alpha$')
    plt.ylabel('$C_m$')
    plt.xlim([-52.0,52.0]) 
#    plt.ylim([-0.1,2.0])  
    plt.minorticks_on()
    plt.grid()
    plt.tight_layout()
    pfpgs.savefig()
    plt.close(fig)


