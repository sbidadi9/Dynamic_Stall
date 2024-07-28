import yaml, json, glob, sys
from pathlib import Path
import numpy as np
import pandas as pd
import os
import shutil

# Specify the location of the directory where the grids reside:
grids_dir = "/scratch/sbidadi/Dynamic_Stall/grids/"

# Specify the run_folder:
runf = "/scratch/sbidadi/Dynamic_Stall/final/nalu_runs"

def gen_static_case(af_name, mesh_file, aoa, term_step_count, frequency, start, time, run_folder='nalu_runs', template="ffa_w3_template.yaml"):

    if ( not Path(template).exists() ):
        print("Template file ", template, " doesn't exist. Please check your inputs")
        sys.exit()

    tfile = yaml.load(open(template),Loader=yaml.UnsafeLoader)

    if ( not Path(mesh_file).exists() ):
        print("Mesh file ", mesh_file, " doesn't exist. Please check your inputs")
        sys.exit()

    Path(run_folder+'/{}/aoa_{}'.format(af_name,aoa)).mkdir(parents=True, exist_ok=True)

    tfile['linear_solvers'][0]['hypre_cfg_file'] = '/scratch/sbidadi/Dynamic_Stall/hypre_file.yaml'
    tfile['linear_solvers'][1]['hypre_cfg_file'] = '/scratch/sbidadi/Dynamic_Stall/hypre_file.yaml'
    tfile['linear_solvers'][2]['hypre_cfg_file'] = '/scratch/sbidadi/Dynamic_Stall/hypre_file.yaml'

    tfile['realms'][0]['mesh'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/restart_background/' + af_name + '_' + str(aoa) + '.rst'
    tfile['realms'][1]['mesh'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/restart_airfoil/' + af_name + '_' + str(aoa) + '.rst'

    tfile['Time_Integrators'][0]['StandardTimeIntegrator']['termination_step_count'] = term_step_count
    tfile['Time_Integrators'][0]['StandardTimeIntegrator']['time_step'] = 0.000266667


    # output:

    tfile['realms'][0]['output']['output_data_base_name'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/results_background/{}_{}.e'.format(af_name, aoa)
    tfile['realms'][1]['output']['output_data_base_name'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/results_airfoil/{}_{}.e'.format(af_name, aoa)

    tfile['realms'][0]['output']['output_frequency'] = frequency
    tfile['realms'][1]['output']['output_frequency'] = frequency

    tfile['realms'][0]['output']['output_start'] = start
    tfile['realms'][1]['output']['output_start'] = start

    tfile['realms'][1]['post_processing'][0]['output_file_name'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/pp_'+af_name+'_'+str(aoa)+'_restart_3.dat'
    tfile['realms'][1]['post_processing'][1]['output_file_name'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/' + af_name+'_'+str(aoa)+'_restart_3.dat'


    # restart:

    tfile['realms'][0]['restart']['restart_data_base_name'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/restart_background/{}_{}.rst'.format(af_name, aoa)
    tfile['realms'][0]['restart']['restart_frequency'] = frequency
    tfile['realms'][0]['restart']['restart_start'] = start
    tfile['realms'][0]['restart']['restart_time'] = time 

    tfile['realms'][1]['restart']['restart_data_base_name'] = run_folder + '/' + af_name + '/aoa_' + str(aoa) + '/restart_airfoil/{}_{}.rst'.format(af_name, aoa)
    tfile['realms'][1]['restart']['restart_frequency'] = frequency
    tfile['realms'][1]['restart']['restart_start'] = start
    tfile['realms'][1]['restart']['restart_time'] = time 

    tfile['realms'][0]['equation_systems']['max_iterations'] = 2

    yaml.dump(tfile, open(run_folder+'/{}/aoa_{}/{}_static_aoa_{}_restart.yaml'.format(af_name, aoa, af_name, aoa),'w'), default_flow_style=False)


def append_to_file(af_name, mesh_file, aoa, term_step_count, run_folder='nalu_runs', template="ffa_w3_restart_template.yaml"):
    atofile = open(run_folder+'/{}/aoa_{}/{}_static_aoa_{}_restart.yaml'.format(af_name, aoa, af_name, aoa),'a')
    bottom = open('bottom_mixed_{}_aoa_{}.yaml'.format(af_name,aoa),'r')
    atofile.write(bottom.read())

def gen_ffa_w3_dynamic_cases(af_name, aoa_range = [32, 50], term_step_count=16000, freq=100, start=100, time=1.0, run_folder='nalu_runs', template="ffa_w3_restart_template.yaml"):

    for i, aoa in enumerate(aoa_range):
        gen_static_case(af_name, grids_dir+af_name+'/'+af_name+'_near_body_aoa_'+str(aoa)+'.exo', aoa, term_step_count, freq, start, time, run_folder, template)
        append_to_file(af_name, grids_dir+af_name+'/'+af_name+'_near_body_aoa_'+str(aoa)+'.exo', aoa, term_step_count, run_folder, template)

if __name__=="__main__":

# ********** Airfoils ************
#    aoa_array = [-50, -45, -40, -35, -30, 30, 35, 40, 45, 50]
    aoa_array = [-50, -40, -30, 30, 40, 50]

#    mesh_files = ["ffa_w3_211", "ffa_w3_241", "ffa_w3_270", "ffa_w3_301", "ffa_w3_330", "ffa_w3_360", "ffa_w3_500"]
#    termination_step_count = [3760, 3080, 2710, 2450, 2100, 1910, 1930]
#    output_restart_freq = [235, 192, 169, 153, 131, 119, 120]

    mesh_files = ["ffa_w3_211"]
    termination_step_count = [9400]
    output_restart_freq = [235]
    output_restart_start = [6815]
    restart_time = [1.81734]

    for i, mfile in enumerate(mesh_files):
        gen_ffa_w3_dynamic_cases(mfile, aoa_range=aoa_array, term_step_count=termination_step_count[i], freq=output_restart_freq[i], start=output_restart_start[i], time=restart_time[i], run_folder=runf)
