import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd

sub_file_name = 'attDMRG.sub'
path_to_repo = '/lustre/fs24/group/cqta/tangelides/OQS'
path_to_project = f'{path_to_repo}/Closed_Schwinger_Staggered'
file_to_run = 'attDMRG.jl'

N_list = [6]
tau_list = [10**(-pow) for pow in [3, 4]] 
tau_previous_list = [0.0] + tau_list[:-1]   
cutoff_list = [1E-9] # this is for compression after gate application
tol_list = [1E-9] # this is when to stop the iTEBD regarding convergence in energy
x_list = [1.0]
l_0_list = [1.0]
mg_list = [1.0]
get_dmrg = 'true'
max_steps_list = [1000] 
measure_every_list = list(np.array(max_steps_list)//10)
project_number = 1

def write_dag():

    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots')
        
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/HDF5'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/HDF5')
        
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Logs/Error'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Logs/Error')
        
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Logs/Output'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Logs/Output')
        
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Logs/Log'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Logs/Log')
        
    with open('attDMRG.dag', 'w') as f:

        f.write(f'CONFIG {path_to_repo}/dagman.config\n') # This will contain DAGMAN_USE_DIRECT_SUBMIT = False to avoid obscure bugs of authentication

        # This defines the Job submit file inline (ie inside the dag submit file)
        counter = 0
        
        for N in N_list:
            
            for tau_idx, tau in enumerate(tau_list):
                
                tau_previous = tau_previous_list[tau_idx]
                
                for cutoff in cutoff_list:
                    
                    for tol in tol_list:
                        
                        for x in x_list:
                            
                            for l_0 in l_0_list:
                                
                                for mg in mg_list:
                                    
                                    for max_steps_idx, max_steps in enumerate(max_steps_list):
                                        
                                        measure_every = measure_every_list[max_steps_idx]
                                        
                                        if tau_idx == 0:
                                            h5_previous_path = 'None'
                                        else:
                                            h5_previous_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau_previous}_x_{x}_l_0_{l_0}_mg_{mg}.h5'
                                        h5_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau}_x_{x}_l_0_{l_0}_mg_{mg}.h5'
                                        
                                        tau_text = str(tau).replace('.', '')
                                        x_text = str(x).replace('.', '')
                                        l_0_text = str(l_0).replace('.', '')
                                        mg_text = str(mg).replace('.', '')
                                        
                                        job_name = f'N_{N}_tau_{tau_text}_x_{x_text}_l_0_{l_0_text}_mg_{mg_text}'
                                        f.write(f'JOB ' + job_name + f' {path_to_project}/{sub_file_name}\n')
                                        f.write(f'VARS ' + job_name + f' N="{N}" tau="{tau}" cutoff="{cutoff}" tol="{tol}" x="{x}" l_0="{l_0}" mg="{mg}" max_steps="{max_steps}" project_number="{project_number}" get_dmrg="{get_dmrg}" h5_path="{h5_path}" measure_every="{measure_every}" h5_previous_path="{h5_previous_path}" path_to_project="{path_to_project}" file_to_run="{file_to_run}"\n')
                                        # f.write('RETRY ' + job_name + ' 3\n')
                                        
                                        if tau_idx != 0:
                                            
                                            tau_previous_text = str(tau_previous).replace('.', '')
                                            previous_job_name = f'N_{N}_tau_{tau_previous_text}_x_{x_text}_l_0_{l_0_text}_mg_{mg_text}'
                                            f.write(f'PARENT ' + previous_job_name + ' CHILD ' + job_name + '\n')

def make_plots():
    
    # columns = ['N', 'dt', 'cutoff', 'tol', 'J', 'g_x', 'g_z', 'avg_M', 'energy', 'dmrg_energy']
    # df = pd.DataFrame(columns = columns)
        
    # for filepath in os.listdir(f'/lustre/fs24/group/cqta/tangelides/OQS/DAGS/{project_number}/HDF5'):
        
    #     file = h5py.File(f'/lustre/fs24/group/cqta/tangelides/OQS/DAGS/{project_number}/HDF5/{filepath}', 'r')

    #     N, dt, cutoff, tol, J, g_x, g_z, _, _, t_taken_for_iTEBD = file['constants_N_dt_cutoff_tol_J_gx_gz_num_ED_duration'][()]
    #     avg_M = file['total_M_vs_iteration'][()][-1]/N
    #     energy = file['energy_vs_iteration'][()][-1]
    #     if 'dmrg_gs_energy' in file.keys():
    #         dmrg_energy = file['dmrg_gs_energy'][()]
    #     else:
    #         dmrg_energy = 0
    #     new_row = pd.DataFrame([[N, dt, cutoff, tol, J, g_x, g_z, avg_M, energy, dmrg_energy]], columns = df.columns)
    #     df = pd.concat([df, new_row], ignore_index=True)
    
    return
    
write_dag()
# make_plots()
