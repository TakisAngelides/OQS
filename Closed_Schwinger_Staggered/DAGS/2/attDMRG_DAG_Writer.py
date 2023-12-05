import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd

sub_file_name = 'attDMRG.sub'
path_to_repo = '/lustre/fs24/group/cqta/tangelides/OQS'
path_to_project = f'{path_to_repo}/Closed_Schwinger_Staggered'
file_to_run = 'attDMRG.jl'

N_list = [4, 6, 8, 10, 20]
tau_list = [10**(-pow) for pow in [3, 4, 5, 6]] 
tau_previous_list = [0.0] + tau_list[:-1]   
cutoff_list = [1E-12] # this is for compression after gate application
tol_list = [1E-12] # this is when to stop the iTEBD regarding convergence in energy
x_list = [0.1, 1.1]
l_0_list = [0.4, 0.8]
mg_list = [0.1, 1.2]
get_dmrg = 'true'
max_steps_list = [50000] 
measure_every_list = list(np.array(max_steps_list)//1000)
project_number = 2

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
                                            h5_previous_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau_previous}_x_{x}_l0_{l_0}_mg_{mg}.h5'
                                        h5_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau}_x_{x}_l0_{l_0}_mg_{mg}.h5'
                                        
                                        tau_text = str(tau).replace('.', '')
                                        x_text = str(x).replace('.', '')
                                        l_0_text = str(l_0).replace('.', '')
                                        mg_text = str(mg).replace('.', '')
                                        
                                        job_name = f'N_{N}_tau_{tau_text}_x_{x_text}_l0_{l_0_text}_mg_{mg_text}'
                                        f.write(f'JOB ' + job_name + f' {path_to_project}/{sub_file_name}\n')
                                        f.write(f'VARS ' + job_name + f' N="{N}" tau="{tau}" cutoff="{cutoff}" tol="{tol}" x="{x}" l_0="{l_0}" mg="{mg}" max_steps="{max_steps}" project_number="{project_number}" get_dmrg="{get_dmrg}" h5_path="{h5_path}" measure_every="{measure_every}" h5_previous_path="{h5_previous_path}" path_to_project="{path_to_project}" file_to_run="{file_to_run}"\n')
                                        # f.write('RETRY ' + job_name + ' 3\n')
                                        
                                        if tau_idx != 0:
                                            
                                            tau_previous_text = str(tau_previous).replace('.', '')
                                            previous_job_name = f'N_{N}_tau_{tau_previous_text}_x_{x_text}_l0_{l_0_text}_mg_{mg_text}'
                                            f.write(f'PARENT ' + previous_job_name + ' CHILD ' + job_name + '\n')

def make_plots():
    
    columns = ['N', 'x', 'mg', 'l_0', 'tau', 'step', 'energy', 'D']
    df = pd.DataFrame(columns = columns)
    
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration')
          
    for filepath in os.listdir(f'{path_to_project}/DAGS/{project_number}/HDF5'):
        
        file = h5py.File(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}', 'r')
        
        N, tau, x, l_0, mg = filepath.strip().split('_')[1::2]
        mg = mg.strip().split('.')[0]
        
        imag = False
        if 'dmrg_energy' in file.keys():
            dmrg_energy = file['dmrg_energy'][()]
            imag = True
            
        energy_list = file['energy_list'][()]
        max_bond_list = file['max_bond_list'][()]
        step_num_list = file['step_num_list'][()]
        
        plt.plot(step_num_list, energy_list, label = f'min: {min(energy_list)}')
        plt.title(f'N = {N}, tau = {tau}, \n x = {x}, l_0 = {l_0}, mg = {mg}')
        if imag:
            plt.hlines(dmrg_energy, 0, max(step_num_list), label = f'dmrg: {dmrg_energy}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_iteration/N_{N}_tau_{tau}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(step_num_list, max_bond_list)
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, mg = {mg}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration/N_{N}_tau_{tau}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
        plt.close()
              
        for i in range(len(step_num_list)):
            
            step, energy, D = step_num_list[i], energy_list[i], max_bond_list[i]
            new_row = pd.DataFrame([[N, x, mg, l_0, tau, step, energy, D]], columns = df.columns)
            df = pd.concat([df, new_row], ignore_index=True)
    
    df.to_csv('results.csv', index = False)

    
# write_dag()
make_plots()
