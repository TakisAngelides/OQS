import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd

sub_file_name = 'thermal_ATD_DMRG.sub'
path_to_repo = '/lustre/fs24/group/cqta/tangelides/OQS'
path_to_project = f'{path_to_repo}/Thermal_MPO_Schwinger_Staggered'
file_to_run = 'thermal_ATD_DMRG.jl'
name_of_dag = 'thermal_ATD_DMRG'

N_list = [4]
tau_list = [10**(-pow) for pow in [3, 4, 5, 6]] 
tau_previous_list = [0.0] + tau_list[:-1]   
cutoff_list = [1E-14] # this is for compression after gate application
tol_list = [1E-12] # this is when to stop the iTEBD regarding convergence in energy
x_list = [0.1, 1.1]
l_0_list = [0.4, 0.8]
mg_list = [0.1, 1.2]
get_dmrg = 'true'
max_steps_list = [20000] 
measure_every_list = [1] # list(np.array(max_steps_list)//1000)
project_number = 2 # <==================================== Important to change according to which project file ================================

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
        
    with open(f'{name_of_dag}.dag', 'w') as f:

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
    
    columns = ['N', 'x', 'mg', 'l_0', 'tau', 'step', 'energy', 'D', 'dmrg_energy']
    df = pd.DataFrame(columns = columns)
    
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_tau'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_tau')    
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_diff_vs_tau'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_diff_vs_tau')  
          
    for filepath in os.listdir(f'{path_to_project}/DAGS/{project_number}/HDF5'):
        
        file = h5py.File(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}', 'r')
        
        row = filepath.strip().split('_')
        N, tau, x, l_0, mg = row[1], row[3], row[5], row[7], row[9]
        N = int(N)
        tau = float(tau)
        x = float(x)
        l_0 = float(l_0)
        mg = mg.strip().split('.')[0]
                
        energy_list = file['energy_list'][()]
        max_bond_list = file['max_bond_list'][()]
        step_num_list = file['step_num_list'][()]
        
        if 'dmrg_energy' in file.keys():
            dmrg_energy = file['dmrg_energy'][()]
            plt.hlines(dmrg_energy, 0, max(step_num_list), label = f'dmrg: {dmrg_energy}')
        else: 
            dmrg_energy = 'None'
            
        plt.plot(step_num_list, energy_list, label = f'min: {min(energy_list)}')
        plt.title(f'N = {N}, tau = {tau}, \n x = {x}, l_0 = {l_0}, mg = {mg}')
        plt.legend()
        plt.ylabel('Energy')
        plt.xlabel('Iteration')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_iteration/N_{N}_tau_{tau}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(step_num_list, max_bond_list)
        plt.ylabel('Maximum bond dimension')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, mg = {mg}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration/N_{N}_tau_{tau}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
        plt.close()
              
        for i in range(len(step_num_list)):
            
            step, energy, D = step_num_list[i], energy_list[i], max_bond_list[i]
            new_row = pd.DataFrame([[N, x, mg, l_0, tau, step, energy, D, dmrg_energy]], columns = df.columns)
            df = pd.concat([df, new_row], ignore_index=True)
    
    df.to_csv('results.csv', index = False)
    
    df = pd.read_csv('results.csv')
    
    # Make plots of energy vs time step size
    for N in df.N.unique():
        for x in df.x.unique():
            for l_0 in df.l_0.unique():
                for mg in df.mg.unique():
                    df_tmp = df[(df.N == N) & (df.x == x) & (df.l_0 == l_0) & (df.mg == mg)]
                    taus = df_tmp.tau.unique()
                    min_energy = []
                    for tau in taus:
                        df_tmp1 = df_tmp[df_tmp.tau == tau]
                        min_energy.append(min(df_tmp1.energy))
                    
                    dmrg_energy = list(df_tmp.dmrg_energy.unique())[0]
                    
                    taus, min_energy = list(zip(*sorted(zip(taus, min_energy), key = lambda x : x[0])))
                    
                    plt.plot(taus, min_energy, '-x', label = f'min: {min(min_energy)}')
                    plt.title(f'N = {N}, x = {x}, \n l_0 = {l_0}, mg = {mg}')
                    if dmrg_energy != 'None':
                        plt.scatter(0, dmrg_energy, label = f'dmrg: {dmrg_energy}')
                    plt.legend()
                    plt.ylabel('Energy')
                    plt.xlabel('tau')
                    plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_tau/N_{N}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
                    plt.close()
                    
                    if dmrg_energy != 'None':
                        plt.plot(taus, abs(np.array(min_energy)-dmrg_energy), '-x', label = f'min: {min(abs(np.array(min_energy)-dmrg_energy))}')
                        plt.title(f'N = {N}, x = {x}, \n l_0 = {l_0}, mg = {mg}')
                        plt.legend()
                        plt.ylabel('Energy diff')
                        plt.xlabel('tau')
                        plt.yscale('log')
                        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_diff_vs_tau/N_{N}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
                        plt.close()
                                        
# write_dag()
make_plots()
