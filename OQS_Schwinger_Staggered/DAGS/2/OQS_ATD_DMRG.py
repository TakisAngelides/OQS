import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd

sub_file_name = 'OQS_ATD_DMRG.sub'
path_to_repo = '/lustre/fs24/group/cqta/tangelides/OQS'
path_to_project = f'{path_to_repo}/OQS_Schwinger_Staggered'
file_to_run = 'OQS_ATD_DMRG.jl'
name_of_dag = 'OQS_ATD_DMRG'

N_list = [8]
tau_list = [10**(-pow) for pow in [4, 5]] # time step size
tau_previous_list = [0] + tau_list[:-1] # this can be left untouched
cutoff_list = [1E-14] # this is for compression after gate application
tol_list = [1E-12] # tol for dmrg stopping condition
x_list = [1/0.8**2]
l_0_list = [0.0]
ma_list = [0.5]
max_steps_list = [10000] # how many ATDDMRG time steps to do 
measure_every_list = [1] # how often to measure observables and store the density matrix
D_list = [1000] # anyway the dmrg will stop from tol
lambd_list = [0.0]
aD_0_list = [1.0]
aT_list = [10.0]
sigma_over_a_list = [3.0]
env_corr_type_list = ["delta"]
max_sweeps_list = [1000] # this is for the dmrg but it will stop from tol
l_0_initial_state = 0.0
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
                                for D in D_list:
                                    for aD_0 in aD_0_list:
                                        for lambd in lambd_list:
                                            for aT in aT_list:
                                                for sigma_over_a in sigma_over_a_list:
                                                    for env_corr_type in env_corr_type_list:
                                                        for max_sweeps in max_sweeps_list:
                                                            sparse_evol = "true" if (max(N_list) < 8) and (tau_idx == 0) else "false"
                                                            for ma in ma_list:
                                                                for max_steps_idx, max_steps in enumerate(max_steps_list):
                                                                    
                                                                    measure_every = measure_every_list[max_steps_idx]
                                                                    
                                                                    if tau_idx == 0:
                                                                        h5_previous_path = 'None'
                                                                    else:
                                                                        h5_previous_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau_previous}_x_{x}_l0_{l_0}_ma_{ma}_env_{env_corr_type}_sig_{sigma_over_a}_aT_{aT}_lam_{lambd}_aD0_{aD_0}_l0init_{l_0_initial_state}.h5'
                                                                    h5_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau}_x_{x}_l0_{l_0}_ma_{ma}_env_{env_corr_type}_sig_{sigma_over_a}_aT_{aT}_lam_{lambd}_aD0_{aD_0}_l0init_{l_0_initial_state}.h5'
                                                                    
                                                                    tau_text = str(tau).replace('.', '')
                                                                    x_text = str(x).replace('.', '')
                                                                    l_0_text = str(l_0).replace('.', '')
                                                                    ma_text = str(ma).replace('.', '')
                                                                    sigma_over_a_text = str(sigma_over_a).replace('.', '')
                                                                    aT_text = str(aT).replace('.', '')
                                                                    lambd_text = str(lambd).replace('.', '')
                                                                    aD_0_text = str(aD_0).replace('.', '')
                                                                    l_0_initial_state_text = str(l_0_initial_state).replace('.', '')
                                                                    
                                                                    job_name = f'N_{N}_tau_{tau_text}_x_{x_text}_l0_{l_0_text}_ma_{ma_text}_env_{env_corr_type}_sig_{sigma_over_a_text}_aT_{aT_text}_lam_{lambd_text}_aD0_{aD_0_text}_l0init_{l_0_initial_state_text}'
                                                                    f.write(f'JOB ' + job_name + f' {path_to_project}/{sub_file_name}\n')
                                                                    f.write(f'VARS ' + job_name + f' N="{N}" tau="{tau}" cutoff="{cutoff}" tol="{tol}" x="{x}" l_0="{l_0}" ma="{ma}" max_steps="{max_steps}" project_number="{project_number}" h5_path="{h5_path}" measure_every="{measure_every}" h5_previous_path="{h5_previous_path}" D="{D}" lambda="{lambd}" aD_0="{aD_0}" aT="{aT}" sigma_over_a="{sigma_over_a}" env_corr_type="{env_corr_type}" max_sweeps="{max_sweeps}" sparse_evol="{sparse_evol}" path_to_project="{path_to_project}" file_to_run="{file_to_run}" l_0_initial_state="{l_0_initial_state}"\n')
                                                                    # f.write('RETRY ' + job_name + ' 3\n')
                                                                    
                                                                    if tau_idx != 0:
                                                                        
                                                                        tau_previous_text = str(tau_previous).replace('.', '')
                                                                        previous_job_name = f'N_{N}_tau_{tau_previous_text}_x_{x_text}_l0_{l_0_text}_ma_{ma_text}_env_{env_corr_type}_sig_{sigma_over_a_text}_aT_{aT_text}_lam_{lambd_text}_aD0_{aD_0_text}_l0init_{l_0_initial_state_text}'
                                                                        f.write(f'PARENT ' + previous_job_name + ' CHILD ' + job_name + '\n')


def make_plots():
    
    columns = ['N', 'x', 'ma', 'l_0_init', 'l_0', 'tau', 'env', 'sig', 'aT', 'lam', 'aD_0', 'step', 'D', 'ee', 'z_mid']
    df = pd.DataFrame(columns = columns)
    df_sparse = pd.DataFrame(columns = columns)
    
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/z_middle_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/z_middle_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration')
                  
    for filepath in os.listdir(f'{path_to_project}/DAGS/{project_number}/HDF5'):
        
        try:
            
            file = h5py.File(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}', 'r')
            
            row = filepath[:-3].strip().split('_')
            N, tau, x, l_0, ma, env, sig, aT, lam, aD_0, l_0_init = row[1], row[3], row[5], row[7], row[9], row[11], row[13], row[15], row[17], row[19], row[21] 
            N = int(N)
            tau = float(tau)
            x = float(x)
            l_0 = float(l_0)
            ma = float(ma)
            sig = float(sig)
            aT = float(aT)
            lam = float(lam)
            aD_0 = float(aD_0)
            l_0_init = float(l_0_init)
                    
            max_bond_list = file['max_bond_list'][()]
            step_num_list = file['step_num_list'][()]
            ee_list = file['ee_list'][()]
            z_middle_list = file['z_middle_list'][()]
            
            # print('-------------------------------------------------------------------------------------------------------------------------')
            # print(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}' + ' was found.')
            # print(filepath)
            print(file.keys())
            # print('-------------------------------------------------------------------------------------------------------------------------')
                        
        except:
            
            # print('*************************************************************************************************************************')
            print(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}' + ' was not found.')
            # print('*************************************************************************************************************************')
            continue
                    
        plt.plot(step_num_list, max_bond_list)
        plt.ylabel('Maximum bond dimension')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(step_num_list, ee_list)
        plt.ylabel('Entanglement Entropy')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(step_num_list, z_middle_list)
        plt.ylabel('Z middle')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/z_middle_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
              
        for i in range(len(step_num_list)):
            
            step, z_mid, D, ee = step_num_list[i], z_middle_list[i], max_bond_list[i], ee_list[i]
            new_row = pd.DataFrame([[N, x, ma, l_0_init, l_0, tau, env, sig, aT, lam, aD_0, step, D, ee, z_mid]], columns = df.columns)
            df = pd.concat([df, new_row], ignore_index = True)
            
        
    
    df.to_csv('results.csv', index = False)
    
    df = pd.read_csv('results.csv')
    
    # Make plots of energy vs time step size
    # for N in df.N.unique():
    #     for x in df.x.unique():
    #         for l_0 in df.l_0.unique():
    #             for mg in df.mg.unique():
    #                 df_tmp = df[(df.N == N) & (df.x == x) & (df.l_0 == l_0) & (df.mg == mg)]
    #                 taus = df_tmp.tau.unique()
    #                 min_energy = []
    #                 for tau in taus:
    #                     df_tmp1 = df_tmp[df_tmp.tau == tau]
    #                     min_energy.append(min(df_tmp1.energy))
                    
    #                 dmrg_energy = list(df_tmp.dmrg_energy.unique())[0]
                    
    #                 taus, min_energy = list(zip(*sorted(zip(taus, min_energy), key = lambda x : x[0])))
                    
    #                 plt.plot(taus, min_energy, '-x', label = f'min: {min(min_energy)}')
    #                 plt.title(f'N = {N}, x = {x}, \n l_0 = {l_0}, mg = {mg}')
    #                 if dmrg_energy != 'None':
    #                     plt.scatter(0, dmrg_energy, label = f'dmrg: {dmrg_energy}')
    #                 plt.legend()
    #                 plt.ylabel('Energy')
    #                 plt.xlabel('tau')
    #                 plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Energy_vs_tau/N_{N}_x_{x}_l_0_{l_0}_mg_{mg}.png', bbox_inches = 'tight')
    #                 plt.close()
                                                            
write_dag()
# make_plots()
