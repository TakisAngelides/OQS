import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler

sub_file_name = 'OQS_ATD_DMRG.sub'
path_to_repo = '/lustre/fs24/group/cqta/tangelides/OQS'
path_to_project = f'{path_to_repo}/OQS_Schwinger_Staggered'
file_to_run = 'OQS_ATD_DMRG.jl'
name_of_dag = 'OQS_ATD_DMRG'

N_list = [12]
tau_list = [0.01, 0.001] # time step size
# tau_previous_list = [0] + tau_list[:-1] # this can be left untouched for continuing from the previous time step size
tau_previous_list = [0 for _ in range(len(tau_list))] # this can be left untouched for having every time step size independent
max_steps_list = [600, 6000] # how many ATDDMRG time steps to do, this needs to be the same size as the tau list
measure_every_list = [1, 10] # how often to measure observables and store the density matrix, this needs to be the same size as the max_steps_list
cutoff_list = [1e-6, 1e-9] # this is for compression after gate application
tol_list = [1E-16] # tol for dmrg stopping condition
x_list = [1] # [np.round(1/0.8**2, 6)]
ma_list = [-0.1]
l_0_list = [-1] 
D_list = [1000] # anyway the dmrg will stop from tol
lambd_list = [0.0]
aD_0_list = [0, 0.5]
aT_list = [5.0, 10.0]
sigma_over_a_list = [3.0]
env_corr_type_list = ["delta"]
max_sweeps_list = [1000] # this is for the dmrg but it will stop from tol
l_0_initial_state = 0.0
dirac_vacuum_initial_state = "false"
max_rho_D = 200
project_number = 16 # <==================================== Important to change according to which project file ================================

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
                max_steps = max_steps_list[tau_idx]
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
                                                            sparse_evol = "true" if (max(N_list) < 8) else "false"
                                                            for ma in ma_list:
                                                                
                                                                measure_every = measure_every_list[tau_idx]
                                                                
                                                                if tau_idx == 0:
                                                                    h5_previous_path = 'None'
                                                                else:
                                                                    h5_previous_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau_previous}_x_{x}_l0_{l_0}_ma_{ma}_env_{env_corr_type}_sig_{sigma_over_a}_aT_{aT}_lam_{lambd}_aD0_{aD_0}_l0init_{l_0_initial_state}_cutoff_{cutoff}.h5'
                                                                h5_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau}_x_{x}_l0_{l_0}_ma_{ma}_env_{env_corr_type}_sig_{sigma_over_a}_aT_{aT}_lam_{lambd}_aD0_{aD_0}_l0init_{l_0_initial_state}_cutoff_{cutoff}.h5'
                                                                
                                                                tau_text = str(tau).replace('.', '')
                                                                x_text = str(x).replace('.', '')
                                                                l_0_text = str(l_0).replace('.', '')
                                                                ma_text = str(ma).replace('.', '')
                                                                sigma_over_a_text = str(sigma_over_a).replace('.', '')
                                                                aT_text = str(aT).replace('.', '')
                                                                lambd_text = str(lambd).replace('.', '')
                                                                aD_0_text = str(aD_0).replace('.', '')
                                                                l_0_initial_state_text = str(l_0_initial_state).replace('.', '')
                                                                
                                                                job_name = f'N_{N}_tau_{tau_text}_x_{x_text}_l0_{l_0_text}_ma_{ma_text}_env_{env_corr_type}_sig_{sigma_over_a_text}_aT_{aT_text}_lam_{lambd_text}_aD0_{aD_0_text}_l0init_{l_0_initial_state_text}_cutoff_{cutoff}'
                                                                f.write(f'JOB ' + job_name + f' {path_to_project}/{sub_file_name}\n')
                                                                f.write(f'VARS ' + job_name + f' N="{N}" tau="{tau}" cutoff="{cutoff}" tol="{tol}" x="{x}" l_0="{l_0}" ma="{ma}" max_steps="{max_steps}" project_number="{project_number}" h5_path="{h5_path}" measure_every="{measure_every}" h5_previous_path="{h5_previous_path}" D="{D}" lambda="{lambd}" aD_0="{aD_0}" aT="{aT}" sigma_over_a="{sigma_over_a}" env_corr_type="{env_corr_type}" max_sweeps="{max_sweeps}" sparse_evol="{sparse_evol}" path_to_project="{path_to_project}" file_to_run="{file_to_run}" l_0_initial_state="{l_0_initial_state}" dirac_vacuum_initial_state="{dirac_vacuum_initial_state}" max_rho_D="{max_rho_D}"\n')
                                                                f.write('RETRY ' + job_name + ' 3\n')
                                                                counter += 1
                                                                
                                                                # if tau_previous != 0:
                                                                    
                                                                #     tau_previous_text = str(tau_previous).replace('.', '')
                                                                #     previous_job_name = f'N_{N}_tau_{tau_previous_text}_x_{x_text}_l0_{l_0_text}_ma_{ma_text}_env_{env_corr_type}_sig_{sigma_over_a_text}_aT_{aT_text}_lam_{lambd_text}_aD0_{aD_0_text}_l0init_{l_0_initial_state_text}_cutoff_{cutoff}'
                                                                #     f.write(f'PARENT ' + previous_job_name + ' CHILD ' + job_name + '\n')
                                                                        
    print(f'Number of jobs: {counter}')

def make_plots():
    
    if max(N_list) <= 8:
        columns = ['N', 'x', 'ma', 'l_0_init', 'l_0', 'tau', 'env', 'sig', 'aT', 'lam', 'aD_0', 'at', 'D', 'ee', 'pn', 'cutoff']
    else:
        columns = ['N', 'x', 'ma', 'l_0_init', 'l_0', 'tau', 'env', 'sig', 'aT', 'lam', 'aD_0', 'at', 'D', 'pn', 'cutoff']
    df = pd.DataFrame(columns = columns)
    
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_D_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_D_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/z_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/z_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_step_time_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_step_time_vs_iteration')
    if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/Particle_Number_Density_vs_iteration'):
        os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/Particle_Number_Density_vs_iteration')
                  
    for filepath in os.listdir(f'{path_to_project}/DAGS/{project_number}/HDF5'):
        
        try:
                        
            file = h5py.File(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}', 'r')
            
            row = filepath[:-3].strip().split('_')
            N, tau, x, l_0, ma, env, sig, aT, lam, aD_0, l_0_init, cutoff = row[1], row[3], row[5], row[7], row[9], row[11], row[13], row[15], row[17], row[19], row[21], row[23]
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
            
            if 'max_bond_list' not in file.keys():
                continue
            max_bond_list = file['max_bond_list'][()]
            at_list = file['step_num_list'][()]*tau
            if max(N_list) <= 8:
                ee_list = file['ee_list'][()]
            avg_bond_list = file['avg_bond_list'][()]
            avg_step_time = file['avg_step_time'][()][1:]
            z_list = [file[f'z_list_{idx}'][()] for idx in range(1, N+1)]
            if N <= 8:
                z_list_sparse = [file[f'z_list_sparse_{idx}'][()] for idx in range(1, N+1)]
            
                # print('-------------------------------------------------------------------------------------------------------------------------')
                # print(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}' + ' was found.')
                # print(filepath)
                # print(file.keys())
                # print('-------------------------------------------------------------------------------------------------------------------------')
                        
        except:
            
            # print('*************************************************************************************************************************')
            print(f'{path_to_project}/DAGS/{project_number}/HDF5/{filepath}' + ' gave an error.')
            # print('*************************************************************************************************************************')
            continue
                    
        plt.plot(at_list, max_bond_list)
        plt.ylabel('Maximum bond dimension')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(at_list, avg_bond_list)
        plt.ylabel('Average bond dimension')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_D_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(np.arange(len(avg_step_time)), avg_step_time)
        plt.ylabel('Average step time')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_step_time_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        particle_number_density = []
        particle_number_density_sparse = []
        for t in range(len(at_list)):
            particle_number_density.append(0.5 + (0.5/N)*sum([z_list[idx][t]*(-1)**(idx) for idx in range(N)]))
            if N <= 8:
                particle_number_density_sparse.append(0.5 + (0.5/N)*sum([z_list_sparse[idx][t]*(-1)**(idx) for idx in range(N)]))
        plt.plot(at_list, particle_number_density)
        if N <= 8:
            plt.plot(at_list, particle_number_density_sparse, '--', label = 'Sparse')
        plt.ylabel('Particle Number Density')
        plt.xlabel('Iteration')
        plt.legend()
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Particle_Number_Density_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        if max(N_list) <= 8:
            plt.plot(at_list, ee_list)
            plt.ylabel('Entanglement Entropy')
            plt.xlabel('Iteration')
            plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
            plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration/{filepath}.png', bbox_inches = 'tight')
            plt.close()
        
        for idx in range(N):
            plt.plot(at_list, z_list[idx], label = f'Site {idx}')
        plt.ylabel('<Z_i>')
        plt.xlabel('Iteration')
        plt.legend()
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/z_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
              
        for i in range(len(at_list)):
            
            if max(N_list) <= 8:
                at, D, ee, pn = at_list[i], max_bond_list[i], ee_list[i], particle_number_density[i]
                new_row = pd.DataFrame([[N, x, ma, l_0_init, l_0, tau, env, sig, aT, lam, aD_0, at, D, ee, pn, cutoff]], columns = df.columns)
            else:
                at, D, pn = at_list[i], max_bond_list[i], particle_number_density[i]
                new_row = pd.DataFrame([[N, x, ma, l_0_init, l_0, tau, env, sig, aT, lam, aD_0, at, D, pn, cutoff]], columns = df.columns)
            df = pd.concat([df, new_row], ignore_index = True)
            
    df.to_csv('results.csv', index = False)
                    
def make_plots_from_df():
    
    df = pd.read_csv('results.csv')
    
    # Make plots of energy vs time at size
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
    
    plt.rc('axes', prop_cycle=(cycler('color', list('rbgk')) + cycler('linestyle', ['-', '--', ':', '-.'])))
    
    for N in df.N.unique():
        for x in df.x.unique():
            for l_0 in df.l_0.unique():
                for ma in df.ma.unique():
                    for aD_0 in df.aD_0.unique():
                        for aT in df.aT.unique():
                            for cutoff in df.cutoff.unique():
                    
                                df_tmp = df[(df.N == N) & (df.x == x) & (df.l_0 == l_0) & (df.ma == ma) & (df.aD_0 == aD_0) & (df.cutoff == cutoff) & (df.aT == aT)]
                            
                                at_list, pn_list = np.array(df_tmp['at']), np.array(df_tmp.pn)
                                
                                try:                                                 
                                    plt.plot(at_list, pn_list, label = f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, max = {max(pn_list)}')
                                except:
                                    print(f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT} gave an error.')
                                plt.title('PND vs at')
                        
                        plt.legend(fontsize = 6)
                        plt.savefig(f'Plots/pnd_vs_at_aD_0_{aD_0}.png')
                        plt.clf()
                                                                          
# write_dag()
make_plots()
make_plots_from_df()
