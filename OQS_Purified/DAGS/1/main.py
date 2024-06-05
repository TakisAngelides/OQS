import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler

# Inputs needed for static field, dirac vacuum initial state ("dirac_vacuum") and delta/constant correlator

# N = inputs["N"]
# x = inputs["x"]
# ma = inputs["ma"]
# lambd = inputs["lambda"]
# aT = inputs["aT"]
# aD = inputs["aD"]
# conserve_qns = inputs["cqns"]
# dissipator_sites = inputs["ds"]
# tau = inputs["tau"]
# number_of_time_steps = inputs["nots"]
# time_varying_applied_field_flag = inputs["tvaff"]
# taylor_expansion_cutoff_1 = inputs["tec_1"]
# taylor_expansion_cutoff_2 = inputs["tec_2"]
# cutoff = inputs["cutoff"]
# maxdim = inputs["md"]
# taylor_expansion_order = inputs["teo"]
# l_0_1 = inputs["l_0_1"]
# which_applied_field = inputs["waf"]
# env_corr_type = inputs["env_corr_type"]
# which_initial_state = inputs["wis"]

# Extra inputs needed when the correlator is gaussian

# sigma_over_a = inputs["sigma_over_a"]

# Extra inputs needed when the applied field is non-static

# l_0_2 = inputs["l_0_2"]
# a_omega = inputs["a_omega"]

# Extra inputs needed when the initial state is not the dirac vacuum but the gs of H ("gs_naive")

# l_0_initial_state = inputs["l_0_initial_state"]
# max_sweeps_dmrg = inputs["msdmrg"]
# maxdim_dmrg = inputs["mddmrg"]
# energy_tol_dmrg = inputs["etdmrg"]
# cutoff_dmrg = inputs["cdmrg"]

# Extra inputs needed when the intial state is not the dirac vacuum but the dirac vacuum with a string ("dirac_vacuum_with_string")

# flip_sites = inputs["fs"]

def write_dag():
    
    # Name of dag
    name_of_dag = 'run.dag'
    
    # Open text file to write the dag instructions
    f_dag = open(name_of_dag, 'w')
    
    # This will contain DAGMAN_USE_DIRECT_SUBMIT = False to avoid obscure bugs of authentication
    f_dag.write(f'CONFIG /lustre/fs24/group/cqta/tangelides/OQS/dagman.config\n')
    
    # The julia file to run with the given inputs
    file_to_run = 'run.jl'
    
    # Path to submission file to run from dag
    path_to_sub = '/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/run.sub'
    
    # Project number of dag
    project_number = os.getcwd().strip().split('/')[-1]
    
    # Where to find the inputs
    path_to_inputs_h5 = f'/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/{project_number}'
    f_h5 = h5py.File(path_to_inputs_h5, 'w')
        
    # Create relevant folders if needed
    if not os.path.exists(f'{path_to_inputs_h5}/Plots'):
        os.makedirs(f'{path_to_inputs_h5}/Plots')
    path_to_HDF5_folder = f'{path_to_inputs_h5}/HDF5' # Where to save the outputs distinguished after by the job id  
    if not os.path.exists(path_to_HDF5_folder):
        os.makedirs(path_to_HDF5_folder)        
    if not os.path.exists(f'{path_to_inputs_h5}/Logs/Error'):
        os.makedirs(f'{path_to_inputs_h5}/Logs/Error')        
    if not os.path.exists(f'{path_to_inputs_h5}/Logs/Output'):
        os.makedirs(f'{path_to_inputs_h5}/Logs/Output')        
    if not os.path.exists(f'{path_to_inputs_h5}/Logs/Log'):
        os.makedirs(f'{path_to_inputs_h5}/Logs/Log')
            
    # This will form the job id
    counter_of_jobs = 1
    
    # Static applied field case and delta correlator
    lambd = 0
    conserve_qns = "true"
    number_of_time_steps_list = [3] # needs same length as tau list
    taylor_expansion_cutoff_1 = 1e-15
    taylor_expansion_cutoff_2 = 1e-12
    maxdim = 300
    which_applied_field = "constant" # options are: "constant", "sauter", "gaussian", "oscillatory"
    time_varying_applied_field_flag = "false" if which_applied_field == "constant" else "true"
    env_corr_type = "delta" # options are: "constant", "delta", "gaussian"
    which_initial_state = "dirac_vacuum_with_string" # options are: "dirac_vacuum", "gs_naive", "dirac_vacuum_with_string"
    
    for N in [4]:
        dissipator_sites = [i for i in range(1, N+1)]
        flip_sites = [N//2, N//2 + 1] # this is for the case when the initial state is the dirac vacuum with a string and specifies where the string should be placed
        for x in [1]:
            for ma in [1]:
                for aD in [1]:
                    for aT in [10]:
                        for tau_idx, tau in enumerate([0.001]):
                            number_of_time_steps = number_of_time_steps_list[tau_idx]
                            for cutoff in [1e-9]:
                                for taylor_expansion_order in [1]:
                                    for l_0_1 in [0]: # this is the constant part of the applied field
                                        
                                        # Memory, CPU and maximum number of days to run
                                        mem, cpu, days = 4, 1, 1
                                        
                                        # Job id for the dag job names and path to h5 for results
                                        job_id = counter_of_jobs
                                        counter_of_jobs += 1 # after assigning the job_id this is incremented for the next job
                                        path_to_save_results_h5 = path_to_HDF5_folder + f'/{job_id}.h5'
                                        
                                        # Write inputs to h5
                                        g = f_h5.create_group(f'{job_id}')       
                                        g.attrs["N"] = N
                                        g.attrs["x"] = x
                                        g.attrs["ma"] = ma
                                        g.attrs["lambda"] = lambd
                                        g.attrs["aT"] = aT
                                        g.attrs["aD"] = aD
                                        g.attrs["cqns"] = conserve_qns
                                        g.attrs["ds"] = dissipator_sites
                                        g.attrs["tau"] = tau
                                        g.attrs["nots"] = number_of_time_steps
                                        g.attrs["tvaff"] = time_varying_applied_field_flag
                                        g.attrs["tec_1"] = taylor_expansion_cutoff_1
                                        g.attrs["tec_2"] = taylor_expansion_cutoff_2
                                        g.attrs["cutoff"] = cutoff
                                        g.attrs["md"] = maxdim
                                        g.attrs["teo"] = taylor_expansion_order
                                        g.attrs["l_0_1"] = l_0_1
                                        g.attrs["waf"] = which_applied_field
                                        g.attrs["env_corr_type"] = env_corr_type
                                        g.attrs["wis"] = which_initial_state
                                        g.attrs["fs"] = flip_sites
            
                                        # Write job to dag
                                        job_name = f'{job_id}'
                                        f_dag.write(f'JOB ' + job_name + f' {path_to_sub}\n')
                                        f_dag.write(f'VARS ' + job_name + f' path_to_inputs_h5="{path_to_inputs_h5}" job_id="{job_id}" path_to_save_results_h5="{path_to_save_results_h5}" file_to_run="{file_to_run}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                                        f_dag.write('RETRY ' + job_name + ' 0\n')
    
    # Close the dag file and the h5 input file
    f_dag.close() 
    f_h5.close()
    print(f'Total number of jobs in the dag is {counter_of_jobs}')                                        

write_dag()

# ---------------------------------------------------------------------------------------------------------------------------------------

# Old stuff that is slowly deprecated

def write_dag_pulse_field():
    
    sub_file_name = 'OQS_ATD_DMRG_pulse_field.sub'
    name_of_dag = 'OQS_ATD_DMRG_pulse_field'

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
                                                            sparse_evol = "true" if (max(N_list) <= 8) else "false"
                                                            for ma in ma_list:
                                                                for l_0_small in l_0_small_list:
                                                                    for omega in omega_list:
                                                                        for type in type_list:
                                                                            
                                                                            mem = 4
                                                                
                                                                            max_steps = max_steps_list[tau_idx]
                                                                            measure_every = measure_every_list[tau_idx]
                                                                            
                                                                            if tau_previous == 0:
                                                                                h5_previous_path = 'None'
                                                                            else:
                                                                                h5_previous_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau_previous}_x_{x}_l0_{l_0}_ma_{ma}_env_{env_corr_type}_sig_{sigma_over_a}_aT_{aT}_lam_{lambd}_aD0_{aD_0}_l0init_{l_0_initial_state}_cutoff_{cutoff}_l0s_{l_0_small}_{type}_om_{omega}.h5'
                                                                            h5_path = f'{path_to_project}/DAGS/{project_number}/HDF5/N_{N}_tau_{tau}_x_{x}_l0_{l_0}_ma_{ma}_env_{env_corr_type}_sig_{sigma_over_a}_aT_{aT}_lam_{lambd}_aD0_{aD_0}_l0init_{l_0_initial_state}_cutoff_{cutoff}_l0s_{l_0_small}_{type}_om_{omega}.h5'
                                                                            
                                                                            tau_text = str(tau).replace('.', '')
                                                                            x_text = str(x).replace('.', '')
                                                                            l_0_text = str(l_0).replace('.', '')
                                                                            l_0_small_text = str(l_0_small).replace('.', '')
                                                                            omega_text = str(omega).replace('.', '')
                                                                            ma_text = str(ma).replace('.', '')
                                                                            sigma_over_a_text = str(sigma_over_a).replace('.', '')
                                                                            aT_text = str(aT).replace('.', '')
                                                                            lambd_text = str(lambd).replace('.', '')
                                                                            aD_0_text = str(aD_0).replace('.', '')
                                                                            l_0_initial_state_text = str(l_0_initial_state).replace('.', '')
                                                                            
                                                                            job_name = f'N_{N}_tau_{tau_text}_x_{x_text}_l0_{l_0_text}_ma_{ma_text}_env_{env_corr_type}_sig_{sigma_over_a_text}_aT_{aT_text}_lam_{lambd_text}_aD0_{aD_0_text}_l0init_{l_0_initial_state_text}_cutoff_{cutoff}_l0s_{l_0_small_text}_{type}_om_{omega_text}'
                                                                            f.write(f'JOB ' + job_name + f' {path_to_project}/{sub_file_name}\n')
                                                                            f.write(f'VARS ' + job_name + f' N="{N}" tau="{tau}" cutoff="{cutoff}" tol="{tol}" x="{x}" l_0="{l_0}" ma="{ma}" max_steps="{max_steps}" project_number="{project_number}" h5_path="{h5_path}" measure_every="{measure_every}" h5_previous_path="{h5_previous_path}" D="{D}" lambda="{lambd}" aD_0="{aD_0}" aT="{aT}" sigma_over_a="{sigma_over_a}" env_corr_type="{env_corr_type}" max_sweeps="{max_sweeps}" sparse_evol="{sparse_evol}" path_to_project="{path_to_project}" file_to_run="{file_to_run}" l_0_initial_state="{l_0_initial_state}" dirac_vacuum_initial_state="{dirac_vacuum_initial_state}" max_rho_D="{max_rho_D}" l_0_small="{l_0_small}" type="{type}" omega="{omega}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                                                                            f.write('RETRY ' + job_name + ' 3\n')
                                                                            counter += 1
                                                                            
                                                                            if tau_previous != 0:
                                                                                                                                                        
                                                                                tau_previous_text = str(tau_previous).replace('.', '')
                                                                                previous_job_name = f'N_{N}_tau_{tau_previous_text}_x_{x_text}_l0_{l_0_text}_ma_{ma_text}_env_{env_corr_type}_sig_{sigma_over_a_text}_aT_{aT_text}_lam_{lambd_text}_aD0_{aD_0_text}_l0init_{l_0_initial_state_text}_cutoff_{cutoff}_l0s_{l_0_small_text}_{type}_om_{omega_text}'
                                                                                f.write(f'PARENT ' + previous_job_name + ' CHILD ' + job_name + '\n')
                                                                        
    print(f'Number of jobs: {counter}')

def make_plots_static_field():
    
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
            if row[-2] == 'om':
                continue
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
            step_num_list = file['step_num_list'][()]
            at_list = step_num_list*tau
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
                    
        plt.plot(step_num_list, max_bond_list)
        plt.ylabel('Maximum bond dimension')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(step_num_list, avg_bond_list)
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
        # plt.legend()
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
        # plt.legend()
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
            
    df.to_csv('results_static.csv', index = False)

def make_plots_from_df_static_field():
    
    df = pd.read_csv('results_static.csv')
    
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
                        plt.savefig(f'Plots/pnd_vs_at_aD_0_{aD_0}_l_0_{l_0}_st.png')
                        plt.clf()

def make_plots_pulse_field():
    
    if max(N_list) <= 8:
        columns = ['N', 'x', 'ma', 'l_0_init', 'l_0', 'tau', 'env', 'sig', 'aT', 'lam', 'aD_0', 'at', 'D', 'ee', 'pn', 'l_0_small', 'omega', 'type', 'cutoff']
    else:
        columns = ['N', 'x', 'ma', 'l_0_init', 'l_0', 'tau', 'env', 'sig', 'aT', 'lam', 'aD_0', 'at', 'D', 'pn', 'l_0_small', 'omega', 'type', 'cutoff']
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
            
            file = h5py.File(f"{path_to_project}/DAGS/{project_number}/HDF5/{filepath}", 'r')
            
            row = filepath[:-3].strip().split('_')
            if row[-2] != 'om':
                continue
            N, tau, x, l_0, ma, env, sig, aT, lam, aD_0, l_0_init, cutoff, l_0_small, type, omega = row[1], row[3], row[5], row[7], row[9], row[11], row[13], row[15], row[17], row[19], row[21], row[23], row[25], row[26], row[28]
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
            l_0_small = float(l_0_small)
            omega = float(omega)
                    
            max_bond_list = file['max_bond_list'][()]
            at_list = file['at_list'][()]
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
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}, omega = {omega}, type = {type}, l_0_small = {l_0_small}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/D_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(at_list, avg_bond_list)
        plt.ylabel('Average bond dimension')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}, omega = {omega}, type = {type}, l_0_small = {l_0_small}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Avg_D_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        plt.plot(np.arange(len(avg_step_time)), avg_step_time)
        plt.ylabel('Average step time')
        plt.xlabel('Iteration')
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}, omega = {omega}, type = {type}, l_0_small = {l_0_small}')
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
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}, omega = {omega}, type = {type}, l_0_small = {l_0_small}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/Particle_Number_Density_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
        
        if max(N_list) <= 8:
            plt.plot(at_list, ee_list)
            plt.ylabel('Entanglement Entropy')
            plt.xlabel('Iteration')
            plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}, omega = {omega}, type = {type}, l_0_small = {l_0_small}')
            plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/ee_vs_iteration/{filepath}.png', bbox_inches = 'tight')
            plt.close()
        
        for idx in range(N):
            plt.plot(at_list, z_list[idx], label = f'Site {idx}')
        plt.ylabel('<Z_i>')
        plt.xlabel('Iteration')
        plt.legend()
        plt.title(f'N = {N}, tau = {tau}, x = {x}, l_0 = {l_0}, ma = {ma}, env = {env}, sig = {sig}, aT = {aT}, lam = {lam}, aD_0 = {aD_0}, omega = {omega}, type = {type}, l_0_small = {l_0_small}')
        plt.savefig(f'{path_to_project}/DAGS/{project_number}/Plots/z_vs_iteration/{filepath}.png', bbox_inches = 'tight')
        plt.close()
              
        for i in range(len(at_list)):
            
            if max(N_list) <= 8:
                at, D, ee, pn = at_list[i], max_bond_list[i], ee_list[i], particle_number_density[i]
                new_row = pd.DataFrame([[N, x, ma, l_0_init, l_0, tau, env, sig, aT, lam, aD_0, at, D, ee, pn, l_0_small, omega, type, cutoff]], columns = df.columns)
            else:
                at, D, pn = at_list[i], max_bond_list[i], particle_number_density[i]
                new_row = pd.DataFrame([[N, x, ma, l_0_init, l_0, tau, env, sig, aT, lam, aD_0, at, D, pn, l_0_small, omega, type, cutoff]], columns = df.columns)
            df = pd.concat([df, new_row], ignore_index = True)
            
    df.to_csv('results_pulse.csv', index = False)
                    
def make_plots_from_df_pulse_field():
    
    # for N in df.N.unique():
    #     for x in df.x.unique():
    #         for l_0 in df.l_0.unique():
    #             for ma in df.ma.unique():
    #                 for aD_0 in df.aD_0.unique():
    #                     for aT in df.aT.unique():
    #                         for cutoff in df.cutoff.unique():
    #                             for tau in df.tau.unique():
    #                                 for l_0_small in df.l_0_small.unique():
    #                                     for type in df.type.unique():
    #                                         for omega in df.omega.unique():
    
    df = pd.read_csv('results_pulse.csv')
    
    line_styles = ['-', '--', '-.', ':', '-', '--', '-.', ':', '--']
    colors = ['b', 'g', 'r', 'orange', 'black', 'cyan', 'purple', 'pink', 'blue']
    plt.rc('axes', prop_cycle=(cycler('color', colors) + cycler('linestyle', line_styles)))
    
    # Everything fixed except from tau in the legend
    def plot_tau():
        
        if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/changing_tau'):
            os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/changing_tau')
        
        for N in df.N.unique():
            for x in df.x.unique():
                for l_0 in df.l_0.unique():
                    for ma in df.ma.unique():
                        for aD_0 in df.aD_0.unique():
                            for aT in df.aT.unique():
                                for cutoff in df.cutoff.unique():
                                    for l_0_small in df.l_0_small.unique():
                                        for type in df.type.unique():
                                            for omega in df.omega.unique():
                                                for tau in df.tau.unique():
                        
                                                    df_tmp = df[(df.N == N) & (df.x == x) & (df.l_0 == l_0) & (df.ma == ma) & (df.aD_0 == aD_0) & (df.cutoff == cutoff) & (df.aT == aT) & (df.tau == tau) & (df.omega == omega) & (df.type == type) & (df.l_0_small == l_0_small)]
                                                    # print(df_tmp)
                                                    at_list, pn_list = np.array(df_tmp['at']), np.array(df_tmp.pn)
                                                    # print(pn_list)
                                                    # sorted_lists = sorted(zip(at_list, pn_list))
                                                    # at_list, pn_list = zip(*sorted_lists)
                                
                                                    try:                                                 
                                                        plt.plot(at_list, pn_list, label = f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, max = {max(pn_list)}, tau = {tau}, l_0_small = {l_0_small}, om = {omega}, type = {type}')
                                                    except:
                                                        print(f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, tau = {tau}, l0s = {l_0_small}, om = {omega}, type = {type} gave an error.')
                                                    plt.title('PND vs at')
                                            
                                                plt.legend(fontsize = 6)
                                                plt.savefig(f'Plots/changing_tau/pnd_vs_at_N_{N}_x_{x}_l_0_{l_0}_ma_{ma}_aD_0_{aD_0}_aT_{aT}_cutoff_{cutoff}_l0s_{l_0_small}_om_{omega}_type_{type}.png')
                                                plt.clf()
                                                plt.close()
                        
    # plot_tau()
    
    # Everything fixed except from cutoff in the legend
    def plot_cutoff():
        
        if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/changing_cutoff'):
            os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/changing_cutoff')
        
        for N in df.N.unique():
            for x in df.x.unique():
                for l_0 in df.l_0.unique():
                    for ma in df.ma.unique():
                        for aD_0 in df.aD_0.unique():
                            for aT in df.aT.unique():
                                for tau in df.tau.unique():
                                    for l_0_small in df.l_0_small.unique():
                                        for type in df.type.unique():
                                            for omega in df.omega.unique():
                                                for cutoff in df.cutoff.unique():
                        
                                                    df_tmp = df[(df.N == N) & (df.x == x) & (df.l_0 == l_0) & (df.ma == ma) & (df.aD_0 == aD_0) & (df.cutoff == cutoff) & (df.aT == aT) & (df.tau == tau) & (df.omega == omega) & (df.type == type) & (df.l_0_small == l_0_small)]
                                                    # print(df_tmp)
                                                    at_list, pn_list = np.array(df_tmp['at']), np.array(df_tmp.pn)
                                                    # print(pn_list)
                                                    # sorted_lists = sorted(zip(at_list, pn_list))
                                                    # at_list, pn_list = zip(*sorted_lists)
                                
                                                    try:                                                 
                                                        plt.plot(at_list, pn_list, label = f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, max = {max(pn_list)}, tau = {tau}, l0s = {l_0_small}, om = {omega}, type = {type}')
                                                    except:
                                                        print(f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, tau = {tau}, l0s = {l_0_small}, om = {omega}, type = {type} gave an error.')
                                                    plt.title('PND vs at')
                                            
                                                plt.legend(fontsize = 6)
                                                plt.savefig(f'Plots/changing_cutoff/pnd_vs_at_N_{N}_x_{x}_l_0_{l_0}_ma_{ma}_aD_0_{aD_0}_aT_{aT}_tau_{tau}_l0s_{l_0_small}_om_{omega}_type_{type}.png')
                                                plt.clf()
                                                plt.close()
                        
    # plot_cutoff()
                    
    def plot_type():
        
        if not os.path.exists(f'{path_to_project}/DAGS/{project_number}/Plots/changing_type'):
            os.makedirs(f'{path_to_project}/DAGS/{project_number}/Plots/changing_type')
            
        for N in df.N.unique():
            for x in df.x.unique():
                for l_0 in df.l_0.unique():
                    for ma in df.ma.unique():
                        for aD_0 in df.aD_0.unique():
                            for aT in df.aT.unique():
                                for tau in df.tau.unique():
                                    for l_0_small in df.l_0_small.unique():
                                        for cutoff in df.cutoff.unique():
                                            for omega in df.omega.unique():
                                                for type_field in df.type.unique():
                                                                    
                                                    df_tmp = df[(df.N == N) & (df.x == x) & (df.l_0 == l_0) & (df.ma == ma) & (df.aD_0 == aD_0) & (df.cutoff == cutoff) & (df.aT == aT) & (df.tau == tau) & (df.omega == omega) & (df.type == type_field) & (df.l_0_small == l_0_small)]
                                                    # print(df_tmp)
                                                    at_list, pn_list = np.array(df_tmp['at']), np.array(df_tmp.pn)
                                                    # print(pn_list)
                                                    # sorted_lists = sorted(zip(at_list, pn_list))
                                                    # at_list, pn_list = zip(*sorted_lists)
                                
                                                    try:                                              
                                                        plt.plot(at_list, pn_list, label = f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, max = {max(pn_list)}, tau = {tau}, l0s = {l_0_small}, om = {omega}, type = {type_field}')
                                                    except:
                                                        print(f'N = {N}, aD_0 = {aD_0}, cutoff = {cutoff}, aT = {aT}, tau = {tau}, l0s = {l_0_small}, om = {omega}, type = {type_field} gave an error.')
                                                    plt.title('PND vs at')
                                            
                                                plt.legend(fontsize = 6)
                                                plt.savefig(f'Plots/changing_type/pnd_vs_at_N_{N}_x_{x}_l_0_{l_0}_ma_{ma}_aD_0_{aD_0}_aT_{aT}_tau_{tau}_l0s_{l_0_small}_om_{omega}.png')
                                                plt.clf()
                                                plt.close()              
    
    plot_type()
                                                                     
# write_dag_static_field()
# write_dag_pulse_field()

# make_plots_static_field()
# make_plots_from_df_static_field()

# make_plots_pulse_field()
# make_plots_from_df_pulse_field()
