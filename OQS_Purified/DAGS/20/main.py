import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd
from cycler import cycler
import seaborn as sns

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

# Project number of dag
project_number = os.getcwd().strip().split('/')[-1]
path_to_project_number = f'/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/{project_number}'

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
    
    # Where to find the inputs
    path_to_inputs_h5 = path_to_project_number + '/inputs.h5'
    f_h5 = h5py.File(path_to_inputs_h5, 'w')
        
    # Create relevant folders if needed
    if not os.path.exists(f'{path_to_project_number}/Plots'):
        os.makedirs(f'{path_to_project_number}/Plots')
    if not os.path.exists(f'{path_to_project_number}/HDF5'):
        os.makedirs(f'{path_to_project_number}/HDF5')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Error'):
        os.makedirs(f'{path_to_project_number}/Logs/Error')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Output'):
        os.makedirs(f'{path_to_project_number}/Logs/Output')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Log'):
        os.makedirs(f'{path_to_project_number}/Logs/Log')
            
    # This will form the job id
    counter_of_jobs = 1
    
    # Static applied field case and delta correlator
    lambd = 0
    number_of_time_steps_list = [1000] # needs same length as tau list
    tau_list = [0.05]*len(number_of_time_steps_list)
    aD_list = [1.0]
    x_list = [7.716]
    ma_list = [0.17]
    taylor_expansion_cutoff_1 = 1e-9
    taylor_expansion_cutoff_2 = 1e-9
    maxdim = 300
    how_many_states_to_save = 10
    which_applied_field = "constant" # options are: "constant", "sauter", "gaussian", "oscillatory"
    time_varying_applied_field_flag = "false" if which_applied_field == "constant" else "true"
    env_corr_type = "delta" # options are: "constant", "delta", "gaussian"
    for N in [12]:
        dissipator_sites = [i for i in range(1, N+1)]
        flip_sites = [N//2-1, N//2 + 2] # this is for the case when the initial state is the dirac vacuum with a string and specifies where the string should be placed
        for aT in [10]:
            for x in x_list:
                for ma in ma_list:
                    for aD in aD_list:
                        for tau_idx, tau in enumerate(tau_list):
                            # Below we define a list of the step numbers at which we want to save the state
                            number_of_time_steps = number_of_time_steps_list[tau_idx]
                            if how_many_states_to_save == 0:
                                which_steps_to_save_state = []
                            else: 
                                step = number_of_time_steps // how_many_states_to_save
                                which_steps_to_save_state = np.arange(0, number_of_time_steps+1, step)
                                which_steps_to_save_state[0] = 1
                            for cutoff in [1e-9]:
                                for taylor_expansion_order in [2]:
                                    for l_0_1 in [-2.5, -1.5, -1.25, -0.99, -0.8, -0.5, 0.2, 0.0, 0.2, 0.5, 0.8, 0.99, 1.25, 1.5, 2.5]: # this is the constant part of the applied field
                                        for conserve_qns in ["true"]:
                                            for which_initial_state in ["dirac_vacuum", "dirac_vacuum_with_string"]: # options are: "dirac_vacuum", "gs_naive", "dirac_vacuum_with_string"
                                        
                                                # Memory, CPU and maximum number of days to run
                                                mem, cpu, days = 16, 8, 6.99
                                                
                                                # Job id for the dag job names and path to h5 for results
                                                job_id = counter_of_jobs
                                                counter_of_jobs += 1 # after assigning the job_id this is incremented for the next job
                                                                                            
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
                                                g.attrs["wstss"] = which_steps_to_save_state
                                                g.attrs["mem"] = mem
                    
                                                # Write job to dag
                                                job_name = f'{job_id}'
                                                f_dag.write(f'JOB ' + job_name + f' {path_to_sub}\n')
                                                f_dag.write(f'VARS ' + job_name + f' job_id="{job_id}" path_to_project_number="{path_to_project_number}" file_to_run="{file_to_run}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                                                f_dag.write('RETRY ' + job_name + ' 2\n')
        
    # Close the dag file and the h5 input file
    f_dag.close() 
    f_h5.close()
    print(f'Total number of jobs in the dag is {counter_of_jobs-1}')                                        

def plot_bond_dimensions():
    
    if not os.path.exists(f'{path_to_project_number}/Plots/D'):
        os.makedirs(f'{path_to_project_number}/Plots/D')
        
    path_to_HDF5 = f'{path_to_project_number}/HDF5'
    path_to_inputs = f'{path_to_project_number}/inputs.h5'
    f_inputs = h5py.File(path_to_inputs, 'r')

    counter = 0

    for file in os.listdir(path_to_HDF5):
        
        try:
        
            g = f_inputs[file[:-3]]
            attributes_dict = {attr_name: attr_value for attr_name, attr_value in g.attrs.items()}
            N = attributes_dict['N']
            l_0_1 = attributes_dict['l_0_1']
            f = h5py.File(f'{path_to_HDF5}/{file}', 'r')
            link_dims = np.asarray(f['link_dims'])
            max_link_dims = np.max(link_dims, axis = 0)
            avg_link_dims = np.average(link_dims, axis = 0)
            
            plt.plot(avg_link_dims, label = "Avg")
            plt.plot(max_link_dims, label = "Max")
            ma, aD, aT, cqns, cutoff, l_0_1, teo, waf, x_val, tau = attributes_dict['ma'], attributes_dict['aD'], attributes_dict['aT'], attributes_dict['cqns'], attributes_dict['cutoff'], attributes_dict['l_0_1'], attributes_dict['teo'], attributes_dict['waf'], np.round(attributes_dict['x'], decimals = 3), attributes_dict["tau"]
            plt.title(f'N_{N}_x_{x_val}_ma_{ma}_aD_{aD}_aT_{aT}_qn_{cqns}\nc_{cutoff}_l01_{l_0_1}_tau_{tau}_taylor_{teo}_waf_{waf}')
            plt.xlabel('Iteration step')
            plt.ylabel('Bond dimension')
            plt.legend()
            plt.savefig(f'Plots/D/{file[:-3]}.png')
            plt.close()
                        
        except:
            
            counter += 1
            # print(file[:-3], attributes_dict)
        
        # staggering = np.array([(-1)**n for n in range(N)])
        # sns.heatmap(link_dims)
        # plt.plot(pn)
        # sns.heatmap(ef_configs)
        # z_configs = np.asarray(f['z_configs'])
        # q_configs = np.array([0.5*(np.real(z_configs[:, i]) + staggering) for i in range(z_configs.shape[1])])
        # ef_configs = np.transpose(np.array([np.array([l_0_1 + sum(q_configs[i][0:j]) for j in range(q_configs.shape[1])]) for i in range(q_configs.shape[0])]))
        # pn = np.array([0.5*N + 0.5*sum(np.real(z_configs[:, i]) * staggering) for i in range(z_configs.shape[1])])

    if counter > 0:
        print('From plotting the bond dimension number of jobs that failed: ', counter)        

def plot_subtracted_observables():
    
    if not os.path.exists(f'{path_to_project_number}/Plots/EF'):
        os.makedirs(f'{path_to_project_number}/Plots/EF')
        
    if not os.path.exists(f'{path_to_project_number}/Plots/EF_Middle'):
        os.makedirs(f'{path_to_project_number}/Plots/EF_Middle')
        
    if not os.path.exists(f'{path_to_project_number}/Plots/PN'):
        os.makedirs(f'{path_to_project_number}/Plots/PN')
        
    if not os.path.exists(f'{path_to_project_number}/Plots/Q'):
        os.makedirs(f'{path_to_project_number}/Plots/Q')
                
    path_to_HDF5 = f'{path_to_project_number}/HDF5'
    
    path_to_inputs = f'{path_to_project_number}/inputs.h5'
    
    f_inputs = h5py.File(path_to_inputs, 'r')

    for file in os.listdir(path_to_HDF5):
        
        try:
        
            g = f_inputs[file[:-3]]
            attributes_dict = {attr_name: attr_value for attr_name, attr_value in g.attrs.items()}
            
            if attributes_dict['wis'] == 'dirac_vacuum_with_string':

                N = attributes_dict['N']
                l_0_1 = attributes_dict['l_0_1']
                tau = attributes_dict['tau']
                staggering = np.array([(-1)**n for n in range(N)])
                f = h5py.File(f'{path_to_HDF5}/{file}', 'r')
                z_configs = np.asarray(f['z_configs'])
                q_configs = np.array([0.5*(np.real(z_configs[:, i]) + staggering) for i in range(z_configs.shape[1])])
                ef_configs = np.transpose(np.array([np.array([l_0_1 + sum(q_configs[i][0:j + 1]) for j in range(q_configs.shape[1] - 1)]) for i in range(q_configs.shape[0])]))
                pn = np.array([0.5*N + 0.5*sum(np.real(z_configs[:, i]) * staggering) for i in range(z_configs.shape[1])])
                f.close()
            
                file_without_string = f'{int(file[:-3])-1}'
                attributes_dict_without_string = {attr_name: attr_value for attr_name, attr_value in f_inputs[file_without_string].attrs.items()}

                N = attributes_dict_without_string['N']
                l_0_1 = attributes_dict_without_string['l_0_1']
                staggering = np.array([(-1)**n for n in range(N)])
                f = h5py.File(f'{path_to_HDF5}/{file_without_string}.h5', 'r')
                z_configs_without_string = np.asarray(f['z_configs'])
                q_configs_without_string = np.array([0.5*(np.real(z_configs_without_string[:, i]) + staggering) for i in range(z_configs_without_string.shape[1])])
                ef_configs_without_string = np.transpose(np.array([np.array([l_0_1 + sum(q_configs_without_string[i][0:j + 1]) for j in range(q_configs_without_string.shape[1] - 1)]) for i in range(q_configs_without_string.shape[0])]))
                pn_without_string = np.array([0.5*N + 0.5*sum(np.real(z_configs_without_string[:, i]) * staggering) for i in range(z_configs_without_string.shape[1])])
                f.close()
                
                z = -ef_configs + ef_configs_without_string
                t_over_a_list = [0] + list(tau*(np.arange(1, z.shape[1])))
                x = np.round(t_over_a_list, decimals = 3)
                y = list(np.arange(1, N))
                
                # x = x[:z.shape[1]//2 + 50]
                # z = z[:, :z.shape[1]//2 + 50]
                  
                sns.heatmap(z, cmap = 'jet', vmin = -1, vmax = 1, yticklabels = y)
                num_xticks_to_display = 10
                step_size = z.shape[1] // num_xticks_to_display
                x_tick_positions = np.arange(0.5, z.shape[1], step_size, dtype = int)
                x_tick_labels = x[x_tick_positions]
                plt.xticks(x_tick_positions, x_tick_labels)
                plt.ylabel('Link')
                plt.xlabel(r'$t/a$')
                ma, aD, aT, cqns, cutoff, l_0_1, teo, waf, x_val = attributes_dict['ma'], attributes_dict['aD'], attributes_dict['aT'], attributes_dict['cqns'], attributes_dict['cutoff'], attributes_dict['l_0_1'], attributes_dict['teo'], attributes_dict['waf'], np.round(attributes_dict['x'], decimals = 3)
                plt.title(f'N_{N}_x_{x_val}_ma_{ma}_aD_{aD}_aT_{aT}_qn_{cqns}\nc_{cutoff}_l01_{l_0_1}_tau_{tau}_taylor_{teo}_waf_{waf}')
                plt.savefig(f'Plots/EF/{file[:-3]}.png')
                plt.close()
                
                for i in range(z.shape[0]):
                    plt.plot(x, z[i, :], label = f'{i}, Max = {max(z[i, :]):.5f}')
                plt.legend()
                plt.title(f'N_{N}_x_{x_val}_ma_{ma}_aD_{aD}_aT_{aT}_qn_{cqns}\nc_{cutoff}_l01_{l_0_1}_tau_{tau}_taylor_{teo}_waf_{waf}')
                plt.ylabel('Middle electric field')
                plt.xlabel(r'$t/a$')
                plt.savefig(f'Plots/EF_Middle/{file[:-3]}.png')
                plt.close()
                
                z = np.transpose(-q_configs + q_configs_without_string)
                t_over_a_list = [0] + list(tau*(np.arange(1, z.shape[1])))
                x = np.round(t_over_a_list, decimals = 3)
                y = list(np.arange(1, N))
                
                # if (file[:-3] == "2") or (file[:-3] == "4"):
                #     x = x[:z.shape[1]//2 - 200]
                #     z = z[:, :z.shape[1]//2 - 200]
                            
                sns.heatmap(z, cmap = 'jet', vmin = -1, vmax = 1, yticklabels = y)
                num_xticks_to_display = 10
                step_size = z.shape[1] // num_xticks_to_display
                x_tick_positions = np.arange(0.5, z.shape[1], step_size, dtype = int)
                x_tick_labels = x[x_tick_positions]
                plt.xticks(x_tick_positions, x_tick_labels)
                plt.ylabel('Site')
                plt.xlabel(r'$t/a$')
                ma, aD, aT, cqns, cutoff, l_0_1, teo, waf, x_val = attributes_dict['ma'], attributes_dict['aD'], attributes_dict['aT'], attributes_dict['cqns'], attributes_dict['cutoff'], attributes_dict['l_0_1'], attributes_dict['teo'], attributes_dict['waf'], np.round(attributes_dict['x'], decimals = 3)
                plt.title(f'N_{N}_x_{x_val}_ma_{ma}_aD_{aD}_aT_{aT}_qn_{cqns}\nc_{cutoff}_l01_{l_0_1}_tau_{tau}_taylor_{teo}_waf_{waf}')
                plt.savefig(f'Plots/Q/{file[:-3]}.png')
                plt.close()
                
                # if (file[:-3] == "2") or (file[:-3] == "4"):
                #     pn = pn[:len(x)]
                #     pn_without_string = pn_without_string[:len(x)]
                
                plt.plot(x, pn - pn[0], label = "string")
                plt.plot(x, pn_without_string - pn_without_string[0], label = "no string")
                plt.plot(x, -pn + pn_without_string + pn[0] - pn_without_string[0], label = f"subtracted, max = {max(-pn + pn_without_string + pn[0] - pn_without_string[0])}")
                plt.title(f'N_{N}_x_{x_val}_ma_{ma}_aD_{aD}_aT_{aT}_qn_{cqns}\nc_{cutoff}_l01_{l_0_1}_tau_{tau}_taylor_{teo}_waf_{waf}')
                plt.legend()
                plt.ylabel('Particle number')
                plt.xlabel(r'$t/a$')
                plt.savefig(f'Plots/PN/{file[:-3]}.png')
                plt.close()
                            
            else:
                
                continue
            
        except:
            
            print(file)
            
# write_dag()

# plot_bond_dimensions()

plot_subtracted_observables()
