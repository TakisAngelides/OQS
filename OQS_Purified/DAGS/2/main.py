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
    path_to_project_number = f'/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/{project_number}'
    path_to_inputs_h5 = path_to_project_number + '/inputs.h5'
    f_h5 = h5py.File(path_to_inputs_h5, 'w')
        
    # Create relevant folders if needed
    if not os.path.exists(f'{path_to_project_number}/Plots'):
        os.makedirs(f'{path_to_project_number}/Plots')
    if not os.path.exists(f'{path_to_project_number}/States'):
        os.makedirs(f'{path_to_project_number}/States')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Error'):
        os.makedirs(f'{path_to_project_number}/Logs/Error')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Output'):
        os.makedirs(f'{path_to_project_number}/Logs/Output')        
    if not os.path.exists(f'{path_to_project_number}/Logs/Log'):
        os.makedirs(f'{path_to_project_number}/Logs/Log')
    if not os.path.exists(f'{path_to_project_number}/Observables'): # here we will store observables in JLD format
        os.makedirs(f'{path_to_project_number}/Observables')
            
    # This will form the job id
    counter_of_jobs = 1
    
    # Static applied field case and delta correlator
    lambd = 0
    number_of_time_steps_list = [800, 1600, 3200] # needs same length as tau list
    tau_list = [0.01, 0.005, 0.0025]
    taylor_expansion_cutoff_1 = 1e-13
    taylor_expansion_cutoff_2 = 1e-13
    maxdim = 500
    how_many_states_to_save = 20
    which_applied_field = "constant" # options are: "constant", "sauter", "gaussian", "oscillatory"
    time_varying_applied_field_flag = "false" if which_applied_field == "constant" else "true"
    env_corr_type = "delta" # options are: "constant", "delta", "gaussian"
    for N in [8, 12]:
        dissipator_sites = [i for i in range(1, N+1)]
        flip_sites = [N//2-1, N//2 + 2] # this is for the case when the initial state is the dirac vacuum with a string and specifies where the string should be placed
        for x in [1/(0.36)**2]:
            for ma in [0.17]:
                for aD in [1]:
                    for aT in [10]:
                        for tau_idx, tau in enumerate(tau_list):
                            # Below we define a list of the step numbers at which we want to save the state
                            number_of_time_steps = number_of_time_steps_list[tau_idx]
                            step = number_of_time_steps // how_many_states_to_save
                            which_steps_to_save_state = np.arange(0, number_of_time_steps+1, step)
                            which_steps_to_save_state[0] = 1
                            for cutoff in [1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]:
                                for taylor_expansion_order in [4]:
                                    for l_0_1 in [0]: # this is the constant part of the applied field
                                        for conserve_qns in ["true", "false"]:
                                            for which_initial_state in ["dirac_vacuum", "dirac_vacuum_with_string"]: # options are: "dirac_vacuum", "gs_naive", "dirac_vacuum_with_string"
                                        
                                                # Memory, CPU and maximum number of days to run
                                                mem, cpu, days = 8, 8, 3
                                                
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
                    
                                                # Write job to dag
                                                job_name = f'{job_id}'
                                                f_dag.write(f'JOB ' + job_name + f' {path_to_sub}\n')
                                                f_dag.write(f'VARS ' + job_name + f' job_id="{job_id}" path_to_project_number="{path_to_project_number}" file_to_run="{file_to_run}" cpu="{cpu}" mem="{mem}" days="{days}"\n')
                                                f_dag.write('RETRY ' + job_name + ' 1\n')
        
    # Close the dag file and the h5 input file
    f_dag.close() 
    f_h5.close()
    print(f'Total number of jobs in the dag is {counter_of_jobs-1}')                                        

write_dag()
