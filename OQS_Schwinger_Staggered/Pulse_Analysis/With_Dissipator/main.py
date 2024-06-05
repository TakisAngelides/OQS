import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "xtick.direction": "in", "ytick.direction": "in"})
xlabel_font = 14
ylabel_font = 14

N_list = [4, 6]
aD_0_list = [0.01, 0.1]
aT_list = [5, 10]
l_0_list = [0, 0.02]
type_field_list = ['static', 'sauter']
omega_list = np.linspace(0, 2, 20)

def write_dag():
    
    if not os.path.exists(f'HDF5'):
        os.makedirs(f'HDF5')
            
    with open(f'main.dag', 'w') as f:

        f.write(f'CONFIG dagman.config\n') # This will contain DAGMAN_USE_DIRECT_SUBMIT = False to avoid obscure bugs of authentication

        # This defines the Job submit file inline (ie inside the dag submit file)
        counter = 0
        
        for N in N_list:
            
            for aD_0 in aD_0_list:
                
                for aT in aT_list:
                    
                    for l_0 in l_0_list:
                        
                        for type_field in type_field_list:
                            
                            if type_field == 'static':
                                omega_list_tmp = [0]
                            else:
                                omega_list_tmp = omega_list
                                
                            for omega in omega_list_tmp:
                                        
                                aD_0_text = str(aD_0).replace('.', '')
                                aT_text = str(aT).replace('.', '')
                                l_0_text = str(l_0).replace('.', '')
                                omega_text = str(omega).replace('.', '')
                                
                                job_name = f'N_{N}_aD0_{aD_0_text}_aT_{aT_text}_l0_{l_0_text}_omega_{omega_text}_type_field_{type_field}'
                                f.write(f'JOB ' + job_name + f' main.sub\n')
                                f.write(f'VARS ' + job_name + f' N="{N}" aD_0="{aD_0}" aT="{aT}" l_0="{l_0}" type_field="{type_field}" omega="{omega}"\n')
                                f.write('RETRY ' + job_name + ' 0\n')
                                counter += 1

                                                                        
    print(f'Number of jobs: {counter}')
    
def collect_data():
    
    cols = ['N', 'aD_0', 'aT', 'l_0', 'type_field', 'omega', 'at', 'spnd']
    df = pd.DataFrame(columns = cols)
    
    max_omega = {}
    
    for file_name in os.listdir('HDF5'):
        
        f = h5py.File('HDF5/' + file_name)
        
        row = file_name.strip().split('_')
                
        N, omega, aD_0, l_0, aT, type_field = row[1], row[3], row[5], row[7], row[9], row[11][:-3]

        at_list, pnd = np.array(f['at_list']), np.array(f['pnd'])
        
        spnd = pnd - pnd[0]
                
        for i in range(len(at_list)):
        
            df_tmp = pd.DataFrame([[N, aD_0, aT, l_0, type_field, omega, at_list[i], spnd[i]]], columns = cols)
                
            df = pd.concat([df, df_tmp], ignore_index = True)
    
    df.to_csv('data.csv', index=False)  # Setting index=False to avoid writing row indices
    
def plot_results():
    
    df = pd.read_csv('data.csv')   
    
    dict_data = {}
    dict_max = {}
        
    for N in df.N.unique():
        for aD_0 in df.aD_0.unique():
            for aT in df.aT.unique():
                for l_0 in df.l_0.unique():
                    for type_field in df.type_field.unique():
                            
                        dict_tmp = {} 
                                                    
                        for omega in df.omega.unique():
                                                                      
                            df_ref = df[(df.N == N) & (df.aD_0 == aD_0) & (df.aT == aT) & (df.l_0 == 0) & (df.type_field == type_field) & (df.omega == 0)]    
                            spnd_ref = list(df_ref['spnd'])
                            
                            df_tmp = df[(df.N == N) & (df.aD_0 == aD_0) & (df.aT == aT) & (df.l_0 == l_0) & (df.type_field == type_field) & (df.omega == omega)]   
                            at, spnd = list(df_tmp['at']), list(df_tmp['spnd'])
                            
                            try:
                                y_data = np.array(spnd)-np.array(spnd_ref)
                            except:
                                continue
                            
                            dict_data[f'N_{N}_aD0_{aD_0}_aT_{aT}_l0_{l_0}_type_{type_field}_omega_{omega}'] = [at, y_data]
                            dict_tmp[f'N_{N}_aD0_{aD_0}_aT_{aT}_l0_{l_0}_type_{type_field}_omega_{omega}'] = max(y_data)
                            
                        dict_max[f'N_{N}_aD0_{aD_0}_aT_{aT}_l0_{l_0}_type_{type_field}'] = max(dict_tmp, key = lambda k : dict_tmp[k]).split('_')[-1]

    fig4, ax4 = plt.subplots(1)
    fig6, ax6 = plt.subplots(1)
    
    for key, value in dict_max.items():
        
        N, aD_0, aT, l_0, type_field = key.strip().split('_')[1::2]
        
        if l_0 == '0.0' or aT == '5.0':
            continue
        
        max_omega = value
        
        dict_data_key = key + '_omega_' + max_omega
    
        if type_field == 'static':
            linestyle = '--'
            zorder = 1
            label = r'$N$' + f' = {N}' + r', $aD$' + f' = {aD_0}'
        else:
            linestyle = '-'
            zorder = 2
            label = r'$N$' + f' = {N}' + r', $aD$' + f' = {aD_0}' + r', $\omega a / ma$' + f' = {float(max_omega):.3f}'
            
        at, y_data = dict_data[dict_data_key]
        
        if aD_0 == '0.01':
            color = 'orange'
        else:
            color = '#1f77b4'
        
        if N == '4':
            ax4.plot(at, y_data, label = label, linestyle = linestyle, zorder = zorder, color = color)
        else:
            ax6.plot(at, y_data, label = label, linestyle = linestyle, zorder = zorder, color = color)
                    
    ax4.legend(fontsize = 5)
    ax6.legend(fontsize = 5)
    ax4.hlines(0, 0, 5, linestyle = ':', color = 'gray', alpha = 0.5)
    ax6.hlines(0, 0, 5, linestyle = ':', color = 'gray', alpha = 0.5)
    ax4.set_ylabel(r'Subtracted Particle Number Density (SPND)', fontsize = ylabel_font)
    ax4.set_xlabel(r'$t/a$', fontsize = xlabel_font)
    ax6.set_ylabel(r'Subtracted Particle Number Density (SPND)', fontsize = ylabel_font)
    ax6.set_xlabel(r'$t/a$', fontsize = xlabel_font)
    ax4.figure.savefig(f'Plots/spnd_vs_at_max_N_4.pdf', dpi = 3000, transparent = True, bbox_inches = 'tight')
    ax6.figure.savefig(f'Plots/spnd_vs_at_max_N_6.pdf', dpi = 3000, transparent = True, bbox_inches = 'tight')
    ax4.clear()
    ax6.clear()
    plt.clf()
    plt.close()
    
    max_y_data_aD_001_N_4 = []
    max_y_data_aD_01_N_4 = []
    max_y_data_aD_001_N_6 = []
    max_y_data_aD_01_N_6 = []
    omega_aD_001_N_4 = []
    omega_aD_01_N_4 = []
    omega_aD_001_N_6 = []
    omega_aD_01_N_6 = []

    # Iterate over the keys of the dictionary
    for key, value in dict_data.items():
        # Extracting parameters from the key
        parameters = key.split('_')
        aD_0 = float(parameters[3])
        N = int(parameters[1])
        l_0 = parameters[-5]
        aT = parameters[-7]
        type_field = parameters[-3]
        
        # Skip if conditions are met
        if l_0 == '0.0' or aT == '5.0' or type_field == 'static':
            continue
        
        # Extracting data
        omega = float(parameters[-1])
        y_data = value[1]
        
        # Update maximum y_data lists accordingly
        if aD_0 == 0.01 and N == 4:
            max_y_data_aD_001_N_4.append(max(y_data))
            omega_aD_001_N_4.append(omega)
        elif aD_0 == 0.1 and N == 4:
            max_y_data_aD_01_N_4.append(max(y_data))
            omega_aD_01_N_4.append(omega)
        elif aD_0 == 0.01 and N == 6:
            max_y_data_aD_001_N_6.append(max(y_data))
            omega_aD_001_N_6.append(omega)
        elif aD_0 == 0.1 and N == 6:
            max_y_data_aD_01_N_6.append(max(y_data))
            omega_aD_01_N_6.append(omega)

    # Sort the lists according to the corresponding omega list
    omega_aD_001_N_4, max_y_data_aD_001_N_4 = zip(*sorted(zip(omega_aD_001_N_4, max_y_data_aD_001_N_4)))
    omega_aD_01_N_4, max_y_data_aD_01_N_4 = zip(*sorted(zip(omega_aD_01_N_4, max_y_data_aD_01_N_4)))
    omega_aD_001_N_6, max_y_data_aD_001_N_6 = zip(*sorted(zip(omega_aD_001_N_6, max_y_data_aD_001_N_6)))
    omega_aD_01_N_6, max_y_data_aD_01_N_6 = zip(*sorted(zip(omega_aD_01_N_6, max_y_data_aD_01_N_6)))

    # Plotting the maximum of y_data vs omega for each condition
        
    # plt.plot(omega_aD_001_N_4, max_y_data_aD_001_N_4, label=r'$N$ = 4, $D$ = 0.01', color = 'r', linestyle = '-')
    # plt.plot(omega_aD_01_N_4, max_y_data_aD_01_N_4, label=r'$N$ = 4, $D$ = 0.1', color = 'g', linestyle = '--')
    plt.plot(omega_aD_001_N_6, max_y_data_aD_001_N_6, label=r'$N$ = 6, $aD$ = 0.01', color = 'orange', linestyle = '--')
    plt.plot(omega_aD_01_N_6, max_y_data_aD_01_N_6, label=r'$N$ = 6, $aD$ = 0.1', color = '#1f77b4', linestyle = '-')
    plt.legend()
    plt.ylabel(r'Maximum of SPND', fontsize = ylabel_font)
    plt.xlabel(r'$\omega a/ma$', fontsize = xlabel_font)
    plt.savefig(f'Plots/max_spnd_vs_omega.pdf', dpi = 3000, transparent = True, bbox_inches = 'tight')
        
    
# write_dag()

# collect_data()

plot_results()
