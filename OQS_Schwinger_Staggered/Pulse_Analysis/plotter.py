import numpy as np
import h5py
import os
import matplotlib.pyplot as plt

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"], "xtick.direction": "in", "ytick.direction": "in"})
xlabel_font = 14
ylabel_font = 14

l_0 = 0.02
l_0_small = 0.5*l_0
ma = 1
e = 1

d = {}
ds = {}

omega_4_list = []
max_4_list = []

omega_6_list = []
max_6_list = []

c_list = {4 : 'r', 6 : 'g', 8 : 'b', 10 : 'k'}

for file in os.listdir('h5_data'):
    
    f = h5py.File('h5_data/'+file, 'r')
    at_list, pnd = np.array(f['at_list']), np.array(f['pnd'])
    spnd = pnd - pnd[0]
    row = file.strip().split('_')
    
    if len(row) == 4:
        
        N, omega = row[1], row[3][:-3]
        
        # plt.plot(at_list, spnd, label = f'N = {N}, omega = {float(omega):.3f}')
        
        d[f'N_{N}_om_{omega}'] = [at_list, spnd, max(spnd)]
        
        if N == '4':
            omega_4_list.append(float(omega))
            max_4_list.append(float(max(spnd)))
        else:
            omega_6_list.append(float(omega))
            max_6_list.append(float(max(spnd)))
        
    else:
        
        N = row[1][:-3]
        
        # plt.plot(at_list, spnd, label = f'N = {N}, max = {max(spnd):.3f}')
        
        ds[f'N_{N}'] = [at_list, spnd, max(spnd)]

max_om_4 = max(d, key = lambda k: d[k][2] if ('N_4_' in k) else 0).split('_')[-1]
max_om_6 = max(d, key = lambda k: d[k][2] if ('N_6_' in k) else 0).split('_')[-1]

for key, value in d.items():
    
    row = key.split('_')
    N, omega = row[1], row[3]
    at_list, spnd, max_val = value
    if (N == '4' and omega == max_om_4) or (N == '6' and omega == max_om_6):
        linestyle = '-'
        plt.plot(at_list, spnd, c = c_list[int(N)], label = r'$N$ = '+ f'{N}' + r', $\frac{\omega a}{ma}$ = ' + f'{float(omega):.3f}', linestyle = linestyle)
        
for key, value in ds.items():
    
    row = key.split('_')
    N = row[1]
    at_list, spnd, max_val = value
    linestyle = '--'
    plt.plot(at_list, spnd, c = c_list[int(N)], label = r'$N$ = ' + f'{N}', linestyle = linestyle)
    
path_tmp = '/Users/takisangelides/Documents/PhD/Project_3_OQS/OQS/Closed_Schwinger_Staggered/HDF5'
for file in os.listdir(path_tmp):
    row = file.strip().split('_')
    # try:
    if len(row) != 12: 
        continue
    N, type_field, tau, cutoff, order, omega = row[1::2]
    omega = omega[:-3]
    if (N != '8') and (N != '10'):
        continue
    if N == '8':
        if (tau != '0.001') or (cutoff != '1.0e-9'):
            continue
    if N == '10': 
        if (tau != '0.0005') or (cutoff != '1.0e-13'):
            continue
    f = h5py.File(f'{path_tmp}/'+file, 'r')
    at_list, pnd = np.array(f['at_list']), np.array(f['pnd'])
    spnd = pnd - pnd[0]
    if type_field == 'sauter':
        linestyle = '-'
        plt.plot(at_list, spnd, c = c_list[int(N)], label = r'$N$ = '+ f'{N}' + r', $\frac{\omega a}{ma}$ = ' + f'{float(omega):.3f}', linestyle = linestyle)
    else:
        linestyle = '--'
        plt.plot(at_list, spnd, c = c_list[int(N)], label = r'$N$ = '+ f'{N}', linestyle = linestyle)
    
plt.hlines(0, 0, max(at_list), color = 'gray', linestyle = ':', alpha = 0.5)
plt.legend(loc = 'lower right', fontsize = 5)
plt.ylabel(r'Subtracted Particle Number Density (SPND)', fontsize = ylabel_font)
plt.xlabel(r'$t/a$', fontsize = xlabel_font)
plt.savefig('SPND_vs_at.pdf', dpi = 3000, transparent = True, bbox_inches = 'tight')
plt.close()

sorted_tuples = sorted(zip(omega_4_list, max_4_list))
omega_4_list, max_4_list = zip(*sorted_tuples)

sorted_tuples = sorted(zip(omega_6_list, max_6_list))
omega_6_list, max_6_list = zip(*sorted_tuples)

plt.plot(omega_6_list, max_6_list, label = r'$N = 6$')
plt.plot(omega_4_list, max_4_list, label = r'$N = 4$', linestyle = '--')
# plt.vlines(omega_6_list[np.argmax(max_6_list)], 0, max(max_6_list), linestyle = ':', alpha = 0.5, color = 'blue')
# plt.vlines(omega_4_list[np.argmax(max_4_list)], 0, max(max_4_list), linestyle = ':', alpha = 0.5, color = 'orange')
plt.legend()
plt.ylabel(r'Maximum of SPND', fontsize = ylabel_font)
plt.xlabel(r'$\omega a/ma$', fontsize = xlabel_font)
plt.savefig('max_vs_omega.pdf', dpi = 3000, transparent = True, bbox_inches = 'tight')
