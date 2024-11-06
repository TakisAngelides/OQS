import pickle
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data_dir = "Plots/data/"
masses = [0.1, 1.0]
observables = ["E", "EFE", "KE", "ME", "Q", "EF"]
l0_values = [0.0, 0.5]
aD_values = [2.0, 5.0]
num_iters = 10000
linestyles = ['-', '--', '-.', ':']

def load_pickle(file_path):
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data

def get_index_of_reduced_value(fraction, l):
    target = l[0] * fraction
    for element_idx, element in enumerate(l):
        if element < target:
            return element_idx
    return -1

project_number = os.getcwd().strip().split('/')[-1]
path_to_project_number = f'/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/{project_number}'
if not os.path.exists(f'{path_to_project_number}/Plots/combined_plots_from_individual_observables'):
    os.makedirs(f'{path_to_project_number}/Plots/combined_plots_from_individual_observables')

linestyle_counter = 0
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

for i, mass in enumerate(masses):
    ax = axes[i]
    inset_ax = inset_axes(ax, width="40%", height="40%", loc="upper right")
    
    # Annotate each main plot
    if i == 0:
        ax.text(0.99, 0.5, r"$(a)$", transform=ax.transAxes, ha='right', va='top', fontsize=16)
        ax.set_ylabel("Subtracted Kinetic Energy (SKE)", fontsize=18)
    else:
        ax.text(0.99, 0.5, r"$(b)$", transform=ax.transAxes, ha='right', va='top', fontsize=16)
    
    for l0 in l0_values:
        for aD in aD_values:
            file_pattern = f"KE_{mass}_{l0}_{aD}.pickle"
            file_path = os.path.join(data_dir, file_pattern)
            
            if os.path.exists(file_path):
                values = load_pickle(file_path)
                time = np.linspace(0, len(values), num_iters) * 0.01
                
                # Plot initial time segment on main axis
                ax.plot(time[:2000], values[:2000], label=fr"$l_0$={l0}, $D$={aD}", linestyle=linestyles[linestyle_counter % len(linestyles)])
                
                # Inset plot showing the full time series
                inset_ax.plot(time, values, linestyle=linestyles[linestyle_counter % len(linestyles)])
                inset_ax.tick_params(axis='both', which='both', direction='in', labelsize=12)
                inset_ax.set_ylabel("SKE", fontsize=16)
                inset_ax.set_xlabel("$t$", fontsize=16)
                
                linestyle_counter += 1
            else:
                print(f"File not found: {file_path}")
                
    ax.tick_params(axis='both', which='both', direction='in', labelsize=16)
    if i == 0:
        ax.legend(loc="upper left", fontsize=16)

plt.tight_layout()
fig.supxlabel(r"$t$", fontsize=20, x=0.53)
fig.subplots_adjust(bottom=0.1) 
output_path = os.path.join("Plots/combined_plots_from_individual_observables/", "KE_mass_0.1_1.0.pdf")
plt.savefig(output_path, dpi=1200)
plt.close(fig)
