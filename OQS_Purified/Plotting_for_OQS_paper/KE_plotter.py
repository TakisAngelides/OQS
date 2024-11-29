import pickle
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('text', usetex=True)

data_dir = "/Users/takisangelides/Downloads/data/"
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

linestyle_counter = 0
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

for i, mass in enumerate(masses):
    ax = axes[i]
    inset_ax = inset_axes(ax, width="60%", height="60%", loc="upper right")  # Increase inset size
    
    # Annotate each main plot
    ax.text(0.01, 0.98, r"$(a)$" if i == 0 else r"$(b)$", transform=ax.transAxes,
            ha='left', va='top', fontsize=16)  # Move label to top left
    
    # Set labels for main and inset plots
    if i == 0:
        ax.set_ylabel("$\Delta K$", fontsize=24)
    ax.set_xlabel("$t$", fontsize=24)
    
    for l0 in l0_values:
        for aD in aD_values:
            file_pattern = f"KE_{mass}_{l0}_{aD}.pickle"
            file_path = os.path.join(data_dir, file_pattern)
            
            if os.path.exists(file_path):
                values = load_pickle(file_path)
                time = np.linspace(0, len(values), num_iters) * 0.01
                
                # Plot initial time segment on main axis
                if i == 0:
                    line, = ax.plot(time[:2000], values[:2000], label=fr"$l_0$={l0}, $D$={aD}", 
                                    linestyle=linestyles[linestyle_counter % len(linestyles)])
                else:
                    line, = ax.plot(time[:2000], values[:2000], label=None, 
                                linestyle=linestyles[linestyle_counter % len(linestyles)])
                
                # Inset plot showing the full time series
                inset_ax.plot(time[2000:], values[2000:], linestyle=linestyles[linestyle_counter % len(linestyles)])
                inset_ax.tick_params(axis='both', which='both', direction='in', labelsize=15)
                inset_ax.set_ylabel("$\Delta K$", fontsize=20)
                inset_ax.set_xlabel("$t$", fontsize=20)
                
                linestyle_counter += 1
            else:
                print(f"File not found: {file_path}")
                
    ax.tick_params(axis='both', which='both', direction='in', labelsize=20)

    # Add legend only for the first subplot
    # if i == 0:
    #     ax.legend(fontsize=16, ncol=2, bbox_to_anchor = (1.5, 1.2))

fig.legend(loc="upper center", fontsize=20, ncol=4, bbox_to_anchor=(0.53, 1))
plt.tight_layout()
# fig.supxlabel(r"$t$", fontsize=20, x=0.53)
fig.subplots_adjust(top=0.88)  # Adjust layout to fit legend
plt.savefig('KE_mass_0.1_1.0.pdf', dpi=1200, bbox_inches='tight')
# plt.savefig('KE_mass_0.1_1.0.png', bbox_inches = 'tight')
plt.close(fig)
