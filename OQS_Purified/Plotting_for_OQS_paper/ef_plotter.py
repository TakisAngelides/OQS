import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc 
rc('text', usetex=True)

# Parameters
base_paths = {
    24: "/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_42",
    12: "/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_47"
}
masses = [0.1, 1.0]
l_0_values = [0.0, 0.5]
aD_values = [2.0, 5.0]
tau = 0.01
steps = 10000
t_over_a_list = [0] + list(tau * np.arange(1, steps))
letters = [
    r'($a$) $m = 0.1, l_0 = 0, D = 2$', r'($b$) $m = 0.1, l_0 = 0, D = 5$', 
    r'($c$) $m = 0.1, l_0 = 0.5, D = 2$', r'($d$) $m = 0.1, l_0 = 0.5, D = 5$', 
    r'($e$) $m = 1.0, l_0 = 0, D = 2$', r'($f$) $m = 1.0, l_0 = 0, D = 5$', 
    r'($g$) $m = 1.0, l_0 = 0.5, D = 2$', r'($h$) $m = 1.0, l_0 = 0.5, D = 5$', 
    r'($i$) $m = 0.1, l_0 = 0, D = 2$', r'($j$) $m = 0.1, l_0 = 0, D = 5$', 
    r'($k$) $m = 0.1, l_0 = 0.5, D = 2$', r'($l$) $m = 0.1, l_0 = 0.5, D = 5$', 
    r'($m$) $m = 1.0, l_0 = 0, D = 2$', r'($n$) $m = 1.0, l_0 = 0, D = 5$', 
    r'($o$) $m = 1.0, l_0 = 0.5, D = 2$', r'($p$) $m = 1.0, l_0 = 0.5, D = 5$'
]

# Collect all values to determine the min and max for color scale
all_values = []

for N in base_paths:
    for mass in masses:
        for l_0 in l_0_values:
            for aD in aD_values:
                file_name = f"EF_{mass}_{l_0}_{aD}.pickle"
                file_path = os.path.join(base_paths[N], file_name)
                
                # Load data
                try:
                    with open(file_path, "rb") as f:
                        values = pickle.load(f)
                        all_values.append(values)
                except FileNotFoundError:
                    print(f"File not found: {file_path}")
                    continue

# Determine global min and max for color scale
global_min = min(np.min(values) for values in all_values)
global_max = max(np.max(values) for values in all_values)

# Plot setup with constrained layout
fig, axes = plt.subplots(nrows=4, ncols=4, figsize=(15, 12), constrained_layout=True)
counter = 0

for i, N in enumerate([12, 24]):
    data_path = base_paths[N]
    for j, mass in enumerate(masses):
        for k, l_0 in enumerate(l_0_values):
            for l, aD in enumerate(aD_values):
                file_name = f"EF_{mass}_{l_0}_{aD}.pickle"
                file_path = os.path.join(data_path, file_name)
                
                try:
                    with open(file_path, "rb") as f:
                        values = pickle.load(f)
                except FileNotFoundError:
                    continue
                
                # Select the correct subplot
                row = 2 * i + j
                col = 2 * k + l
                ax = axes[row, col]
                
                # Plot heatmap with fixed color scale
                im = ax.imshow(values, aspect='auto', origin='lower',
                               extent=[t_over_a_list[0], t_over_a_list[-1], 0, values.shape[0]],
                               vmin=global_min, vmax=global_max, cmap='jet')
                
                # White text at bottom right for subplot label
                ax.text(0.005, 0.02, f"{letters[counter]}", color="white", ha="left", va="bottom", transform=ax.transAxes, fontsize=19.5)
                counter += 1

                # Configure axis tick labels only for first column and last row
                if col == 0:
                    ax.set_ylabel(r'$n$', fontsize=24)
                    ax.tick_params(axis='y', labelsize=24)
                    # if row < 2:  # No decimals for first two rows
                    #     ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: rf"${int(x)}$"))
                else:
                    ax.set_yticklabels([])
                    
                if row == 3:
                    ax.set_xlabel(r"$t$", fontsize=24)
                    ax.tick_params(axis='x', labelsize=24)
                    # if col < 2:  # No decimals for first two rows
                    #     ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: rf"${int(x)}$"))
                else:
                    ax.set_xticklabels([])
                
                for spine in ax.spines.values():
                    spine.set_visible(False)

# Add colorbar with fixed scale
cbar = fig.colorbar(im, ax=axes, orientation='horizontal', location='top', fraction=0.03)
cbar.ax.tick_params(labelsize=24)
cbar.set_label(r'$\Delta F(n)$', fontsize=24)

# Hide colorbar borders
for spine in cbar.ax.spines.values():
    spine.set_visible(False)

plt.savefig('ef_N_12_24.pdf', bbox_inches='tight', dpi=1200)
# plt.savefig('ef_N_12_24.png')
# plt.show()
