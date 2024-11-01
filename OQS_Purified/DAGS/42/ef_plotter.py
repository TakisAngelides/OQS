import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# Parameters
base_paths = {
    24: "/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/42/Plots/data",
    12: "/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/47/Plots/data"
}
masses = [0.1, 1.0]
l_0_values = [0.0, 0.5]
aD_values = [2.0, 5.0]
tau = 0.01
steps = 10000
t_over_a_list = [0] + list(tau * np.arange(1, steps))
letters = [r'$a$', r'$b$', r'$c$', r'$d$', r'$e$', r'$f$', r'$g$', r'$h$', r'$i$', r'$j$', r'$k$', r'$l$', r'$m$', r'$n$', r'$o$', r'$p$']

# Collect all values to determine the min and max for color scale
all_values = []

for N in base_paths:
    for mass in masses:
        for l_0 in l_0_values:
            for aD in aD_values:
                # Construct file path
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
                # Construct file path
                file_name = f"EF_{mass}_{l_0}_{aD}.pickle"
                file_path = os.path.join(data_path, file_name)
                
                # Load data
                try:
                    with open(file_path, "rb") as f:
                        values = pickle.load(f)
                except FileNotFoundError:
                    continue
                
                # Select the correct subplot
                ax = axes[2 * i + j, 2 * k + l]
                
                # Plot heatmap with fixed color scale
                im = ax.imshow(values, aspect='auto', origin='lower',
                               extent=[t_over_a_list[0], t_over_a_list[-1], 0, values.shape[0]],
                               vmin=global_min, vmax=global_max, cmap='jet')
                
                # Set subplot title including N
                ax.set_title(f"({letters[counter]}): $l_0$ = {l_0}, $D$ = {aD}", fontsize=16)
                if k == 0 and l == 0:
                    ax.set_ylabel(f'N = {N}, m = {mass}', fontsize = 16)
                counter += 1

                # Keep ticks for readability
                ax.tick_params(labelsize=16)
                
                for spine in ax.spines.values():
                    spine.set_visible(False)

# Set common x and y labels
fig.text(0.515, -0.02, r'$t$', ha='center', fontsize=16)
fig.text(-0.02, 0.45, 'Link', va='center', rotation='vertical', fontsize=16)

# Add colorbar with fixed scale
# fig.subplots_adjust(top=0.9) 
cbar = fig.colorbar(im, ax=axes, orientation='horizontal', location='top', fraction=0.03)
cbar.ax.tick_params(labelsize=16)  # Change this value as needed
cbar.set_label("SEF", fontsize=16)
# Remove borders from the colorbar
for spine in cbar.ax.spines.values():
    spine.set_visible(False)

# Adjust layout
# fig.tight_layout(rect=[0.05, 0.05, 1, 0.95])

plt.savefig('Plots/ef_N_12_24.pdf', bbox_inches = 'tight', dpi = 1200)
# plt.savefig('Plots/ef_N_12_24.png')
