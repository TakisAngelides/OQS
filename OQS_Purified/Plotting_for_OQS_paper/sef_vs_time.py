import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import rc 

rc('text', usetex=True)
plt.tight_layout()

# Parameters
base_paths = {
    12: "/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_47"
}
masses = [0.1, 1.0]
l_0_values = [0.0, 0.5]
aD_values = [2.0, 5.0]
tau = 0.01
steps = 10000
t_over_a_list = [0] + list(tau * np.arange(1, steps))

# Define colors and linestyles
colors = ['#1f77b4', 'g', 'r', 'c', 'm', 'y', 'k', 'orange']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':']

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

# Plot setup
plt.figure(figsize=(15, 10))
main_ax = plt.gca()

line_counter = 0

for i, N in enumerate([12]):
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
                
                # Plot with color and linestyle cycling
                main_ax.plot(
                    t_over_a_list, 
                    values[5], 
                    label=rf'$N = {N}, D = {aD}, l_0 = {l_0}, m = {mass}$',
                    color=colors[line_counter % len(colors)],
                    linestyle=linestyles[line_counter % len(linestyles)]
                )
                line_counter += 1

# Add inset for the first 100 points
inset_ax = inset_axes(main_ax, width="50%", height="50%", loc="upper right", borderpad=2.5)
line_counter = 0  # Reset counter for inset

for i, N in enumerate([12]):
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
                
                # Plot inset with same color and linestyle
                inset_ax.plot(
                    t_over_a_list[:10], 
                    values[5][:10], 
                    color=colors[line_counter % len(colors)],
                    linestyle=linestyles[line_counter % len(linestyles)]
                )
                line_counter += 1

# Customize ticks, labels, and legend
main_ax.tick_params(axis='both', direction='in', labelsize=24)
inset_ax.tick_params(axis='both', direction='in', labelsize=24)

main_ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, fontsize=18)
main_ax.set_ylabel(r'$\Delta F(n = 5)$', fontsize=24)
main_ax.set_xlabel(r'$t$', fontsize=24)
inset_ax.set_ylabel(r'$\Delta F(n = 5)$', fontsize=20)
inset_ax.set_xlabel(r'$t$', fontsize=20)
inset_ax.yaxis.get_offset_text().set_fontsize(20)  # Set to desired font size

# Save the plot
plt.savefig('sef_vs_time.pdf', bbox_inches='tight', dpi=1200)
# plt.show()
