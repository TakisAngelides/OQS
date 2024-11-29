import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from collections import defaultdict

rc('text', usetex=True)

# Define the directories and file names
thermalization_dir = 'Thermalization_Heatmap'
file_names = ['/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_32/0.1.pickle', '/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_32/1.0.pickle']
labels = [r'm=0.1', r'm=1.0']

# Load the data from the pickle files for the first two subplots
plots = []
for file_name in file_names:
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
        plots.append(data)

# Dictionary to store thermalization times for the third subplot
data_dict = defaultdict(dict)

# Load thermalization time vs ma data for third subplot
for filename in os.listdir(thermalization_dir):
    if filename.endswith('.pickle'):
        parts = filename.split('_')
        ma = float(parts[0])
        l_0_1 = float(parts[1])
        aD = float(parts[2].replace('.pickle', ''))
        
        filepath = os.path.join(thermalization_dir, filename)
        with open(filepath, 'rb') as f:
            thermalization_time = pickle.load(f)
        
        data_dict[(l_0_1, aD)][ma] = thermalization_time

# Define a colorblind-friendly color palette and linestyles
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
linestyles = ['-'] #, '--', '-.', ':']
marker_list = ['o', 's', '^', 'D', '>', 'v', '*', '<', 'p', 'h']
line_combinations = [(color, linestyle) for color in colors for linestyle in linestyles]

# Create the subplots with 3 columns
fig, axes = plt.subplots(1, 3, figsize=(12, 6))
plt.tight_layout()

# Collect y-axis limits across all plots for uniform ticks
all_thermalization_times = []

# First subplot: thermalization time vs D for each mass and l_0 combination
counter = 0
for i, df in enumerate(plots):
    for j, l0 in enumerate(plots[i].columns.astype(float)):
        D_values = df.index
        thermalization_times = df[l0]
        all_thermalization_times.extend(thermalization_times)  # Gather y-values for scaling

        if j != 0 and j != 19 and j != 10:
            continue
        
        axes[0].errorbar(D_values, thermalization_times, yerr=0.1, fmt=marker_list[counter],
                         linestyle=linestyles[counter % len(linestyles)], color=colors[counter % len(colors)],
                         label=rf'${labels[i]}, l_0 = {l0:.2f}$', markerfacecolor='none')        
        counter += 1

axes[0].set_xlabel(r'$D$', fontsize=18)
axes[0].set_ylabel(r'$\mathcal{T}$', fontsize=18)
axes[0].legend(fontsize = 10, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
axes[0].tick_params(labelsize=18)

# Second subplot: thermalization time vs l_0 for D = 2, 5 and mass = 0.1, 1.0
counter = 0
for i, df in enumerate(plots):
    for k, D in enumerate(plots[i].index.astype(float)):
        l0_values = df.columns
        thermalization_times = df.loc[D, :]
        all_thermalization_times.extend(thermalization_times)  # Gather y-values for scaling

        if k != 0 and k != 19 and k != 10:
            continue
        
        axes[1].errorbar(l0_values, thermalization_times, yerr=0.1, fmt=marker_list[counter],
                         linestyle=linestyles[counter % len(linestyles)], color=colors[counter % len(colors)],
                         label=f'${labels[i]}, D = {D:.2f}$', markerfacecolor='none')
        counter += 1

axes[1].set_xlabel(r'$l_0$', fontsize=18)
axes[1].legend(fontsize = 10, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2)
axes[1].tick_params(labelsize=18)
axes[1].set_yticklabels([])  # Remove y-axis labels for second subplot
# axes[1].set_yscale('log')
# axes[1].set_xscale('log')

# Determine common y-axis limits and ticks based on data
y_min, y_max = min(all_thermalization_times)-1, max(all_thermalization_times)+1
y_ticks = np.linspace(y_min, y_max, num=6)

# Apply the determined limits and ticks to both y-axes
for ax in axes:
    ax.set_ylim(y_min, y_max)
    ax.set_yticks(y_ticks)

# Third subplot: thermalization time vs ma for each (l_0, aD) pair, sorted by l_0 and then aD
# Sort data_dict items by l_0 and then aD
sorted_data_dict = sorted(data_dict.items(), key=lambda x: (x[0][0], x[0][1]))

# Plot lines in sorted order with distinct color and linestyle
colors = ['#1f77b4', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'k']
for idx, ((l_0_1, aD), ma_dict) in enumerate(sorted_data_dict):
    ma_values = sorted(ma_dict.keys())
    thermalization_times = [ma_dict[ma] for ma in ma_values]
    color, linestyle = line_combinations[idx % len(line_combinations)]
    label = f"$l_0={l_0_1:.2f}, aD={aD:.2f}$"
    axes[2].errorbar(ma_values, thermalization_times, yerr=0.1, fmt=marker_list[idx],
                     label=label, color=colors[idx], linestyle=linestyle, markerfacecolor='none')

axes[2].set_xlabel(r"$m$", fontsize=18)
axes[2].legend(loc="upper center", bbox_to_anchor=(0.5, 1.2215), ncol=2, fontsize = 10)
axes[2].set_yticklabels([])  # Remove y-axis labels for third subplot
axes[2].tick_params(labelsize=18)

annotations = [r'$(a)$', r'$(b)$', r'$(c)$']
for ax, annotation in zip(axes, annotations):
    ax.text(0.97, 0.08, annotation, transform=ax.transAxes, 
            fontsize=20, fontweight='bold', va='top', ha='right')

fig.subplots_adjust(wspace=0.05)  # Decrease wspace for tighter horizontal spacing

# Adjust layout and show/save the plot
plt.savefig('/Users/takisangelides/Downloads/Plotting_for_OQS_paper/thermalization_vs_individual.pdf', bbox_inches='tight', dpi = 1200)
# plt.show()