import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from collections import defaultdict
from matplotlib import rc 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('text', usetex=True)

# Set the directory where files are located
directory_path = '/Users/takisangelides/Downloads/Plotting_for_OQS_paper/mutual_info_data'

# Initialize a dictionary to store data, sorted by D, l_0, and m
data_dict = defaultdict(lambda: defaultdict(dict))

# Time array
time = np.linspace(1*0.05, 2000*0.05, 2000)

# Define colors and linestyles to use
colors = ['#1f77b4', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':']

# Iterate over all files in the specified directory
for filename in os.listdir(directory_path):
    
    if not "2_sites" in filename:
        continue
    
    # Extract D, m, and l_0 values from filename
    D, m, l_0, _, _ = filename.strip().split('_')
    # l_0 = l_0[:-3]

    # Construct full file path
    file_path = os.path.join(directory_path, filename)
    
    # Open the HDF5 file and retrieve the mutual_info_list data
    with h5py.File(file_path, 'r') as file:
        if 'mutual_info_list' in file:
            mutual_info_list = file['mutual_info_list'][:]
            data_dict[D][l_0][m] = mutual_info_list
        else:
            print(f"'mutual_info_list' not found in {filename}")

# Plotting
fig, ax = plt.subplots(figsize=(15, 10))

# Create inset
ax_inset = inset_axes(ax, width="70%", height="70%", loc="upper right", bbox_to_anchor=(-0.01, -0.01, 1, 1), bbox_transform=ax.transAxes)
# ax_inset.set_xlim(time[0], time[24])
# ax_inset.set_ylim(-0.001, 0.003)

# Sort and plot data with unique color and linestyle combinations
line_index = 0
for D in sorted(data_dict):
    for l_0 in sorted(data_dict[D]):
        for m in sorted(data_dict[D][l_0]):
            mutual_info_list = data_dict[D][l_0][m]
            label = f'$D$ = {D}, $m$ = {m}, $l_0$ = {l_0}'
            
            # Use unique color and linestyle for each line
            color = colors[line_index % len(colors)]
            linestyle = linestyles[line_index % len(linestyles)]
            ax.plot(time, mutual_info_list, label=label, color=color, linestyle=linestyle)
            ax_inset.plot(time[30:250], mutual_info_list[30:250], color=color, linestyle=linestyle)
            ax_inset.set_xlabel(r"$t$", fontsize = 24)
            ax_inset.set_ylabel(r"$\Delta I$", fontsize = 24)
            
            line_index += 1


# plt.tight_layout()

# Set labels
ax.set_xlabel(r"$t$", fontsize = 24)
ax.set_ylabel(r"$\Delta I$", fontsize = 24)

ax.tick_params(axis='both', direction='in', labelsize=24)
ax_inset.tick_params(axis='both', direction='in', labelsize=24)

# Customize legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3, fontsize = 19)

# Save the plot
plt.savefig('mutual_info.pdf', bbox_inches = 'tight', dpi = 1200)

# plt.savefig('mutual_info_2_sites.png', bbox_inches = 'tight')
plt.close()
