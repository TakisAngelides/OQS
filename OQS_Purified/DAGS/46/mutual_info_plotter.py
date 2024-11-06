import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from collections import defaultdict

# Set the directory where files are located
directory_path = '/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/46/Plots/mutual_info_data'

# Initialize a dictionary to store data, sorted by D, l_0, and m
data_dict = defaultdict(lambda: defaultdict(dict))

# Time array
time = np.linspace(1*0.05, 2000*0.05, 2000)

# Define colors and linestyles to use
colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray']
linestyles = ['-', '--', '-.', ':', '-', '--', '-.', ':']

# Iterate over all files in the specified directory
for filename in os.listdir(directory_path):
    # Extract D, m, and l_0 values from filename
    D, m, l_0 = filename.strip().split('_')
    l_0 = l_0[:-3]

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
fig, ax = plt.subplots(figsize=(10, 6))

# Create inset
ax_inset = ax.inset_axes([0.55, 0.45, 0.4, 0.5])
# ax_inset.set_xlim(time[0], time[24])
# ax_inset.set_ylim(-0.001, 0.003)

# Sort and plot data with unique color and linestyle combinations
line_index = 0
for D in sorted(data_dict):
    for l_0 in sorted(data_dict[D]):
        for m in sorted(data_dict[D][l_0]):
            mutual_info_list = data_dict[D][l_0][m]
            label = f'D = {D}, m = {m}, $l_0$ = {l_0}'
            
            # Use unique color and linestyle for each line
            color = colors[line_index % len(colors)]
            linestyle = linestyles[line_index % len(linestyles)]
            ax.plot(time, mutual_info_list, label=label, color=color, linestyle=linestyle)
            ax_inset.plot(time[:25], mutual_info_list[:25], color=color, linestyle=linestyle)
            ax_inset.set_xlabel("$t$", fontsize = 18)
            ax_inset.set_ylabel("SMI", fontsize = 16)
            
            line_index += 1

# Set labels
ax.set_xlabel("$t$", fontsize = 18)
ax.set_ylabel("Subtracted Mutual Information (SMI)", fontsize = 16)

ax.tick_params(axis='both', direction='in', labelsize=16)
ax_inset.tick_params(axis='both', direction='in', labelsize=16)

# Customize legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=3, fontsize = 13)

# Save the plot
plt.savefig('/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/46/Plots/mutual_info.pdf', bbox_inches = 'tight', dpi = 1200)
