import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from collections import defaultdict
from matplotlib import rc 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
rc('text', usetex=True)

# Set the directory where files are located
directory_path = '/Users/takisangelides/Downloads/Plotting_for_OQS_paper/mutual_info_data_50'

# Initialize a dictionary to store data, sorted by D, l_0, m, and T
data_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))

# Time array
time = np.linspace(1 * 0.05, 2000 * 0.05, 2000)

# Define colors and linestyles to use
colors = ['orange', 'purple', 'green', 'red', '#1f77b4']
linestyles = ['-', '--', ':', '-.', '-', '--', '-.', ':']

D_filter = "5"
m_filter = "1.0"
l_0_filter = "0.5"

# Iterate over all files in the specified directory
for filename in os.listdir(directory_path):
    
    if not "2_sites" in filename:
        continue
    
    # Extract D, m, l_0, and T values from filename
    D, m, l_0, T, _, _ = filename.strip().split('_')
    
    if D != D_filter or m != m_filter or l_0 != l_0_filter:
        continue

    # Construct full file path
    file_path = os.path.join(directory_path, filename)
    
    # Open the HDF5 file and retrieve the mutual_info_list data
    with h5py.File(file_path, 'r') as file:
        if 'mutual_info_list' in file:
            mutual_info_list = file['mutual_info_list'][:]
            data_dict[float(D)][float(l_0)][float(m)][float(T)] = mutual_info_list
        else:
            print(f"'mutual_info_list' not found in {filename}")

# Plotting
line_index = 0
fig, ax = plt.subplots(figsize=(15, 10))

# Create inset
ax_inset = inset_axes(ax, width="37%", height="37%", loc="lower right", bbox_to_anchor=(-0.05, 0.15, 1, 1), bbox_transform=ax.transAxes)

ax_inset2 = inset_axes(ax, width="37%", height="37%", loc="upper right", bbox_to_anchor=(-0.05, 0, 1, 1), bbox_transform=ax.transAxes)

# Sort and plot data with unique color and linestyle combinations

for D in sorted(data_dict):
    for l_0 in sorted(data_dict[D]):
        for m in sorted(data_dict[D][l_0]):
            for T in sorted(data_dict[D][l_0][m]):
                mutual_info_list = data_dict[D][l_0][m][T]
                label = f'$D$ = {D}, $m$ = {m}, $l_0$ = {l_0}, $T$ = {T:.1f}'
                
                # Use unique color and linestyle for each line
                color = colors[line_index % len(colors)]
                linestyle = linestyles[line_index % len(linestyles)]
                ax.plot(time, mutual_info_list, label=label, color=color, linestyle=linestyle)
                ax_inset2.plot(time[95:150], mutual_info_list[95:150], color=color, linestyle=linestyle)
                ax_inset.plot(time[600:625], mutual_info_list[600:625], color=color, linestyle=linestyle)
                ax_inset.set_xlabel(r"$t$", fontsize=24)
                ax_inset2.set_xlabel(r"$t$", fontsize=24)
                ax_inset2.set_ylabel(r"$\Delta I$", fontsize=24)
                ax_inset.set_ylabel(r"$\Delta I$", fontsize=24)
                
                line_index += 1

# Set labels
ax.set_xlabel(r"$t$", fontsize=24)
ax.set_ylabel(r"$\Delta I$", fontsize=24)

ax.tick_params(axis='both', direction='in', labelsize=24)
ax_inset.tick_params(axis='both', direction='in', labelsize=24)
ax_inset2.tick_params(axis='both', direction='in', labelsize=24)

ax_inset2.text(0.99, 0.95, r"$(a)$", transform=ax_inset2.transAxes, fontsize=24,
               verticalalignment='top', horizontalalignment='right')

ax_inset.text(0.99, 0.95, r"$(b)$", transform=ax_inset.transAxes, fontsize=24,
              verticalalignment='top', horizontalalignment='right')

# Customize legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=2, fontsize=19)

# Save the plot
plt.savefig('mutual_info_with_T.pdf', bbox_inches='tight', dpi=1200)
plt.close()
