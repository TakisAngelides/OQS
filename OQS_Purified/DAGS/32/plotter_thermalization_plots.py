import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.pyplot import xticks, yticks
from matplotlib.ticker import FuncFormatter

def round_ticks(value, pos):
    return f'{round(value, 2)}'

# Define the directory and file names
file_names = [
    'Plots/Thermalization_Heatmap/data/0.1.pickle',
    'Plots/Thermalization_Heatmap/data/0.5.pickle',
    'Plots/Thermalization_Heatmap/data/0.75.pickle',
    'Plots/Thermalization_Heatmap/data/1.0.pickle'
]
letters = [r'$(a)$: $m = 0.1$', r'$(b)$: $m = 0.5$ ', r'$(c)$: $m = 0.75$', r'$(d)$: $m = 1.0$']

# Load the data from the pickle files
plots = []
for file_name in file_names:
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
        plots.append(data)

# Determine global min and max for color scale across all plots
vmin = min(np.min(np.array(data)) for data in plots)
vmax = max(np.max(np.array(data)) for data in plots)

# Create a figure and axes with specific positions for a nicer layout
fig = plt.figure(figsize=(12, 10))  # Adjust figure size

# Create subplots in a 2x2 grid
ax1 = fig.add_subplot(221)  # First plot in the top left
ax2 = fig.add_subplot(222)  # Second plot in the top right
ax3 = fig.add_subplot(223)  # Third plot in the bottom left
ax4 = fig.add_subplot(224)  # Fourth plot in the bottom right

# List of axes for convenience
axes = [ax1, ax2, ax3, ax4]

# Plot all heatmaps using imshow in each subplot
for i, ax in enumerate(axes):

    # Get l_0 and aD values from the data directly, assuming DataFrame
    l_0_values = plots[i].columns.astype(float)  # assuming columns are l_0 values
    aD_values = plots[i].index.astype(float)     # assuming index represents aD values
    
    # Use linspace for evenly distributed ticks
    xticks_positions = np.linspace(0, len(l_0_values)-1, num=5)
    yticks_positions = np.linspace(0, len(aD_values)-1, num=5)
    
    # Set the tick labels at the calculated positions
    ax.set_xticks(xticks_positions)
    ax.set_yticks(yticks_positions)
    ax.set_xticklabels(np.round(np.linspace(l_0_values[0], l_0_values[-1], num=5), 2))
    ax.set_yticklabels(np.round(np.linspace(aD_values[0], aD_values[-1], num=5), 2))
    
    # Display the heatmap
    im = ax.imshow(plots[i], vmin=vmin, vmax=vmax, cmap='magma')

    # Set axis labels and ticks for readability
    # ax.set_xlabel(r'$l_0$', fontsize=16)
    # ax.set_ylabel(r'$D$', fontsize=16)
    ax.tick_params(labelsize=16)

    # Add subplot labels (a), (b), etc.
    ax.set_title(letters[i], fontsize=16)
    # ax.text(1.15, 0.05, letters[i], fontsize=16, ha='right', va='bottom', transform=ax.transAxes)
    
    for spine in ax.spines.values():
        spine.set_visible(False)


fig.text(0.5, 0.04, r'$l_0$', ha='center', fontsize=16)
fig.text(0.13, 0.48, r'$D$', va='center', rotation='vertical', fontsize=16)

# Add colorbar with consistent scale across subplots
plt.subplots_adjust(hspace=0.3, wspace = -0.1)
cbar = fig.colorbar(im, ax=axes, orientation='horizontal', location='top', fraction=0.03)
cbar.ax.tick_params(labelsize=16)
cbar.set_label("Thermalization Time", fontsize=16)
# Remove borders from the colorbar
for spine in cbar.ax.spines.values():
    spine.set_visible(False)

# Adjust layout to ensure everything is properly contained
# plt.tight_layout()

# plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.1, wspace=-0.25, hspace=0.2)

# Save the figure as a high-resolution PDF (DPI = 1200)
plt.savefig('Plots/Thermalization_Heatmap/thermalization_times_ma_combined.pdf', bbox_inches = 'tight', dpi = 1200)
plt.close()
