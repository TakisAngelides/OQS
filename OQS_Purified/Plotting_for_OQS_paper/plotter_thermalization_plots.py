import pickle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from matplotlib import rc
rc('text', usetex=True)

# Define the directory and file names
file_names = [
    '/Users/takisangelides/Downloads/data_32/0.1.pickle',
    '/Users/takisangelides/Downloads/data_32/0.5.pickle',
    '/Users/takisangelides/Downloads/data_32/0.75.pickle',
    '/Users/takisangelides/Downloads/data_32/1.0.pickle'
]
labels = [r'$(a)$', r'$(b)$', r'$(c)$', r'$(d)$']

# Load the data from the pickle files
plots = []
for file_name in file_names:
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
        plots.append(data)

# Determine global min and max for color scale across all plots
vmin = min(np.min(np.array(data)) for data in plots)
vmax = max(np.max(np.array(data)) for data in plots)

# Define contour levels for specific thermalization times
contour_levels = [30, 40, 50, 60, 70]

# Create a figure and axes with specific positions for a nicer layout
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

# Plot all heatmaps using imshow in each subplot
for i, ax in enumerate(axes):
    l_0_values = plots[i].columns.astype(float)  # assuming columns are l_0 values
    aD_values = plots[i].index.astype(float)     # assuming index represents aD values
    
    # Use linspace for evenly distributed ticks
    xticks_positions = np.linspace(0, len(l_0_values)-1, num=5)
    yticks_positions = np.linspace(0, len(aD_values)-1, num=5)
    
    ax.set_xticks(xticks_positions)
    ax.set_yticks(yticks_positions)
    ax.set_xticklabels(np.round(np.linspace(l_0_values[0], l_0_values[-1], num=5), 2))
    ax.set_yticklabels(np.round(np.linspace(aD_values[0], aD_values[-1], num=5), 2))

    # Display the heatmap
    im = ax.imshow(plots[i], vmin=vmin, vmax=vmax, cmap='magma', aspect='auto')

    # Add contour lines for specified thermalization times with labels
    if i == 0:
        contours = ax.contour(plots[i], levels=contour_levels[:-2], colors='cyan', linewidths=0.7)
    else:
        contours = ax.contour(plots[i], levels=contour_levels, colors='cyan', linewidths=0.7)
    ax.clabel(contours, fmt='%d', colors='cyan', fontsize=20)  # Label the contours

    # Set only bottom row x-label and first column y-label, remove ticks for others
    if i >= 2:
        ax.set_xlabel(r'$l_0$', fontsize=24)
    else:
        ax.set_xticklabels([])  # Remove x-axis tick labels if not bottom row

    if i % 2 == 0:
        ax.set_ylabel(r'$D$', fontsize=24)
    else:
        ax.set_yticklabels([])  # Remove y-axis tick labels if not first column
    
    ax.tick_params(labelsize=20)
    ax.text(0.1, 0.02, labels[i], fontsize=20, ha='right', va='bottom', transform=ax.transAxes)

    for spine in ax.spines.values():
        spine.set_visible(False)

plt.tight_layout()

# Add colorbar with consistent scale across subplots
cbar = fig.colorbar(im, ax=axes, orientation='horizontal', location='top', fraction=0.03)
cbar.ax.tick_params(labelsize=20)
cbar.set_label(r"$\mathcal{T}$", fontsize=24)

for spine in cbar.ax.spines.values():
    spine.set_visible(False)

# Adjust layout and save the figure
plt.savefig('thermalization_times_ma_combined.pdf', bbox_inches='tight', dpi=1200)
# plt.savefig('thermalization_times_ma_combined.png', bbox_inches='tight')
plt.close()
