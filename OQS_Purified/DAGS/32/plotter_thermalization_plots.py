import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from matplotlib.pyplot import xticks, yticks

# Define the directory and file names
file_names = ['Plots/Thermalization_Heatmap/data/0.1.pickle', 'Plots/Thermalization_Heatmap/data/0.5.pickle', 'Plots/Thermalization_Heatmap/data/0.75.pickle', 'Plots/Thermalization_Heatmap/data/1.0.pickle']
letters = [r'$(a)$',r'$(b)$',r'$(c)$',r'$(d)$',r'$(e)$',r'$(f)$']

# Load the data from the pickle files
plots = []
for file_name in file_names:
    with open(file_name, 'rb') as f:
        data = pickle.load(f)
        plots.append(data)
        
# Create a figure and axes with specific positions for a nicer layout
fig = plt.figure(figsize=(12, 10))  # Adjust figure size

# Create subplots in a 2x2 grid
ax1 = fig.add_subplot(221)  # First plot in the top left
ax2 = fig.add_subplot(222)  # Second plot in the top right
ax3 = fig.add_subplot(223)  # Third plot in the bottom left
ax4 = fig.add_subplot(224)  # Fourth plot in the bottom right

# List of axes for convenience
axes = [ax1, ax2, ax3, ax4]

vmin = min([np.min(np.array(data)) for data in plots])
vmax = max([np.max(np.array(data)) for data in plots])

# Plot all heatmaps in subplots
for i, ax in enumerate(axes):  # We only have 5 plots
    
    sns.heatmap(plots[i], ax=ax, linewidths=0, vmin = vmin, vmax = vmax)  # Use the jet colormap
    
    labels = [item.get_text() for item in ax.get_xticklabels()]
    ax.set_xticklabels([str(round(float(label), 2)) for label in labels])
            
    ax.set_xlabel(r'$l_0$', fontsize=14)
    ax.set_ylabel(r'$D$', fontsize=14)
        
    # Set tick parameters for better readability
    ax.tick_params(labelsize=12)

    # Add subplot labels (a), (b), etc. using LaTeX formatting
    ax.text(1.1, -0.1, letters[i], fontsize=14, ha='center', transform=ax.transAxes)

# Adjust layout to ensure everything is properly contained
plt.tight_layout()

# Save the figure as a high-resolution PDF (DPI = 1200)
plt.savefig('Plots/Thermalization_Heatmap/thermalization_times_ma_combined.pdf', dpi=1200)
plt.close()