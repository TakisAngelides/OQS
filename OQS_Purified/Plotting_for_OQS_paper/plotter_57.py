import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc 

# Increase font size globally
plt.rcParams.update({
    'font.size': 16,  # Base font size
    'axes.titlesize': 18,  # Title font size
    'axes.labelsize': 18,  # Axes label font size
    'xtick.labelsize': 16,  # X-tick font size
    'ytick.labelsize': 16,  # Y-tick font size
    'legend.fontsize': 16,  # Legend font size
    'figure.titlesize': 20  # Figure title font size
})

rc('text', usetex=True)

# Define the directory path
directory_path = "/Users/takisangelides/Downloads/Plotting_for_OQS_paper/data_57"  

# Initialize an empty dictionary to store the data
data_dict = {}

# Iterate over all files in the directory
for filename in os.listdir(directory_path):
    # Construct the full path to the file
    file_path = os.path.join(directory_path, filename)
    
    # Ensure it's a file (not a directory)
    if os.path.isfile(file_path) and filename.endswith('.pickle'):
        key_name = filename.replace('.pickle', '')
        
        # Load the pickle file and store its content in the dictionary
        with open(file_path, 'rb') as f:
            data_dict[key_name] = np.real(pickle.load(f))

N = 14
D_list = [2.0, 2.75, 3.5, 4.25, 5.0]
t = np.linspace(0, 2000*0.05, 2001)
fraction = 0.3

def get_index_of_reduced_value(fraction, l):
    target = l[0]*fraction
    for element_idx, element in enumerate(l):
        if element < target:
            return element_idx
    return -1

# Define the values of D for which we want to create subplots
D_panel = [2.0, 3.5, 5.0]

# Create a figure with 2 rows and 3 columns of subplots
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Flatten the axes array for easier indexing
axes = axes.ravel()

keep = 150  # Limit the time range for visualization

# Define linestyles for lines
linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]

text_list = ['$(a)$', '$(b)$', '$(c)$']

# Iterate through specified values of D for EF panels
for idx, D in enumerate(D_panel):
    
    EF = data_dict[f'EF_{D}']
    axes[idx].set_xlabel('$t$')
    axes[idx].set_ylabel('$n$')
    axes[idx].text(6.7, 11.7, rf'${text_list[idx]}$', fontsize=20)
    
    im = axes[idx].imshow(
        EF[:, :keep],
        aspect='auto',
        origin='lower',
        extent=[t[0], t[keep-1], 0, EF[:, :keep].shape[0]],
        cmap='jet'  # Ensure colormap is specified
    )
    
    for spine in axes[idx].spines.values():
        spine.set_visible(False)
        
# Add a single colorbar for the EF panels with global vmin and vmax
cbar_ax = fig.add_axes([0.31, 0.9, 0.4, 0.03])  # Define the location and size of the colorbar
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', location = 'top')
cbar.set_label('$\Delta F_{B}(n)$')

for spine in cbar_ax.spines.values():
    spine.set_visible(False)

# Adjust subplot layout
plt.subplots_adjust(hspace=0.4, wspace=0.25)  # Increase vertical and horizontal spacing

# Subplot 3: Thermalization times with error bars
thermalization_times = []
for D in D_list:
    ef_middle = data_dict[f'EF_{D}'][N//2-1, :]
    thermalization_times.append(t[get_index_of_reduced_value(fraction, ef_middle)])
axes[3].errorbar(
    D_list, 
    thermalization_times, 
    markerfacecolor='none',
    yerr=0.1, 
    fmt='-o', 
    label='$\mathcal{T}_{\Delta F_B}$',
    markersize = 15
)

thermalization_times = []
for D in D_list:
    E = data_dict[f'E_{D}']
    thermalization_times.append(t[get_index_of_reduced_value(fraction, E)])
axes[3].errorbar(
    D_list, 
    thermalization_times, 
    yerr=0.1, 
    fmt='d--', 
    markerfacecolor='none', 
    label='$\mathcal{T}_{\Delta E_B}$',
    markersize = 15
)
axes[3].set_xlabel('$D$')
axes[3].set_ylabel('$\mathcal{T}_{\Delta F_B}$, $\mathcal{T}_{\Delta E_B}$')
# axes[3].text(4.8, 2.5, r'$(d)$', fontsize=20)
axes[3].text(0.97, 0.13, r'$(d)$', transform=axes[3].transAxes, fontsize=20, va='bottom', ha='right')
axes[3].legend()

# Subplots 4 and 5: Middle EF and Energy with linestyles and insets
handles = []
labels = []
start = 0
for i, D in enumerate(D_list):
    ef_middle = data_dict[f'EF_{D}'][N//2-1, :]
    line, = axes[4].plot(
        t[start:], ef_middle[start:], linestyle=linestyles[i], label=f'$D = {D}$'
    )
    handles.append(line)
    labels.append(f'$D = {D}$')

    E = data_dict[f'E_{D}']
    line, = axes[5].plot(
        t[start:], E[start:], linestyle=linestyles[i], label=f'$D = {D}$'
    )
    handles.append(line)

# Add insets for full evolution
start_inset = 1000
for idx, (ax, data_key) in enumerate(zip([axes[4], axes[5]], ['EF', 'E'])):
    inset = ax.inset_axes([0.3, 0.3, 0.65, 0.65])
    for i, D in enumerate(D_list):
        inset.plot(
            t[start_inset:], data_dict[f'{data_key}_{D}'][N//2-1, start_inset:] if data_key == 'EF' else data_dict[f'{data_key}_{D}'][start_inset:],
            linestyle=linestyles[i]
        )
    inset.set_xlabel('$t$', fontsize=14)
    inset.set_ylabel('$\Delta F_{B}(n = 6)$' if data_key == 'EF' else '$\Delta E_{B}$', fontsize=14)
    inset.tick_params(axis='both', which='both', labelsize=12, direction = 'in')

# Add a common legend for Subplots 4 and 5
axes[4].legend(
    loc='upper center', bbox_to_anchor=(1.12, 1.2), ncol=5, frameon=True
)

axes[4].set_xlabel('$t$')
axes[4].set_ylabel('$\Delta F_{B}(n = 6)$')
# axes[4].text(94.0, 0.01, r'$(e)$', fontsize=20)
axes[4].text(0.97, 0.13, r'$(e)$', transform=axes[4].transAxes, fontsize=20, va='bottom', ha='right')

axes[5].set_xlabel('$t$')
axes[5].set_ylabel('$\Delta E_{B}$')
# axes[5].text(0.9, 0.1, r'$(f)$', fontsize=20)
axes[5].text(0.97, 0.13, r'$(f)$', transform=axes[5].transAxes, fontsize=20, va='bottom', ha='right')


for i in [3, 4, 5]:
    axes[i].tick_params(axis = "both", direction = "in")

# Save the figure with manual colorbar positioning
plt.savefig('plot_57.pdf', bbox_inches='tight', dpi=1200)
plt.close()

# Create a new figure for the PN and KE plots in one row
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5))

# Plot PN observable for all D values as 'a'
for D in D_list:
    PN = data_dict[f'PN_{D}']
    ax1.plot(t, PN, linestyle='-', label=f'$D = {D}$')

ax1.set_xlabel('$t$')
ax1.set_ylabel('$\Delta \mathcal{P}_B$')
ax1.text(0.95, 0.965, r'$(a)$', transform=ax1.transAxes, fontsize=20, va='top', ha='left')

# Plot KE observable for all D values as 'b'
for D in D_list:
    KE = data_dict[f'KE_{D}']
    ax2.plot(t, KE, linestyle='-')

ax2.set_xlabel('$t$')
ax2.set_ylabel('$\Delta K_B$')
ax2.text(0.95, 0.965, r'$(b)$', transform=ax2.transAxes, fontsize=20, va='top', ha='left')


# Add a common legend for both subplots with 5 columns in a single row
fig.legend(
    loc='upper center', 
    bbox_to_anchor=(0.5, 1.1), 
    ncol=5, 
    frameon=True
)

# Add insets similar to the previous plot for PN and KE observables
start_inset = 800
for ax, data_key in zip([ax1, ax2], ['PN', 'KE']):
    ax.tick_params(axis = 'both', direction = 'in')
    inset = ax.inset_axes([0.25, 0.3, 0.65, 0.65])
    for i, D in enumerate(D_list):
        data = data_dict[f'{data_key}_{D}'][start_inset:] if data_key == 'PN' else data_dict[f'{data_key}_{D}'][start_inset:]
        inset.plot(t[start_inset:], data, linestyle=linestyles[i])
    inset.set_xlabel('$t$', fontsize=14)
    if data_key == 'PN':
        data_key_label = '$\Delta \mathcal{P}_B$'
    else:
        data_key_label = '$\Delta K_B$'
    inset.set_ylabel(f'{data_key_label}', fontsize=14)
    inset.tick_params(axis='both', which='both', labelsize=14, direction = 'in')
    

# Adjust layout for better visibility and add space for legends
plt.tight_layout()
plt.savefig('plot_57_pn_ke.pdf', bbox_inches='tight', dpi=1200)


