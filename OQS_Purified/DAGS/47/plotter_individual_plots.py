import pickle
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import seaborn as sn
from matplotlib import rc

# Define the directory and observables
data_dir = "Plots/data/"
# masses = [0.1, 0.25, 0.5, 0.75, 1.0]
masses = [0.1, 1.0]
observables = ["E", "EFE", "KE", "ME", "Q", "EF"]  # Added "Q" and "EF"
# l0_values = [0.0, 0.26315789, 0.5]
# aD_values = [2.0, 3.57894737, 5.0]
l0_values = [0.0, 0.5]
aD_values = [2.0, 5.0]
num_iters = 10000

# List of linestyles to cycle through for line plots
linestyles = ['-', '--', '-.', ':']

# Function to load pickle data
def load_pickle(file_path):
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    return data

def get_index_of_reduced_value(fraction, l):
        
    target = l[0]*fraction
    for element_idx, element in enumerate(l):
        if element < target:
            return element_idx
    return -1

project_number = os.getcwd().strip().split('/')[-1]
path_to_project_number = f'/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/{project_number}'
if not os.path.exists(f'{path_to_project_number}/Plots/combined_plots_from_individual_observables'):
    os.makedirs(f'{path_to_project_number}/Plots/combined_plots_from_individual_observables')

# Iterate over each mass and observable to plot the data
# for mass in masses:
#     for observable in observables:
        
#         if observable in ["Q", "EF"]:  # Handle heatmap cases with subplots
#             fig, axes = plt.subplots(4, 4, figsize=(12, 10))  # 2x2 grid of subplots
            
#             # Flatten axes for easier iteration
#             axes = axes.flatten()
            
#             subplot_counter = 0  # Counter for subplot index
            
#             # Iterate over all combinations of l0 and aD
#             for l0 in l0_values:
#                 for aD in aD_values:
#                     # Construct the file pattern
#                     file_pattern = f"{observable}_{mass}_{l0}_{aD}.pickle"
#                     file_path = os.path.join(data_dir, file_pattern)
                    
#                     if os.path.exists(file_path):
#                         # Load the data (assuming it's a 2D array: sites vs time)
#                         values = load_pickle(file_path)
                        
#                         fraction = 0.3
#                         N = 12
#                         tau = 0.01
#                         steps = 10000
#                         t_over_a_list = [0] + list(tau*(np.arange(1, steps)))
#                         thermalization_time = t_over_a_list[get_index_of_reduced_value(fraction, values[N//2-1,:])]
                        
#                         # if aD == 2.0 and mass == 1.0 and observable == 'EF':
#                         #     print(mass, aD, l0, thermalization_time)
                        
#                         if observable == 'EF':
#                             print(f'D = {aD}, l_0 = {l0}, m = {mass}, t = {thermalization_time}')

#                         # Plot the heatmap
#                         im = axes[subplot_counter].imshow(values, aspect='auto', origin='lower')
#                         axes[subplot_counter].set_title(f"l0={l0}, aD={aD}")
#                         axes[subplot_counter].set_xlabel('Time')
#                         axes[subplot_counter].set_ylabel('Site')
#                         subplot_counter += 1
#                     else:
#                         print(f"File not found: {file_path}")
            
#             # Adjust layout and add colorbar
#             fig.suptitle(f"{observable} Heatmap for Mass {mass}", fontsize=16)
#             fig.tight_layout(rect=[0, 0, 0.9, 0.96])
#             cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Position for colorbar
#             fig.colorbar(im, cax=cbar_ax)
            
#             # Save the figure
#             output_path = os.path.join("Plots/combined_plots_from_individual_observables/", f"{observable}_mass_{mass}_heatmap.png")
#             plt.savefig(output_path)
#             plt.close()  # Close the figure to avoid display if running in non-interactive environments
        
#         else:  # Line plots for other observables
#             plt.figure(figsize=(10, 6))  # Create a new figure for each mass and observable
            
#             # Initialize a counter to apply different linestyles
#             linestyle_counter = 0
            
#             # Iterate over all combinations of l0 and aD
#             for l0 in l0_values:
#                 for aD in aD_values:
#                     # Construct the file pattern
#                     file_pattern = f"{observable}_{mass}_{l0}_{aD}.pickle"
#                     file_path = os.path.join(data_dir, file_pattern)
                    
#                     if os.path.exists(file_path):
#                         # Load the data
#                         values = load_pickle(file_path)
                        
#                         time = np.linspace(0, len(values), 10000)
                        
#                         # Plot the data with a different linestyle
#                         plt.plot(time, values, label=f"l0={l0}, aD={aD}", linestyle=linestyles[linestyle_counter % len(linestyles)])
                        
#                         # Increment the counter
#                         linestyle_counter += 1
#                     else:
#                         print(f"File not found: {file_path}")
            
#             # Customize the plot
#             plt.title(f"{observable} vs Time for Mass {mass}")
#             plt.xlabel("Time")
#             plt.ylabel(f"{observable}")
#             plt.legend()
#             plt.grid(True)
            
#             # Save the figure
#             output_path = os.path.join("Plots/combined_plots_from_individual_observables/", f"{observable}_mass_{mass}.png")
#             plt.savefig(output_path)
#             plt.close()  # Close the figure to avoid display if running in non-interactive environments
            
#     print('---------------------------------------------------')


# Create a figure with two subplots in a row
linestyle_counter = 0
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

for i, mass in enumerate(masses):
    ax = axes[i]
    
    for l0 in l0_values:
        for aD in aD_values:
            
            # Construct the file pattern
            file_pattern = f"KE_{mass}_{l0}_{aD}.pickle"
            file_path = os.path.join(data_dir, file_pattern)
            
            if os.path.exists(file_path):
                # Load the data
                values = load_pickle(file_path)
                time = np.linspace(0, len(values), num_iters)*0.01
                
                # print(mass, l0, aD, max(values))
                # print('----------------')
                
                # Plot on the specific subplot
                ax.plot(time[:2000], values[:2000], label=fr"$l_0$={l0}, $D$={aD}", linestyle=linestyles[linestyle_counter % len(linestyles)])
                
                # Increment the counter
                linestyle_counter += 1
            else:
                print(f"File not found: {file_path}")
                
    # Customize the subplot
    # ax.set_title(f"KE vs Time for Mass {mass}")
    # ax.set_xlabel("t")
    
    if i == 0:
        ax.text(0.95, 0.95, r"$(a)$", transform=ax.transAxes, ha='right', va='top', fontsize=16)
    else:
        ax.text(0.95, 0.95, r"$(b)$", transform=ax.transAxes, ha='right', va='top', fontsize=16)

    ax.tick_params(axis='both', which='both', direction='in', labelsize=10)  # Tick bars inwards, larger font
    if i == 0:
        ax.set_ylabel("SKE", fontsize=16)
        ax.legend(loc="upper left", fontsize=16)
    # ax.grid(True)

# Adjust layout and save the figure
plt.tight_layout()
fig.supxlabel(r"$t$", fontsize=16)
fig.subplots_adjust(bottom=0.1) 
output_path = os.path.join("Plots/combined_plots_from_individual_observables/", "KE_mass_0.1_1.0.pdf")
plt.savefig(output_path, dpi=1200)
plt.close(fig)


# ---------------------------------------------------
# D = 2.0, l_0 = 0.0, m = 0.1, t = 23.32
# D = 3.57894737, l_0 = 0.0, m = 0.1, t = 40.89
# D = 5.0, l_0 = 0.0, m = 0.1, t = 56.56
# D = 2.0, l_0 = 0.26315789, m = 0.1, t = 26.150000000000002
# D = 3.57894737, l_0 = 0.26315789, m = 0.1, t = 42.9
# D = 5.0, l_0 = 0.26315789, m = 0.1, t = 58.2
# D = 2.0, l_0 = 0.5, m = 0.1, t = 29.72
# D = 3.57894737, l_0 = 0.5, m = 0.1, t = 45.1
# D = 5.0, l_0 = 0.5, m = 0.1, t = 59.72
# ---------------------------------------------------
# D = 2.0, l_0 = 0.0, m = 1.0, t = 43.660000000000004
# D = 3.57894737, l_0 = 0.0, m = 1.0, t = 54.14
# D = 5.0, l_0 = 0.0, m = 1.0, t = 66.98
# D = 2.0, l_0 = 0.26315789, m = 1.0, t = 48.68
# D = 3.57894737, l_0 = 0.26315789, m = 1.0, t = 56.75
# D = 5.0, l_0 = 0.26315789, m = 1.0, t = 68.86
# D = 2.0, l_0 = 0.5, m = 1.0, t = 53.980000000000004
# D = 3.57894737, l_0 = 0.5, m = 1.0, t = 59.99
# D = 5.0, l_0 = 0.5, m = 1.0, t = 71.12
# ---------------------------------------------------
