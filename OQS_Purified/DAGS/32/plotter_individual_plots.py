import pickle
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import seaborn as sn

# Define the directory and observables
data_dir = "Plots/data/"
masses = [0.1, 0.25, 0.5, 0.75, 1.0]
observables = ["E", "EFE", "KE", "ME", "Q", "EF"]  # Added "Q" and "EF"
l0_values = [0.0, 0.5]
aD_values = [2.0, 5.0]
num_iters = 2000

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

# Iterate over each mass and observable to plot the data
for mass in masses:
    for observable in observables:
        
        if observable in ["Q", "EF"]:  # Handle heatmap cases with subplots
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))  # 2x2 grid of subplots
            
            # Flatten axes for easier iteration
            axes = axes.flatten()
            
            subplot_counter = 0  # Counter for subplot index
            
            # Iterate over all combinations of l0 and aD
            for l0 in l0_values:
                for aD in aD_values:
                    # Construct the file pattern
                    file_pattern = f"{observable}_{mass}_{l0}_{aD}.pickle"
                    file_path = os.path.join(data_dir, file_pattern)
                    
                    if os.path.exists(file_path):
                        # Load the data (assuming it's a 2D array: sites vs time)
                        values = load_pickle(file_path)
                        
                        fraction = 0.3
                        N = 12
                        tau = 0.05
                        steps = 2000
                        t_over_a_list = [0] + list(tau*(np.arange(1, steps+1)))
                        thermalization_time = t_over_a_list[get_index_of_reduced_value(fraction, values[N//2-1,:])]
                        
                        if aD == 2.0 and mass == 1.0 and observable == 'EF':
                            print(mass, aD, l0, thermalization_time)

                        # Plot the heatmap
                        im = axes[subplot_counter].imshow(values, aspect='auto', origin='lower')
                        axes[subplot_counter].set_title(f"l0={l0}, aD={aD}")
                        axes[subplot_counter].set_xlabel('Time')
                        axes[subplot_counter].set_ylabel('Site')
                        subplot_counter += 1
                    else:
                        print(f"File not found: {file_path}")
            
            # Adjust layout and add colorbar
            fig.suptitle(f"{observable} Heatmap for Mass {mass}", fontsize=16)
            fig.tight_layout(rect=[0, 0, 0.9, 0.96])
            cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Position for colorbar
            fig.colorbar(im, cax=cbar_ax)
            
            # Save the figure
            output_path = os.path.join("Plots/combined_plots_from_individual_observables/", f"{observable}_mass_{mass}_heatmap.png")
            plt.savefig(output_path)
            plt.close()  # Close the figure to avoid display if running in non-interactive environments
        
        else:  # Line plots for other observables
            plt.figure(figsize=(10, 6))  # Create a new figure for each mass and observable
            
            # Initialize a counter to apply different linestyles
            linestyle_counter = 0
            
            # Iterate over all combinations of l0 and aD
            for l0 in l0_values:
                for aD in aD_values:
                    # Construct the file pattern
                    file_pattern = f"{observable}_{mass}_{l0}_{aD}.pickle"
                    file_path = os.path.join(data_dir, file_pattern)
                    
                    if os.path.exists(file_path):
                        # Load the data
                        values = load_pickle(file_path)
                        
                        time = np.linspace(0, len(values), 2000)
                        
                        # Plot the data with a different linestyle
                        plt.plot(time, values, label=f"l0={l0}, aD={aD}", linestyle=linestyles[linestyle_counter % len(linestyles)])
                        
                        # Increment the counter
                        linestyle_counter += 1
                    else:
                        print(f"File not found: {file_path}")
            
            # Customize the plot
            plt.title(f"{observable} vs Time for Mass {mass}")
            plt.xlabel("Time")
            plt.ylabel(f"{observable}")
            plt.legend()
            plt.grid(True)
            
            # Save the figure
            output_path = os.path.join("Plots/combined_plots_from_individual_observables/", f"{observable}_mass_{mass}.png")
            plt.savefig(output_path)
            plt.close()  # Close the figure to avoid display if running in non-interactive environments
